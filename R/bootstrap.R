# ---- Bootstrap pathway uncertainty ----
#
# Methodologically built on Saqr, Tikka & López-Pernas (2025), `tna`
# package: extends the bootstrap-and-stability framework from edges to
# variable-depth pathways.
#
# Design rules:
# - The bootstrap operates on raw counts only. No `nmin` filter, no
#   smoothing inside the resampling loop. Pathways tracked are
#   exactly those in the original tree; their values in each
#   resample are read directly from the per-depth count matrices.
# - Stability: p_stability = bootstrap tail probability that the
#   chosen `stat` (default `count`) falls outside
#   [cr[1] x observed, cr[2] x observed], with a +1 correction.
#   stability_rate = 1 - uncorrected instability fraction is retained
#   as a descriptive companion. count = 0 (pathway absent in a
#   resample) automatically lands outside the band; no special
#   "appearance rate" concept needed.
# - Informativeness: informative_rate = fraction of resamples where
#   the pathway's empirical likelihood-ratio statistic G^2 against
#   its parent context exceeds the chi-square critical value at
#   level `alpha_g2` (df = |alphabet| - 1). Parallels the
#   parametric `prune_pathtree(criterion = "G2")` decision, but reproducibly
#   across resamples.
# - Performance: per-sequence (context, next_state) pair-index
#   vectors are precomputed once per depth. Each bootstrap
#   iteration is one `unlist`+`tabulate` per depth followed by a
#   per-pathway lookup; no `data.frame` construction inside the
#   loop.

#' @noRd
.pt_pathway_depth <- function(p) {
  is_root <- p == .ROOT_LABEL
  out <- integer(length(p))
  if (any(!is_root)) {
    parts <- strsplit(p[!is_root], " -> ", fixed = TRUE)
    out[!is_root] <- lengths(parts)
  }
  out
}

#' @noRd
.pt_precompute_for_boot <- function(trajs, max_depth, alphabet) {
  ## For each (sequence s, depth d) precompute the integer-encoded
  ## (context, next_state) pair vector. Keys are pair_idx =
  ## (ctx_idx - 1) * na + next_idx, where ctx_idx indexes the
  ## global context dictionary at depth d. Bootstrap aggregation is
  ## then a single `unlist(pair[idx])` + `tabulate` per depth —
  ## Nestimate's transition fast-path lifted from edges to
  ## pathways.
  na <- length(alphabet)
  pair_per_depth     <- vector("list", max_depth + 1L)
  ctx_dict_per_depth <- vector("list", max_depth + 1L)
  for (d in seq.int(0L, max_depth)) {
    if (d == 0L) {
      ctx_per_seq <- lapply(trajs,
                            function(traj) rep(.ROOT, length(traj)))
      nxt_per_seq <- lapply(trajs,
                            function(traj) match(traj, alphabet))
    } else {
      ctx_per_seq <- lapply(trajs, function(traj) {
        L <- length(traj)
        if (L < d + 1L) return(character(0))
        starts <- seq_len(L - d)
        vapply(starts, function(t)
          paste(traj[t:(t + d - 1L)], collapse = " -> "),
          character(1))
      })
      nxt_per_seq <- lapply(trajs, function(traj) {
        L <- length(traj)
        if (L < d + 1L) return(integer(0))
        starts <- seq_len(L - d)
        match(traj[starts + d], alphabet)
      })
    }
    all_ctx <- unique(unlist(ctx_per_seq, use.names = FALSE))
    ctx_dict_per_depth[[d + 1L]] <- all_ctx

    pair_per_seq <- mapply(function(ctx_v, nxt_v) {
      keep <- !is.na(nxt_v) & !is.na(ctx_v)
      if (!any(keep)) return(integer(0))
      ctx_v <- ctx_v[keep]; nxt_v <- nxt_v[keep]
      ctx_idx <- match(ctx_v, all_ctx)
      (ctx_idx - 1L) * na + nxt_v
    }, ctx_per_seq, nxt_per_seq, SIMPLIFY = FALSE)
    pair_per_depth[[d + 1L]] <- pair_per_seq
  }
  list(pair      = pair_per_depth,
       ctx_dict  = ctx_dict_per_depth,
       alphabet  = alphabet,
       max_depth = max_depth,
       na        = na)
}

#' @noRd
.pt_resolve_pathway_layout <- function(pathways, precomp) {
  ## Resolve each original-tree pathway to (depth, ctx_idx in dict)
  ## and its parent's (depth, ctx_idx). Returns NA for the root's
  ## parent. Done once before the bootstrap loop.
  na <- precomp$na
  P <- length(pathways)
  pw_depth        <- integer(P)
  pw_ctx_idx      <- integer(P)
  parent_depth    <- integer(P)
  parent_ctx_idx  <- integer(P)
  for (j in seq_len(P)) {
    p <- pathways[[j]]
    if (identical(p, .ROOT)) {
      pw_depth[j]       <- 0L
      pw_ctx_idx[j]     <- match(.ROOT, precomp$ctx_dict[[1L]])
      parent_depth[j]   <- NA_integer_
      parent_ctx_idx[j] <- NA_integer_
    } else {
      d <- length(strsplit(p, " -> ", fixed = TRUE)[[1L]])
      pw_depth[j]       <- d
      pw_ctx_idx[j]     <- match(p, precomp$ctx_dict[[d + 1L]])
      parent_depth[j]   <- d - 1L
      parent_ctx_idx[j] <- match(.pt_parent_ctx(p),
                                 precomp$ctx_dict[[d]])
    }
  }
  list(pw_depth = pw_depth, pw_ctx_idx = pw_ctx_idx,
       parent_depth = parent_depth,
       parent_ctx_idx = parent_ctx_idx)
}

#' @noRd
.pt_resample_count_matrices <- function(precomp, idx) {
  ## Return a list of n_ctx_d x na count matrices, one per depth.
  na <- precomp$na
  max_depth <- precomp$max_depth
  out <- vector("list", max_depth + 1L)
  for (d in seq.int(0L, max_depth)) {
    n_ctx_d <- length(precomp$ctx_dict[[d + 1L]])
    if (n_ctx_d == 0L) {
      out[[d + 1L]] <- matrix(0L, 0L, na)
      next
    }
    pooled <- unlist(precomp$pair[[d + 1L]][idx],
                     use.names = FALSE)
    flat <- tabulate(pooled, nbins = n_ctx_d * na)
    out[[d + 1L]] <- matrix(flat, n_ctx_d, na, byrow = TRUE)
  }
  out
}

#' @noRd
.pt_pathway_stats_from_counts <- function(cnt_by_depth, layout) {
  ## Empirical pathway stats from per-depth count matrices. No
  ## smoothing. Returns parallel numeric vectors over the layout's
  ## pathways.
  P <- length(layout$pw_depth)
  count_total <- integer(P)
  prob_next   <- rep(NA_real_,    P)
  modal       <- rep(NA_integer_, P)
  KL          <- rep(NA_real_,    P)
  G2          <- rep(NA_real_,    P)
  flips       <- rep(NA_integer_, P)

  for (j in seq_len(P)) {
    cnt_vec <- cnt_by_depth[[layout$pw_depth[j] + 1L]][
      layout$pw_ctx_idx[j], ]
    tot <- sum(cnt_vec)
    count_total[j] <- as.integer(tot)
    if (tot == 0L) next
    m <- which.max(cnt_vec)[1L]
    modal[j]     <- m
    prob_next[j] <- cnt_vec[m] / tot

    pd  <- layout$parent_depth[j]
    pix <- layout$parent_ctx_idx[j]
    if (is.na(pd) || is.na(pix)) next
    pcnt <- cnt_by_depth[[pd + 1L]][pix, ]
    ptot <- sum(pcnt)
    if (ptot == 0L) next
    pprob <- pcnt / ptot
    flips[j] <- as.integer(m != which.max(pcnt)[1L])

    p <- cnt_vec / tot
    msk <- p > 0
    if (any(p > 0 & pprob == 0)) {
      KL[j] <- Inf
    } else {
      KL[j] <- sum(p[msk] * log(p[msk] / pprob[msk], base = 2))
    }

    G2[j] <- .ct_g2(cnt_vec, pprob)
  }
  list(count = count_total, prob_next = prob_next, modal = modal,
       KL = KL, G2 = G2, flips = flips)
}

#' Bootstrap Pathway Stability and Informativeness
#'
#' @description
#' Non-parametric sequence bootstrap for a fitted \code{pathtree}.
#' Methodologically built on Saqr, Tikka & López-Pernas (2025),
#' \code{tna}, extending the edge-level bootstrap framework to
#' variable-depth pathways.
#'
#' The bootstrap tracks every pathway in the original tree. Each
#' iteration resamples whole sequences with replacement, aggregates
#' raw counts per depth, and reads each pathway's count vector
#' directly from the resample. \strong{No smoothing, no \code{nmin}
#' filter, no extra parameters} inside the loop — the bootstrap
#' operates on counts the same way \code{tna::bootstrap} operates on
#' edge weights.
#'
#' Two complementary measures are reported per pathway:
#' \describe{
#'   \item{\code{p_stability}}{Bootstrap-estimated probability that
#'     the chosen \code{stat} (default \code{count}) falls outside
#'     \code{[cr[1] * observed, cr[2] * observed]}, with a +1
#'     correction. This is a stability p-value: small values mean the
#'     pathway rarely fails the chosen reproducibility criterion under
#'     sequence-level resampling.}
#'   \item{\code{stability_rate}}{Uncorrected descriptive companion:
#'     the fraction of resamples where the chosen \code{stat} lies
#'     inside the consistency band. A pathway whose count drops to zero
#'     in a resample fails the band test automatically.}
#'   \item{\code{informative_rate}}{Fraction of resamples where the
#'     pathway's empirical \eqn{G^2} likelihood-ratio statistic
#'     against its parent context exceeds the chi-square critical
#'     value at level \code{alpha_g2} (df = \code{|alphabet| - 1}).
#'     Tests \emph{reproducibly significant divergence from the
#'     shorter-history baseline}.}
#' }
#'
#' Read \code{stable} and \code{informative} together:
#' \itemize{
#'   \item \code{stable && informative}: reproducible and predictively
#'     distinctive pathway.
#'   \item \code{stable && !informative}: reproducible pathway count /
#'     statistic, but not predictively distinctive from its parent.
#'   \item \code{!stable && informative}: sharp or divergent pathway
#'     carried by an unstable subset of sequences.
#'   \item \code{!stable && !informative}: weak or sample-fragile
#'     pathway.
#' }
#'
#' @param tree A fitted \code{pathtree} carrying \code{tree$data}.
#' @param iter Integer. Number of bootstrap iterations.
#'   Default \code{1000}.
#' @param stat Character. Pathway statistic on which
#'   \code{p_stability} is measured. One of \code{"count"}
#'   (default; tna's edge-weight analogue), \code{"next_probability"}
#'   (the most-likely next-state probability), or \code{"divergence"}.
#' @param consistency_range Numeric vector of length 2.
#'   Multiplicative tolerance band around the observed value.
#'   Default \code{c(0.75, 1.25)}.
#' @param stability_threshold Numeric in \eqn{(0, 1)}. Backward-
#'   compatible stability-rate threshold. A pathway is
#'   \code{stable = TRUE} when
#'   \code{p_stability < 1 - stability_threshold}. Default
#'   \code{0.75}, equivalent to a 0.25 instability tolerance with the
#'   bootstrap +1 correction.
#' @param informative_threshold Numeric in \eqn{(0, 1)}. A pathway
#'   is \code{informative = TRUE} when \code{informative_rate >=
#'   informative_threshold}. Default \code{0.95}: at the standard
#'   \eqn{G^2} significance level (\code{alpha = 0.05}) the
#'   chi-square test has a 5\% Type-I rate, so requiring 95\% of
#'   resamples to clear the critical value rules out finite-sample
#'   chance deviations from the parent distribution.
#' @param alpha Numeric in \eqn{(0, 1)}. Significance level for
#'   the \eqn{G^2} test against parent. Default \code{0.05}.
#' @param ci_level Numeric in \eqn{(0, 1)}. Tail probability for the
#'   bootstrap CIs on \code{count}, \code{next_probability},
#'   \code{divergence}, \code{G2}. Default \code{0.05} (95\% CI).
#' @param seed Integer or \code{NULL}. RNG seed.
#'   Default \code{1L}.
#' @param keep_resamples Logical. If \code{TRUE} (default), the
#'   per-iteration resample matrices \code{M_count},
#'   \code{M_next_probability}, \code{M_divergence}, \code{M_G2},
#'   \code{M_changes_prediction} are retained on the
#'   returned object. Set to \code{FALSE} to drop them (each is
#'   \code{iter x n_pathways}) when memory matters; the summary table
#'   is unaffected.
#' @param progress Logical. Show a progress bar.
#'   Default \code{FALSE}.
#'
#' @return A \code{pathtree_bootstrap} object: a list with
#'   \describe{
#'     \item{summary}{Per-pathway data.frame, sorted so that
#'       \code{stable & informative} pathways come first then by
#'       \code{stability_rate} descending.}
#'     \item{pathways_orig}{Empirical original-pathway statistics
#'       (no smoothing) as a tidy data.frame.}
#'     \item{M_count, M_next_probability, M_divergence, M_G2,
#'       M_changes_prediction}{Raw resample matrices:
#'       \code{iter x n_pathways}, columns named by pathway.}
#'     \item{iter, stat, consistency_range, stability_threshold,
#'       informative_threshold, alpha_g2, level, seed,
#'       g2_critical_value}{Configuration.}
#'   }
#'
#' @references
#' Saqr, M., Tikka, S., & López-Pernas, S. (2025). Transition
#' Network Analysis. \emph{LAK '25},
#' doi:10.1145/3706468.3706513.
#' @export
bootstrap_pathways <- function(tree,
                               iter                  = 1000L,
                               stat                  = c("count",
                                                         "next_probability",
                                                         "divergence"),
                               consistency_range     = c(0.75, 1.25),
                               stability_threshold   = 0.75,
                               informative_threshold = 0.95,
                               alpha                 = 0.05,
                               ci_level              = 0.05,
                               seed                  = 1L,
                               keep_resamples        = TRUE,
                               progress              = FALSE) {
  stopifnot(inherits(tree, "pathtree"))
  if (is.null(tree$data))
    stop("tree$data is missing; bootstrap requires the original ",
         "sequences. Refit with context_tree() (>=0.1.1 keeps data).",
         call. = FALSE)
  iter <- as.integer(iter)
  stopifnot(iter >= 2L,
            ci_level > 0, ci_level < 1,
            alpha > 0, alpha < 1,
            stability_threshold > 0, stability_threshold < 1,
            informative_threshold > 0, informative_threshold < 1,
            length(consistency_range) == 2L,
            consistency_range[1] > 0, consistency_range[2] > 0,
            consistency_range[1] < consistency_range[2])
  stat <- match.arg(stat)
  if (!is.null(seed)) set.seed(as.integer(seed))

  trajs    <- tree$data
  n        <- length(trajs)
  alphabet <- tree$alphabet
  na       <- length(alphabet)
  max_depth <- tree$max_depth
  g2_crit   <- stats::qchisq(1 - alpha, df = na - 1L)

  precomp <- .pt_precompute_for_boot(trajs, max_depth, alphabet)

  ## Pathways from the original tree (raw — no resample dependency).
  orig_path_keys <- names(tree$nodes)
  if (length(orig_path_keys) == 0L)
    stop("Original tree has no pathways. Did fitting succeed?",
         call. = FALSE)
  layout <- .pt_resolve_pathway_layout(orig_path_keys, precomp)
  P <- length(orig_path_keys)
  pathway_label <- ifelse(orig_path_keys == .ROOT, .ROOT_LABEL,
                          orig_path_keys)

  ## Original pathway stats — recompute empirically (no smoothing)
  ## so original and bootstrap values are on the same scale.
  full_idx  <- seq_len(n)
  cnt_full  <- .pt_resample_count_matrices(precomp, full_idx)
  orig_stat <- .pt_pathway_stats_from_counts(cnt_full, layout)
  orig_modal_next <- ifelse(is.na(orig_stat$modal), NA_character_,
                            alphabet[orig_stat$modal])
  orig_flips      <- ifelse(is.na(orig_stat$flips), NA, as.logical(orig_stat$flips))

  ## Bootstrap loop.
  M_count <- matrix(0L,           iter, P,
                    dimnames = list(NULL, pathway_label))
  M_prob  <- matrix(NA_real_,     iter, P,
                    dimnames = list(NULL, pathway_label))
  M_KL    <- matrix(NA_real_,     iter, P,
                    dimnames = list(NULL, pathway_label))
  M_G2    <- matrix(NA_real_,     iter, P,
                    dimnames = list(NULL, pathway_label))
  M_flips <- matrix(NA_integer_,  iter, P,
                    dimnames = list(NULL, pathway_label))

  if (progress) pb <- utils::txtProgressBar(min = 0, max = iter,
                                            style = 3)
  for (b in seq_len(iter)) {
    idx <- sample.int(n, n, replace = TRUE)
    cnt_b <- .pt_resample_count_matrices(precomp, idx)
    s     <- .pt_pathway_stats_from_counts(cnt_b, layout)
    M_count[b, ] <- s$count
    M_prob[b, ]  <- s$prob_next
    M_KL[b, ]    <- s$KL
    M_G2[b, ]    <- s$G2
    M_flips[b, ] <- s$flips
    if (progress) utils::setTxtProgressBar(pb, b)
  }
  if (progress) close(pb)

  ## CIs (ci_level = 0.05 -> 95%).
  q_probs    <- c(ci_level / 2, 1 - ci_level / 2)
  col_quants <- function(M) apply(M, 2L, stats::quantile,
                                  probs = q_probs, na.rm = TRUE)
  ci_count <- col_quants(M_count)
  ci_prob  <- col_quants(M_prob)
  ci_KL    <- col_quants(M_KL)
  ci_G2    <- col_quants(M_G2)

  mean_count <- colMeans(M_count)
  sd_count   <- apply(M_count,    2L, stats::sd)
  mean_prob  <- colMeans(M_prob,  na.rm = TRUE)
  sd_prob    <- apply(M_prob,     2L, stats::sd, na.rm = TRUE)
  mean_KL    <- colMeans(M_KL,    na.rm = TRUE)
  sd_KL      <- apply(M_KL,       2L, stats::sd, na.rm = TRUE)
  mean_G2    <- colMeans(M_G2,    na.rm = TRUE)
  sd_G2      <- apply(M_G2,       2L, stats::sd, na.rm = TRUE)

  ## Stability rate on the chosen stat (band test).
  M_stat_arr  <- switch(stat,
                        count            = M_count,
                        next_probability = M_prob,
                        divergence       = M_KL)
  orig_stat_v <- switch(stat,
                        count            = orig_stat$count,
                        next_probability = orig_stat$prob_next,
                        divergence       = orig_stat$KL)
  cr_lo <- pmin(orig_stat_v * consistency_range[1],
                orig_stat_v * consistency_range[2])
  cr_hi <- pmax(orig_stat_v * consistency_range[1],
                orig_stat_v * consistency_range[2])
  n_in <- vapply(seq_len(P), function(j) {
    if (is.na(orig_stat_v[j])) return(NA_integer_)
    sum(M_stat_arr[, j] >= cr_lo[j] & M_stat_arr[, j] <= cr_hi[j],
        na.rm = TRUE)
  }, integer(1L))
  stability_rate <- n_in / iter
  n_out <- iter - n_in
  p_stability <- (n_out + 1) / (iter + 1)
  p_stability[is.na(n_in)] <- NA_real_
  stable <- !is.na(p_stability) &
            p_stability < (1 - stability_threshold)

  ## Informative rate: fraction of resamples where empirical G^2
  ## exceeds the chi-square critical value (df = na - 1).
  n_inf <- vapply(seq_len(P), function(j) {
    g <- M_G2[, j]
    valid <- !is.na(g) & is.finite(g)
    if (!any(valid)) return(NA_integer_)
    sum(g[valid] > g2_crit)
  }, integer(1L))
  informative_rate <- n_inf / iter
  informative_rate[is.na(orig_stat$G2)] <- NA_real_
  informative <- !is.na(informative_rate) &
                 informative_rate >= informative_threshold

  ## Modal-flip consistency vs the original empirical modal-flip
  ## flag. Descriptive companion to stability/informativeness.
  flip_consistency <- vapply(seq_len(P), function(j) {
    of <- orig_stat$flips[j]
    if (is.na(of)) return(NA_real_)
    fb <- M_flips[, j]
    valid <- !is.na(fb)
    if (!any(valid)) return(NA_real_)
    mean(fb[valid] == as.integer(of))
  }, numeric(1L))

  summary_df <- data.frame(
    pathway              = pathway_label,
    depth                = layout$pw_depth,
    count                = orig_stat$count,
    likely_next          = orig_modal_next,
    next_probability     = orig_stat$prob_next,
    divergence           = orig_stat$KL,
    changes_prediction   = orig_flips,
    G2                   = orig_stat$G2,
    p_stability          = p_stability,
    stability_rate       = stability_rate,
    stable               = stable,
    informative_rate     = informative_rate,
    informative          = informative,
    flip_consistency     = flip_consistency,
    mean_count           = mean_count,
    sd_count             = sd_count,
    ci_count_lo          = ci_count[1L, ],
    ci_count_hi          = ci_count[2L, ],
    mean_next_probability   = mean_prob,
    sd_next_probability     = sd_prob,
    ci_next_probability_lo  = ci_prob[1L, ],
    ci_next_probability_hi  = ci_prob[2L, ],
    mean_divergence      = mean_KL,
    sd_divergence        = sd_KL,
    ci_divergence_lo     = ci_KL[1L, ],
    ci_divergence_hi     = ci_KL[2L, ],
    mean_G2              = mean_G2,
    sd_G2                = sd_G2,
    ci_G2_lo             = ci_G2[1L, ],
    ci_G2_hi             = ci_G2[2L, ],
    stringsAsFactors     = FALSE
  )
  ## Sort key: 2 (both flags), 1 (one flag), 0 (neither).
  ## Pathways with both stable + informative come first; ties
  ## break on stability_rate then count.
  trust_score <- as.integer(summary_df$stable) +
                 as.integer(summary_df$informative)
  trust_score[is.na(trust_score)] <- 0L
  ord <- order(-trust_score,
               -ifelse(is.na(summary_df$stability_rate), 0,
                       summary_df$stability_rate),
               -ifelse(is.na(summary_df$count), 0,
                       summary_df$count))
  summary_df <- summary_df[ord, , drop = FALSE]
  rownames(summary_df) <- NULL

  pathways_orig <- data.frame(
    pathway            = pathway_label,
    depth              = layout$pw_depth,
    count              = orig_stat$count,
    likely_next        = orig_modal_next,
    next_probability   = orig_stat$prob_next,
    divergence         = orig_stat$KL,
    changes_prediction = orig_flips,
    G2                 = orig_stat$G2,
    stringsAsFactors   = FALSE
  )

  if (!isTRUE(keep_resamples)) {
    M_count <- NULL; M_prob <- NULL; M_KL <- NULL
    M_G2    <- NULL; M_flips <- NULL
  }

  structure(
    list(
      pathways_orig         = pathways_orig,
      summary               = summary_df,
      M_count               = M_count,
      M_next_probability    = M_prob,
      M_divergence          = M_KL,
      M_G2                  = M_G2,
      M_changes_prediction  = M_flips,
      iter                  = iter,
      stat                  = stat,
      consistency_range     = consistency_range,
      stability_threshold   = stability_threshold,
      informative_threshold = informative_threshold,
      alpha_g2              = alpha,
      ci_level              = ci_level,
      g2_critical_value     = g2_crit,
      seed                  = seed
    ),
    class = "pathtree_bootstrap"
  )
}

#' @export
print.pathtree_bootstrap <- function(x, n = 10L, digits = 3L, ...) {
  cat(sprintf(
    "<pathtree_bootstrap>  %d resamples\n", x$iter))
  cat(sprintf(
    "  stability  : %s in [%.2f, %.2f] x observed, p < %.2f\n",
    x$stat, x$consistency_range[1L], x$consistency_range[2L],
    1 - x$stability_threshold))
  cat(sprintf(
    "  informative: G^2 > qchisq(%.2f, df=k-1) = %.2f, threshold %.2f\n",
    1 - x$alpha_g2, x$g2_critical_value,
    x$informative_threshold))
  s <- x$summary
  n_total <- nrow(s)
  n_stable <- sum(s$stable, na.rm = TRUE)
  n_info   <- sum(s$informative, na.rm = TRUE)
  n_both   <- sum(s$stable & s$informative, na.rm = TRUE)
  cat(sprintf(
    "  pathways   : %d total, %d stable, %d informative, %d both\n",
    n_total, n_stable, n_info, n_both))
  cat("\ntop pathways (stable + informative first):\n")
  show <- utils::head(s, n)
  show_disp <- show[, c("pathway", "depth", "count",
                        "p_stability", "stability_rate", "stable",
                        "informative_rate", "informative",
                        "mean_G2", "ci_G2_lo", "ci_G2_hi")]
  num_cols <- vapply(show_disp, is.numeric, logical(1))
  show_disp[, num_cols] <- lapply(show_disp[, num_cols], round,
                                  digits = digits)
  print.data.frame(show_disp, row.names = FALSE)
  if (n_total > n) {
    cat(sprintf(
      "# ... %d more pathways (use summary(x) for full table)\n",
      n_total - n))
  }
  invisible(x)
}

#' @export
summary.pathtree_bootstrap <- function(object, ...) {
  object$summary
}

#' Coerce a Pathtree Bootstrap to a Tidy Data Frame
#'
#' @description
#' Uniform tidy-extract: returns the per-pathway summary table
#' (\code{object$summary}), so \code{as.data.frame(boot)} and
#' \code{summary(boot)} are interchangeable extractors.
#'
#' @param x A \code{pathtree_bootstrap}.
#' @param row.names,optional Ignored.
#' @param ... Ignored.
#' @return A data.frame; see \code{\link{bootstrap_pathways}} for the
#'   full column vocabulary.
#' @export
as.data.frame.pathtree_bootstrap <- function(x, row.names = NULL,
                                              optional = FALSE, ...) {
  x$summary
}

#' Plot a Pathtree Bootstrap
#'
#' @description
#' Forest plot of per-pathway \eqn{G^2} (likelihood-ratio against
#' parent) with bootstrap 95\% CI bars. Pathways are ordered with
#' \emph{stable & informative} ones first, then by stability rate.
#' The chi-square critical value at \code{alpha_g2} is shown as a
#' dashed reference line: pathways whose CI lies entirely above it
#' are reproducibly informative.
#'
#' @param x A \code{pathtree_bootstrap} object.
#' @param top Integer. Maximum pathways to show. Default 25.
#' @param min_stability Numeric. Minimum stability_rate to display.
#'   Default \code{NULL} (use \code{x$stability_threshold}).
#' @param ... Ignored.
#' @return A ggplot object.
#' @export
plot.pathtree_bootstrap <- function(x, top = 25L,
                                    min_stability = NULL, ...) {
  if (is.null(min_stability)) min_stability <- x$stability_threshold
  s <- x$summary
  s <- s[!is.na(s$stability_rate) &
         s$stability_rate >= min_stability &
         !is.na(s$mean_G2), , drop = FALSE]
  trust <- as.integer(s$stable) + as.integer(s$informative)
  trust[is.na(trust)] <- 0L
  s <- s[order(-trust, -s$stability_rate, -s$count), ,
         drop = FALSE]
  s <- utils::head(s, top)
  if (nrow(s) == 0L)
    stop("No pathways meet min_stability = ", min_stability, ".",
         call. = FALSE)
  s$pathway <- factor(s$pathway, levels = rev(s$pathway))
  s$status <- factor(
    ifelse(s$stable & s$informative,        "stable + informative",
    ifelse(s$stable & !s$informative,        "stable only",
    ifelse(!s$stable & s$informative,        "informative only",
                                              "neither"))),
    levels = c("stable + informative", "stable only",
               "informative only", "neither"))
  ggplot2::ggplot(s,
    ggplot2::aes(x = .data$mean_G2, y = .data$pathway,
                 colour = .data$status)) +
    ggplot2::geom_vline(xintercept = x$g2_critical_value,
                        linetype = "dashed", colour = "grey50") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$ci_G2_lo, xmax = .data$ci_G2_hi),
      height = 0.25, linewidth = 0.6) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::scale_colour_manual(
      values = c(`stable + informative` = "#0072B2",
                 `stable only`          = "#888888",
                 `informative only`     = "#D55E00",
                 `neither`              = "#cccccc"),
      drop = FALSE, name = NULL) +
    ggplot2::labs(
      x = sprintf("bootstrap G^2 (mean and %d%% CI)",
                  round((1 - x$ci_level) * 100)),
      y = NULL,
      title = sprintf("Pathway bootstrap (%d resamples)", x$iter),
      subtitle = sprintf(
        "showing %d / %d pathways with stability_rate >= %.2f  |  dashed line: G^2 critical value at alpha = %.2f",
        nrow(s), nrow(x$summary), min_stability,
        x$alpha_g2)) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40",
                                             size = 9),
      legend.position = "bottom",
      panel.grid.major.y = ggplot2::element_blank())
}

#' Plot Bootstrap Resample Distributions per Pathway
#'
#' @description
#' Faceted histogram of the bootstrap resample values for a chosen
#' pathway statistic, one panel per pathway.
#'
#' @param x A \code{pathtree_bootstrap} object.
#' @param pathways Character vector of pathway names. \code{NULL}
#'   (default) picks the top \code{top} pathways from the summary.
#' @param stat Character. One of \code{"count"} (default),
#'   \code{"next_probability"}, \code{"divergence"}, \code{"G2"}.
#' @param top Integer. Default 6.
#' @param bins Integer. Histogram bins. Default 30.
#' @return A ggplot object.
#' @export
plot_pathway_resamples <- function(x, pathways = NULL,
                                   stat = c("count", "next_probability",
                                            "divergence", "G2"),
                                   top = 6L, bins = 30L) {
  stopifnot(inherits(x, "pathtree_bootstrap"))
  stat <- match.arg(stat)
  M <- switch(stat,
              count            = x$M_count,
              next_probability = x$M_next_probability,
              divergence       = x$M_divergence,
              G2               = x$M_G2)
  if (is.null(M))
    stop("Resample matrices were not retained on this bootstrap ",
         "object. Re-run with keep_resamples = TRUE.", call. = FALSE)
  if (is.null(pathways)) {
    pathways <- utils::head(x$summary$pathway, top)
  }
  missing_pw <- setdiff(pathways, colnames(M))
  if (length(missing_pw) > 0L)
    stop("Unknown pathway(s): ",
         paste(missing_pw, collapse = ", "), call. = FALSE)
  rows <- lapply(pathways, function(p) {
    data.frame(pathway = p, value = M[, p],
               stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  df <- df[!is.na(df$value), , drop = FALSE]
  df$pathway <- factor(df$pathway, levels = pathways)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_histogram(bins = bins, fill = "#0072B2",
                            colour = "white", alpha = 0.85) +
    ggplot2::facet_wrap(~ pathway, scales = "free_y") +
    ggplot2::labs(
      x = sprintf("bootstrap %s", stat), y = "resamples",
      title = sprintf("Bootstrap distribution of %s per pathway",
                       stat),
      subtitle = sprintf("%d resamples; %d pathways shown",
                         x$iter, length(pathways))) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40",
                                             size = 9),
      strip.text    = ggplot2::element_text(face = "bold"))
}
