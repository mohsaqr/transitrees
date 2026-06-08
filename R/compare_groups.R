# ---- Multi-group comparison ----

#' @noRd
.cg_depth <- function(ctx) {
  if (identical(ctx, .ROOT)) return(0L)
  length(strsplit(ctx, " -> ", fixed = TRUE)[[1L]])
}

#' Shannon entropy in bits of a probability vector.
#' @noRd
.cg_entropy_bits <- function(p) {
  p <- p[p > 0]
  if (length(p) == 0L) return(0)
  -sum(p * log2(p))
}

#' Count-weighted Jensen-Shannon divergence (bits) across the rows of a
#' K x |alphabet| count matrix. Groups with zero total at this context
#' are dropped; returns 0 when fewer than two groups are present.
#' @noRd
.cg_jsd_bits <- function(counts_mat) {
  tot  <- rowSums(counts_mat)
  keep <- tot > 0
  if (sum(keep) < 2L) return(0)
  cm <- counts_mat[keep, , drop = FALSE]
  tt <- tot[keep]
  w  <- tt / sum(tt)
  p  <- cm / tt                                   # row-wise MLE
  m  <- colSums(w * p)                            # mixture
  h_each <- sum(w * apply(p, 1L, .cg_entropy_bits))
  max(.cg_entropy_bits(m) - h_each, 0)
}

#' Expected counts under independence for a contingency table (the
#' margin-product / grand-total). Shared by the usage G^2 and the
#' Pearson-residual difference map so they cannot drift apart.
#' @noRd
.cg_expected <- function(tab) outer(rowSums(tab), colSums(tab)) / sum(tab)

#' G^2 homogeneity statistic for a pathway's prevalence across groups:
#' the K x 2 table [count_c,g ; opportunities_g - count_c,g].
#' @noRd
.cg_usage_g2 <- function(count_c, n_g) {
  tab <- cbind(count_c, pmax(n_g - count_c, 0))
  if (sum(tab) == 0) return(0)
  e  <- .cg_expected(tab)
  nz <- tab > 0 & e > 0
  2 * sum(tab[nz] * log(tab[nz] / e[nz]))
}

#' Per-context usage G^2 for one group assignment. The denominator for a
#' context is the number of windows at \emph{its own depth} (the sum of
#' all same-depth contexts' counts per group), not the total token count
#' --- so prevalence is "share of the depth-d windows" and is comparable
#' across groups regardless of mean sequence length.
#' @noRd
.cg_usage_vec <- function(cc, depths) {
  K        <- nrow(cc[[1L]])
  d_levels <- sort(unique(depths))
  tot_by_d <- lapply(d_levels, function(d)
    Reduce(`+`, lapply(cc[depths == d], rowSums), numeric(K)))
  names(tot_by_d) <- as.character(d_levels)
  vapply(seq_along(cc), function(i)
    .cg_usage_g2(rowSums(cc[[i]]), tot_by_d[[as.character(depths[i])]]),
    numeric(1))
}

#' Per-context K x |alphabet| count matrices for one group assignment.
#' @noRd
.cg_context_counts <- function(trajs_by_grp, grp_names, contexts, depths,
                               alphabet) {
  k      <- length(alphabet)
  uniq_d <- sort(unique(depths))
  ## one count table per group per needed depth (the same routine that
  ## builds the fitted node counts, so observed and permuted agree)
  per_grp <- lapply(grp_names, function(g) {
    tg <- trajs_by_grp[[g]]
    if (is.null(tg) || length(tg) == 0L) tg <- list(character(0))
    setNames(lapply(uniq_d, function(d)
      .ct_count_table(tg, depth = d, alphabet = alphabet)),
      as.character(uniq_d))
  })
  names(per_grp) <- grp_names
  out <- lapply(seq_along(contexts), function(i) {
    ctx <- contexts[i]; d <- as.character(depths[i])
    mat <- t(vapply(grp_names, function(g) {
      cv <- per_grp[[g]][[d]][[ctx]]
      if (is.null(cv)) rep(0, k) else cv
    }, numeric(k)))
    rownames(mat) <- grp_names
    mat
  })
  names(out) <- contexts
  out
}

#' Per-context statistics (JSD, usage G^2) and the two omnibus totals
#' for one group assignment.
#' @noRd
.cg_stats <- function(labels, pooled, grp_names, contexts, depths,
                      alphabet, root_ctx) {
  tbg <- split(pooled, factor(labels, levels = grp_names))
  cc  <- .cg_context_counts(tbg, grp_names, contexts, depths, alphabet)
  jsd <- vapply(cc, .cg_jsd_bits, numeric(1))
  use <- .cg_usage_vec(cc, depths)
  n_c <- vapply(cc, sum, numeric(1))
  list(jsd = jsd, usage = use,
       omni_beh = sum(n_c * jsd), omni_use = sum(use))
}

#' Compare Groups of Sequences for Structural Differences
#'
#' @description
#' Test how a set of fitted group trees (a \code{transitiontrees_group} from
#' \code{context_tree(..., group =)}) differ, on two complementary axes:
#' \describe{
#'   \item{behavioral}{given the \emph{same} context, do the groups
#'     predict a different next state? Measured per context by the
#'     count-weighted Jensen-Shannon divergence (bits) across the
#'     groups' next-state distributions.}
#'   \item{usage}{do the groups \emph{reach} a context at different
#'     rates? Measured per context by a \eqn{G^2} homogeneity statistic
#'     on its prevalence (its share of each group's positions).}
#' }
#' Significance is assessed by a label-permutation null throughout
#' (per-pathway and omnibus), with Benjamini-Hochberg FDR on the
#' per-pathway p-values.
#'
#' @details
#' The permutation pools every sequence, shuffles the group labels
#' (preserving group sizes), and recomputes the statistics from raw
#' counts using the same counting routine as the fit. The tested
#' context set is the union of the contexts the groups' trees actually
#' represent. For two groups the behavioral measure is JSD, which is
#' \strong{not} the symmetric-KL distance used by
#' \code{\link{compare_trees}()}; the \code{distance_matrix} component
#' does use \code{\link{tree_distance}()} (symmetric KL) for
#' consistency with the pairwise function.
#'
#' @param group A \code{transitiontrees_group}.
#' @param iter Integer. Number of label permutations. Default 999.
#' @param min_count Integer. Drop contexts whose total count across all
#'   groups is below this. Default 1.
#' @param seed Integer or \code{NULL}. RNG seed. Default 1.
#' @param block Block ids for a \strong{stratified} permutation: group
#'   labels are shuffled only \emph{within} each block, so the null
#'   respects nested / repeated-measures structure (e.g. several
#'   sequences from one subject) and holds any between-block difference
#'   fixed. Normally you do not pass this --- fit with
#'   \code{context_tree(..., block = )} and it is carried on the object
#'   and used automatically. Passing a vector here (one id per sequence,
#'   in pooled group-then-row order) overrides that. \code{NULL} with no
#'   stored block shuffles labels freely.
#'
#' @return A \code{transitiontrees_group_comparison}: a list with
#'   \describe{
#'     \item{pathways}{Per-context data.frame sorted by \code{jsd_bits}
#'       descending, with columns \code{pathway}, \code{depth},
#'       \code{count_total}, one \code{count_<group>} and one
#'       \code{modal_<group>} column per group (most likely next state,
#'       ties broken by alphabet order), \code{flips} (do the groups'
#'       modal next states disagree?), \code{jsd_bits}, \code{jsd_p},
#'       \code{jsd_padj}, \code{usage_g2}, \code{usage_p},
#'       \code{usage_padj}. \code{usage_*} is \code{NA} for the root,
#'       which has no prevalence test.}
#'     \item{omnibus}{Two-row data.frame: the behavioral and usage
#'       global statistics with permutation p-values.}
#'     \item{distance_matrix}{K x K symmetric-KL distance matrix between
#'       the groups (from \code{tree_distance()}).}
#'     \item{groups, iter, seed, n_contexts}{Configuration.}
#'   }
#'
#' @examples
#' \donttest{
#' gx <- replicate(40, sample(c("A","B","C"), 8, replace = TRUE,
#'                            prob = c(.2,.6,.2)), simplify = FALSE)
#' gy <- replicate(40, sample(c("A","B","C"), 8, replace = TRUE,
#'                            prob = c(.2,.2,.6)), simplify = FALSE)
#' grp <- context_tree(c(gx, gy), group = rep(c("x","y"), each = 40),
#'                     max_depth = 1L)
#' cmp <- compare_groups(grp, iter = 199L)
#' cmp
#' }
#' @seealso \code{\link{compare_trees}} for the pairwise permutation test.
#' @export
compare_groups <- function(group, iter = 999L, min_count = 1L, seed = 1L,
                           block = NULL) {
  stopifnot(inherits(group, "transitiontrees_group"))
  grp_names <- names(group)
  K <- length(grp_names)
  if (K < 2L) stop("'group' must contain at least two groups.", call. = FALSE)
  for (m in group) .pt_assert_unweighted(m, "compare_groups()")
  alphabet <- group[[1L]]$alphabet
  if (!is.null(seed)) set.seed(seed)

  trajs_by_grp <- lapply(group, function(m) .ct_traj(m$data))
  names(trajs_by_grp) <- grp_names
  pooled <- unlist(trajs_by_grp, recursive = FALSE, use.names = FALSE)
  labels <- rep(grp_names, vapply(trajs_by_grp, length, integer(1)))

  ## A block id carried on each member tree (from context_tree(block = ))
  ## is used by default. Filter it through the same survivor index the
  ## fit used (.ct_traj drops empty rows) so it stays aligned to `pooled`.
  if (is.null(block)) {
    mb <- lapply(seq_along(group), function(i) {
      b <- attr(group[[i]], "block")
      if (is.null(b)) NULL else b[attr(trajs_by_grp[[i]], "idx")]
    })
    if (!any(vapply(mb, is.null, logical(1))))
      block <- unlist(mb, use.names = FALSE)
  }

  ## Permutation scheme. Without 'block', labels are shuffled freely
  ## (sequences assumed exchangeable). With 'block' (one id per pooled
  ## sequence, in group-then-row order), the shuffle is stratified
  ## *within* each block, so nested/repeated-measures structure — e.g.
  ## several prompts from one student — is respected and any
  ## between-block confound is held fixed.
  if (is.null(block)) {
    perm_labels <- function() sample(labels)
  } else {
    if (length(block) != length(labels))
      stop("'block' must have one entry per sequence (length ",
           length(labels), "), aligned to the group-then-row order.",
           call. = FALSE)
    block <- as.character(block)
    blk_idx <- split(seq_along(labels), block)
    ## A block only contributes to the null if it holds >1 distinct
    ## label; if none do (e.g. a block id that is unique per sequence),
    ## no label can ever move and the p-values are a degenerate ~1.
    if (!any(vapply(blk_idx,
                    function(idx) length(unique(labels[idx])) > 1L,
                    logical(1))))
      warning("Stratified permutation is degenerate: no block contains ",
              "more than one group, so labels never move and p-values ",
              "collapse to ~1. Is 'block' unique per sequence?",
              call. = FALSE)
    perm_labels <- function() {
      out <- labels
      ## permute the label *values* within each block; safe for
      ## single-sequence blocks (sample() of a length-1 character vector
      ## returns it, unlike sample() of a single index).
      for (idx in blk_idx) out[idx] <- sample(labels[idx])
      out
    }
  }

  contexts <- unique(unlist(lapply(group, function(m) names(m$nodes)),
                            use.names = FALSE))
  if (!.ROOT %in% contexts) contexts <- c(.ROOT, contexts)
  depths   <- vapply(contexts, .cg_depth, integer(1))

  ## ---- observed ----
  obs_cc <- .cg_context_counts(
    split(pooled, factor(labels, levels = grp_names)),
    grp_names, contexts, depths, alphabet)
  obs <- list(
    jsd   = vapply(obs_cc, .cg_jsd_bits, numeric(1)),
    usage = .cg_usage_vec(obs_cc, depths))
  obs$omni_beh <- sum(vapply(obs_cc, sum, numeric(1)) * obs$jsd)
  obs$omni_use <- sum(obs$usage)

  ## ---- permutation null (label shuffle) ----
  ge_jsd  <- integer(length(contexts))
  ge_use  <- integer(length(contexts))
  ge_obeh <- 0L; ge_ouse <- 0L
  for (b in seq_len(iter)) {
    s <- .cg_stats(perm_labels(), pooled, grp_names, contexts, depths,
                   alphabet, .ROOT)
    ge_jsd  <- ge_jsd  + (s$jsd   >= obs$jsd   - 1e-12)
    ge_use  <- ge_use  + (s$usage >= obs$usage - 1e-12)
    ge_obeh <- ge_obeh + (s$omni_beh >= obs$omni_beh - 1e-12)
    ge_ouse <- ge_ouse + (s$omni_use >= obs$omni_use - 1e-12)
  }
  jsd_p   <- (1 + ge_jsd) / (iter + 1)
  usage_p <- (1 + ge_use) / (iter + 1)

  ## ---- per-group count / modal columns ----
  count_g <- vapply(grp_names, function(g)
    vapply(obs_cc, function(m) sum(m[g, ]), numeric(1)), numeric(length(contexts)))
  modal_g <- vapply(grp_names, function(g)
    vapply(obs_cc, function(m) {
      if (sum(m[g, ]) == 0) NA_character_ else alphabet[which.max(m[g, ])]
    }, character(1)), character(length(contexts)))
  if (length(contexts) == 1L) {                   # vapply drops to vector
    count_g <- matrix(count_g, nrow = 1L)
    modal_g <- matrix(modal_g, nrow = 1L)
  }
  colnames(count_g) <- grp_names
  colnames(modal_g) <- grp_names
  flips <- apply(modal_g, 1L, function(r) length(unique(r[!is.na(r)])) > 1L)

  pathway <- ifelse(contexts == .ROOT, .ROOT_LABEL, contexts)
  tab <- data.frame(pathway = pathway, depth = depths,
                    count_total = rowSums(count_g),
                    stringsAsFactors = FALSE)
  for (g in grp_names) tab[[paste0("count_", g)]] <- count_g[, g]
  for (g in grp_names) tab[[paste0("modal_", g)]] <- modal_g[, g]
  tab$flips      <- flips
  tab$jsd_bits   <- obs$jsd
  tab$jsd_p      <- jsd_p
  tab$jsd_padj   <- stats::p.adjust(jsd_p, method = "BH")
  ## The root has no prevalence test (it occupies 100% of positions by
  ## definition), so its usage statistic is structurally 0/p=1; report it
  ## as NA and keep it out of the usage FDR family so it does not pad the
  ## multiple-testing count.
  is_root        <- contexts == .ROOT
  tab$usage_g2   <- ifelse(is_root, NA_real_, obs$usage)
  tab$usage_p    <- ifelse(is_root, NA_real_, usage_p)
  tab$usage_padj <- NA_real_
  tab$usage_padj[!is_root] <- stats::p.adjust(usage_p[!is_root],
                                              method = "BH")

  tab <- tab[tab$count_total >= min_count, , drop = FALSE]
  tab <- tab[order(-tab$jsd_bits, -tab$usage_g2), , drop = FALSE]
  rownames(tab) <- NULL

  omnibus <- data.frame(
    axis      = c("behavioral", "usage"),
    statistic = c("count-weighted JSD (bits)", "sum G^2"),
    value     = c(obs$omni_beh, obs$omni_use),
    p_value   = c((1 + ge_obeh) / (iter + 1), (1 + ge_ouse) / (iter + 1)),
    stringsAsFactors = FALSE)

  ## ---- pairwise symmetric-KL distance matrix ----
  dm <- matrix(0, K, K, dimnames = list(grp_names, grp_names))
  if (K >= 2L) {
    pr <- utils::combn(K, 2L)
    for (j in seq_len(ncol(pr))) {
      d <- tree_distance(group[[pr[1, j]]], group[[pr[2, j]]])
      dm[pr[1, j], pr[2, j]] <- dm[pr[2, j], pr[1, j]] <- d
    }
  }

  structure(list(
    pathways = tab, omnibus = omnibus, distance_matrix = dm,
    groups = grp_names, iter = iter, seed = seed,
    n_contexts = nrow(tab), stratified = !is.null(block)),
    class = "transitiontrees_group_comparison")
}

#' Print a Group Comparison
#'
#' @description
#' Print method for a \code{transitiontrees_group_comparison}: a header with
#' the groups and permutation setup, the omnibus behavioral / usage
#' statistics, and the top pathways ranked by behavioral divergence.
#'
#' @param x A \code{transitiontrees_group_comparison} object.
#' @param n Integer. Number of top pathways to print. Default 10.
#' @param digits Integer. Numeric digits for the printed tables.
#'   Default 3.
#' @param ... Ignored.
#' @return \code{x} invisibly.
#' @export
print.transitiontrees_group_comparison <- function(x, n = 10L, digits = 3L, ...) {
  cat(sprintf("<transitiontrees_group_comparison>  %d groups, %d %spermutations\n",
              length(x$groups), x$iter,
              if (isTRUE(x$stratified)) "stratified " else ""))
  cat(sprintf("  groups   : %s\n", paste(x$groups, collapse = ", ")))
  cat(sprintf("  contexts : %d tested\n", x$n_contexts))
  cat("\nomnibus (permutation):\n")
  om <- x$omnibus
  om$value   <- round(om$value, digits)
  om$p_value <- round(om$p_value, digits)
  print.data.frame(om, row.names = FALSE)
  cat("\ntop pathways by behavioral divergence (JSD):\n")
  keep <- c("pathway", "depth", "count_total",
            grep("^modal_", names(x$pathways), value = TRUE),
            "flips", "jsd_bits", "jsd_padj", "usage_padj")
  show <- utils::head(x$pathways[, keep, drop = FALSE], n)
  show$jsd_bits   <- round(show$jsd_bits, digits)
  show$jsd_padj   <- round(show$jsd_padj, digits)
  show$usage_padj <- round(show$usage_padj, digits)
  print.data.frame(show, row.names = FALSE)
  if (nrow(x$pathways) > n)
    cat(sprintf("# ... %d more (use as.data.frame(x) for the full table)\n",
                nrow(x$pathways) - n))
  invisible(x)
}

#' Summarise a Group Comparison
#'
#' @description
#' Prints a compact verdict for a \code{transitiontrees_group_comparison}:
#' the omnibus behavioral and usage permutation p-values, and how many
#' pathways pass the FDR cutoff on each axis or flip their modal next
#' state between groups. Returns the per-pathway table invisibly.
#'
#' @param object A \code{transitiontrees_group_comparison} object.
#' @param alpha Numeric. FDR cutoff used when counting significant
#'   pathways. Default 0.05.
#' @param ... Ignored.
#' @return Invisibly, the per-pathway data.frame (\code{object$pathways});
#'   see \code{\link{compare_groups}} for the column vocabulary.
#' @export
summary.transitiontrees_group_comparison <- function(object, alpha = 0.05, ...) {
  p <- object$pathways
  cat(sprintf("Group comparison: %s\n", paste(object$groups, collapse = " vs ")))
  cat(sprintf("  behavioral omnibus p = %.3f | usage omnibus p = %.3f\n",
              object$omnibus$p_value[1], object$omnibus$p_value[2]))
  cat(sprintf("  pathways with behavioral difference (jsd_padj < %.2f): %d\n",
              alpha, sum(p$jsd_padj < alpha, na.rm = TRUE)))
  cat(sprintf("  pathways with usage difference      (usage_padj < %.2f): %d\n",
              alpha, sum(p$usage_padj < alpha, na.rm = TRUE)))
  cat(sprintf("  pathways whose modal next state flips between groups: %d\n",
              sum(p$flips)))
  invisible(object$pathways)
}

#' Coerce a Group Comparison to a Tidy Data Frame
#'
#' @description
#' Uniform tidy-extract: returns the per-pathway comparison table
#' (\code{object$pathways}) so \code{as.data.frame(cmp)} yields the full
#' divergence / usage breakdown as a base \code{data.frame}.
#'
#' @param x A \code{transitiontrees_group_comparison}.
#' @param row.names,optional Ignored.
#' @param ... Ignored.
#' @return A data.frame of per-pathway results; see
#'   \code{\link{compare_groups}} for the column vocabulary.
#' @export
as.data.frame.transitiontrees_group_comparison <- function(x, row.names = NULL,
                                                       optional = FALSE, ...) {
  x$pathways
}

#' Plot a Group Comparison
#'
#' @param x A \code{transitiontrees_group_comparison}.
#' @param style One of \code{"divergence"} (default; top pathways ranked
#'   by behavioral JSD, significant ones highlighted) or \code{"matrix"}
#'   (heatmap of the between-group symmetric-KL distance matrix).
#' @param top Integer. Pathways to show in \code{"divergence"}. Default 15.
#' @param alpha Numeric. FDR cutoff for highlighting. Default 0.05.
#' @param ... Ignored.
#' @return A ggplot object.
#' @export
plot.transitiontrees_group_comparison <- function(x, style = c("divergence",
                                                            "matrix"),
                                               top = 15L, alpha = 0.05, ...) {
  style <- match.arg(style)
  if (style == "matrix") {
    dm <- x$distance_matrix
    df <- expand.grid(a = rownames(dm), b = colnames(dm),
                      stringsAsFactors = FALSE)
    df$dist <- as.vector(dm)
    return(
      ggplot2::ggplot(df, ggplot2::aes(.data$a, .data$b, fill = .data$dist)) +
        ggplot2::geom_tile(colour = "white") +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", .data$dist)),
                           size = 3.2) +
        ggplot2::scale_fill_gradient(low = "#deebf7", high = "#08519c",
                                     name = "sym KL") +
        ggplot2::labs(x = NULL, y = NULL,
                      title = "Between-group distance (symmetric KL)") +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(plot.title =
                         ggplot2::element_text(face = "bold")))
  }

  d <- utils::head(x$pathways[x$pathways$jsd_bits > 0, , drop = FALSE], top)
  if (nrow(d) == 0L)
    stop("No pathways with positive behavioral divergence to plot.",
         call. = FALSE)
  d$pathway <- factor(d$pathway, levels = rev(d$pathway))
  d$sig <- ifelse(d$jsd_padj < alpha, sprintf("FDR < %.2f", alpha), "n.s.")
  ggplot2::ggplot(d, ggplot2::aes(x = .data$jsd_bits, y = .data$pathway,
                                  colour = .data$sig)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$jsd_bits,
                                       yend = .data$pathway),
                          colour = "grey75") +
    ggplot2::geom_point(ggplot2::aes(size = .data$count_total)) +
    ggplot2::scale_colour_manual(
      values = stats::setNames(c("#D55E00", "#bbbbbb"),
                               c(sprintf("FDR < %.2f", alpha), "n.s.")),
      name = NULL) +
    ggplot2::scale_size_continuous(name = "count", range = c(2, 7)) +
    ggplot2::labs(
      x = "behavioral divergence between groups (JSD, bits)", y = NULL,
      title = sprintf("How %s differ by context",
                      paste(x$groups, collapse = " vs ")),
      subtitle = "given the same history, how differently do the groups continue?") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40", size = 9),
      panel.grid.major.y = ggplot2::element_blank())
}
