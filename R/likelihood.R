# ---- Likelihood, perplexity, and stats S3 generics ----
#
# Implements:
#   logLik.transitrees   - in-sample if newdata = NULL, held-out otherwise
#   nobs.transitrees
#   perplexity()      - exp(- ll / n)
#   score_sequences() - per-sequence log-lik
#   score_positions() - per-position log-lik

#' @noRd
.pt_unique_counts <- function(tree) {
  ## Each token is attributed to exactly one node (its deepest matching
  ## context); summing unique counts across all nodes equals n_obs.
  children_of <- .pt_children_of(tree)
  lapply(setNames(names(tree$nodes), names(tree$nodes)), function(ctx) {
    info     <- tree$nodes[[ctx]]
    children <- children_of[[ctx]]
    if (is.null(children) || length(children) == 0L) return(info$counts)
    child_sum <- Reduce("+", lapply(children,
                                    function(c) tree$nodes[[c]]$counts))
    info$counts - child_sum
  })
}

#' @noRd
.pt_loglik_in_sample <- function(tree) {
  uc <- .pt_unique_counts(tree)
  contribs <- vapply(names(tree$nodes), function(ctx) {
    counts <- uc[[ctx]]
    p      <- tree$nodes[[ctx]]$prob
    msk    <- counts > 0 & p > 0
    if (!any(msk)) return(0)
    sum(counts[msk] * log(p[msk]))
  }, numeric(1))
  list(ll = sum(contribs), n = sum(vapply(uc, sum, numeric(1))))
}

#' @noRd
.pt_score_walk <- function(tree, newdata) {
  trajs <- .ct_traj(.ct_coerce(newdata))
  if (length(trajs) == 0L)
    stop("No usable held-out sequences after coercion.", call. = FALSE)
  alpha <- tree$alphabet

  per_seq <- lapply(seq_along(trajs), function(i) {
    traj <- trajs[[i]]
    L    <- length(traj)
    if (L < 1L) return(NULL)
    ctx_per_pos <- vapply(seq_len(L), function(t) {
      hist <- if (t == 1L) character(0) else traj[seq_len(t - 1L)]
      .ct_match_context(tree, hist)
    }, character(1))
    obs    <- traj
    idx    <- match(obs, alpha)
    keep   <- !is.na(idx)
    if (!any(keep)) return(NULL)
    ctx_k  <- ctx_per_pos[keep]
    pos_k  <- seq_len(L)[keep]
    obs_k  <- obs[keep]
    idx_k  <- idx[keep]
    p_k    <- vapply(seq_along(ctx_k),
                     function(j) tree$nodes[[ctx_k[j]]]$prob[idx_k[j]],
                     numeric(1))
    list(sequence_id     = rep.int(i, length(pos_k)),
         position        = pos_k,
         matched_context = ifelse(ctx_k == .ROOT, .ROOT_LABEL, ctx_k),
         observed        = obs_k,
         predicted_prob  = p_k,
         log_lik         = ifelse(p_k > 0, log(p_k), -Inf))
  })
  per_seq <- per_seq[!vapply(per_seq, is.null, logical(1))]

  if (length(per_seq) == 0L)
    return(data.frame(sequence_id = integer(0), position = integer(0),
                      matched_context = character(0),
                      observed = character(0),
                      predicted_prob = numeric(0),
                      log_lik = numeric(0),
                      stringsAsFactors = FALSE))

  bind <- function(field, fn) fn(lapply(per_seq, `[[`, field))
  data.frame(
    sequence_id     = bind("sequence_id",     function(x) unlist(x, use.names = FALSE)),
    position        = bind("position",        function(x) unlist(x, use.names = FALSE)),
    matched_context = bind("matched_context", function(x) unlist(x, use.names = FALSE)),
    observed        = bind("observed",        function(x) unlist(x, use.names = FALSE)),
    predicted_prob  = bind("predicted_prob",  function(x) unlist(x, use.names = FALSE)),
    log_lik         = bind("log_lik",         function(x) unlist(x, use.names = FALSE)),
    stringsAsFactors = FALSE
  )
}

#' Log-Likelihood of a Pathtree
#'
#' @description
#' Returns a \code{logLik} object compatible with \code{stats::AIC()},
#' \code{stats::BIC()}, and the rest of the model-comparison toolchain.
#' If \code{newdata} is \code{NULL}, returns the in-sample log-likelihood
#' computed from the fitted node counts. Otherwise returns the held-out
#' log-likelihood scoring \code{newdata} under the fitted tree.
#'
#' @param object A \code{transitrees}.
#' @param newdata Optional. Sequence data in any format accepted by
#'   \code{context_tree()}. \code{NULL} (default) returns in-sample
#'   log-likelihood.
#' @param ... Ignored.
#'
#' @return A \code{logLik} object with attributes \code{nobs} and
#'   \code{df} (number of free parameters in the fitted tree).
#'
#' @examples
#' \donttest{
#' tree <- context_tree(matrix(sample(c("A","B","C"), 200, TRUE), 20),
#'                      max_depth = 2, min_count = 2)
#' logLik(tree)
#' AIC(tree); BIC(tree)
#' }
#' @export
logLik.transitrees <- function(object, newdata = NULL, ...) {
  alpha_size <- length(object$alphabet)
  df <- length(object$nodes) * (alpha_size - 1L)
  if (is.null(newdata)) {
    res <- .pt_loglik_in_sample(object)
    val <- res$ll; n <- res$n
  } else {
    pos <- .pt_score_walk(object, newdata)
    val <- sum(pos$log_lik)
    n   <- nrow(pos)
  }
  structure(val, nobs = as.integer(n), df = as.integer(df),
            class = "logLik")
}

#' Number of Observations Used to Fit a Pathtree
#'
#' @param object A \code{transitrees}.
#' @param ... Ignored.
#' @return Integer. Number of state observations used to fit the tree.
#' @export
nobs.transitrees <- function(object, ...) {
  as.integer(object$n_obs)
}

#' Perplexity of a Pathtree
#'
#' @description
#' \code{exp(-mean log-likelihood per observation)}, the standard
#' language-modelling evaluation metric. Lower is better. A perplexity
#' of \eqn{k} on an alphabet of size \eqn{|S|} means the model is as
#' predictive as a uniform distribution over \eqn{k} symbols.
#' \eqn{k = |S|} is the uniform baseline; \eqn{k = 1} is perfect
#' deterministic prediction.
#'
#' @param tree A \code{transitrees}.
#' @param newdata Sequence data; \code{NULL} (default) returns
#'   in-sample perplexity.
#' @return Numeric scalar.
#' @examples
#' \donttest{
#' tree <- context_tree(matrix(sample(c("A","B","C"), 200, TRUE), 20),
#'                      max_depth = 2, min_count = 2)
#' perplexity(tree)
#' }
#' @export
perplexity <- function(tree, newdata = NULL) {
  stopifnot(inherits(tree, "transitrees"))
  ll <- logLik.transitrees(tree, newdata = newdata)
  n  <- attr(ll, "nobs")
  if (n == 0L) return(NA_real_)
  exp(-as.numeric(ll) / n)
}

#' Per-Sequence Scoring
#'
#' @description
#' Returns one row per held-out sequence with its log-likelihood,
#' number of scored positions, and per-position perplexity.
#'
#' @param tree A \code{transitrees}.
#' @param newdata Sequence data.
#' @return A data.frame with columns \code{sequence_id},
#'   \code{n_scored}, \code{log_lik}, \code{perplexity}.
#' @export
score_sequences <- function(tree, newdata) {
  stopifnot(inherits(tree, "transitrees"))
  pos <- .pt_score_walk(tree, newdata)
  if (nrow(pos) == 0L)
    return(data.frame(sequence_id = integer(0), n_scored = integer(0),
                      log_lik = numeric(0), perplexity = numeric(0),
                      stringsAsFactors = FALSE))
  seq_ids <- sort(unique(pos$sequence_id))
  ll  <- vapply(seq_ids,
                function(i) sum(pos$log_lik[pos$sequence_id == i]),
                numeric(1))
  n   <- vapply(seq_ids,
                function(i) sum(pos$sequence_id == i), integer(1))
  pp  <- exp(-ll / n)
  out <- data.frame(sequence_id = seq_ids, n_scored = n,
                    log_lik = ll, perplexity = pp,
                    stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}

#' Per-Position Scoring
#'
#' @description
#' Returns one row per held-out (sequence, position) with the matched
#' context, predicted probability of the observed next state, and
#' log-likelihood contribution. Useful for diagnostic plots showing
#' where the model is confident vs. surprised.
#'
#' @param tree A \code{transitrees}.
#' @param newdata Sequence data.
#' @param worst Integer or \code{NULL}. If given, return only the
#'   \code{worst} positions — those with the lowest
#'   \code{predicted_prob} (the moves the model was most surprised by).
#'   Default \code{NULL} (all positions, in sequence order).
#' @return A data.frame with columns \code{sequence_id},
#'   \code{position}, \code{matched_context}, \code{observed},
#'   \code{predicted_prob}, \code{log_lik}.
#' @export
score_positions <- function(tree, newdata, worst = NULL) {
  stopifnot(inherits(tree, "transitrees"))
  out <- .pt_score_walk(tree, newdata)
  if (!is.null(worst) && nrow(out) > 0L) {
    out <- out[order(out$predicted_prob), , drop = FALSE]
    out <- utils::head(out, n = as.integer(worst))
    rownames(out) <- NULL
  }
  out
}

#' Model-Fit Scalars in One Call
#'
#' @description
#' Bundles the standard goodness-of-fit scalars for a fitted tree into
#' one tidy row: \code{logLik}, the parameter count \code{df}, the
#' observation count \code{nobs}, \code{AIC}, \code{BIC}, and
#' \code{perplexity}. A one-call replacement for
#' \code{logLik(); nobs(); AIC(); BIC(); perplexity()}.
#'
#' @details
#' With \code{newdata}, every scalar is computed \strong{out-of-sample}
#' (\code{AIC}/\code{BIC} use the held-out deviance with the model's
#' training \code{df}). A \code{transitrees_group} returns one row per
#' group, tagged with a leading \code{group} column.
#'
#' @param tree A \code{transitrees} or \code{transitrees_group}.
#' @param newdata Optional sequence data. If supplied, the scalars are
#'   evaluated on it (held-out); if \code{NULL} (default), in-sample.
#'
#' @return A one-row \code{data.frame} (one row per group for a
#'   \code{transitrees_group}) with columns \code{logLik}, \code{df},
#'   \code{nobs}, \code{AIC}, \code{BIC}, \code{perplexity}.
#'
#' @examples
#' \donttest{
#' seqs <- replicate(60, sample(c("A","B","C"), 12, replace = TRUE),
#'                   simplify = FALSE)
#' tree <- context_tree(seqs, max_depth = 2L, min_count = 3L)
#' model_fit(tree)
#' }
#'
#' @seealso \code{\link{perplexity}}, \code{\link{logLik.transitrees}}.
#' @export
model_fit <- function(tree, newdata = NULL) {
  if (inherits(tree, "transitrees_group")) {
    parts <- lapply(names(tree), function(nm)
      cbind(group = nm, model_fit(tree[[nm]], newdata = newdata),
            stringsAsFactors = FALSE))
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    return(out)
  }
  stopifnot(inherits(tree, "transitrees"))
  ll  <- logLik(tree, newdata = newdata)
  df  <- attr(ll, "df")
  n   <- attr(ll, "nobs")
  llv <- as.numeric(ll)
  data.frame(
    logLik     = llv,
    df         = df,
    nobs       = n,
    AIC        = -2 * llv + 2 * df,
    BIC        = -2 * llv + log(n) * df,
    perplexity = perplexity(tree, newdata = newdata),
    row.names  = NULL
  )
}

#' Number of Contexts (Nodes) in a Tree
#'
#' @description
#' The count of contexts the tree represents — an intuitive accessor for
#' \code{length(tree$nodes)} (the number printed in the tree banner).
#'
#' @param tree A \code{transitrees} or \code{transitrees_group}.
#' @return An integer. For a \code{transitrees_group}, a named integer
#'   vector with one count per group.
#'
#' @examples
#' \donttest{
#' seqs <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
#'                   simplify = FALSE)
#' n_nodes(context_tree(seqs, max_depth = 3L))
#' }
#' @export
n_nodes <- function(tree) {
  if (inherits(tree, "transitrees_group"))
    return(vapply(tree, function(t) length(t$nodes), integer(1)))
  stopifnot(inherits(tree, "transitrees"))
  length(tree$nodes)
}
