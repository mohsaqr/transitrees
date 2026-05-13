# ---- Likelihood, perplexity, and stats S3 generics ----
#
# Implements:
#   logLik.pathtree   - in-sample if newdata = NULL, held-out otherwise
#   nobs.pathtree
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
         matched_context = ifelse(ctx_k == .ROOT, "(root)", ctx_k),
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
#' @param object A \code{pathtree}.
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
#'                      max_depth = 2, nmin = 2)
#' logLik(tree)
#' AIC(tree); BIC(tree)
#' }
#' @export
logLik.pathtree <- function(object, newdata = NULL, ...) {
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
#' @param object A \code{pathtree}.
#' @param ... Ignored.
#' @return Integer. Number of state observations used to fit the tree.
#' @export
nobs.pathtree <- function(object, ...) {
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
#' @param tree A \code{pathtree}.
#' @param newdata Sequence data; \code{NULL} (default) returns
#'   in-sample perplexity.
#' @return Numeric scalar.
#' @examples
#' \donttest{
#' tree <- context_tree(matrix(sample(c("A","B","C"), 200, TRUE), 20),
#'                      max_depth = 2, nmin = 2)
#' perplexity(tree)
#' }
#' @export
perplexity <- function(tree, newdata = NULL) {
  stopifnot(inherits(tree, "pathtree"))
  ll <- logLik.pathtree(tree, newdata = newdata)
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
#' @param tree A \code{pathtree}.
#' @param newdata Sequence data.
#' @return A data.frame with columns \code{sequence_id},
#'   \code{n_scored}, \code{log_lik}, \code{perplexity}.
#' @export
score_sequences <- function(tree, newdata) {
  stopifnot(inherits(tree, "pathtree"))
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
#' @param tree A \code{pathtree}.
#' @param newdata Sequence data.
#' @return A data.frame with columns \code{sequence_id},
#'   \code{position}, \code{matched_context}, \code{observed},
#'   \code{predicted_prob}, \code{log_lik}.
#' @export
score_positions <- function(tree, newdata) {
  stopifnot(inherits(tree, "pathtree"))
  .pt_score_walk(tree, newdata)
}
