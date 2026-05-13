# ---- Pruning ----

#' Prune a Context Tree
#'
#' @description
#' Removes nodes that do not earn their depth under the chosen
#' criterion. Pruning is applied bottom-up: a node is dropped when
#' extending its parent's prediction at this context produces less
#' information / likelihood than the depth penalty allows.
#'
#' Renamed from \code{prune()} to avoid collision with
#' \code{tna::prune}, \code{rpart::prune}, and \code{tree::prune}.
#'
#' @param tree A \code{pathtree}.
#' @param criterion One of \code{"G2"} (likelihood-ratio test against
#'   parent; default), \code{"KL"} (per-context Kullback-Leibler
#'   against parent), \code{"AIC"} (Akaike penalty), \code{"BIC"}
#'   (Bayesian penalty). Case-sensitive.
#' @param alpha Numeric in (0, 1). Significance level for \code{"G2"};
#'   ignored otherwise. Default 0.05.
#' @param threshold Numeric. Minimum information gain in nats for
#'   \code{"KL"}; ignored otherwise. Default 0.005.
#'
#' @return A pruned \code{pathtree} with \code{tree$pruned = TRUE} and
#'   \code{tree$pruning} carrying the criterion + threshold settings.
#'
#' @details
#' For each leaf, compute the criterion against its parent. If the
#' criterion does not exceed its threshold, drop the leaf and revisit
#' the parent. Repeat until stable. The root is never dropped.
#' Surviving nodes keep their original smoothed \code{prob} vector
#' (whatever smoothing scheme was applied at fit time).
#'
#' @examples
#' \donttest{
#' seqs   <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
#'                     simplify = FALSE)
#' tree   <- context_tree(seqs, max_depth = 4)
#' pruned <- prune_pathtree(tree, criterion = "G2", alpha = 0.05)
#' }
#'
#' @references
#' Ron, D., Singer, Y., Tishby, N. (1996). The power of amnesia.
#' \emph{Machine Learning}, 25, 117-149.
#'
#' @export
prune_pathtree <- function(tree,
                           criterion = c("G2", "KL", "AIC", "BIC"),
                           alpha = 0.05, threshold = 0.005) {
  stopifnot(inherits(tree, "pathtree"))
  criterion <- match.arg(criterion)

  ## Critical value for G2: chi-square at (k-1) df
  k_states <- length(tree$alphabet)
  g2_cv    <- stats::qchisq(1 - alpha, df = k_states - 1L)

  nodes <- tree$nodes
  edges <- tree$edges
  alpha_size <- k_states

  ## Repeat until no further removals
  repeat {
    leaves <- setdiff(names(nodes), edges$parent)
    leaves <- leaves[leaves != .ROOT]  # never the root
    drop   <- character(0)

    for (leaf in leaves) {
      parts  <- strsplit(leaf, " -> ", fixed = TRUE)[[1L]]
      parent <- if (length(parts) == 1L) .ROOT else
        paste(parts[-1L], collapse = " -> ")
      if (!parent %in% names(nodes)) next

      child  <- nodes[[leaf]]
      par    <- nodes[[parent]]
      keep <- switch(
        criterion,
        G2  = .ct_g2(child$counts, par$prob) > g2_cv,
        KL  = .ct_kl(child$prob, par$prob) > threshold,
        AIC = {
          ll_child  <- sum(child$counts * log(pmax(child$prob,
                                                   .Machine$double.eps)))
          ll_parent <- sum(child$counts * log(pmax(par$prob,
                                                   .Machine$double.eps)))
          aic_child  <- 2 * (alpha_size - 1L) - 2 * ll_child
          aic_parent <- -2 * ll_parent
          aic_child < aic_parent
        },
        BIC = {
          ll_child  <- sum(child$counts * log(pmax(child$prob,
                                                   .Machine$double.eps)))
          ll_parent <- sum(child$counts * log(pmax(par$prob,
                                                   .Machine$double.eps)))
          bic_child  <- log(child$n) * (alpha_size - 1L) - 2 * ll_child
          bic_parent <- -2 * ll_parent
          bic_child < bic_parent
        }
      )
      if (!keep) drop <- c(drop, leaf)
    }
    if (length(drop) == 0L) break
    nodes <- nodes[!names(nodes) %in% drop]
    edges <- edges[!edges$child %in% drop, , drop = FALSE]
  }

  tree$nodes  <- nodes
  tree$edges  <- edges
  tree$pruned <- TRUE
  tree$pruning <- list(criterion = criterion, alpha = alpha,
                       threshold = threshold)
  tree
}
