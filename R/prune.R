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
#' @param tree A \code{transitrees}, or a \code{transitrees_group}, in
#'   which case each member is pruned and the group wrapper is
#'   preserved.
#' @param criterion One of \code{"G2"} (likelihood-ratio test against
#'   parent; default), \code{"KL"} (per-context Kullback-Leibler
#'   against parent), \code{"AIC"} (Akaike penalty), \code{"BIC"}
#'   (Bayesian penalty). Case-sensitive.
#' @param alpha Numeric in (0, 1). Significance level for \code{"G2"};
#'   ignored otherwise. Default 0.05.
#' @param threshold Numeric. Minimum information gain in nats for
#'   \code{"KL"}; ignored otherwise. Default 0.005.
#'
#' @return A pruned \code{transitrees} with \code{tree$pruned = TRUE} and
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
#' pruned <- prune_tree(tree, criterion = "G2", alpha = 0.05)
#' }
#'
#' @references
#' Ron, D., Singer, Y., Tishby, N. (1996). The power of amnesia.
#' \emph{Machine Learning}, 25, 117-149.
#'
#' @export
prune_tree <- function(tree,
                           criterion = c("G2", "KL", "AIC", "BIC"),
                           alpha = 0.05, threshold = 0.005) {
  ## A transitrees_group prunes each member, preserving the group wrapper.
  if (inherits(tree, "transitrees_group")) {
    out <- lapply(tree, prune_tree, criterion = criterion,
                  alpha = alpha, threshold = threshold)
    return(structure(out, class = class(tree),
                     group = attr(tree, "group")))
  }
  stopifnot(inherits(tree, "transitrees"))
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

#' Compare Pruning Criteria on One Tree
#'
#' @description
#' Prunes a fitted tree under several criteria — holding \code{alpha}
#' and \code{threshold} fixed — and returns a tidy one-row-per-criterion
#' summary of how aggressively each trims the tree. A convenience
#' wrapper over repeated \code{\link{prune_tree}} calls that
#' collapses the usual \code{vapply()} criterion loop into one call.
#'
#' @param tree A \code{transitrees} (typically unpruned).
#' @param criterion Character vector of criteria to compare. Defaults
#'   to all four: \code{"G2"}, \code{"KL"}, \code{"AIC"}, \code{"BIC"}.
#' @param alpha Significance level for \code{"G2"} (and the AIC/BIC
#'   penalties' chi-square cutoff). Default 0.05.
#' @param threshold Minimum information gain in nats for \code{"KL"}.
#'   Default 0.005.
#'
#' @return A \code{data.frame} with one row per criterion (in the order
#'   given by \code{criterion}) and columns \code{criterion},
#'   \code{n_nodes} (post-prune size) and \code{reduction_pct} (percent
#'   of the original nodes removed).
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' seqs <- replicate(80, sample(c("A", "B", "C"), 14, replace = TRUE),
#'                   simplify = FALSE)
#' tree <- context_tree(seqs, max_depth = 4L, min_count = 3L)
#' compare_pruning(tree)
#' compare_pruning(tree, criterion = c("G2", "BIC"), alpha = 0.01)
#' }
#'
#' @seealso \code{\link{prune_tree}} to apply one criterion,
#'   \code{\link{tune_tree}} for cross-validated selection.
#' @export
compare_pruning <- function(tree,
                            criterion = c("G2", "KL", "AIC", "BIC"),
                            alpha = 0.05, threshold = 0.005) {
  stopifnot(inherits(tree, "transitrees"))
  if (!is.character(criterion) || length(criterion) < 1L)
    stop("'criterion' must be a non-empty character vector of ",
         "criterion names.", call. = FALSE)
  n0 <- length(tree$nodes)
  n  <- vapply(criterion, function(cr)
    length(prune_tree(tree, criterion = cr, alpha = alpha,
                          threshold = threshold)$nodes), integer(1))
  data.frame(
    criterion     = criterion,
    n_nodes       = n,
    reduction_pct = round(100 * (1 - n / n0), 1),
    stringsAsFactors = FALSE,
    row.names     = NULL)
}
