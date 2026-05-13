# ---- pathtree_dependence(): per-context KL diagnostic ----

#' @noRd
.pt_empty_dependence_df <- function() {
  data.frame(
    pathway      = character(0),
    depth        = integer(0),
    count        = numeric(0),
    KL           = numeric(0),
    H_node       = numeric(0),
    H_parent     = numeric(0),
    H_drop       = numeric(0),
    modal_next   = character(0),
    modal_parent = character(0),
    flips        = logical(0),
    stringsAsFactors = FALSE
  )
}

#' Per-Context Path Dependence of a Context Tree
#'
#' @description
#' For each non-root node in the tree, reports the Kullback-Leibler
#' divergence of its conditional next-state distribution against its
#' parent's. Large values flag contexts where extending memory by one
#' more step changes the prediction; the \code{flips} column flags
#' contexts where the modal next state changes between the node and its
#' parent.
#'
#' Renamed from \code{path_dependence()} to avoid collision with
#' \code{Nestimate::path_dependence}.
#'
#' @param tree A \code{pathtree}.
#' @param base Numeric. Logarithm base for the KL divergence. Default
#'   2 (bits). Use \code{exp(1)} for nats or \code{10} for hartleys.
#'
#' @return A data.frame with one row per non-root pathway, sorted by
#'   KL descending. Columns: \code{pathway}, \code{depth}, \code{count},
#'   \code{KL}, \code{H_node} (entropy of the pathway's next-state
#'   distribution), \code{H_parent} (entropy of the parent's
#'   distribution), \code{H_drop} (\code{H_parent - H_node}),
#'   \code{modal_next} (this node's modal next state),
#'   \code{modal_parent} (the parent context's modal next state),
#'   \code{flips} (\code{modal_next != modal_parent}). The empty case
#'   returns a 0-row data.frame with the same schema.
#'
#' @details
#' This is the \emph{diagnostic} that the tree's pruning rule (under
#' \code{criterion = "KL"}) is comparing against its threshold. It
#' answers the substantive question: for which contexts does this
#' tree disagree with a memoryless / shorter-memory model, and where
#' does that disagreement actually flip the prediction?
#'
#' The mean of \code{n * KL} across rows recovers, up to constants,
#' the chain-level mutual-information gain from the variable-depth
#' model over the order-1 model.
#'
#' @examples
#' \donttest{
#' seqs   <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
#'                     simplify = FALSE)
#' tree   <- context_tree(seqs, max_depth = 3)
#' pruned <- prune_pathtree(tree, criterion = "G2")
#' pathtree_dependence(pruned)
#' }
#'
#' @references
#' Cover, T.M. & Thomas, J.A. (2006). \emph{Elements of Information
#' Theory}, 2nd ed. Wiley.
#'
#' @export
pathtree_dependence <- function(tree, base = 2) {
  stopifnot(inherits(tree, "pathtree"))
  stopifnot(is.numeric(base), length(base) == 1L, base > 0, base != 1)

  H <- function(p) {
    msk <- p > 0
    -sum(p[msk] * log(p[msk], base = base))
  }
  KL <- function(p, q) {
    msk <- p > 0
    if (any(p[msk] > 0 & q[msk] == 0)) return(Inf)
    sum(p[msk] * log(p[msk] / q[msk], base = base))
  }

  rows <- lapply(setdiff(names(tree$nodes), .ROOT), function(ctx) {
    parts  <- strsplit(ctx, " -> ", fixed = TRUE)[[1L]]
    parent <- if (length(parts) == 1L) .ROOT else
      paste(parts[-1L], collapse = " -> ")
    if (!parent %in% names(tree$nodes)) return(NULL)

    node <- tree$nodes[[ctx]]
    par  <- tree$nodes[[parent]]
    p    <- node$prob
    q    <- par$prob

    modal_next   <- tree$alphabet[which.max(p)]
    modal_parent <- tree$alphabet[which.max(q)]

    data.frame(
      pathway      = ctx,
      depth        = node$depth,
      count        = node$n,
      KL           = KL(p, q),
      H_node       = H(p),
      H_parent     = H(q),
      H_drop       = H(q) - H(p),
      modal_next   = modal_next,
      modal_parent = modal_parent,
      flips        = modal_next != modal_parent,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(out) || nrow(out) == 0L) return(.pt_empty_dependence_df())
  out <- out[order(-out$KL), , drop = FALSE]
  rownames(out) <- NULL
  out
}
