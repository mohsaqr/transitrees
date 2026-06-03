# ---- tree_dependence(): per-context KL diagnostic ----

#' @noRd
.pt_empty_dependence_df <- function() {
  data.frame(
    pathway            = character(0),
    depth              = integer(0),
    count              = numeric(0),
    divergence         = numeric(0),
    entropy            = numeric(0),
    entropy_before     = numeric(0),
    entropy_drop       = numeric(0),
    likely_next        = character(0),
    likely_before      = character(0),
    changes_prediction = logical(0),
    stringsAsFactors   = FALSE
  )
}

#' Per-Context Path Dependence of a Context Tree
#'
#' @description
#' For each non-root node in the tree, reports the Kullback-Leibler
#' divergence of its conditional next-state distribution against its
#' parent's. Large values flag contexts where extending memory by one
#' more step changes the prediction; the \code{changes_prediction}
#' column flags contexts where the most likely next state changes
#' between the node and its parent.
#'
#' Renamed from \code{path_dependence()} to avoid collision with
#' \code{Nestimate::path_dependence}.
#'
#' @param tree A \code{transitrees}.
#' @param base Numeric. Logarithm base for the KL divergence. Default
#'   2 (bits). Use \code{exp(1)} for nats or \code{10} for hartleys.
#' @param sort_by Character. Column to sort by, descending. One of
#'   \code{"divergence"} (default), \code{"entropy_drop"},
#'   \code{"entropy"}, \code{"count"}, \code{"depth"}.
#' @param top Integer or \code{NULL}. If given, keep only the top
#'   \code{top} rows after sorting. Default \code{NULL} (all rows).
#'
#' @return A data.frame with one row per non-root pathway, sorted by
#'   \code{divergence} descending. Columns: \code{pathway},
#'   \code{depth}, \code{count}, \code{divergence} (Kullback-Leibler
#'   divergence from the parent's prediction), \code{entropy} (Shannon
#'   entropy of this pathway's next-state distribution),
#'   \code{entropy_before} (entropy of the parent's distribution),
#'   \code{entropy_drop} (\code{entropy_before - entropy}, the
#'   uncertainty this step of history removes), \code{likely_next}
#'   (this node's most likely next state), \code{likely_before} (the
#'   parent context's most likely next state), \code{changes_prediction}
#'   (\code{likely_next != likely_before}). The empty case returns a
#'   0-row data.frame with the same schema.
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
#' pruned <- prune_tree(tree, criterion = "G2")
#' tree_dependence(pruned)
#' }
#'
#' @references
#' Cover, T.M. & Thomas, J.A. (2006). \emph{Elements of Information
#' Theory}, 2nd ed. Wiley.
#'
#' @export
tree_dependence <- function(tree, base = 2,
                                sort_by = c("divergence", "entropy_drop",
                                            "entropy", "count", "depth"),
                                top = NULL) {
  stopifnot(inherits(tree, "transitrees"))
  stopifnot(is.numeric(base), length(base) == 1L, base > 0, base != 1)
  sort_by <- match.arg(sort_by)

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
      pathway            = ctx,
      depth              = node$depth,
      count              = node$n,
      divergence         = KL(p, q),
      entropy            = H(p),
      entropy_before     = H(q),
      entropy_drop       = H(q) - H(p),
      likely_next        = modal_next,
      likely_before      = modal_parent,
      changes_prediction = modal_next != modal_parent,
      stringsAsFactors   = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(out) || nrow(out) == 0L) return(.pt_empty_dependence_df())
  out <- out[order(-out[[sort_by]]), , drop = FALSE]
  rownames(out) <- NULL
  if (!is.null(top)) out <- utils::head(out, n = as.integer(top))
  out
}
