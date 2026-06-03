# ---- Pathway-centric API ----

#' @noRd
.pt_empty_pathways_df <- function() {
  data.frame(
    pathway            = character(0),
    depth              = integer(0),
    count              = numeric(0),
    likely_next        = character(0),
    next_probability   = numeric(0),
    divergence         = numeric(0),
    changes_prediction = logical(0),
    stringsAsFactors   = FALSE
  )
}

#' Pathways from a Fitted Tree
#'
#' @description
#' Returns a tidy data.frame with one row per pathway (= context) in the
#' tree. The pathway is the sequence of states the tree conditions on;
#' each row reports the count, depth, modal next state, and how
#' surprising the next-state distribution is relative to a shorter
#' history. This is the substantive consumer-facing API of a fitted
#' tree -- a ranked list of trajectories that the data actually
#' supports, with a consistent interpretive frame.
#'
#' Renamed from \code{pathways()} to avoid collision with
#' \code{Nestimate::pathways}.
#'
#' @param tree A \code{transitrees} object.
#' @param min_count Integer. Drop pathways with fewer than this many
#'   occurrences. Default 1.
#' @param sort_by Character. One of \code{"count"} (default),
#'   \code{"divergence"}, or \code{"depth"}. Sorts the returned
#'   data.frame.
#' @param decreasing Logical. Default \code{TRUE}.
#' @param ... Ignored.
#'
#' @return A data.frame with columns \code{pathway} (arrow notation,
#'   e.g. \code{"A -> B -> C"}; the root is reported as \code{"(start)"}),
#'   \code{depth} (history length), \code{count}, \code{likely_next}
#'   (the most likely next state), \code{next_probability} (its
#'   probability), \code{divergence} (Kullback-Leibler divergence from
#'   the parent context's prediction, in bits; \code{NA} for the root),
#'   and \code{changes_prediction} (logical, did the most likely next
#'   state change vs the parent context?). The empty case returns a
#'   0-row data.frame with the same schema.
#'
#' @details
#' Each row is a pathway -- a (possibly empty) sequence of states ending
#' at the point where a prediction is made. The \code{divergence} column
#' quantifies how much more information the pathway carries than its
#' (k-1)-suffix in bits. \code{changes_prediction = TRUE} marks pathways
#' where the longer history changes which next state is most likely.
#'
#' @examples
#' \donttest{
#' seqs <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
#'                   simplify = FALSE)
#' tree <- context_tree(seqs, max_depth = 3)
#' tree_pathways(tree)
#' }
#'
#' @export
tree_pathways <- function(tree, min_count = 1L,
                              sort_by = c("count", "divergence", "depth"),
                              decreasing = TRUE, ...) {
  stopifnot(inherits(tree, "transitrees"))
  sort_by   <- match.arg(sort_by)
  min_count <- as.integer(min_count)

  rows <- lapply(names(tree$nodes), function(ctx) {
    info <- tree$nodes[[ctx]]
    if (info$n < min_count) return(NULL)

    pathway <- if (identical(ctx, .ROOT)) .ROOT_LABEL else ctx
    parts <- if (identical(ctx, .ROOT)) character(0) else
      strsplit(ctx, " -> ", fixed = TRUE)[[1L]]

    modal_idx  <- which.max(info$prob)
    modal_next <- tree$alphabet[modal_idx]
    prob_next  <- info$prob[modal_idx]

    ## KL vs parent (NA for root)
    if (identical(ctx, .ROOT)) {
      KL    <- NA_real_
      flips <- NA
    } else {
      parent <- if (length(parts) == 1L) .ROOT else
        paste(parts[-1L], collapse = " -> ")
      par_info <- tree$nodes[[parent]]
      if (is.null(par_info)) {
        KL    <- NA_real_
        flips <- NA
      } else {
        p <- info$prob; q <- par_info$prob
        msk <- p > 0
        KL <- if (any(p[msk] > 0 & q[msk] == 0)) Inf else
          sum(p[msk] * log(p[msk] / q[msk], base = 2))
        modal_parent <- tree$alphabet[which.max(par_info$prob)]
        flips <- modal_next != modal_parent
      }
    }

    data.frame(
      pathway            = pathway,
      depth              = info$depth,
      count              = info$n,
      likely_next        = modal_next,
      next_probability   = prob_next,
      divergence         = KL,
      changes_prediction = flips,
      stringsAsFactors   = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(out) || nrow(out) == 0L) return(.pt_empty_pathways_df())

  key <- switch(sort_by,
                count      = out$count,
                divergence = ifelse(is.na(out$divergence), -Inf,
                                    out$divergence),
                depth      = out$depth)
  ord <- order(key, decreasing = decreasing)
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a


# ---- High-level convenience functions -----------------------------------

#' Most Common Pathways in a Fitted Tree
#'
#' @description
#' Returns the top \code{n} pathways by occurrence count -- the
#' trajectories the data actually contains many copies of.
#'
#' @param tree A \code{transitrees}.
#' @param top Integer. Number of pathways to return. Default 10.
#' @param depth Integer or NULL. Restrict to pathways of this exact
#'   depth. \code{NULL} (default) keeps all depths.
#' @param min_count Integer. Minimum count cut-off. Default 1.
#'
#' @return A data.frame, same columns as \code{\link{tree_pathways}},
#'   sorted by count descending.
#'
#' @examples
#' \donttest{
#' seqs <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
#'                   simplify = FALSE)
#' tree <- context_tree(seqs, max_depth = 3)
#' common_pathways(tree, top = 8)
#' common_pathways(tree, top = 8, depth = 3L)   # restrict to depth-3
#' }
#'
#' @export
common_pathways <- function(tree, top = 10L, depth = NULL, min_count = 1L) {
  stopifnot(inherits(tree, "transitrees"))
  out <- tree_pathways(tree, min_count = min_count, sort_by = "count",
                           decreasing = TRUE)
  if (!is.null(depth))
    out <- out[out$depth == as.integer(depth), , drop = FALSE]
  utils::head(out, n = as.integer(top))
}

#' Most Predictively Divergent Pathways
#'
#' @description
#' Returns the top \code{n} pathways by Kullback-Leibler divergence
#' from their (k-1)-suffix. These are the pathways whose extended
#' history adds the most predictive information over the shorter one.
#' Pathways whose most likely next state actually flips between orders
#' are marked in the \code{changes_prediction} column.
#'
#' @param tree A \code{transitrees}.
#' @param top Integer. Number of pathways to return. Default 10.
#' @param min_count Integer. Drop pathways with fewer than this many
#'   occurrences. Default 1.
#' @param flips_only Logical. If \code{TRUE}, return only pathways
#'   that flip the most likely next state between orders. Default
#'   \code{FALSE}.
#'
#' @return A data.frame, same columns as \code{\link{tree_pathways}},
#'   sorted by \code{divergence} descending.
#'
#' @examples
#' \donttest{
#' seqs <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
#'                   simplify = FALSE)
#' tree <- context_tree(seqs, max_depth = 3)
#' divergent_pathways(tree, top = 6)
#' divergent_pathways(tree, flips_only = TRUE)
#' }
#'
#' @export
divergent_pathways <- function(tree, top = 10L, min_count = 1L,
                               flips_only = FALSE) {
  stopifnot(inherits(tree, "transitrees"))
  out <- tree_pathways(tree, min_count = min_count,
                           sort_by = "divergence", decreasing = TRUE)
  out <- out[!is.na(out$divergence), , drop = FALSE]
  if (isTRUE(flips_only))
    out <- out[!is.na(out$changes_prediction) & out$changes_prediction,
               , drop = FALSE]
  utils::head(out, n = as.integer(top))
}

#' Sharpest Pathways
#'
#' @description
#' Returns the top \code{n} pathways by predictive sharpness -- the
#' probability mass on their modal next state. High values indicate
#' strongly deterministic continuations; low values indicate
#' ambiguous next-state distributions.
#'
#' @param tree A \code{transitrees}.
#' @param top Integer. Number of pathways to return. Default 10.
#' @param min_count Integer. Drop pathways with fewer than this many
#'   occurrences. Default 1.
#'
#' @return A data.frame, same columns as \code{\link{tree_pathways}},
#'   sorted by \code{next_probability} descending.
#'
#' @examples
#' \donttest{
#' seqs <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
#'                   simplify = FALSE)
#' tree <- context_tree(seqs, max_depth = 3)
#' sharp_pathways(tree, top = 5)
#' }
#'
#' @export
sharp_pathways <- function(tree, top = 10L, min_count = 1L) {
  stopifnot(inherits(tree, "transitrees"))
  out <- tree_pathways(tree, min_count = min_count, sort_by = "count")
  out <- out[order(-out$next_probability), , drop = FALSE]
  rownames(out) <- NULL
  utils::head(out, n = as.integer(top))
}
