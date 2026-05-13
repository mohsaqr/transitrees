# ---- Pathway-centric API ----

#' @noRd
.pt_empty_pathways_df <- function() {
  data.frame(
    pathway    = character(0),
    depth      = integer(0),
    count      = numeric(0),
    modal_next = character(0),
    prob_next  = numeric(0),
    KL         = numeric(0),
    flips      = logical(0),
    stringsAsFactors = FALSE
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
#' @param tree A \code{pathtree} object.
#' @param min_count Integer. Drop pathways with fewer than this many
#'   occurrences. Default 1.
#' @param sort_by Character. One of \code{"count"} (default),
#'   \code{"KL"}, or \code{"depth"}. Sorts the returned data.frame.
#' @param decreasing Logical. Default \code{TRUE}.
#' @param ... Ignored.
#'
#' @return A data.frame with columns \code{pathway} (arrow notation,
#'   e.g. \code{"A -> B -> C"}; the root is reported as \code{"(root)"}),
#'   \code{depth}, \code{count}, \code{modal_next}, \code{prob_next}
#'   (probability of the modal next state), \code{KL} (divergence from
#'   the parent context's prediction, in bits; \code{NA} for the root),
#'   \code{flips} (logical, did the modal next state change vs the
#'   parent context?). The empty case returns a 0-row data.frame with
#'   the same schema.
#'
#' @details
#' Each row is a pathway -- a (possibly empty) sequence of states ending
#' at the point where a prediction is made. The KL column quantifies
#' how much more information the pathway carries than its
#' (k-1)-suffix in bits. \code{flips = TRUE} marks pathways where the
#' longer history changes which next state is most likely.
#'
#' @examples
#' \donttest{
#' seqs <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
#'                   simplify = FALSE)
#' tree <- context_tree(seqs, max_depth = 3)
#' pathtree_pathways(tree)
#' }
#'
#' @export
pathtree_pathways <- function(tree, min_count = 1L,
                              sort_by = c("count", "KL", "depth"),
                              decreasing = TRUE, ...) {
  stopifnot(inherits(tree, "pathtree"))
  sort_by   <- match.arg(sort_by)
  min_count <- as.integer(min_count)

  rows <- lapply(names(tree$nodes), function(ctx) {
    info <- tree$nodes[[ctx]]
    if (info$n < min_count) return(NULL)

    pathway <- if (identical(ctx, .ROOT)) "(root)" else ctx
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
      pathway    = pathway,
      depth      = info$depth,
      count      = info$n,
      modal_next = modal_next,
      prob_next  = prob_next,
      KL         = KL,
      flips      = flips,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(out) || nrow(out) == 0L) return(.pt_empty_pathways_df())

  key <- switch(sort_by,
                count = out$count,
                KL    = ifelse(is.na(out$KL), -Inf, out$KL),
                depth = out$depth)
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
#' @param tree A \code{pathtree}.
#' @param n Integer. Number of pathways to return. Default 10.
#' @param depth Integer or NULL. Restrict to pathways of this exact
#'   depth. \code{NULL} (default) keeps all depths.
#' @param min_count Integer. Minimum count cut-off. Default 1.
#'
#' @return A data.frame, same columns as \code{\link{pathtree_pathways}},
#'   sorted by count descending.
#'
#' @examples
#' \donttest{
#' seqs <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
#'                   simplify = FALSE)
#' tree <- context_tree(seqs, max_depth = 3)
#' common_pathways(tree, n = 8)
#' common_pathways(tree, n = 8, depth = 3L)   # restrict to depth-3
#' }
#'
#' @export
common_pathways <- function(tree, n = 10L, depth = NULL, min_count = 1L) {
  stopifnot(inherits(tree, "pathtree"))
  out <- pathtree_pathways(tree, min_count = min_count, sort_by = "count",
                           decreasing = TRUE)
  if (!is.null(depth))
    out <- out[out$depth == as.integer(depth), , drop = FALSE]
  utils::head(out, n = as.integer(n))
}

#' Most Predictively Divergent Pathways
#'
#' @description
#' Returns the top \code{n} pathways by Kullback-Leibler divergence
#' from their (k-1)-suffix. These are the pathways whose extended
#' history adds the most predictive information over the shorter one.
#' Pathways whose modal next state actually flips between orders are
#' marked in the \code{flips} column.
#'
#' @param tree A \code{pathtree}.
#' @param n Integer. Number of pathways to return. Default 10.
#' @param min_count Integer. Drop pathways with fewer than this many
#'   occurrences. Default 1.
#' @param flips_only Logical. If \code{TRUE}, return only pathways
#'   that flip the modal next state between orders. Default \code{FALSE}.
#'
#' @return A data.frame, same columns as \code{\link{pathtree_pathways}},
#'   sorted by KL descending.
#'
#' @examples
#' \donttest{
#' seqs <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
#'                   simplify = FALSE)
#' tree <- context_tree(seqs, max_depth = 3)
#' divergent_pathways(tree, n = 6)
#' divergent_pathways(tree, flips_only = TRUE)
#' }
#'
#' @export
divergent_pathways <- function(tree, n = 10L, min_count = 1L,
                               flips_only = FALSE) {
  stopifnot(inherits(tree, "pathtree"))
  out <- pathtree_pathways(tree, min_count = min_count, sort_by = "KL",
                           decreasing = TRUE)
  out <- out[!is.na(out$KL), , drop = FALSE]
  if (isTRUE(flips_only))
    out <- out[!is.na(out$flips) & out$flips, , drop = FALSE]
  utils::head(out, n = as.integer(n))
}

#' Sharpest Pathways
#'
#' @description
#' Returns the top \code{n} pathways by predictive sharpness -- the
#' probability mass on their modal next state. High values indicate
#' strongly deterministic continuations; low values indicate
#' ambiguous next-state distributions.
#'
#' @param tree A \code{pathtree}.
#' @param n Integer. Number of pathways to return. Default 10.
#' @param min_count Integer. Drop pathways with fewer than this many
#'   occurrences. Default 1.
#'
#' @return A data.frame, same columns as \code{\link{pathtree_pathways}},
#'   sorted by \code{prob_next} descending.
#'
#' @examples
#' \donttest{
#' seqs <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
#'                   simplify = FALSE)
#' tree <- context_tree(seqs, max_depth = 3)
#' sharp_pathways(tree, n = 5)
#' }
#'
#' @export
sharp_pathways <- function(tree, n = 10L, min_count = 1L) {
  stopifnot(inherits(tree, "pathtree"))
  out <- pathtree_pathways(tree, min_count = min_count, sort_by = "count")
  out <- out[order(-out$prob_next), , drop = FALSE]
  rownames(out) <- NULL
  utils::head(out, n = as.integer(n))
}
