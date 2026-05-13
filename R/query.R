# ---- query_pathway() / subtree() / pathway_exists(): introspection ----

#' @noRd
.pt_normalise_pathway <- function(pathway) {
  ## Accept either a character vector of states or a single arrow-
  ## notation string. Returns a single canonical " -> " string. The
  ## special token "(root)" or "<root>" maps to .ROOT.
  if (is.null(pathway) || length(pathway) == 0L) return(.ROOT)
  if (length(pathway) == 1L && is.character(pathway)) {
    p <- pathway
    if (identical(p, "(root)") || identical(p, .ROOT)) return(.ROOT)
    return(trimws(p))
  }
  paste(as.character(pathway), collapse = " -> ")
}

#' Query the Probability of a Specific Pathway -> Next State
#'
#' @description
#' Returns the probability the fitted tree assigns to a given pathway
#' / next-state pair. Two lookup modes:
#' \itemize{
#'   \item \code{exact = TRUE}: the pathway must appear as a node;
#'     otherwise returns \code{NA}.
#'   \item \code{exact = FALSE} (default): if the pathway is missing,
#'     falls back to the longest matching suffix that *is* in the tree
#'     (mirrors \code{predict.pathtree()}).
#' }
#'
#' @param tree A \code{pathtree}.
#' @param pathway Character. The conditioning pathway, either as a
#'   single arrow-notation string ("A -> B -> C") or as a character
#'   vector of states (\code{c("A","B","C")}).
#' @param next_state Character. Next-state symbol to query, or
#'   \code{NULL} (default) to return the full conditional distribution.
#' @param exact Logical. Default \code{FALSE} — fall back to longest
#'   matching suffix.
#'
#' @return If \code{next_state} is supplied, a numeric scalar. Otherwise
#'   a named numeric vector indexed by alphabet.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' m <- matrix(sample(c("A","B","C"), 200, TRUE), 20)
#' tr <- context_tree(m, max_depth = 2L, nmin = 3L)
#' query_pathway(tr, c("A","B"))
#' query_pathway(tr, "A -> B", next_state = "C")
#' }
#' @export
query_pathway <- function(tree, pathway, next_state = NULL,
                          exact = FALSE) {
  stopifnot(inherits(tree, "pathtree"))
  key <- .pt_normalise_pathway(pathway)

  if (key %in% names(tree$nodes)) {
    p <- tree$nodes[[key]]$prob
  } else if (isTRUE(exact)) {
    p <- rep(NA_real_, length(tree$alphabet))
  } else {
    ## Use longest-matching-suffix fall-back.
    parts <- if (identical(key, .ROOT)) character(0) else
      strsplit(key, " -> ", fixed = TRUE)[[1L]]
    ctx <- .ct_match_context(tree, parts)
    p   <- tree$nodes[[ctx]]$prob
  }
  names(p) <- tree$alphabet

  if (is.null(next_state)) return(p)
  idx <- match(as.character(next_state), tree$alphabet)
  if (is.na(idx))
    stop("'next_state' (", next_state,
         ") is not in the tree's alphabet.", call. = FALSE)
  p[[idx]]
}

#' Test Whether a Pathway Exists in the Tree
#'
#' @param tree A \code{pathtree}.
#' @param pathway Character. Pathway as arrow-notation string or
#'   character vector.
#'
#' @return Logical scalar.
#'
#' @examples
#' \donttest{
#' tr <- context_tree(matrix(sample(c("A","B"), 50, TRUE), 5),
#'                    max_depth = 2L, nmin = 1L)
#' pathway_exists(tr, "A")
#' }
#' @export
pathway_exists <- function(tree, pathway) {
  stopifnot(inherits(tree, "pathtree"))
  key <- .pt_normalise_pathway(pathway)
  key %in% names(tree$nodes)
}

#' Extract the Subtree Rooted at a Pathway
#'
#' @description
#' Returns a new \code{pathtree} containing only the queried node and
#' its descendants. Node names are kept absolute (so the original
#' pathway is preserved as a key), but the returned object has a
#' \code{local_root} attribute pointing at the queried pathway.
#'
#' @param tree A \code{pathtree}.
#' @param pathway Character. The root pathway (arrow-notation or
#'   character vector). Must exist in the tree.
#'
#' @return A new \code{pathtree} whose nodes and edges are restricted
#'   to descendants of \code{pathway}. The alphabet, smoothing, and
#'   other hyperparameters are copied unchanged. \code{attr(., "local_root")}
#'   carries the queried pathway in canonical form.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' m <- matrix(sample(c("A","B","C"), 200, TRUE), 20)
#' tr <- context_tree(m, max_depth = 2L, nmin = 3L)
#' sub <- subtree(tr, "A")
#' attr(sub, "local_root")
#' }
#' @export
subtree <- function(tree, pathway) {
  stopifnot(inherits(tree, "pathtree"))
  key <- .pt_normalise_pathway(pathway)
  if (!key %in% names(tree$nodes))
    stop("Pathway '", key, "' is not a node in the tree.", call. = FALSE)

  edges_by_parent <- .pt_children_of(tree)
  collect <- function(ctx) {
    children <- edges_by_parent[[ctx]]
    if (is.null(children) || length(children) == 0L) return(ctx)
    c(ctx, unlist(lapply(children, collect), use.names = FALSE))
  }
  keep <- collect(key)

  new_nodes <- tree$nodes[keep]
  new_edges <- tree$edges[tree$edges$parent %in% keep &
                           tree$edges$child  %in% keep, , drop = FALSE]
  rownames(new_edges) <- NULL

  new_tree <- tree
  new_tree$nodes <- new_nodes
  new_tree$edges <- new_edges
  attr(new_tree, "local_root") <- key
  new_tree
}
