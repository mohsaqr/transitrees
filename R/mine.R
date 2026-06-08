# ---- Mining contexts and sequences ----

#' Mine Contexts by Next-State Probability
#'
#' @description
#' Scan every context in a fitted tree for a chosen next state and
#' return those whose predicted probability for that state falls in a
#' requested range. A tidy context-mining table:
#' "in which histories is the next move \code{state} unusually likely or unlikely?"
#'
#' @param tree A \code{transitiontrees}.
#' @param state Character. The next state to score, one of the tree's
#'   alphabet.
#' @param min_prob,max_prob Numeric in \eqn{[0, 1]} or \code{NULL}.
#'   Keep contexts whose \code{P(state | context)} is at least
#'   \code{min_prob} and/or at most \code{max_prob}. \code{NULL}
#'   (default) leaves that side unbounded.
#' @param min_count Integer. Drop contexts with fewer than this many
#'   occurrences. Default 1.
#'
#' @return A data.frame with columns \code{pathway}, \code{depth},
#'   \code{count}, \code{state}, \code{prob} (\code{P(state | context)}),
#'   and \code{is_modal} (whether \code{state} is the context's most
#'   likely next state; ties broken by alphabet order), sorted by
#'   \code{prob} descending. The empty case returns a 0-row data.frame
#'   with the same schema.
#'
#' @examples
#' seqs <- replicate(60, sample(c("A", "B", "C"), 10, replace = TRUE),
#'                   simplify = FALSE)
#' tree <- context_tree(seqs, max_depth = 2L)
#' mine_contexts(tree, state = "A", min_prob = 0.4)
#' @export
mine_contexts <- function(tree, state, min_prob = NULL, max_prob = NULL,
                          min_count = 1L) {
  stopifnot(inherits(tree, "transitiontrees"))
  alpha <- tree$alphabet
  if (missing(state) || length(state) != 1L || !state %in% alpha)
    stop("'state' must be a single state from the alphabet: ",
         paste(alpha, collapse = ", "), call. = FALSE)
  sidx      <- match(state, alpha)
  min_count <- as.integer(min_count)

  empty <- data.frame(pathway = character(0), depth = integer(0),
                      count = numeric(0), state = character(0),
                      prob = numeric(0), is_modal = logical(0),
                      stringsAsFactors = FALSE)

  rows <- lapply(names(tree$nodes), function(ctx) {
    nd <- tree$nodes[[ctx]]
    if (nd$n < min_count) return(NULL)
    data.frame(
      pathway  = if (identical(ctx, .ROOT)) .ROOT_LABEL else ctx,
      depth    = nd$depth,
      count    = as.numeric(nd$n),   # numeric, matching the empty template
      state    = state,
      prob     = nd$prob[sidx],
      is_modal = which.max(nd$prob) == sidx,
      stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0L) return(empty)

  out <- do.call(rbind, rows)
  if (!is.null(min_prob)) out <- out[out$prob >= min_prob, , drop = FALSE]
  if (!is.null(max_prob)) out <- out[out$prob <= max_prob, , drop = FALSE]
  out <- out[order(-out$prob, -out$count), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Mine Sequences by Predictive Surprise
#'
#' @description
#' Rank held-out sequences by how well the fitted tree predicts them.
#' A tidy pattern-mining table: surface the
#' subsequences the model finds most \emph{surprising} (poor fit, high
#' perplexity) or most \emph{expected} (good fit, low perplexity).
#'
#' @param tree A \code{transitiontrees}.
#' @param newdata Sequence data in any format accepted by
#'   \code{context_tree()}.
#' @param n Integer. Number of sequences to return. Default 10.
#' @param which One of \code{"surprising"} (default; highest perplexity
#'   first) or \code{"expected"} (lowest perplexity first).
#'
#' @return A data.frame with the \code{\link{score_sequences}} columns
#'   (\code{sequence_id}, \code{n_scored}, \code{log_lik},
#'   \code{perplexity}), the top \code{n} by the chosen direction.
#'
#' @examples
#' fit  <- replicate(60, sample(c("A", "B", "C"), 10, replace = TRUE),
#'                   simplify = FALSE)
#' tree <- context_tree(fit, max_depth = 2L)
#' new  <- replicate(20, sample(c("A", "B", "C"), 10, replace = TRUE),
#'                   simplify = FALSE)
#' mine_sequences(tree, new, n = 5, which = "surprising")
#' @export
mine_sequences <- function(tree, newdata, n = 10L,
                           which = c("surprising", "expected")) {
  stopifnot(inherits(tree, "transitiontrees"))
  which <- match.arg(which)
  sc <- score_sequences(tree, newdata)
  if (nrow(sc) == 0L) return(sc)
  ord <- if (which == "surprising") order(-sc$perplexity) else
    order(sc$perplexity)
  out <- utils::head(sc[ord, , drop = FALSE], as.integer(n))
  rownames(out) <- NULL
  out
}
