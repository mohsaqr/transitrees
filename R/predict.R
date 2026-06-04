# ---- Prediction and generation ----

#' @noRd
.ct_match_context <- function(tree, history) {
  ## Walks history backwards, returning the deepest surviving context
  ## that matches the tail of `history`. Always falls back to root ("").
  history <- as.character(history)
  L <- length(history)
  if (L == 0L) return(.ROOT)
  for (k in seq.int(min(tree$max_depth, L), 1L)) {
    ctx <- paste(history[(L - k + 1L):L], collapse = " -> ")
    if (ctx %in% names(tree$nodes)) return(ctx)
  }
  .ROOT
}

#' Predict Next-State Probabilities from a Context Tree
#'
#' @param object A \code{transitiontrees}.
#' @param newdata Either (i) a list of character vectors (each is the
#'   "history" leading up to the prediction point), (ii) a wide
#'   data.frame / matrix whose rows are histories, or (iii) a single
#'   character vector treated as one history.
#' @param type One of \code{"prob"} (default; named numeric matrix of
#'   next-state probabilities) or \code{"class"} (character vector of
#'   modal predictions).
#' @param ... Ignored.
#'
#' @return If \code{type = "prob"}: a matrix with one row per history
#'   and one column per state. A list/data.frame/matrix \code{newdata}
#'   always returns a matrix (1 x k for a single-history container);
#'   a bare character vector returns a named vector for interactive
#'   convenience. If \code{type = "class"}: a character vector of
#'   modal predictions.
#'
#' @examples
#' \donttest{
#' seqs <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
#'                   simplify = FALSE)
#' tree <- context_tree(seqs, max_depth = 3)
#' predict(tree, newdata = list(c("A","B"), c("C","C","B")))
#' predict(tree, newdata = list(c("A","B")), type = "class")
#' predict(tree, newdata = c("A","B"))   # bare vector → named vector
#' }
#'
#' @export
predict.transitiontrees <- function(object, newdata, type = c("prob", "class"),
                            ...) {
  type <- match.arg(type)
  ## A bare character vector is the interactive shortcut; a character
  ## *matrix* is a container (it has a dim) and must yield a matrix, so
  ## exclude it here even though is.character() is TRUE for it.
  bare_vec <- is.character(newdata) && !is.list(newdata) &&
    is.null(dim(newdata))

  histories <- if (is.list(newdata) && !is.data.frame(newdata)) {
    lapply(newdata, as.character)
  } else if (is.data.frame(newdata) || is.matrix(newdata)) {
    m <- as.matrix(newdata); storage.mode(m) <- "character"
    lapply(seq_len(nrow(m)), function(i) {
      r <- m[i, ]
      r[!is.na(r) & nzchar(r)]
    })
  } else if (bare_vec) {
    list(newdata)
  } else {
    stop("'newdata' must be a list, data.frame, matrix, or character vector.",
         call. = FALSE)
  }

  probs <- t(vapply(histories, function(h) {
    ctx  <- .ct_match_context(object, h)
    info <- object$nodes[[ctx]]
    if (is.null(info)) {
      ## Can only happen if root is missing — return uniform
      rep(1 / length(object$alphabet), length(object$alphabet))
    } else info$prob
  }, numeric(length(object$alphabet))))
  colnames(probs) <- object$alphabet

  if (type == "class") {
    return(object$alphabet[max.col(probs, ties.method = "first")])
  }
  ## Schema-stable: container inputs always yield a matrix. Only a
  ## bare character vector — an interactive shortcut — collapses to
  ## a named vector.
  if (bare_vec) return(setNames(as.vector(probs), object$alphabet))
  probs
}

#' Sample Sequences from a Fitted Context Tree
#'
#' @param tree A \code{transitiontrees}.
#' @param n Integer. Number of sequences to sample.
#' @param length Integer. Length of each sampled sequence.
#' @param start NULL or character vector. If NULL (default), each
#'   sequence starts from the root marginal; otherwise must have
#'   length \code{n} and supply the first state of each sequence.
#'
#' @return A character matrix of dimension \code{n x length}.
#'
#' @examples
#' \donttest{
#' seqs <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
#'                   simplify = FALSE)
#' tree <- context_tree(seqs, max_depth = 3)
#' generate_sequences(tree, n = 5, length = 10)
#' }
#'
#' @export
generate_sequences <- function(tree, n = 5L, length = 10L,
                               start = NULL) {
  stopifnot(inherits(tree, "transitiontrees"))
  ## base::length() because the `length` argument shadows the function.
  if (!is.numeric(n) || base::length(n) != 1L || is.na(n) || n < 1)
    stop("'n' must be a single positive integer.", call. = FALSE)
  if (!is.numeric(length) || base::length(length) != 1L ||
      is.na(length) || length < 1)
    stop("'length' must be a single positive integer.", call. = FALSE)
  n <- as.integer(n); length <- as.integer(length)
  alpha <- tree$alphabet

  draw <- function(p) sample(alpha, 1L, prob = p)

  if (is.null(start)) {
    start <- replicate(n, draw(tree$nodes[[.ROOT]]$prob))
  } else {
    start <- as.character(start)
    if (length(start) != n)
      stop("'start' must have length 'n'.", call. = FALSE)
  }

  out <- matrix(NA_character_, n, length)
  out[, 1L] <- start
  ## length == 1 is the start column alone; seq.int(2L, 1L) counts *down*,
  ## so only iterate when there is a second position to fill.
  if (length >= 2L) {
    for (t in seq.int(2L, length)) {
      for (i in seq_len(n)) {
        hist <- out[i, seq_len(t - 1L)]
        ctx  <- .ct_match_context(tree, hist)
        info <- tree$nodes[[ctx]]
        out[i, t] <- draw(info$prob)
      }
    }
  }
  out
}

#' Simulate Sequences from a Fitted Pathtree
#'
#' @description
#' S3 \code{simulate()} method for \code{transitiontrees} objects. Wraps
#' \code{\link{generate_sequences}} with the standard \code{nsim}
#' argument name and an optional \code{seed} (set via
#' \code{set.seed()} when supplied).
#'
#' @param object A \code{transitiontrees}.
#' @param nsim Integer. Number of sequences to simulate. Default 5.
#' @param seed Integer or \code{NULL}. Optional RNG seed.
#' @param length Integer. Length of each simulated sequence.
#' @param start NULL or character vector. Optional first state for each
#'   sequence; see \code{\link{generate_sequences}}.
#' @param ... Ignored.
#'
#' @return A character matrix of dimension \code{nsim x length}.
#' @export
simulate.transitiontrees <- function(object, nsim = 5L, seed = NULL,
                              length = 10L, start = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  generate_sequences(object, n = nsim, length = length, start = start)
}
