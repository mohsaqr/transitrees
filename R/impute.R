# ---- Imputation of missing states ----

#' Impute Missing States in Sequences
#'
#' @description
#' Fill the gaps in incomplete sequences using a fitted context tree:
#' each missing state is predicted from the longest matching context of
#' the states that precede it. Filling proceeds left to right, so a
#' just-imputed state becomes part of the context for later gaps.
#'
#' @details
#' Only \strong{internal} gaps are imputed. A run of trailing \code{NA}
#' / \code{""} cells (end-of-sequence padding in a wide frame) is left
#' untouched, since there is no observed state after it to mark the
#' sequence as continuing. A sequence that is entirely missing is
#' returned unchanged (there is nothing to condition on).
#'
#' @param tree A \code{transitiontrees}.
#' @param newdata Sequences with gaps: a list of character vectors, a
#'   character matrix / data.frame (one row per sequence, \code{NA} or
#'   \code{""} marking a gap), or a single character vector.
#' @param method One of \code{"modal"} (default; fill with the most
#'   likely state) or \code{"prob"} (sample from the predicted
#'   distribution).
#' @param seed Integer or \code{NULL}. Optional RNG seed, used only when
#'   \code{method = "prob"}.
#'
#' @return The same container shape as \code{newdata} (list, matrix,
#'   data.frame, or character vector) with internal gaps filled.
#'
#' @examples
#' seqs <- replicate(60, sample(c("A", "B", "C"), 8, replace = TRUE),
#'                   simplify = FALSE)
#' tree <- context_tree(seqs, max_depth = 2L)
#' gappy <- list(c("A", NA, "C"), c("B", "B", NA, "A"))
#' impute_sequences(tree, gappy)
#' @export
impute_sequences <- function(tree, newdata,
                             method = c("modal", "prob"), seed = NULL) {
  stopifnot(inherits(tree, "transitiontrees"))
  method <- match.arg(method)
  if (!is.null(seed)) {
    ## seed locally without clobbering the caller's RNG stream
    if (exists(".Random.seed", envir = globalenv())) {
      old_seed <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()))
    } else {
      on.exit(rm(".Random.seed", envir = globalenv()))
    }
    set.seed(seed)
  }
  alpha <- tree$alphabet

  fill_one <- function(seq_chr) {
    seq_chr <- as.character(seq_chr)
    present <- !is.na(seq_chr) & nzchar(seq_chr)
    if (all(present)) return(seq_chr)
    last_obs <- if (any(present)) max(which(present)) else 0L
    if (last_obs == 0L) return(seq_chr)   # nothing to condition on
    ## Left-to-right with a true sequential dependency: a filled state
    ## becomes context for the gaps that follow it.
    for (t in seq_len(last_obs)) {
      if (present[t]) next
      hist <- seq_chr[seq_len(t - 1L)]
      hist <- hist[!is.na(hist) & nzchar(hist)]
      p    <- tree$nodes[[.ct_match_context(tree, hist)]]$prob
      seq_chr[t] <- if (method == "modal")
        alpha[which.max(p)] else sample(alpha, 1L, prob = p)
      present[t] <- TRUE
    }
    seq_chr
  }

  if (is.list(newdata) && !is.data.frame(newdata)) {
    return(lapply(newdata, fill_one))
  }
  if (is.matrix(newdata) || is.data.frame(newdata)) {
    m <- as.matrix(newdata); storage.mode(m) <- "character"
    filled <- t(vapply(seq_len(nrow(m)), function(i) fill_one(m[i, ]),
                       character(ncol(m))))
    dimnames(filled) <- dimnames(m)
    if (is.data.frame(newdata))
      return(as.data.frame(filled, stringsAsFactors = FALSE))
    return(filled)
  }
  if (is.character(newdata) && is.null(dim(newdata))) {
    return(fill_one(newdata))
  }
  stop("'newdata' must be a list, matrix, data.frame, or character vector.",
       call. = FALSE)
}
