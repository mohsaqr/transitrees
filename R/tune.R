# ---- tune_pathtree(): k-fold cross-validated hyperparameter selection ----

#' @noRd
.pt_make_folds <- function(n, k, seed = 1L) {
  set.seed(seed)
  ord  <- sample.int(n)
  fold <- ((seq_len(n) - 1L) %% k) + 1L
  split(ord, fold[seq_len(n)])
}

#' @noRd
.pt_fit_score_fold <- function(train_seqs, test_seqs, alphabet, cfg) {
  tr <- context_tree(train_seqs,
                     max_depth = cfg$max_depth,
                     nmin      = cfg$nmin,
                     smoothing = cfg$smoothing,
                     alphabet  = alphabet)
  if (isTRUE(cfg$prune))
    tr <- prune_pathtree(tr, criterion = "G2", alpha = cfg$prune_alpha)
  ll <- logLik.pathtree(tr, newdata = test_seqs)
  list(ll = as.numeric(ll), n = attr(ll, "nobs"),
       n_nodes = length(tr$nodes))
}

#' @noRd
.pt_smoothing_grid <- function(spec) {
  ## Normalise a user-supplied smoothing grid to a list of resolved
  ## smoothing specs. Accepts a character vector of method names, or a
  ## list of method-name-or-list elements.
  if (is.character(spec)) return(lapply(spec, .pt_resolve_smoothing))
  if (is.list(spec)) {
    if (!is.null(spec$method) || (length(spec) >= 1L &&
        is.character(spec[[1L]]) && length(spec[[1L]]) == 1L &&
        spec[[1L]] %in% names(.pt_smoothing_defaults)))
      return(list(.pt_resolve_smoothing(spec)))
    return(lapply(spec, .pt_resolve_smoothing))
  }
  stop("'smoothing' must be a character vector or a list of specs.",
       call. = FALSE)
}

#' @noRd
.pt_smoothing_label <- function(sm) {
  extras <- sm[setdiff(names(sm), "method")]
  if (length(extras) == 0L) return(sm$method)
  paste0(sm$method, "(",
         paste(names(extras), unlist(extras), sep = "=",
               collapse = ", "),
         ")")
}

#' Cross-Validated Hyperparameter Tuning for Pathtrees
#'
#' @description
#' Runs k-fold cross-validation over a grid of fitting and pruning
#' hyperparameters. Returns a data.frame ranked by held-out perplexity.
#' The configuration with minimum perplexity is exposed via
#' \code{attr(result, "best")}.
#'
#' Folds are at the sequence level (each fold holds out whole
#' sequences, not positions within sequences).
#'
#' @param data Sequence data; format accepted by \code{context_tree()}.
#' @param max_depth Integer vector. Grid values for tree depth.
#'   Default \code{2:5}.
#' @param nmin Integer vector. Grid for minimum-count threshold.
#'   Default \code{c(3L, 5L, 10L)}.
#' @param smoothing Smoothing grid. A character vector of method names
#'   (e.g. \code{c("floor", "kneser_ney")}) — each method is tried with
#'   its default hyperparameters — or a list of explicit specs (e.g.
#'   \code{list(list("floor", ymin = 0.001), list("floor", ymin = 0.005))}
#'   for a hyperparameter sweep within one method).
#' @param prune Logical vector. Whether to apply G^2 pruning.
#'   Default \code{c(FALSE, TRUE)}.
#' @param prune_alpha Numeric. Significance level for G^2 pruning when
#'   \code{prune = TRUE}. Default \code{0.05}.
#' @param k Integer. Number of CV folds. Default 5.
#' @param seed Integer. RNG seed for reproducible folds. Default 1.
#'
#' @return A \code{pathtree_tune} object: a data.frame with one row per
#'   grid point and columns \code{max_depth}, \code{nmin},
#'   \code{smoothing}, \code{prune}, \code{logLik}, \code{n_scored},
#'   \code{perplexity}, \code{n_nodes_avg}, sorted by \code{perplexity}
#'   ascending. \code{attr(result, "best")} carries the
#'   minimum-perplexity row.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' m <- matrix(sample(c("A","B","C"), 30 * 12, replace = TRUE), 30, 12)
#' tune_pathtree(m, max_depth = 1:3,
#'               smoothing = c("floor", "kneser_ney"),
#'               prune = FALSE, k = 4)
#' }
#' @export
tune_pathtree <- function(data,
                          max_depth   = 2L:5L,
                          nmin        = c(3L, 5L, 10L),
                          smoothing   = "floor",
                          prune       = c(FALSE, TRUE),
                          prune_alpha = 0.05,
                          k           = 5L,
                          seed        = 1L) {
  trajs <- .ct_traj(.ct_coerce(data))
  if (length(trajs) < k)
    stop("Not enough sequences (", length(trajs), ") for ", k, " folds.",
         call. = FALSE)
  alphabet <- .ct_alphabet(trajs)
  folds    <- .pt_make_folds(length(trajs), k, seed = seed)

  smoothing_grid <- .pt_smoothing_grid(smoothing)
  sm_labels      <- vapply(smoothing_grid, .pt_smoothing_label,
                            character(1))

  base <- expand.grid(max_depth   = as.integer(max_depth),
                      nmin        = as.integer(nmin),
                      smoothing_i = seq_along(smoothing_grid),
                      prune       = as.logical(prune),
                      stringsAsFactors = FALSE)

  scored <- lapply(seq_len(nrow(base)), function(gi) {
    sm <- smoothing_grid[[base$smoothing_i[gi]]]
    cfg <- list(max_depth   = base$max_depth[gi],
                nmin        = base$nmin[gi],
                smoothing   = sm,
                prune       = base$prune[gi],
                prune_alpha = prune_alpha)
    fold_results <- lapply(folds, function(test_idx) {
      train_idx <- setdiff(seq_along(trajs), test_idx)
      tryCatch(
        .pt_fit_score_fold(train_seqs = trajs[train_idx],
                            test_seqs = trajs[test_idx],
                            alphabet  = alphabet,
                            cfg       = cfg),
        error = function(e) list(ll = NA_real_, n = 0L,
                                  n_nodes = NA_integer_)
      )
    })
    ll_total    <- sum(vapply(fold_results, `[[`, numeric(1), "ll"),
                       na.rm = TRUE)
    n_total     <- sum(vapply(fold_results, `[[`, numeric(1), "n"),
                       na.rm = TRUE)
    n_nodes_avg <- mean(vapply(fold_results, `[[`, numeric(1), "n_nodes"),
                        na.rm = TRUE)
    pp <- if (n_total > 0) exp(-ll_total / n_total) else NA_real_
    data.frame(logLik      = ll_total,
               n_scored    = as.integer(n_total),
               perplexity  = pp,
               n_nodes_avg = n_nodes_avg,
               stringsAsFactors = FALSE)
  })
  metrics <- do.call(rbind, scored)

  out <- data.frame(
    max_depth  = base$max_depth,
    nmin       = base$nmin,
    smoothing  = sm_labels[base$smoothing_i],
    prune      = base$prune,
    metrics,
    stringsAsFactors = FALSE
  )
  ord <- order(out$perplexity, na.last = TRUE)
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL

  best <- out[1L, , drop = FALSE]
  structure(out, best = best, class = c("pathtree_tune", "data.frame"))
}

#' Plot a Pathtree CV Grid
#'
#' @description
#' Visualises the held-out perplexity surface returned by
#' \code{tune_pathtree()}. Lines track perplexity vs. \code{max_depth};
#' facets split by smoothing scheme and \code{prune}; colour encodes
#' \code{nmin}. The minimum-perplexity configuration is highlighted
#' with a star.
#'
#' @param x A \code{pathtree_tune} object.
#' @param ... Ignored.
#' @return A ggplot object.
#' @export
plot.pathtree_tune <- function(x, ...) {
  best <- attr(x, "best")
  df   <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  df$prune_label <- ifelse(df$prune, "pruned", "unpruned")
  df$nmin <- factor(df$nmin)

  facet <- if (length(unique(df$prune_label)) > 1L)
    ggplot2::facet_grid(prune_label ~ smoothing) else
    ggplot2::facet_wrap(~ smoothing)

  ggplot2::ggplot(df,
    ggplot2::aes(x = .data$max_depth, y = .data$perplexity,
                 colour = .data$nmin, group = .data$nmin)) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::geom_point(data = best,
                        ggplot2::aes(x = .data$max_depth,
                                     y = .data$perplexity),
                        shape = 8, size = 4, colour = "#D55E00",
                        inherit.aes = FALSE) +
    facet +
    ggplot2::labs(
      x = "max_depth", y = "held-out perplexity",
      colour = "nmin",
      title = "Pathtree CV grid",
      subtitle = "lower is better; star = minimum-perplexity configuration"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40", size = 9),
      legend.position = "bottom",
      strip.text = ggplot2::element_text(face = "bold")
    )
}

#' @export
print.pathtree_tune <- function(x, n = 10L, ...) {
  cat(sprintf("<pathtree_tune>  %d configurations\n", nrow(x)))
  print.data.frame(utils::head(x, n), row.names = FALSE)
  best <- attr(x, "best")
  if (!is.null(best)) {
    cat("\nbest (min perplexity):\n")
    print.data.frame(best, row.names = FALSE)
  }
  invisible(x)
}
