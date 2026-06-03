# ---- compare_trees(): two-tree divergence + permutation test ----

#' @noRd
.pt_lookup_prob <- function(tree, ctx) {
  ## Lookup prob at ctx; fall back to longest-matching-suffix.
  p <- if (ctx %in% names(tree$nodes)) {
    tree$nodes[[ctx]]$prob
  } else {
    parts <- if (identical(ctx, .ROOT)) character(0) else
      strsplit(ctx, " -> ", fixed = TRUE)[[1L]]
    matched <- .ct_match_context(tree, parts)
    tree$nodes[[matched]]$prob
  }
  stats::setNames(as.numeric(p), tree$alphabet)
}

#' @noRd
.pt_lookup_count <- function(tree, ctx) {
  if (ctx %in% names(tree$nodes)) return(tree$nodes[[ctx]]$n)
  0L
}

#' @noRd
.pt_distance_breakdown <- function(tree_a, tree_b) {
  ctx_union <- unique(c(names(tree_a$nodes), names(tree_b$nodes)))
  alphabet <- tree_a$alphabet
  p_a <- lapply(ctx_union, function(c) .pt_lookup_prob(tree_a, c)[alphabet])
  p_b <- lapply(ctx_union, function(c) .pt_lookup_prob(tree_b, c)[alphabet])
  n_a <- vapply(ctx_union, function(c) .pt_lookup_count(tree_a, c),
                numeric(1))
  n_b <- vapply(ctx_union, function(c) .pt_lookup_count(tree_b, c),
                numeric(1))
  KL_ab <- mapply(.ct_kl, p_a, p_b, USE.NAMES = FALSE)
  KL_ba <- mapply(.ct_kl, p_b, p_a, USE.NAMES = FALSE)
  data.frame(
    pathway        = ifelse(ctx_union == .ROOT, .ROOT_LABEL, ctx_union),
    count_a        = n_a, count_b = n_b,
    divergence_ab  = KL_ab,
    divergence_ba  = KL_ba,
    divergence_sym = 0.5 * (KL_ab + KL_ba),
    stringsAsFactors = FALSE
  )
}

#' @noRd
.pt_check_alphabet_compatible <- function(tree_a, tree_b) {
  if (!setequal(tree_a$alphabet, tree_b$alphabet))
    stop("Trees have incompatible alphabets.", call. = FALSE)
}

#' @noRd
.pt_reapply_pruning <- function(tree, template) {
  if (!isTRUE(template$pruned)) return(tree)
  pr <- template$pruning
  if (is.null(pr) || is.null(pr$criterion))
    stop("Pruned comparison tree is missing pruning metadata.",
         call. = FALSE)
  prune_tree(tree, criterion = pr$criterion,
                 alpha = pr$alpha %||% 0.05,
                 threshold = pr$threshold %||% 0.005)
}

#' Symmetric KL Distance Between Two Pathtrees
#'
#' @description
#' Bare-metal scalar: a count-weighted average of per-context symmetric
#' Kullback-Leibler divergence over the union of pathways present in
#' either tree. No null distribution.
#'
#' @param tree_a,tree_b Pathtrees with matching alphabets.
#' @param symmetric Logical. \code{TRUE} (default) returns
#'   \eqn{0.5(D_{KL}(A\|B) + D_{KL}(B\|A))}; \code{FALSE} returns
#'   \eqn{D_{KL}(A\|B)} only.
#' @return Numeric scalar.
#' @examples
#' \donttest{
#' set.seed(1)
#' m1 <- matrix(sample(c("A","B","C"), 200, TRUE), 20)
#' m2 <- matrix(sample(c("A","B","C"), 200, TRUE), 20)
#' tr1 <- context_tree(m1, max_depth = 2L, min_count = 3L)
#' tr2 <- context_tree(m2, max_depth = 2L, min_count = 3L)
#' tree_distance(tr1, tr2)
#' }
#' @export
tree_distance <- function(tree_a, tree_b, symmetric = TRUE) {
  stopifnot(inherits(tree_a, "transitrees"),
            inherits(tree_b, "transitrees"))
  .pt_check_alphabet_compatible(tree_a, tree_b)

  bd <- .pt_distance_breakdown(tree_a, tree_b)
  per <- if (isTRUE(symmetric)) bd$divergence_sym else bd$divergence_ab
  w   <- bd$count_a + bd$count_b
  total <- sum(w)
  if (total == 0) return(0)
  ## Drop any infinite contributions from zero weights to avoid NaN.
  msk <- w > 0
  sum(per[msk] * w[msk]) / total
}

#' Compare Two Pathtrees by Symmetric Divergence with Permutation Test
#'
#' @description
#' Computes the count-weighted symmetric Kullback-Leibler divergence
#' between two fitted transitreess, then provides a reference distribution
#' by permuting sequence-to-tree assignments.
#'
#' Use this to ask: do two cohorts (group A vs. group B, baseline vs.
#' intervention) generate significantly different pathway dynamics?
#'
#' @param tree_a,tree_b Pathtrees fit on data subsets A and B.
#'   Alternatively, pass a two-element \code{transitrees_group} (from
#'   \code{context_tree(..., group =)}) as \code{tree_a} and leave
#'   \code{tree_b = NULL}; its two trees are compared in key order.
#' @param iter Integer. Number of permutations. Default 200.
#' @param seed Integer. RNG seed for reproducibility. Default 1.
#' @param symmetric Logical. Default \code{TRUE}.
#'
#' @return A \code{transitrees_comparison} object with components:
#' \describe{
#'   \item{pdist}{observed scalar distance}
#'   \item{null_dist}{numeric vector, length \code{iter}}
#'   \item{p_value}{one-sided p-value (proportion of null at least as
#'     extreme as observed)}
#'   \item{pathways}{per-pathway breakdown data.frame}
#' }
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' m1 <- matrix(sample(c("A","B","C"), 200, TRUE), 20)
#' m2 <- matrix(sample(c("A","B","C"), 200, TRUE), 20)
#' tr1 <- context_tree(m1, max_depth = 2L, min_count = 3L)
#' tr2 <- context_tree(m2, max_depth = 2L, min_count = 3L)
#' compare_trees(tr1, tr2, iter = 50)
#' }
#' @export
compare_trees <- function(tree_a, tree_b = NULL, iter = 200L,
                              seed = 1L, symmetric = TRUE) {
  ## A 2-group transitrees_group as the sole argument -> compare its pair.
  if (inherits(tree_a, "transitrees_group")) {
    if (!is.null(tree_b))
      stop("Pass a 2-group 'transitrees_group' as the only tree argument, ",
           "not a group plus 'tree_b'.", call. = FALSE)
    if (length(tree_a) != 2L)
      stop("compare_trees() is pairwise; the transitrees_group must ",
           "have exactly 2 groups (got ", length(tree_a), ").",
           call. = FALSE)
    grp    <- tree_a
    tree_a <- grp[[1L]]
    tree_b <- grp[[2L]]
  }
  if (is.null(tree_b))
    stop("'tree_b' is required (or pass a 2-group 'transitrees_group' as ",
         "the first argument).", call. = FALSE)
  stopifnot(inherits(tree_a, "transitrees"),
            inherits(tree_b, "transitrees"))
  if (!is.numeric(iter) || length(iter) != 1L ||
      is.na(iter) || iter < 1)
    stop("'iter' must be a single positive integer.", call. = FALSE)
  .pt_check_alphabet_compatible(tree_a, tree_b)

  pdist_obs <- tree_distance(tree_a, tree_b, symmetric = symmetric)
  bd        <- .pt_distance_breakdown(tree_a, tree_b)

  pooled  <- c(tree_a$data, tree_b$data)
  n_a     <- length(tree_a$data)
  n_total <- length(pooled)
  alphabet <- tree_a$alphabet

  ## Refit using each observed tree's hyperparameters and pruning state.
  set.seed(seed)
  null_dist <- vapply(seq_len(iter), function(b) {
    perm  <- sample.int(n_total)
    a_idx <- perm[seq_len(n_a)]
    b_idx <- perm[(n_a + 1L):n_total]
    tr_a  <- context_tree(pooled[a_idx],
                          max_depth = tree_a$max_depth,
                          min_count = tree_a$nmin,
                          smoothing = tree_a$smoothing,
                          alphabet  = alphabet)
    tr_b  <- context_tree(pooled[b_idx],
                          max_depth = tree_b$max_depth,
                          min_count = tree_b$nmin,
                          smoothing = tree_b$smoothing,
                          alphabet  = alphabet)
    tr_a <- .pt_reapply_pruning(tr_a, tree_a)
    tr_b <- .pt_reapply_pruning(tr_b, tree_b)
    tree_distance(tr_a, tr_b, symmetric = symmetric)
  }, numeric(1))

  p_value <- (1 + sum(null_dist >= pdist_obs)) / (iter + 1)

  structure(
    list(pdist     = pdist_obs,
         null_dist = null_dist,
         p_value   = p_value,
         pathways  = bd[order(-bd$divergence_sym), , drop = FALSE]),
    class = "transitrees_comparison"
  )
}

#' Plot a Pathtree Comparison
#'
#' @description
#' Visualises the permutation-test result. Histogram of the null
#' distribution of transitrees distances under shuffled group labels;
#' a vertical line marks the observed distance; the panel header
#' carries the p-value.
#'
#' @param x A \code{transitrees_comparison} object.
#' @param bins Integer. Histogram bins. Default 30.
#' @param ... Ignored.
#' @return A ggplot object.
#' @export
plot.transitrees_comparison <- function(x, bins = 30L, ...) {
  null_df <- data.frame(d = x$null_dist)
  ggplot2::ggplot(null_df, ggplot2::aes(x = .data$d)) +
    ggplot2::geom_histogram(bins = bins, fill = "#0072B2",
                             colour = "white", alpha = 0.85) +
    ggplot2::geom_vline(xintercept = x$pdist, colour = "#D55E00",
                        linewidth = 1.1) +
    ggplot2::annotate("text", x = x$pdist, y = Inf,
                       label = sprintf("observed = %.3f",
                                        x$pdist),
                       hjust = -0.05, vjust = 1.4, colour = "#D55E00",
                       size = 3.6, fontface = "bold") +
    ggplot2::labs(
      x = "symmetric KL between trees", y = "permutations",
      title = "Pathtree comparison: observed vs. null distribution",
      subtitle = sprintf("permutation p-value = %.3f  (n = %d)",
                          x$p_value, length(x$null_dist))) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40",
                                             size = 9),
      panel.grid.minor = ggplot2::element_blank()
    )
}

#' Coerce a Pathtree Comparison to a Tidy Data Frame
#'
#' @description
#' Uniform tidy-extract: returns the per-pathway divergence breakdown
#' (\code{object$pathways}) — the consumer-facing detail behind the
#' scalar \code{pdist} and the permutation \code{p_value}.
#'
#' @param x A \code{transitrees_comparison}.
#' @param row.names,optional Ignored.
#' @param ... Ignored.
#' @return A data.frame with columns \code{pathway}, \code{count_a},
#'   \code{count_b}, \code{divergence_ab}, \code{divergence_ba},
#'   \code{divergence_sym}.
#' @export
as.data.frame.transitrees_comparison <- function(x, row.names = NULL,
                                               optional = FALSE, ...) {
  x$pathways
}

#' @export
print.transitrees_comparison <- function(x, digits = 3L, n = 6L, ...) {
  cat(sprintf("<transitrees_comparison>  iter = %d\n",
              length(x$null_dist)))
  cat(sprintf("  observed distance : %s\n",
              format(x$pdist, digits = digits)))
  cat(sprintf("  null mean         : %s\n",
              format(mean(x$null_dist), digits = digits)))
  cat(sprintf("  p-value           : %s\n",
              format(x$p_value, digits = digits)))
  cat("\ntop divergent pathways:\n")
  bd <- utils::head(x$pathways, n)
  bd$divergence_ab  <- round(bd$divergence_ab, digits)
  bd$divergence_ba  <- round(bd$divergence_ba, digits)
  bd$divergence_sym <- round(bd$divergence_sym, digits)
  print.data.frame(bd, row.names = FALSE)
  invisible(x)
}
