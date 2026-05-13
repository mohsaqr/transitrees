# ---- Pathway-centric visualisations ----

#' Plot Pathways as a Probability Heatmap
#'
#' @description
#' Heatmap visualisation of the pathway table: rows are pathways
#' (sorted), columns are next-state probabilities under the fitted
#' tree. The modal next state of each row is annotated in bold; rows
#' whose modal next state \emph{flips} relative to their parent
#' pathway are flagged in the row labels (with a leading caret
#' \code{>}). A side strip on the left encodes the pathway count on a
#' log scale.
#'
#' This is the natural pathway-focused visualisation: one glance shows
#' which pathways are common, which are sharp (high mass on a single
#' next state), which are diffuse (mass spread evenly), and which
#' carry trajectory-specific structure that order-1 misses.
#'
#' @param tree A \code{pathtree}.
#' @param top Integer. Maximum number of pathways to show. Default 20.
#' @param sort_by Character. One of \code{"count"} (default), \code{"KL"},
#'   or \code{"depth"}.
#' @param min_count Integer. Drop pathways with fewer than this many
#'   occurrences. Default 5.
#' @param show_flips Logical. Mark modal-flip pathways with a leading
#'   caret in the label. Default \code{TRUE}.
#' @param title Character. Plot title.
#' @param ... Ignored.
#' @return A ggplot object.
#'
#' @examples
#' \donttest{
#' seqs <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
#'                   simplify = FALSE)
#' tree <- context_tree(seqs, max_depth = 3)
#' plot_pathways(tree)
#' plot_pathways(tree, sort_by = "KL", top = 12)
#' }
#'
#' @export
plot_pathways <- function(tree,
                          top = 20L,
                          sort_by = c("count", "KL", "depth"),
                          min_count = 5L,
                          show_flips = TRUE,
                          title = NULL, ...) {
  stopifnot(inherits(tree, "pathtree"))
  sort_by <- match.arg(sort_by)

  pw <- pathtree_pathways(tree, min_count = min_count, sort_by = sort_by,
                          decreasing = TRUE)
  pw <- pw[pw$pathway != "(root)", , drop = FALSE]
  pw <- utils::head(pw, n = as.integer(top))
  if (nrow(pw) == 0L)
    stop("No pathways meet the threshold.", call. = FALSE)

  ## Build the long form: one row per (pathway, next_state)
  states <- tree$alphabet
  long <- do.call(rbind, lapply(seq_len(nrow(pw)), function(i) {
    ctx <- pw$pathway[i]
    info <- tree$nodes[[ctx]]
    data.frame(
      pathway    = ctx,
      next_state = factor(states, levels = states),
      prob       = info$prob,
      stringsAsFactors = FALSE
    )
  }))

  ## Display label: prefix flippers with a caret
  pw$label <- ifelse(isTRUE(show_flips) & !is.na(pw$flips) & pw$flips,
                     sprintf("> %s", pw$pathway),
                     pw$pathway)
  ## Reverse ordering so highest-ranked appears at top of the y axis
  pw <- pw[seq.int(nrow(pw), 1L), , drop = FALSE]
  long$label <- factor(pw$label[match(long$pathway, pw$pathway)],
                        levels = pw$label)

  ## Highlight modal cells
  long$is_modal <- mapply(function(pa, ns) {
    info <- tree$nodes[[pa]]
    states[which.max(info$prob)] == as.character(ns)
  }, long$pathway, long$next_state)

  if (is.null(title)) {
    sort_label <- switch(sort_by,
                         count = "by frequency",
                         KL    = "by KL from shorter history",
                         depth = "by length (deepest first)")
    title <- sprintf("Top %d pathways %s", nrow(pw), sort_label)
  }

  ggplot2::ggplot(long, ggplot2::aes(x = .data$next_state,
                                     y = .data$label,
                                     fill = .data$prob)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::geom_text(
      data = long[long$is_modal, , drop = FALSE],
      ggplot2::aes(label = sprintf("%.2f", .data$prob)),
      colour = "white", fontface = "bold", size = 3.2) +
    ggplot2::geom_text(
      data = long[!long$is_modal & long$prob >= 0.05, , drop = FALSE],
      ggplot2::aes(label = sprintf("%.2f", .data$prob)),
      colour = "grey20", size = 3.0) +
    ggplot2::scale_fill_gradient(
      low = "#deebf7", high = "#08519c",
      name = "P(next | pathway)",
      limits = c(0, 1)) +
    ggplot2::labs(
      x = "Next state", y = NULL,
      title = title,
      subtitle = sprintf(
        "%s; n_min = %d; %d flips marked with leading >",
        switch(sort_by,
               count = "ranked by pathway frequency",
               KL    = "ranked by KL from (k-1)-suffix",
               depth = "ranked by depth"),
        min_count, sum(pw$flips, na.rm = TRUE))) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(family = "mono"),
      plot.title = ggplot2::element_text(face = "bold")
    )
}


#' Lollipop Chart of Pathway Divergence
#'
#' @description
#' Lollipop chart of per-pathway Kullback-Leibler divergence from the
#' (k-1)-suffix. Point size is proportional to pathway count; orange
#' points mark pathways whose modal next state flips between orders.
#' Annotates each flip with the prediction change, e.g. "Disengaged
#' -> Active".
#'
#' @param tree A \code{pathtree}.
#' @param top Integer. Number of pathways to show. Default 15.
#' @param min_count Integer. Drop pathways below this count. Default 5.
#' @param title Character. Plot title.
#' @param ... Ignored.
#' @return A ggplot object.
#'
#' @export
plot_divergence <- function(tree, top = 15L, min_count = 5L,
                            title = NULL, ...) {
  stopifnot(inherits(tree, "pathtree"))
  pw <- divergent_pathways(tree, n = top, min_count = min_count)
  if (nrow(pw) == 0L) stop("No pathways meet the threshold.", call. = FALSE)
  pw$pathway <- factor(pw$pathway, levels = rev(pw$pathway))

  if (is.null(title))
    title <- sprintf("Top %d pathways by KL divergence", nrow(pw))

  ggplot2::ggplot(pw, ggplot2::aes(x = .data$KL, y = .data$pathway)) +
    ggplot2::geom_segment(ggplot2::aes(xend = 0, yend = .data$pathway),
                          colour = "grey70") +
    ggplot2::geom_point(ggplot2::aes(size = .data$count,
                                     colour = .data$flips)) +
    ggplot2::geom_text(
      data = pw[!is.na(pw$flips) & pw$flips, ],
      ggplot2::aes(label = sprintf("%.0f%% %s",
                                   100 * .data$prob_next,
                                   .data$modal_next)),
      hjust = -0.15, size = 3.2, colour = "#D55E00") +
    ggplot2::scale_colour_manual(
      values = c(`FALSE` = "#0072B2", `TRUE` = "#D55E00"),
      labels = c("modal next-state same",
                 "modal next-state flips"),
      name = NULL,
      na.value = "grey60") +
    ggplot2::scale_size_continuous(name = "pathway count") +
    ggplot2::labs(
      x = "KL divergence (bits)  vs (k-1)-suffix",
      y = NULL,
      title = title) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )
}
