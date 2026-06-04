# ---- Distribution and predictive-diagnostic plots ----

#' Per-Context Next-State Distributions
#'
#' @description
#' Small-multiples bar chart of the full next-state distribution
#' \code{P(next | context)}, one panel per context. The bar-chart
#' analogue of the per-node probability display in \code{PST::plot()} /
#' \code{cplot()}: where \code{\link{plot_pathways}()} renders the same
#' numbers as a heatmap, this shows each context's distribution as its
#' own panel, with the modal bar highlighted.
#'
#' @param tree A \code{transitiontrees}.
#' @param contexts Character vector of pathway strings to show, or
#'   \code{NULL} (default) to take the \code{top} most frequent
#'   contexts.
#' @param top Integer. Number of contexts to show when \code{contexts}
#'   is \code{NULL}. Default 12.
#' @param min_count Integer. Drop contexts below this count. Default 1.
#'
#' @return A ggplot object.
#'
#' @examples
#' \donttest{
#' seqs <- replicate(60, sample(c("A", "B", "C"), 10, replace = TRUE),
#'                   simplify = FALSE)
#' tree <- context_tree(seqs, max_depth = 2L)
#' plot_distributions(tree, top = 9)
#' }
#' @export
plot_distributions <- function(tree, contexts = NULL, top = 12L,
                               min_count = 1L) {
  stopifnot(inherits(tree, "transitiontrees"))
  alpha <- tree$alphabet
  pw <- tree_pathways(tree, min_count = min_count, sort_by = "count")
  pw <- if (is.null(contexts))
    utils::head(pw, as.integer(top)) else
    pw[pw$pathway %in% as.character(contexts), , drop = FALSE]
  if (nrow(pw) == 0L)
    stop("No contexts to plot.", call. = FALSE)

  long <- do.call(rbind, lapply(pw$pathway, function(p) {
    ctx  <- if (identical(p, .ROOT_LABEL)) .ROOT else p
    info <- tree$nodes[[ctx]]
    data.frame(pathway  = p,
               state    = factor(alpha, levels = alpha),
               prob     = info$prob,
               is_modal = seq_along(alpha) == which.max(info$prob),
               stringsAsFactors = FALSE)
  }))
  long$pathway <- factor(long$pathway, levels = pw$pathway)

  ggplot2::ggplot(long, ggplot2::aes(x = .data$state, y = .data$prob,
                                     fill = .data$is_modal)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::facet_wrap(~ pathway) +
    ggplot2::scale_fill_manual(
      values = c(`TRUE` = "#08519c", `FALSE` = "#9ecae1"),
      guide = "none") +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "next state", y = "P(next | context)",
      title = "Per-context next-state distributions",
      subtitle = sprintf("%d contexts; modal next state highlighted",
                         nrow(pw))) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40", size = 9),
      axis.text.x   = ggplot2::element_text(angle = 45, hjust = 1),
      strip.text    = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank())
}

#' Predictive Diagnostics for Held-Out Scoring
#'
#' @description
#' Visual diagnostics for how a fitted tree scores held-out sequences,
#' built on \code{\link{score_positions}()}. The analogues of
#' \code{PST::ppplot()} / \code{pqplot()}:
#' \describe{
#'   \item{\code{type = "position"}}{predicted probability of the
#'     observed next state against position in the sequence — where the
#'     model is confident vs. surprised as a sequence unfolds.}
#'   \item{\code{type = "ecdf"}}{the empirical cumulative distribution of
#'     those predicted probabilities — a calibration-style view of how
#'     often the model assigns high vs. low probability to what actually
#'     happened.}
#' }
#'
#' @param tree A \code{transitiontrees}.
#' @param newdata Held-out sequence data in any format accepted by
#'   \code{context_tree()}.
#' @param type One of \code{"position"} (default) or \code{"ecdf"}.
#'
#' @return A ggplot object.
#'
#' @examples
#' \donttest{
#' fit  <- replicate(60, sample(c("A", "B", "C"), 10, replace = TRUE),
#'                   simplify = FALSE)
#' tree <- context_tree(fit, max_depth = 2L)
#' new  <- replicate(15, sample(c("A", "B", "C"), 10, replace = TRUE),
#'                   simplify = FALSE)
#' plot_predictive(tree, new, type = "position")
#' plot_predictive(tree, new, type = "ecdf")
#' }
#' @export
plot_predictive <- function(tree, newdata, type = c("position", "ecdf")) {
  stopifnot(inherits(tree, "transitiontrees"))
  type <- match.arg(type)
  sp <- score_positions(tree, newdata)
  if (nrow(sp) == 0L)
    stop("No held-out positions to score.", call. = FALSE)

  base_theme <- ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40", size = 9),
      panel.grid.minor = ggplot2::element_blank())

  if (type == "position") {
    sp$sequence_id <- factor(sp$sequence_id)
    return(
      ggplot2::ggplot(sp, ggplot2::aes(x = .data$position,
                                       y = .data$predicted_prob)) +
        ggplot2::geom_line(ggplot2::aes(group = .data$sequence_id),
                           colour = "grey70", alpha = 0.5) +
        ggplot2::geom_point(ggplot2::aes(colour = .data$predicted_prob),
                            size = 2) +
        ggplot2::scale_colour_gradient(low = "#D55E00", high = "#0072B2",
                                       limits = c(0, 1), name = "P(obs)") +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(
          x = "position in sequence",
          y = "predicted probability of observed move",
          title = "Predictive confidence by position",
          subtitle = "low (orange) = the model was surprised here") +
        base_theme)
  }

  ggplot2::ggplot(sp, ggplot2::aes(x = .data$predicted_prob)) +
    ggplot2::stat_ecdf(geom = "step", colour = "#0072B2", linewidth = 1) +
    ggplot2::xlim(0, 1) +
    ggplot2::labs(
      x = "predicted probability of observed move",
      y = "cumulative fraction of positions",
      title = "Predictive probability ECDF",
      subtitle = sprintf("%d held-out positions scored", nrow(sp))) +
    base_theme
}

#' Difference (Subtraction) Map Between Two Groups
#'
#' @description
#' A per-context map of how two groups differ in their next-state
#' predictions: one row per shared context, one column per next state.
#' By default each cell is the Pearson standardized residual of the first
#' group against the no-difference null (red = more than expected, blue =
#' less; \code{|r| > 2} notable), which is support-aware and decomposes
#' the per-context \eqn{G^2}; \code{measure = "probability"} shows the raw
#' \code{P(group1) - P(group2)} instead.
#'
#' @param group A \code{transitiontrees_group}.
#' @param groups Optional length-2 character vector naming the two groups
#'   to subtract (\code{group1 - group2}). Defaults to the first two; it
#'   is required when the object has more than two groups.
#' @param depth Integer or \code{NULL}. Restrict to contexts of this
#'   depth (\code{depth = 1} gives the order-1 transition-matrix
#'   difference). \code{NULL} (default) uses all shared contexts.
#' @param min_count Integer. Drop contexts whose count in \emph{either}
#'   group is below this. Default 1.
#' @param comparison Optional \code{transitiontrees_group_comparison} (from
#'   \code{\link{compare_groups}()}). When supplied, contexts whose
#'   behavioral difference is significant (\code{jsd_padj < alpha}) are
#'   starred, so the map shows which differences survived the
#'   permutation test.
#' @param alpha Numeric. FDR threshold for the significance stars when
#'   \code{comparison} is given. Default 0.05.
#' @param annotate Logical. Print the signed difference in each cell
#'   (\code{layout = "tile"} only). Default \code{TRUE}.
#' @param layout One of \code{"tile"} (default; per-context heatmap over
#'   all shared contexts) or \code{"tree"} (the horizontal context-tree
#'   phylogram drawn on a pooled backbone, with each node and branch
#'   coloured by which group reaches that context more — red for the first
#'   group, blue for the second; node size = pooled count).
#' @param measure For \code{layout = "tile"}, what each cell encodes:
#'   \code{"residual"} (default) is the Pearson standardized residual of
#'   the first group against the no-group null (observed vs expected from
#'   the context's margins) — support-aware and decomposing the per-context
#'   \eqn{G^2}; \code{"probability"} is the raw next-state probability
#'   difference \code{P(group1) - P(group2)}. Ignored by \code{"tree"}.
#'
#' @return A ggplot object.
#'
#' @examples
#' \donttest{
#' gx <- replicate(60, sample(c("A","B","C"), 8, replace = TRUE,
#'                            prob = c(.2,.6,.2)), simplify = FALSE)
#' gy <- replicate(60, sample(c("A","B","C"), 8, replace = TRUE,
#'                            prob = c(.2,.2,.6)), simplify = FALSE)
#' grp <- context_tree(c(gx, gy), group = rep(c("x","y"), each = 60),
#'                     max_depth = 2L)
#' plot_difference(grp)                       # residual heatmap
#' plot_difference(grp, layout = "tree")      # difference on the tree map
#' }
#' @seealso \code{\link{compare_groups}} for the significance test.
#' @export
plot_difference <- function(group, groups = NULL, depth = NULL,
                            min_count = 1L, comparison = NULL,
                            alpha = 0.05, annotate = TRUE,
                            layout = c("tile", "tree"),
                            measure = c("residual", "probability")) {
  stopifnot(inherits(group, "transitiontrees_group"))
  layout  <- match.arg(layout)
  measure <- match.arg(measure)
  gn <- names(group)
  if (is.null(groups)) {
    if (length(gn) != 2L)
      stop("The object has ", length(gn), " groups; name the two to ",
           "subtract via 'groups = c(\"a\", \"b\")'.", call. = FALSE)
    groups <- gn[1:2]
  }
  if (length(groups) != 2L || !all(groups %in% gn))
    stop("'groups' must name two groups of the object.", call. = FALSE)
  a <- group[[groups[1L]]]; b <- group[[groups[2L]]]
  alphabet <- a$alphabet

  if (layout == "tree")
    return(.diff_tree(a, b, groups, alphabet))

  ctx <- intersect(names(a$nodes), names(b$nodes))
  if (!is.null(depth))
    ctx <- ctx[vapply(ctx, .cg_depth, integer(1)) == as.integer(depth)]
  ctx <- ctx[vapply(ctx, function(cc)
    a$nodes[[cc]]$n >= min_count && b$nodes[[cc]]$n >= min_count,
    logical(1))]
  if (length(ctx) == 0L)
    stop("No shared contexts to compare (try a lower 'min_count' or ",
         "another 'depth').", call. = FALSE)

  ## significance stars from the permutation comparison, if supplied
  sig_label <- function(p) p
  if (!is.null(comparison)) {
    stopifnot(inherits(comparison, "transitiontrees_group_comparison"))
    cp <- comparison$pathways
    sig <- cp$pathway[cp$jsd_padj < alpha]
    sig_label <- function(p) ifelse(p %in% sig, paste0(p, " *"), p)
  }

  pretty <- function(cc) {
    lab <- if (identical(cc, .ROOT)) .ROOT_LABEL else cc
    sig_label(lab)
  }
  ## per-cell value: Pearson residual of group 1 against the no-group
  ## null (observed vs expected from the context's margins), or the raw
  ## probability difference.
  cell_vals <- function(cc) {
    if (measure == "probability")
      return(a$nodes[[cc]]$prob - b$nodes[[cc]]$prob)
    o <- rbind(a$nodes[[cc]]$counts, b$nodes[[cc]]$counts)
    e <- .cg_expected(o)
    r <- (o[1L, ] - e[1L, ]) / sqrt(e[1L, ])
    r[!is.finite(r)] <- 0
    r
  }
  long <- do.call(rbind, lapply(ctx, function(cc) {
    v <- cell_vals(cc)
    data.frame(pathway    = pretty(cc),
               next_state = factor(alphabet, levels = alphabet),
               val        = v, abs_tot = sum(abs(v)),
               stringsAsFactors = FALSE)
  }))
  ## most-different contexts at the top
  row_order <- unique(long$pathway[order(long$abs_tot)])
  long$pathway <- factor(long$pathway, levels = row_order)

  fmt   <- if (measure == "probability") "%+.2f" else "%+.1f"
  thresh <- if (measure == "probability") 0.01 else 0.5
  fill_name <- if (measure == "probability")
    sprintf("P(%s) - P(%s)", groups[1L], groups[2L]) else
    sprintf("Pearson\nresidual\n(%s)", groups[1L])
  M <- max(abs(long$val)); if (M == 0) M <- 1

  p <- ggplot2::ggplot(long, ggplot2::aes(.data$next_state, .data$pathway,
                                          fill = .data$val)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4)
  if (isTRUE(annotate))
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = ifelse(abs(.data$val) >= thresh,
                                  sprintf(fmt, .data$val), "")),
      size = 3)
  sub <- if (measure == "probability")
    sprintf("red = more likely in %s, blue = more likely in %s",
            groups[1L], groups[2L]) else
    sprintf(paste0("observed vs expected (no-group null): red = %s does ",
                   "the move MORE than expected, blue = less; |r| > 2 notable"),
            groups[1L])
  p +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      limits = c(-M, M), name = fill_name) +
    ggplot2::labs(
      x = "next state", y = NULL,
      title = sprintf("Difference map: %s vs %s", groups[1L], groups[2L]),
      subtitle = sprintf("%s%s", sub,
        if (is.null(comparison)) "" else
          sprintf("  |  * = FDR < %.2f", alpha))) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40", size = 9),
      axis.text.y   = ggplot2::element_text(family = "mono"),
      panel.grid    = ggplot2::element_blank())
}

#' Difference rendered on the horizontal context-tree layout: a pooled
#' backbone whose nodes and edges are coloured by which group reaches a
#' context more (prevalence early - late).
#' @noRd
.diff_tree <- function(a, b, groups, alphabet) {
  trajs_a <- .ct_traj(a$data); trajs_b <- .ct_traj(b$data)
  pooled <- context_tree(c(trajs_a, trajs_b),
                         max_depth = a$max_depth, min_count = a$nmin)
  if (length(pooled$nodes) == 0L)
    stop("Pooled tree is empty; nothing to plot.", call. = FALSE)

  contexts <- names(pooled$nodes)
  depths   <- vapply(contexts, .cg_depth, integer(1))
  cc <- .cg_context_counts(stats::setNames(list(trajs_a, trajs_b), groups),
                           groups, contexts, depths, alphabet)
  Ng <- vapply(groups, function(g) sum(cc[[.ROOT]][g, ]), numeric(1))
  delta <- vapply(contexts, function(ctx)
    sum(cc[[ctx]][groups[1L], ]) / Ng[1L] -
    sum(cc[[ctx]][groups[2L], ]) / Ng[2L], numeric(1))
  names(delta) <- contexts

  layout  <- .ct_horizontal_layout(pooled)
  e_paths <- .ct_horizontal_edge_paths(layout, pooled$edges)
  layout$delta  <- delta[layout$context]
  e_paths$delta <- delta[pooled$edges$child][e_paths$edge_id]

  is_root <- layout$context == .ROOT
  body_df <- layout[!is_root, , drop = FALSE]
  root_df <- layout[ is_root, , drop = FALSE]
  is_leaf <- !(layout$context %in% unique(pooled$edges$parent))
  leaf_df <- layout[is_leaf & !is_root, , drop = FALSE]
  leaf_df$label <- leaf_df$context
  internal_df <- layout[!is_leaf & !is_root, , drop = FALSE]
  internal_df$abbr <- vapply(internal_df$context, .pt_last_state,
                             character(1))

  M <- max(abs(layout$delta), na.rm = TRUE); if (!is.finite(M) || M == 0) M <- 1
  div_fill <- ggplot2::scale_fill_gradient2(
    low = "#2166AC", mid = "grey92", high = "#B2182B", midpoint = 0,
    limits = c(-M, M),
    name = sprintf("prevalence\n%s - %s", groups[1L], groups[2L]))
  div_col <- ggplot2::scale_colour_gradient2(
    low = "#2166AC", mid = "grey75", high = "#B2182B", midpoint = 0,
    limits = c(-M, M), guide = "none")

  ggplot2::ggplot() +
    ggplot2::geom_path(
      data = e_paths,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$edge_id,
                   linewidth = .data$edge_weight, colour = .data$delta),
      lineend = "round", linejoin = "round") +
    ggplot2::geom_point(
      data = body_df,
      ggplot2::aes(x = .data$x, y = .data$y, size = .data$n,
                   fill = .data$delta),
      shape = 21, colour = "grey25", stroke = 0.2, na.rm = TRUE) +
    ggplot2::geom_text(
      data = internal_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$abbr),
      vjust = 1, nudge_y = -0.05, size = 2.7, colour = "grey20",
      fontface = "bold", na.rm = TRUE) +
    ggplot2::geom_text(
      data = leaf_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      hjust = 0, nudge_x = 0.10, size = 3, colour = "grey20",
      family = "mono") +
    ggplot2::geom_point(
      data = root_df, ggplot2::aes(x = .data$x, y = .data$y),
      shape = 16, size = 12, colour = "black") +
    ggplot2::geom_text(
      data = root_df, ggplot2::aes(x = .data$x, y = .data$y),
      label = .ROOT_LABEL, hjust = 1, nudge_x = -0.14, size = 3,
      colour = "grey25", fontface = "bold") +
    div_fill + div_col +
    ggplot2::scale_size_continuous(range = c(4, 14), name = "context count") +
    ggplot2::scale_linewidth_continuous(range = c(0.3, 2.5),
                                        name = "edge flow") +
    ggplot2::scale_x_continuous(
      breaks = seq.int(0L, pooled$max_depth, by = 1L),
      expand = ggplot2::expansion(mult = c(0.17, 0.32))) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40", size = 9),
      legend.position = "right",
      axis.title.y = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid   = ggplot2::element_blank()) +
    ggplot2::labs(
      x = "depth",
      title = sprintf("Difference tree: %s vs %s", groups[1L], groups[2L]),
      subtitle = sprintf(
        "node & edge colour = which phase reaches the context more (red = %s, blue = %s); size = pooled count",
        groups[1L], groups[2L]))
}
