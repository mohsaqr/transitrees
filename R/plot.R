# ---- plot.pathtree() dispatcher + dendrogram + horizontal phylogram ----

#' @noRd
.pt_state_palette <- function(state_levels) {
  ## Base-R qualitative palettes from grDevices::hcl.colors().
  ## "Set 2" is the colour-blind-safe 8-class default; falls through
  ## to higher-cardinality palettes for richer alphabets.
  n <- length(state_levels)
  pal_name <- if (n <= 8L) "Set 2"
              else if (n <= 12L) "Set 3"
              else "Dark 3"
  pal <- grDevices::hcl.colors(max(2L, n), palette = pal_name)
  setNames(pal[seq_len(n)], state_levels)
}

#' @noRd
.ct_radial_layout <- function(tree) {
  ## Radial coordinates for each node:
  ##   r     = depth
  ##   theta = leaf-rank for leaves (evenly spaced 0..2pi),
  ##           mean of children's theta for internal nodes.
  ## Result includes (x, y) = (r * cos theta, r * sin theta).
  if (length(tree$nodes) == 0L) return(NULL)

  children_by_parent <- .pt_children_of(tree)
  children_of <- function(parent)
    children_by_parent[[parent]] %||% character(0)

  ## DFS to collect leaves left-to-right.
  leaves <- character(0)
  walk <- function(ctx) {
    ch <- children_of(ctx)
    if (length(ch) == 0L) {
      leaves[[length(leaves) + 1L]] <<- ctx
    } else {
      lapply(ch, walk)
    }
    invisible(NULL)
  }
  walk(.ROOT)

  theta_env <- new.env(hash = TRUE, parent = emptyenv())

  if (length(leaves) <= 1L) {
    ## Degenerate trees (root only, or root + single leaf) get a fixed angle.
    Map(function(nm) assign(nm, 0, envir = theta_env), names(tree$nodes))
  } else {
    leaf_theta <- 2 * pi * (seq_along(leaves) - 0.5) / length(leaves)
    Map(function(nm, t) assign(nm, t, envir = theta_env),
        leaves, leaf_theta)

    assign_theta <- function(ctx) {
      ch <- children_of(ctx)
      if (length(ch) == 0L) return(theta_env[[ctx]])
      child_thetas <- vapply(ch, assign_theta, numeric(1))
      t <- mean(child_thetas)
      assign(ctx, t, envir = theta_env)
      t
    }
    assign_theta(.ROOT)
  }

  do.call(rbind, lapply(names(tree$nodes), function(ctx) {
    d  <- tree$nodes[[ctx]]$depth
    th <- theta_env[[ctx]]
    data.frame(context = ctx, depth = d, r = d, theta = th,
               x = d * cos(th), y = d * sin(th),
               n = tree$nodes[[ctx]]$n,
               stringsAsFactors = FALSE)
  }))
}

#' @noRd
.ct_radial_edge_paths <- function(layout, edges, n_arc = 24L) {
  ## Phylogram-style elbow edges with an edge-weight column equal to
  ## the child's count, so callers can map it to linewidth.
  if (nrow(edges) == 0L)
    return(data.frame(edge_id = integer(0), order = integer(0),
                      x = numeric(0), y = numeric(0),
                      edge_weight = numeric(0)))

  by_ctx <- setNames(seq_len(nrow(layout)), layout$context)
  paths  <- Map(function(p_ctx, c_ctx, eid) {
    rp <- layout$r[[by_ctx[[p_ctx]]]];  tp <- layout$theta[[by_ctx[[p_ctx]]]]
    rc <- layout$r[[by_ctx[[c_ctx]]]];  tc <- layout$theta[[by_ctx[[c_ctx]]]]
    nc <- layout$n[[by_ctx[[c_ctx]]]]
    arc_t <- seq.int(tp, tc, length.out = n_arc)
    data.frame(
      edge_id     = eid,
      order       = seq_len(n_arc + 1L),
      x           = c(rp * cos(arc_t), rc * cos(tc)),
      y           = c(rp * sin(arc_t), rc * sin(tc)),
      edge_weight = nc
    )
  }, edges$parent, edges$child, seq_len(nrow(edges)))
  do.call(rbind, paths)
}

#' @noRd
.ct_horizontal_layout <- function(tree) {
  ## Cartesian coordinates for a left-to-right phylogram:
  ##   x = depth (root at x=0, leaves at x=max_depth)
  ##   y = leaf-rank for leaves (evenly spaced in [0, 1]),
  ##       mean of children's y for internal nodes.
  if (length(tree$nodes) == 0L) return(NULL)

  children_by_parent <- .pt_children_of(tree)
  children_of <- function(parent)
    children_by_parent[[parent]] %||% character(0)

  ## DFS to collect leaves top-to-bottom.
  leaves <- character(0)
  walk <- function(ctx) {
    ch <- children_of(ctx)
    if (length(ch) == 0L) {
      leaves[[length(leaves) + 1L]] <<- ctx
    } else {
      lapply(ch, walk)
    }
    invisible(NULL)
  }
  walk(.ROOT)

  y_env <- new.env(hash = TRUE, parent = emptyenv())

  if (length(leaves) <= 1L) {
    Map(function(nm) assign(nm, 0.5, envir = y_env), names(tree$nodes))
  } else {
    leaf_y <- seq.int(0, 1, length.out = length(leaves))
    ## Reverse so the first leaf is at the top of the plot canvas.
    leaf_y <- rev(leaf_y)
    Map(function(nm, t) assign(nm, t, envir = y_env), leaves, leaf_y)

    assign_y <- function(ctx) {
      ch <- children_of(ctx)
      if (length(ch) == 0L) return(y_env[[ctx]])
      child_ys <- vapply(ch, assign_y, numeric(1))
      t <- mean(child_ys)
      assign(ctx, t, envir = y_env)
      t
    }
    assign_y(.ROOT)
  }

  do.call(rbind, lapply(names(tree$nodes), function(ctx) {
    d <- tree$nodes[[ctx]]$depth
    y <- y_env[[ctx]]
    data.frame(context = ctx, depth = d, x = d, y = y,
               n = tree$nodes[[ctx]]$n,
               stringsAsFactors = FALSE)
  }))
}

#' @noRd
.ct_horizontal_edge_paths <- function(layout, edges) {
  ## Right-angle elbow edges:
  ##   1. vertical from (x_parent, y_parent) -> (x_parent, y_child)
  ##   2. horizontal from (x_parent, y_child) -> (x_child, y_child)
  ## Each edge is a 3-vertex polyline; edge_weight = child's count.
  if (nrow(edges) == 0L)
    return(data.frame(edge_id = integer(0), order = integer(0),
                      x = numeric(0), y = numeric(0),
                      edge_weight = numeric(0)))

  by_ctx <- setNames(seq_len(nrow(layout)), layout$context)
  paths <- Map(function(p_ctx, c_ctx, eid) {
    xp <- layout$x[[by_ctx[[p_ctx]]]]; yp <- layout$y[[by_ctx[[p_ctx]]]]
    xc <- layout$x[[by_ctx[[c_ctx]]]]; yc <- layout$y[[by_ctx[[c_ctx]]]]
    nc <- layout$n[[by_ctx[[c_ctx]]]]
    data.frame(
      edge_id     = eid,
      order       = 1:3,
      x           = c(xp, xp, xc),
      y           = c(yp, yc, yc),
      edge_weight = nc
    )
  }, edges$parent, edges$child, seq_len(nrow(edges)))
  do.call(rbind, paths)
}

#' @noRd
.plot_dendrogram <- function(x,
                              point_size_range = c(5, 16),
                              edge_size_range  = c(0.3, 2.5),
                              ...) {
  if (length(x$nodes) == 0L)
    stop("Cannot plot an empty tree.", call. = FALSE)

  layout  <- .ct_radial_layout(x)
  e_paths <- .ct_radial_edge_paths(layout, x$edges)

  state_levels   <- x$alphabet
  last_state_str <- vapply(layout$context, .pt_last_state, character(1))
  layout$state   <- factor(last_state_str, levels = state_levels)
  layout$abbr    <- ifelse(is.na(layout$state), NA_character_,
                            substr(as.character(layout$state), 1L, 3L))
  pal <- .pt_state_palette(state_levels)

  ## Root is rendered as a separate layer (last, on top, fixed size,
  ## fixed black fill) so it never gets buried under depth-1 children
  ## drawn on top of the (0, 0) origin.
  is_root  <- layout$context == .ROOT
  body_df  <- layout[!is_root, , drop = FALSE]
  root_df  <- layout[ is_root, , drop = FALSE]

  ggplot2::ggplot() +
    ggplot2::geom_path(
      data = e_paths,
      ggplot2::aes(x = .data$x, y = .data$y,
                   group = .data$edge_id,
                   linewidth = .data$edge_weight),
      colour = "grey60", lineend = "round"
    ) +
    ggplot2::geom_point(
      data = body_df,
      ggplot2::aes(x = .data$x, y = .data$y,
                   size = .data$n, fill = .data$state),
      shape = 21, colour = "grey25", stroke = 0.15, na.rm = TRUE
    ) +
    ggplot2::geom_text(
      data = body_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$abbr),
      size = 2.4, colour = "grey15", fontface = "bold", na.rm = TRUE
    ) +
    ggplot2::geom_point(
      data = root_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      shape = 16, size = max(point_size_range) * 0.8,
      colour = "black"
    ) +
    ggplot2::scale_size_continuous(range = point_size_range,
                                   name = "context count") +
    ggplot2::scale_linewidth_continuous(range = edge_size_range,
                                         name = "edge flow") +
    ggplot2::scale_fill_manual(values = pal, name = "State",
                               na.translate = FALSE,
                               drop = FALSE) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40", size = 9),
      legend.position = "right"
    ) +
    ggplot2::labs(
      title    = "Context tree (radial)",
      subtitle = sprintf(
        "%d nodes, depth <= %d;  node size = count;  edge thickness = flow",
        length(x$nodes), x$max_depth)
    )
}

#' @noRd
.plot_horizontal <- function(x,
                              point_size_range = c(4, 14),
                              edge_size_range  = c(0.3, 2.5),
                              ...) {
  if (length(x$nodes) == 0L)
    stop("Cannot plot an empty tree.", call. = FALSE)

  layout  <- .ct_horizontal_layout(x)
  e_paths <- .ct_horizontal_edge_paths(layout, x$edges)

  state_levels   <- x$alphabet
  last_state_str <- vapply(layout$context, .pt_last_state, character(1))
  layout$state   <- factor(last_state_str, levels = state_levels)
  layout$abbr    <- ifelse(is.na(layout$state), NA_character_,
                            substr(as.character(layout$state), 1L, 3L))
  pal <- .pt_state_palette(state_levels)

  is_root  <- layout$context == .ROOT
  body_df  <- layout[!is_root, , drop = FALSE]
  root_df  <- layout[ is_root, , drop = FALSE]
  ## Leaf labels: the full pathway in arrow notation, drawn to the
  ## right of each leaf node.
  is_leaf  <- !(layout$context %in% unique(x$edges$parent))
  leaf_df  <- layout[is_leaf & !is_root, , drop = FALSE]
  leaf_df$label <- ifelse(leaf_df$context == .ROOT, "(root)",
                          leaf_df$context)

  ggplot2::ggplot() +
    ggplot2::geom_path(
      data = e_paths,
      ggplot2::aes(x = .data$x, y = .data$y,
                   group = .data$edge_id,
                   linewidth = .data$edge_weight),
      colour = "grey60", lineend = "round"
    ) +
    ggplot2::geom_point(
      data = body_df,
      ggplot2::aes(x = .data$x, y = .data$y,
                   size = .data$n, fill = .data$state),
      shape = 21, colour = "grey25", stroke = 0.2, na.rm = TRUE
    ) +
    ggplot2::geom_text(
      data = body_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$abbr),
      size = 2.6, colour = "grey15", fontface = "bold", na.rm = TRUE
    ) +
    ggplot2::geom_text(
      data = leaf_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      hjust = -0.15, size = 2.8, colour = "grey25",
      family = "mono"
    ) +
    ggplot2::geom_point(
      data = root_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      shape = 16, size = max(point_size_range) * 0.85,
      colour = "black"
    ) +
    ggplot2::geom_text(
      data = root_df,
      ggplot2::aes(x = .data$x, y = .data$y), label = "(root)",
      hjust = 1.25, size = 3, colour = "grey25", fontface = "bold"
    ) +
    ggplot2::scale_size_continuous(range = point_size_range,
                                   name = "context count") +
    ggplot2::scale_linewidth_continuous(range = edge_size_range,
                                         name = "edge flow") +
    ggplot2::scale_fill_manual(values = pal, name = "State",
                               na.translate = FALSE,
                               drop = FALSE) +
    ggplot2::scale_x_continuous(
      breaks = seq.int(0L, x$max_depth, by = 1L),
      expand = ggplot2::expansion(mult = c(0.06, 0.22))
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold"),
      plot.subtitle    = ggplot2::element_text(colour = "grey40", size = 9),
      legend.position  = "right",
      axis.title.y     = ggplot2::element_blank(),
      axis.text.y      = ggplot2::element_blank(),
      axis.ticks.y     = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      x        = "depth",
      title    = "Context tree (horizontal)",
      subtitle = sprintf(
        "%d nodes, depth <= %d;  node size = count;  edge thickness = flow",
        length(x$nodes), x$max_depth)
    )
}

#' Plot a Context Tree
#'
#' @description
#' Renders a fitted pathtree in one of four styles:
#' \itemize{
#'   \item \code{"dendrogram"} (default) — pure-ggplot2 radial tree:
#'     root at the centre, leaves on the outer ring.
#'   \item \code{"horizontal"} — pure-ggplot2 left-to-right
#'     phylogram: root on the left, leaves fanned out vertically on
#'     the right with full arrow-notation leaf labels.
#'   \item \code{"icicle"} — circular partition / sunburst via
#'     \code{ggraph} (Suggests); inner ring carries full state
#'     names, outer rings carry 3-letter abbreviations.
#'   \item \code{"interactive"} — \code{collapsibleTree} htmlwidget
#'     (Suggests). Click a node to collapse / expand its subtree;
#'     hover for a tooltip with the full pathway, count, modal next
#'     state, and the complete next-state distribution.
#' }
#'
#' In the three static styles: node fill = last state of the
#' pathway, node size = context count, edge thickness = child's
#' count (\dQuote{flow}). In the interactive style, fill encodes the
#' last state and the subtree depth is the visual cue for count.
#'
#' @param x A \code{pathtree}.
#' @param style One of \code{"dendrogram"} (default),
#'   \code{"horizontal"}, \code{"icicle"}, or \code{"interactive"}.
#' @param point_size_range Numeric length-2 vector controlling the
#'   minimum and maximum node-point size. Default \code{c(5, 16)} for
#'   \code{"dendrogram"}, \code{c(4, 14)} for \code{"horizontal"}.
#'   Ignored by \code{"icicle"} and \code{"interactive"}.
#' @param edge_size_range Numeric length-2 vector for edge linewidth.
#'   Default \code{c(0.3, 2.5)}. Ignored by \code{"icicle"} and
#'   \code{"interactive"}.
#' @param ... Passed to the chosen backend (e.g. \code{width} /
#'   \code{height} for \code{"interactive"}).
#'
#' @return A ggplot object for the three static styles; an
#'   \code{htmlwidget} for \code{"interactive"}.
#'
#' @details
#' \code{"icicle"} requires \code{ggraph} + \code{tidygraph};
#' \code{"interactive"} requires \code{collapsibleTree}. The
#' dispatcher errors informatively if a needed Suggests dependency
#' is missing.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' m <- matrix(sample(c("A","B","C"), 200, TRUE), 20)
#' tr <- context_tree(m, max_depth = 2L, nmin = 3L)
#' plot(tr)                           # radial dendrogram (default)
#' plot(tr, style = "horizontal")     # left-to-right phylogram
#' plot(tr, style = "icicle")         # sunburst (needs Suggests)
#' plot(tr, style = "interactive")    # collapsibleTree (needs Suggests)
#' }
#' @export
plot.pathtree <- function(x,
                          style = c("dendrogram", "horizontal",
                                     "icicle", "interactive"),
                          point_size_range = NULL,
                          edge_size_range  = c(0.3, 2.5),
                          ...) {
  style <- match.arg(style)
  if (style == "icicle" &&
      (!requireNamespace("ggraph", quietly = TRUE) ||
       !requireNamespace("tidygraph", quietly = TRUE))) {
    stop("plot(style = 'icicle') requires 'ggraph' and 'tidygraph'. ",
         "Install both, or use style = 'dendrogram' / 'horizontal'.",
         call. = FALSE)
  }
  if (style == "interactive" &&
      !requireNamespace("collapsibleTree", quietly = TRUE)) {
    stop("plot(style = 'interactive') requires 'collapsibleTree'. ",
         "Install it, or use style = 'dendrogram' / 'horizontal' / ",
         "'icicle'.", call. = FALSE)
  }
  switch(style,
    icicle      = .plot_icicle(x, ...),
    interactive = .plot_interactive(x, ...),
    dendrogram  = .plot_dendrogram(x,
                   point_size_range = point_size_range %||% c(5, 16),
                   edge_size_range  = edge_size_range,
                   ...),
    horizontal  = .plot_horizontal(x,
                   point_size_range = point_size_range %||% c(4, 14),
                   edge_size_range  = edge_size_range,
                   ...))
}
