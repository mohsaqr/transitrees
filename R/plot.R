# ---- plot.transitiontrees() dispatcher + dendrogram + horizontal phylogram ----

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
.ct_horizontal_edge_paths <- function(layout, edges, n_pt = 40L) {
  ## Smooth S-curve edges: a cosine smoothstep carries each branch from
  ## (x_parent, y_parent) to (x_child, y_child), leaving the parent and
  ## arriving at the child horizontally, so branches read as flowing
  ## curves rather than mechanical right-angle elbows. Each edge is an
  ## n_pt-vertex polyline; edge_weight = child's count.
  if (nrow(edges) == 0L)
    return(data.frame(edge_id = integer(0), order = integer(0),
                      x = numeric(0), y = numeric(0),
                      edge_weight = numeric(0)))

  by_ctx <- setNames(seq_len(nrow(layout)), layout$context)
  t  <- seq(0, 1, length.out = n_pt)
  sm <- (1 - cos(pi * t)) / 2          # 0 -> 1 smoothstep
  paths <- Map(function(p_ctx, c_ctx, eid) {
    xp <- layout$x[[by_ctx[[p_ctx]]]]; yp <- layout$y[[by_ctx[[p_ctx]]]]
    xc <- layout$x[[by_ctx[[c_ctx]]]]; yc <- layout$y[[by_ctx[[c_ctx]]]]
    nc <- layout$n[[by_ctx[[c_ctx]]]]
    data.frame(
      edge_id     = eid,
      order       = seq_len(n_pt),
      x           = xp + t  * (xc - xp),
      y           = yp + sm * (yc - yp),
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
  layout$abbr    <- as.character(layout$state)
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
                              show_prediction  = TRUE,
                              ...) {
  if (length(x$nodes) == 0L)
    stop("Cannot plot an empty tree.", call. = FALSE)

  layout  <- .ct_horizontal_layout(x)
  e_paths <- .ct_horizontal_edge_paths(layout, x$edges)

  ## Node colour = the MOST-RECENT move in the context (the rightmost
  ## token), i.e. "what just happened". This is the intuitive reading and
  ## makes each branch off a depth-1 hub share one colour (every "X ->
  ## Specify" ends in Specify, so the whole Specify branch is one colour).
  state_levels     <- x$alphabet
  recent_state_str <- vapply(layout$context, function(ctx) {
    if (identical(ctx, .ROOT)) return(NA_character_)
    parts <- strsplit(ctx, " -> ", fixed = TRUE)[[1L]]
    parts[[length(parts)]]
  }, character(1))
  layout$state <- factor(recent_state_str, levels = state_levels)
  pal <- .pt_state_palette(state_levels)

  is_root  <- layout$context == .ROOT
  body_df  <- layout[!is_root, , drop = FALSE]
  root_df  <- layout[ is_root, , drop = FALSE]
  ## Labels: full pathway in arrow notation drawn UNDER every non-root
  ## node. Simple by default (context only); set show_prediction = TRUE to
  ## add the modal prediction "(state pct%)" on a second line.
  body_df$label <- if (isTRUE(show_prediction)) {
    vapply(body_df$context, function(ctx) {
      pr   <- x$nodes[[ctx]]$prob
      best <- which.max(pr)
      sprintf("%s\n(%s %.0f%%)", ctx, x$alphabet[best], 100 * pr[best])
    }, character(1))
  } else {
    body_df$context
  }

  ggplot2::ggplot() +
    ggplot2::geom_path(
      data = e_paths,
      ggplot2::aes(x = .data$x, y = .data$y,
                   group = .data$edge_id,
                   linewidth = .data$edge_weight),
      colour = "grey60", lineend = "round", linejoin = "round"
    ) +
    ggplot2::geom_point(
      data = body_df,
      ggplot2::aes(x = .data$x, y = .data$y,
                   size = .data$n, fill = .data$state),
      shape = 21, colour = "grey25", stroke = 0.2, na.rm = TRUE
    ) +
    ggplot2::geom_text(
      data = body_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      hjust = 0.5, vjust = 1, nudge_y = -0.035, size = 2.5,
      lineheight = 0.9, colour = "grey20", family = "mono", na.rm = TRUE
    ) +
    ggplot2::geom_point(
      data = root_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      shape = 16, size = max(point_size_range) * 0.85,
      colour = "black"
    ) +
    ggplot2::geom_text(
      data = root_df,
      ggplot2::aes(x = .data$x, y = .data$y), label = .ROOT_LABEL,
      hjust = 1, nudge_x = -0.14, size = 3,
      colour = "grey25", fontface = "bold"
    ) +
    ggplot2::scale_size_continuous(range = point_size_range,
                                   name = "context count") +
    ggplot2::scale_linewidth_continuous(range = edge_size_range,
                                         name = "edge flow") +
    ggplot2::scale_fill_manual(values = pal, name = "Most recent move",
                               na.translate = FALSE,
                               drop = FALSE) +
    ggplot2::scale_x_continuous(
      breaks = seq.int(0L, x$max_depth, by = 1L),
      expand = ggplot2::expansion(mult = c(0.17, 0.32))
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold"),
      plot.subtitle    = ggplot2::element_text(colour = "grey40", size = 9),
      legend.position  = "right",
      axis.title.y     = ggplot2::element_blank(),
      axis.text.y      = ggplot2::element_blank(),
      axis.ticks.y     = ggplot2::element_blank(),
      panel.grid       = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      x        = "depth",
      title    = "Context tree",
      subtitle = "colour = most recent move; size = count; ( ) = predicted next move"
    )
}

#' Plot a Context Tree
#'
#' @description
#' Renders a fitted transitiontrees in one of four styles:
#' \itemize{
#'   \item \code{"horizontal"} (default) — pure-ggplot2 left-to-right
#'     phylogram: root on the left, contexts fanned out vertically to
#'     the right, each labelled beneath with its full arrow-notation
#'     context and (by default) its modal prediction
#'     \code{"(state pct\%)"} on a second line. Node fill is the
#'     \emph{most recent} move (rightmost token), so each branch off a
#'     depth-1 hub shares a colour. Set \code{show_prediction = FALSE}
#'     for context-only labels.
#'   \item \code{"dendrogram"} — pure-ggplot2 radial tree: root at the
#'     centre, leaves on the outer ring.
#'   \item \code{"icicle"} — circular partition / sunburst via
#'     \code{ggraph} (Suggests); inner ring carries full state
#'     names, outer rings carry 3-letter abbreviations.
#'   \item \code{"interactive"} — \code{visNetwork} htmlwidget
#'     (Suggests). A draggable, zoomable hierarchical tree; node size
#'     = context count and edge width = child's count (\dQuote{flow}),
#'     the same encoding as the static styles. Hover for a tooltip
#'     with the full pathway, count, modal next state, and the
#'     complete next-state distribution.
#' }
#'
#' Common encoding: node size = context count, edge thickness = child's
#' count (\dQuote{flow}). Node fill is the most recent move (rightmost
#' token of the context) in the \code{"horizontal"} style; the
#' \code{"dendrogram"}, \code{"icicle"}, and \code{"interactive"} styles
#' colour by the branching (oldest) token.
#'
#' @param x A \code{transitiontrees}.
#' @param style One of \code{"horizontal"} (default),
#'   \code{"dendrogram"}, \code{"icicle"}, or \code{"interactive"}.
#' @param point_size_range Numeric length-2 vector controlling the
#'   minimum and maximum node-point size. Default \code{c(5, 16)} for
#'   \code{"dendrogram"}, \code{c(4, 14)} for \code{"horizontal"},
#'   \code{c(10, 45)} (pixels) for \code{"interactive"}. Ignored by
#'   \code{"icicle"}.
#' @param edge_size_range Numeric length-2 vector for edge width.
#'   Default \code{c(0.3, 2.5)} for the static styles, \code{c(1, 10)}
#'   (pixels) for \code{"interactive"}. Ignored by \code{"icicle"}.
#' @param ... Passed to the chosen backend. For
#'   \code{style = "horizontal"}, \code{show_prediction} (logical, default
#'   \code{TRUE}) toggles each node's modal prediction
#'   \code{"(state pct\%)"} on a second label line; set \code{FALSE} for
#'   context-only labels. For \code{style = "interactive"}, \code{width} /
#'   \code{height} size the htmlwidget.
#'
#' @return A ggplot object for the three static styles; an
#'   \code{htmlwidget} for \code{"interactive"}.
#'
#' @details
#' \code{"icicle"} requires \code{ggraph} + \code{tidygraph};
#' \code{"interactive"} requires \code{visNetwork}. The dispatcher
#' errors informatively if a needed Suggests dependency is missing.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' m <- matrix(sample(c("A","B","C"), 200, TRUE), 20)
#' tr <- context_tree(m, max_depth = 2L, min_count = 3L)
#' plot(tr)                           # left-to-right phylogram (default)
#' plot(tr, style = "dendrogram")     # radial dendrogram
#' if (requireNamespace("ggraph", quietly = TRUE) &&
#'     requireNamespace("tidygraph", quietly = TRUE))
#'   plot(tr, style = "icicle")       # sunburst (needs ggraph + tidygraph)
#' if (requireNamespace("visNetwork", quietly = TRUE))
#'   plot(tr, style = "interactive")  # visNetwork htmlwidget (Suggests)
#' }
#' @export
plot.transitiontrees <- function(x,
                          style = c("horizontal", "dendrogram",
                                     "icicle", "interactive"),
                          point_size_range = NULL,
                          edge_size_range  = NULL,
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
      !requireNamespace("visNetwork", quietly = TRUE)) {
    stop("plot(style = 'interactive') requires 'visNetwork'. ",
         "Install it, or use style = 'dendrogram' / 'horizontal' / ",
         "'icicle'.", call. = FALSE)
  }
  switch(style,
    icicle      = .plot_icicle(x, ...),
    interactive = .plot_interactive(x,
                   point_size_range = point_size_range %||% c(10, 45),
                   edge_size_range  = edge_size_range  %||% c(1, 10),
                   ...),
    dendrogram  = .plot_dendrogram(x,
                   point_size_range = point_size_range %||% c(5, 16),
                   edge_size_range  = edge_size_range  %||% c(0.3, 2.5),
                   ...),
    horizontal  = .plot_horizontal(x,
                   point_size_range = point_size_range %||% c(4, 14),
                   edge_size_range  = edge_size_range  %||% c(0.3, 2.5),
                   ...))
}

#' Plot Each Tree in a context tree Group
#'
#' @description
#' Draw every member of a \code{transitiontrees_group} in turn via
#' \code{\link{plot.transitiontrees}}. Each member's plot is printed (so the
#' call produces one figure per group, e.g. in an R Markdown chunk),
#' captioned with its group name; the named list of plot objects is
#' returned invisibly for further use.
#'
#' @param x A \code{transitiontrees_group}.
#' @param ... Passed to \code{\link{plot.transitiontrees}} (e.g.
#'   \code{style}).
#'
#' @return Invisibly, a named list (one entry per group, in group order)
#'   of the per-member plot objects.
#'
#' @examples
#' \donttest{
#' m   <- matrix(sample(c("A","B","C"), 200, replace = TRUE), 40, 5)
#' grp <- context_tree(m, group = rep(c("x","y"), each = 20),
#'                     max_depth = 2L)
#' plot(grp)                         # one tree per group
#' }
#' @export
plot.transitiontrees_group <- function(x, ...) {
  plots <- lapply(names(x), function(nm) {
    p <- plot(x[[nm]], ...)
    if (inherits(p, "ggplot"))
      p <- p + ggplot2::labs(caption = paste0("group: ", nm))
    print(p)
    p
  })
  invisible(setNames(plots, names(x)))
}
