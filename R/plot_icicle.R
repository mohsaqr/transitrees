# ---- Static icicle / partition visualisation ----

#' @noRd
.pt_arc_text_angle <- function(x, y) {
  ## Tangent angle (in degrees) for text running along the arc of a
  ## circular sunburst, with an upside-down flip so text stays
  ## right-side-up. (x, y) are the layout-assigned node centroids.
  ang <- atan2(y, x) * 180 / pi - 90
  ifelse(ang < -90, ang + 180,
    ifelse(ang > 90, ang - 180, ang))
}

#' @noRd
.plot_icicle <- function(tree,
                          label_min_fraction    = 0.03,
                          label_abbrev_fraction = 0.08,
                          gap_size              = 5,
                          gap_size_within       = 0.3,
                          abbrev_len            = 2L,
                          abbrev_fun            = NULL,
                          jitter                = 0.04,
                          jitter_seed           = 1L,
                          ...) {
  ## Two-tier (default) vs three-tier label policy:
  ## - If `label_abbrev_fraction` is NULL (default), the icicle is
  ##   binary: any slice >= label_min_fraction gets the full state
  ##   name; anything smaller gets no label.
  ## - If `label_abbrev_fraction` is supplied AND > label_min_fraction,
  ##   slices in [label_min_fraction, label_abbrev_fraction) get an
  ##   abbreviation (via `abbrev_fun` or `substr` to `abbrev_len`);
  ##   slices >= label_abbrev_fraction still get the full name.
  if (is.null(label_abbrev_fraction))
    label_abbrev_fraction <- label_min_fraction
  if (!requireNamespace("ggraph", quietly = TRUE) ||
      !requireNamespace("tidygraph", quietly = TRUE))
    stop("plot(style = 'icicle') requires 'ggraph' and 'tidygraph'.",
         call. = FALSE)
  stopifnot(
    is.numeric(label_min_fraction),    label_min_fraction    >= 0,
    is.numeric(label_abbrev_fraction), label_abbrev_fraction >= 0,
    is.numeric(gap_size),              gap_size        >= 0,
    is.numeric(gap_size_within),       gap_size_within >= 0,
    is.numeric(abbrev_len),            abbrev_len >= 1L,
    is.null(abbrev_fun) || is.function(abbrev_fun),
    is.numeric(jitter),                jitter >= 0
  )
  ## Abbreviation function: by default a simple substr to abbrev_len
  ## characters. Pass abbrev_fun = base::abbreviate (vowel-dropping
  ## smart shortening with uniqueness guarantees) or any custom
  ## function taking a character vector and returning short codes.
  .abbrev <- if (!is.null(abbrev_fun)) abbrev_fun
             else function(s) substr(s, 1L, as.integer(abbrev_len))

  node_names <- names(tree$nodes)
  is_leaf <- !(node_names %in% unique(tree$edges$parent))

  ## Each tile's fill encodes the *state* it represents — the last
  ## token of the pathway. Root has no last state.
  last_state <- vapply(node_names, .pt_last_state, character(1))
  state_levels <- tree$alphabet
  state_factor <- factor(last_state, levels = state_levels)

  pt_depth <- vapply(tree$nodes, function(n) n$depth, integer(1))
  pt_count <- vapply(tree$nodes, function(n) as.numeric(n$n), numeric(1))

  ## Label policy (size-based, not depth-based):
  ##   arc_fraction < label_min_fraction    → no label  (too thin to read)
  ##   arc_fraction < label_abbrev_fraction → abbreviation (abbrev_len chars)
  ##   otherwise                            → full state name
  ##
  ## arc_fraction is the *actual* visual arc fraction in the ggraph
  ## partition: each leaf contributes its count as its weight; each
  ## internal node's visual arc is the sum of its descendant leaves'
  ## counts. Using the node's own count would overstate the fraction
  ## for internal nodes (parent counts > sum of leaf-descendant
  ## counts whenever some observations end before reaching the leaf).
  total_leaf_count <- sum(pt_count[is_leaf])
  children_by_parent <- if (nrow(tree$edges) > 0L)
                          split(tree$edges$child, tree$edges$parent)
                        else list()
  .subtree_leaf <- function(ctx) {
    ch <- children_by_parent[[ctx]]
    if (is.null(ch) || length(ch) == 0L) return(pt_count[[ctx]])
    sum(vapply(ch, .subtree_leaf, numeric(1)))
  }
  subtree_leaf_w <- vapply(node_names, .subtree_leaf, numeric(1))
  arc_fraction <- subtree_leaf_w / max(total_leaf_count, 1)
  ## Pre-compute abbreviations for every alphabet state once so the
  ## abbrev_fun sees the whole vocabulary at once (matters for
  ## base::abbreviate which needs the full set to enforce
  ## uniqueness). Then index by state.
  state_abbrev <- as.character(.abbrev(state_levels))
  abbrev_lookup <- setNames(state_abbrev, state_levels)
  abbreviated   <- ifelse(is.na(last_state) | last_state == .ROOT_LABEL,
                          NA_character_,
                          abbrev_lookup[last_state])
  label_text <- ifelse(
    pt_depth == 0L | last_state == .ROOT_LABEL | is.na(last_state) |
      arc_fraction < label_min_fraction,
    NA_character_,
    ifelse(arc_fraction < label_abbrev_fraction,
           abbreviated,
           last_state)
  )

  nodes <- data.frame(
    name         = node_names,
    state        = state_factor,
    count        = pt_count,
    pt_depth     = pt_depth,
    label_text   = label_text,
    arc_fraction = arc_fraction,
    stringsAsFactors = FALSE
  )
  ## ggraph's partition layout derives parent tile sizes from children;
  ## non-zero internal-node weights trigger a "Non-leaf weights ignored"
  ## message. Pass leaf-only weights to communicate that explicitly.
  nodes$leaf_weight <- ifelse(is_leaf, nodes$count, 0)

  edges <- if (nrow(tree$edges) > 0L)
    data.frame(from = tree$edges$parent,
               to   = tree$edges$child,
               stringsAsFactors = FALSE) else
    data.frame(from = character(0), to = character(0),
               stringsAsFactors = FALSE)

  tg <- tidygraph::tbl_graph(nodes = nodes, edges = edges,
                             directed = TRUE)

  pal <- .pt_state_palette(state_levels)

  ## Build the layout explicitly so we can pull out the angular
  ## boundaries of the depth-1 slices.
  lay <- ggraph::create_layout(tg, layout = "partition",
                                circular = TRUE,
                                weight   = leaf_weight)

  ## Optional radial jitter on leaf tiles: each leaf's outer radius
  ## is perturbed by up to `jitter * ring_width`, giving the outer
  ## edge of the icicle an organic, slightly-irregular look without
  ## breaking the within-ring layout. RNG state is saved & restored
  ## so we don't clobber the caller's seed.
  if (jitter > 0 && any(lay$leaf)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv))
                  get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (!is.null(old_seed))
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(as.integer(jitter_seed))
    ring_width <- mean(diff(sort(unique(lay$r))), na.rm = TRUE)
    if (!is.finite(ring_width) || ring_width <= 0)
      ring_width <- 0.2
    leaf_idx <- which(lay$leaf)
    delta <- stats::runif(length(leaf_idx),
                          -jitter * ring_width,
                          +jitter * ring_width)
    lay$r[leaf_idx] <- lay$r[leaf_idx] + delta
  }

  ## Compute the radial separator segments — one per depth-1 boundary
  ## angle, extending from the inner edge of the depth-1 ring all the
  ## way to the outer edge of the deepest ring. These read as wedge
  ## dividers in a pie chart: the separation is defined at the first
  ## level, then "stays" out through every deeper ring.
  ##
  ## ggraph's circular partition uses the convention
  ##   theta_cartesian = pi/2 - theta_partition
  ## so (x, y) = (r * sin(theta), r * cos(theta)).
  d1 <- lay[lay$pt_depth == 1L, , drop = FALSE]
  if (nrow(d1) > 1L) {
    angles <- sort(unique(c(d1$start, d1$end) %% (2 * pi)))
    r_in   <- min(d1$r0)
    r_out  <- max(lay$r, na.rm = TRUE)
    sep_df <- data.frame(
      x    = r_in  * sin(angles),
      y    = r_in  * cos(angles),
      xend = r_out * sin(angles),
      yend = r_out * cos(angles)
    )
  } else {
    sep_df <- data.frame(x = numeric(0), y = numeric(0),
                         xend = numeric(0), yend = numeric(0))
  }

  ggraph::ggraph(lay) +
    ggraph::geom_node_arc_bar(
      ggplot2::aes(fill = .data$state),
      colour    = "white",
      linewidth = gap_size_within,
      na.rm     = TRUE
    ) +
    ggplot2::geom_segment(
      data = sep_df,
      ggplot2::aes(x = .data$x, y = .data$y,
                   xend = .data$xend, yend = .data$yend),
      colour    = "white",
      linewidth = gap_size,
      lineend   = "butt"
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = .data$label_text,
                   angle = .pt_arc_text_angle(.data$x, .data$y)),
      size = 2.8, colour = "grey15", fontface = "bold",
      na.rm = TRUE
    ) +
    ggplot2::scale_fill_manual(values = pal, name = "State",
                               na.value = "grey92",
                               na.translate = FALSE,
                               drop = FALSE) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40"),
      legend.position = "right"
    ) +
    ggplot2::labs(
      title    = "Context tree (sunburst)",
      subtitle = sprintf(
        "%d nodes, depth <= %d;  arc width = count;  fill = state;  rings = depth",
        length(tree$nodes), tree$max_depth)
    )
}
