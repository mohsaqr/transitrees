# Tests for Phase G: plot.transitiontrees style dispatcher (ggplot only)

mk_plot_tree <- function() {
  set.seed(1)
  m <- matrix(sample(c("A","B","C"), 30 * 12, replace = TRUE), 30, 12)
  context_tree(m, max_depth = 2L, min_count = 3L)
}

test_that("plot.transitiontrees default style is horizontal and returns a ggplot", {
  tr <- mk_plot_tree()
  p  <- plot(tr)
  expect_s3_class(p, "ggplot")
  ## default == horizontal: x-axis is depth, same as style = "horizontal"
  expect_identical(p$labels$x, plot(tr, style = "horizontal")$labels$x)
})

test_that("plot.transitiontrees style='dendrogram' returns a ggplot", {
  tr <- mk_plot_tree()
  p  <- plot(tr, style = "dendrogram")
  expect_s3_class(p, "ggplot")
})

test_that("plot.transitiontrees style='dendrogram' accepts custom point_size_range", {
  tr <- mk_plot_tree()
  p  <- plot(tr, style = "dendrogram", point_size_range = c(3, 12))
  expect_s3_class(p, "ggplot")
})

test_that("plot.transitiontrees style='dendrogram' edges carry per-edge linewidth", {
  tr <- mk_plot_tree()
  p  <- plot(tr, style = "dendrogram")
  ## The first layer is the edge geom_path; its data should carry an
  ## edge_weight column (= child's count) for the linewidth aes.
  edge_layer <- p$layers[[1L]]
  expect_true("edge_weight" %in% names(edge_layer$data))
  expect_true(all(edge_layer$data$edge_weight > 0))
})

test_that("plot.transitiontrees style='horizontal' returns a ggplot with x = depth", {
  tr <- mk_plot_tree()
  p  <- plot(tr, style = "horizontal")
  expect_s3_class(p, "ggplot")
  ## Edge layer carries edge_weight column for linewidth aes.
  edge_layer <- p$layers[[1L]]
  expect_true("edge_weight" %in% names(edge_layer$data))
  ## The body-node layer's x is the integer depth.
  body_layer <- p$layers[[2L]]
  expect_true(all(body_layer$data$x == as.integer(body_layer$data$x)))
})

test_that("plot.transitiontrees style='horizontal' accepts custom size ranges", {
  tr <- mk_plot_tree()
  p  <- plot(tr, style = "horizontal",
             point_size_range = c(3, 12),
             edge_size_range  = c(0.5, 4))
  expect_s3_class(p, "ggplot")
})

test_that("plot.transitiontrees style='icicle' returns a ggplot when ggraph + tidygraph are installed", {
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  tr <- mk_plot_tree()
  p  <- plot(tr, style = "icicle")
  expect_s3_class(p, "ggplot")
})

test_that("plot.transitiontrees style='interactive' returns an htmlwidget with sized nodes and edges", {
  skip_if_not_installed("visNetwork")
  tr <- mk_plot_tree()
  w  <- plot(tr, style = "interactive")
  expect_s3_class(w, "htmlwidget")
  expect_s3_class(w, "visNetwork")
  ## node size encodes count, edge width encodes child-count flow
  expect_true("size"  %in% names(w$x$nodes))
  expect_true("width" %in% names(w$x$edges))
  expect_true(all(w$x$nodes$size  >= 10 & w$x$nodes$size  <= 45))
  expect_true(all(w$x$edges$width >= 1  & w$x$edges$width <= 10))
})

test_that("plot.transitiontrees style='icicle' uses size-based abbreviation", {
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  set.seed(7)
  states <- c("Active", "Average", "Disengaged")
  P <- matrix(c(0.7, 0.2, 0.1,
                0.1, 0.7, 0.2,
                0.2, 0.2, 0.6),
              nrow = 3, byrow = TRUE,
              dimnames = list(states, states))
  seqs <- lapply(seq_len(80L), function(i) {
    s <- character(15L); s[1L] <- sample(states, 1L)
    for (t in 2:15) s[t] <- sample(states, 1L, prob = P[s[t - 1L], ])
    s
  })
  tr <- context_tree(seqs, max_depth = 3L, min_count = 5L)
  ## Binary mode (label_abbrev_fraction == label_min_fraction): any slice
  ## >= label_min_fraction gets the full state name; the rest no label.
  ## (The DEFAULT is now three-tier, so we request binary explicitly.)
  p  <- plot(tr, style = "icicle",
             label_min_fraction = 0.10, label_abbrev_fraction = 0.10)
  graph_df <- p$data
  expect_true(all(c("label_text","arc_fraction") %in% names(graph_df)))
  labelled <- graph_df[!is.na(graph_df$label_text), ]
  if (nrow(labelled) > 0L)
    expect_true(all(labelled$label_text %in% states))

  ## Three-tier mode kicks in when label_abbrev_fraction is supplied
  ## and > label_min_fraction. Below that, slices get the
  ## abbreviation.
  p3 <- plot(tr, style = "icicle",
             label_min_fraction    = 0.05,
             label_abbrev_fraction = 0.20)
  d3 <- p3$data
  small <- d3[!is.na(d3$label_text) & d3$arc_fraction < 0.20, ]
  big   <- d3[!is.na(d3$label_text) & d3$arc_fraction >= 0.20, ]
  if (nrow(big) > 0L)
    expect_true(all(big$label_text %in% states))
  if (nrow(small) > 0L)
    expect_true(all(nchar(small$label_text) <= 3L))
})

test_that("plot.transitiontrees style='icicle' respects override args", {
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  set.seed(7)
  states <- c("Active", "Average", "Disengaged")
  P <- matrix(c(0.7, 0.2, 0.1,
                0.1, 0.7, 0.2,
                0.2, 0.2, 0.6),
              nrow = 3, byrow = TRUE,
              dimnames = list(states, states))
  seqs <- lapply(seq_len(80L), function(i) {
    s <- character(15L); s[1L] <- sample(states, 1L)
    for (t in 2:15) s[t] <- sample(states, 1L, prob = P[s[t - 1L], ])
    s
  })
  tr <- context_tree(seqs, max_depth = 3L, min_count = 5L)
  ## abbrev_fraction = 0 → never abbreviate (every labelled tile shows
  ## the full state name).
  p_full <- plot(tr, style = "icicle", label_abbrev_fraction = 0)
  d_full <- p_full$data
  d_full <- d_full[!is.na(d_full$label_text), ]
  expect_true(all(d_full$label_text %in% states))
  ## abbrev_fraction = 1 → always abbreviate (every labelled tile is
  ## <= 3 chars).
  p_abbr <- plot(tr, style = "icicle", label_abbrev_fraction = 1)
  d_abbr <- p_abbr$data
  d_abbr <- d_abbr[!is.na(d_abbr$label_text), ]
  if (nrow(d_abbr) > 0L)
    expect_true(all(nchar(d_abbr$label_text) <= 3L))
  ## gap_size override flows through to the arc-bar linewidth.
  p_gap <- plot(tr, style = "icicle", gap_size = 4)
  expect_s3_class(p_gap, "ggplot")
})

test_that("plot.transitiontrees errors on unrecognised style", {
  tr <- mk_plot_tree()
  expect_error(plot(tr, style = "sankey"))
  expect_error(plot(tr, style = "treemap"))
})

# ---- plot_pathways() (pathway x next-move heatmap) ----------------------

test_that("plot_pathways returns a ggplot for each sort_by", {
  tr <- mk_plot_tree()
  for (sb in c("count", "divergence", "depth"))
    expect_s3_class(plot_pathways(tr, sort_by = sb, min_count = 1L), "ggplot")
})

test_that("plot_pathways honours top, title, and show_flips", {
  tr <- mk_plot_tree()
  p  <- plot_pathways(tr, top = 4L, min_count = 1L, title = "custom")
  expect_s3_class(p, "ggplot")
  ## at most `top` distinct pathways on the y axis
  expect_lte(length(unique(p$data$pathway)), 4L)
  expect_identical(p$labels$title, "custom")
  ## show_flips = FALSE drops the caret prefix from every label
  p_noflip <- plot_pathways(tr, min_count = 1L, show_flips = FALSE)
  expect_false(any(grepl("^> ", as.character(p_noflip$data$label))))
})

test_that("plot_pathways errors when no pathway clears the threshold", {
  tr <- mk_plot_tree()
  expect_error(plot_pathways(tr, min_count = 1e6),
               "No pathways meet the threshold")
})

# ---- plot_divergence() (KL lollipop) ------------------------------------

test_that("plot_divergence returns a ggplot and honours a custom title", {
  tr <- mk_plot_tree()
  p  <- plot_divergence(tr, top = 6L, min_count = 1L, title = "div")
  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$title, "div")
})

test_that("plot_divergence errors when no pathway clears the threshold", {
  tr <- mk_plot_tree()
  expect_error(plot_divergence(tr, min_count = 1e6),
               "No pathways meet the threshold")
})
