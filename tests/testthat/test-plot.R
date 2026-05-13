# Tests for Phase G: plot.pathtree style dispatcher (ggplot only)

mk_plot_tree <- function() {
  set.seed(1)
  m <- matrix(sample(c("A","B","C"), 30 * 12, replace = TRUE), 30, 12)
  context_tree(m, max_depth = 2L, nmin = 3L)
}

test_that("plot.pathtree default style is dendrogram and returns a ggplot", {
  tr <- mk_plot_tree()
  p  <- plot(tr)
  expect_s3_class(p, "ggplot")
})

test_that("plot.pathtree style='dendrogram' returns a ggplot", {
  tr <- mk_plot_tree()
  p  <- plot(tr, style = "dendrogram")
  expect_s3_class(p, "ggplot")
})

test_that("plot.pathtree style='dendrogram' accepts custom point_size_range", {
  tr <- mk_plot_tree()
  p  <- plot(tr, style = "dendrogram", point_size_range = c(3, 12))
  expect_s3_class(p, "ggplot")
})

test_that("plot.pathtree style='dendrogram' edges carry per-edge linewidth", {
  tr <- mk_plot_tree()
  p  <- plot(tr, style = "dendrogram")
  ## The first layer is the edge geom_path; its data should carry an
  ## edge_weight column (= child's count) for the linewidth aes.
  edge_layer <- p$layers[[1L]]
  expect_true("edge_weight" %in% names(edge_layer$data))
  expect_true(all(edge_layer$data$edge_weight > 0))
})

test_that("plot.pathtree style='horizontal' returns a ggplot with x = depth", {
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

test_that("plot.pathtree style='horizontal' accepts custom size ranges", {
  tr <- mk_plot_tree()
  p  <- plot(tr, style = "horizontal",
             point_size_range = c(3, 12),
             edge_size_range  = c(0.5, 4))
  expect_s3_class(p, "ggplot")
})

test_that("plot.pathtree style='icicle' returns a ggplot when ggraph + tidygraph are installed", {
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  tr <- mk_plot_tree()
  p  <- plot(tr, style = "icicle")
  expect_s3_class(p, "ggplot")
})

test_that("plot.pathtree style='interactive' returns an htmlwidget when collapsibleTree is installed", {
  skip_if_not_installed("collapsibleTree")
  tr <- mk_plot_tree()
  w  <- plot(tr, style = "interactive")
  expect_s3_class(w, "htmlwidget")
})

test_that("plot.pathtree style='icicle' uses size-based abbreviation", {
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
  tr <- context_tree(seqs, max_depth = 3L, nmin = 5L)
  ## Default mode (label_abbrev_fraction = NULL) is binary:
  ## any slice >= label_min_fraction gets the full state name; the
  ## rest get no label.
  p  <- plot(tr, style = "icicle")
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

test_that("plot.pathtree style='icicle' respects override args", {
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
  tr <- context_tree(seqs, max_depth = 3L, nmin = 5L)
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

test_that("plot.pathtree errors on unrecognised style", {
  tr <- mk_plot_tree()
  expect_error(plot(tr, style = "sankey"))
  expect_error(plot(tr, style = "treemap"))
})
