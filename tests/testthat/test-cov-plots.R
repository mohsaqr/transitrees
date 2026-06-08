# Coverage-targeted tests for the plotting source files.
# These hit branches not exercised by the existing plot tests:
# degenerate (root-only / single-leaf / single-depth-1) trees, the
# zero-difference paths, the many-state palette fallback, and the
# show_prediction / NULL-abbrev option branches. Drawing code is forced
# with pdf(NULL); print(p) so aes-time helpers run too.

# ---- shared fixtures -------------------------------------------------

# Root-only tree (no depth-1 context survives) -> 0 edges, 1 leaf.
mk_root_only <- function() {
  set.seed(1)
  m <- matrix(sample(c("A", "B", "C"), 30 * 5, replace = TRUE), 30, 5)
  context_tree(m, max_depth = 2L, min_count = 1000L)
}

# Tree with a single depth-1 context: only "A" ever precedes a move.
mk_single_d1 <- function() {
  set.seed(3)
  seqs <- replicate(60, c(rep("A", sample(3:6, 1)), sample(c("B", "C"), 1)),
                    simplify = FALSE)
  context_tree(seqs, max_depth = 2L, min_count = 2L)
}

draw <- function(p) {
  pdf(NULL)
  on.exit(grDevices::dev.off())
  suppressWarnings(print(p))
  invisible(TRUE)
}

# ---- plot.R: degenerate layouts (lines 46, 78-80, 127, 162-164) ------

test_that("horizontal/dendrogram render a root-only tree (empty-edge + single-leaf branches)", {
  tr <- mk_root_only()
  expect_identical(nrow(tr$edges), 0L)
  ph <- plot(tr, style = "horizontal")
  pd <- plot(tr, style = "dendrogram")
  expect_s3_class(ph, "ggplot")
  expect_s3_class(pd, "ggplot")
  expect_true(draw(ph))
  expect_true(draw(pd))
})

# ---- plot.R: show_prediction = FALSE (line 293) ----------------------

test_that("horizontal style with show_prediction = FALSE uses context-only labels", {
  set.seed(1)
  m <- matrix(sample(c("A", "B", "C"), 30 * 12, replace = TRUE), 30, 12)
  tr <- context_tree(m, max_depth = 2L, min_count = 3L)
  p <- plot(tr, style = "horizontal", show_prediction = FALSE)
  expect_s3_class(p, "ggplot")
  ## context-only labels carry no "(state pct%)" prediction line
  body_layer <- p$layers[[2L]]
  expect_false(any(grepl("%", body_layer$data$label, fixed = TRUE)))
  expect_true(draw(p))
})

# ---- plot_interactive.R: degenerate scale + empty edges (9, 65-66) ---

test_that("interactive style renders a root-only tree (single-count scale, empty edges)", {
  skip_if_not_installed("visNetwork")
  tr <- mk_root_only()
  w <- plot(tr, style = "interactive")
  expect_s3_class(w, "htmlwidget")
  ## single node -> degenerate range -> midpoint size; no edges
  expect_identical(nrow(w$x$edges), 0L)
  expect_true(all(w$x$nodes$size > 0))
})

# ---- plot_icicle.R: NULL abbrev (33), single-d1 sep_df else (182-183),
#       depth-1 ring_width fallback (153), arc-text angle (8-10) --------

test_that("icicle renders with label_abbrev_fraction = NULL", {
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  set.seed(1)
  m <- matrix(sample(c("A", "B", "C"), 30 * 12, replace = TRUE), 30, 12)
  tr <- context_tree(m, max_depth = 2L, min_count = 3L)
  p <- plot(tr, style = "icicle", label_abbrev_fraction = NULL)
  expect_s3_class(p, "ggplot")
  expect_true(draw(p))
})

test_that("icicle renders a single-depth-1 tree (separator else-branch)", {
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  tr <- mk_single_d1()
  d1 <- Filter(function(k) k != "<root>" &&
                 length(strsplit(k, " -> ", fixed = TRUE)[[1L]]) == 1L,
               names(tr$nodes))
  expect_length(d1, 1L)
  p <- plot(tr, style = "icicle")
  expect_s3_class(p, "ggplot")
  expect_true(draw(p))
})

test_that("icicle renders a depth-1-only tree (ring-width fallback)", {
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  set.seed(3)
  seqs <- replicate(80, sample(c("A", "B", "C"), 5, replace = TRUE),
                    simplify = FALSE)
  tr <- context_tree(seqs, max_depth = 1L, min_count = 3L)
  p <- plot(tr, style = "icicle")
  expect_s3_class(p, "ggplot")
  expect_true(draw(p))
})

# ---- plot_diagnostics.R: okabe-ito many-state fallback (430) ----------

test_that("plot_pruning uses the hcl fallback palette for > 9 states", {
  set.seed(2)
  states <- paste0("S", 1:11)
  seqs <- replicate(200, sample(states, 6, replace = TRUE), simplify = FALSE)
  tr <- context_tree(seqs, max_depth = 2L, min_count = 2L)
  expect_gt(length(tr$alphabet), 9L)
  keys <- setdiff(names(tr$nodes), "<root>")
  deep <- Filter(function(k)
    length(strsplit(k, " -> ", fixed = TRUE)[[1L]]) == 2L, keys)
  expect_gt(length(deep), 0L)
  p <- plot_pruning(tr, deep[[1L]])
  expect_s3_class(p, "ggplot")
  expect_true(draw(p))
})

# ---- plot_diagnostics.R: plot_difference branches (505, 519-520,
#       563 tile M==0, 633 tree M==0) ------------------------------------

test_that("plot_difference errors on bad/length-wrong groups (line 505)", {
  set.seed(1)
  base <- replicate(40, sample(c("A", "B", "C"), 8, replace = TRUE),
                    simplify = FALSE)
  grp <- context_tree(c(base, base), group = rep(c("x", "y"), each = 40),
                      max_depth = 2L, min_count = 3L)
  expect_error(plot_difference(grp, groups = c("x", "z")),
               "must name two groups")
  expect_error(plot_difference(grp, groups = "x"),
               "must name two groups")
})

test_that("plot_difference errors when no contexts are shared (lines 519-520)", {
  set.seed(1)
  base <- replicate(40, sample(c("A", "B", "C"), 8, replace = TRUE),
                    simplify = FALSE)
  grp <- context_tree(c(base, base), group = rep(c("x", "y"), each = 40),
                      max_depth = 2L, min_count = 3L)
  expect_error(plot_difference(grp, min_count = 1e6L),
               "No shared contexts")
})

test_that("plot_difference handles a zero-difference object (tile + tree, M==0)", {
  ## Identical sequences in both groups -> all residuals/deltas are 0,
  ## driving the M==0 guards in both layouts.
  set.seed(1)
  base <- replicate(40, sample(c("A", "B", "C"), 8, replace = TRUE),
                    simplify = FALSE)
  grp <- context_tree(c(base, base), group = rep(c("x", "y"), each = 40),
                      max_depth = 2L, min_count = 3L)
  pt <- plot_difference(grp, measure = "probability")
  pr <- plot_difference(grp, layout = "tree")
  expect_s3_class(pt, "ggplot")
  expect_s3_class(pr, "ggplot")
  expect_true(all(pt$data$val == 0))
  expect_true(draw(pt))
  expect_true(draw(pr))
})

# ---- plot_pathways.R: plot_divergence default title (line 164) --------

test_that("plot_divergence builds a default title when none is given", {
  set.seed(1)
  m <- matrix(sample(c("A", "B", "C"), 30 * 12, replace = TRUE), 30, 12)
  tr <- context_tree(m, max_depth = 2L, min_count = 3L)
  p <- plot_divergence(tr, min_count = 1L)
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$title, "by KL divergence")
})
