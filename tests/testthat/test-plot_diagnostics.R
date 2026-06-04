# Tests for plot_distributions() and plot_predictive()

mk_diag_tree <- function(seed = 1L) {
  set.seed(seed)
  seqs <- replicate(100, sample(c("A", "B", "C"), 10, replace = TRUE),
                    simplify = FALSE)
  context_tree(seqs, max_depth = 2L)
}

test_that("plot_distributions returns a ggplot", {
  tree <- mk_diag_tree()
  p <- plot_distributions(tree, top = 6)
  expect_s3_class(p, "ggplot")
})

test_that("plot_distributions accepts an explicit context set", {
  tree <- mk_diag_tree()
  p <- plot_distributions(tree, contexts = c("(start)", "A"))
  expect_s3_class(p, "ggplot")
})

test_that("plot_distributions errors when no contexts qualify", {
  tree <- mk_diag_tree()
  expect_error(plot_distributions(tree, min_count = 1e6), "No contexts")
})

test_that("plot_predictive returns a ggplot for both types", {
  tree <- mk_diag_tree()
  set.seed(3)
  new <- replicate(12, sample(c("A", "B", "C"), 10, replace = TRUE),
                   simplify = FALSE)
  expect_s3_class(plot_predictive(tree, new, type = "position"), "ggplot")
  expect_s3_class(plot_predictive(tree, new, type = "ecdf"), "ggplot")
})

test_that("plot_predictive errors when nothing scores", {
  tree <- mk_diag_tree()
  # all-out-of-alphabet held-out data -> zero scored positions
  expect_error(plot_predictive(tree, list(c("Z", "Z", "Z"))),
               "No held-out positions")
})

test_that("plot_difference returns a ggplot for a two-group object", {
  set.seed(1)
  gx <- replicate(40, sample(c("A","B","C"), 8, replace = TRUE,
                             prob = c(.2,.6,.2)), simplify = FALSE)
  gy <- replicate(40, sample(c("A","B","C"), 8, replace = TRUE,
                             prob = c(.2,.2,.6)), simplify = FALSE)
  grp <- context_tree(c(gx, gy), group = rep(c("x","y"), each = 40),
                      max_depth = 2L, min_count = 3L)
  expect_s3_class(plot_difference(grp), "ggplot")
  expect_s3_class(plot_difference(grp, depth = 1), "ggplot")      # order-1
  # significance overlay from a comparison
  cmp <- compare_groups(grp, iter = 99L)
  expect_s3_class(plot_difference(grp, comparison = cmp), "ggplot")
})

test_that("plot_difference requires naming two groups when there are >2", {
  set.seed(2)
  g <- context_tree(
    replicate(90, sample(c("A","B","C"), 8, replace = TRUE), simplify = FALSE),
    group = rep(c("p","q","r"), each = 30), max_depth = 1L, min_count = 3L)
  expect_error(plot_difference(g), "name the two")
  expect_s3_class(plot_difference(g, groups = c("p", "r")), "ggplot")
})

test_that("plot_difference layout = tree colours the phylogram by group difference", {
  set.seed(1)
  gx <- replicate(50, sample(c("A","B","C"), 8, replace = TRUE,
                             prob = c(.2,.6,.2)), simplify = FALSE)
  gy <- replicate(50, sample(c("A","B","C"), 8, replace = TRUE,
                             prob = c(.2,.2,.6)), simplify = FALSE)
  grp <- context_tree(c(gx, gy), group = rep(c("x","y"), each = 50),
                      max_depth = 2L, min_count = 3L)
  p <- plot_difference(grp, layout = "tree")
  expect_s3_class(p, "ggplot")
  pdf(NULL); on.exit(dev.off()); expect_silent(print(p))
})

test_that("plot_difference measure = residual/probability both render", {
  set.seed(1)
  gx <- replicate(40, sample(c("A","B","C"), 8, replace = TRUE,
                             prob = c(.2,.6,.2)), simplify = FALSE)
  gy <- replicate(40, sample(c("A","B","C"), 8, replace = TRUE,
                             prob = c(.2,.2,.6)), simplify = FALSE)
  grp <- context_tree(c(gx, gy), group = rep(c("x","y"), each = 40),
                      max_depth = 2L, min_count = 3L)
  expect_s3_class(plot_difference(grp, measure = "residual"), "ggplot")
  expect_s3_class(plot_difference(grp, measure = "probability"), "ggplot")
})
