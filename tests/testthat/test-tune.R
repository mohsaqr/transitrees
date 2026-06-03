# Tests for tune_tree

mk_iid_data <- function(n = 60, len = 14, seed = 1,
                        states = c("A","B","C","D")) {
  set.seed(seed)
  matrix(sample(states, n * len, replace = TRUE), n, len)
}

test_that("tune_tree returns a sorted transitrees_tune data.frame", {
  m <- mk_iid_data()
  tg <- tune_tree(m, max_depth = 1L:2L, min_count = c(3L, 5L),
                      smoothing = "floor",
                      prune = c(FALSE, TRUE), folds = 3L)
  expect_s3_class(tg, "transitrees_tune")
  expect_s3_class(tg, "data.frame")
  expect_named(tg, c("max_depth", "nmin", "smoothing", "prune",
                     "logLik", "n_scored", "perplexity", "n_nodes_avg"))
  expect_equal(nrow(tg), 2L * 2L * 1L * 2L)
  expect_equal(order(tg$perplexity, na.last = TRUE), seq_len(nrow(tg)))
  expect_false(is.null(attr(tg, "best")))
})

test_that("tune_tree compares smoothing methods cleanly", {
  m <- mk_iid_data()
  tg <- tune_tree(m, max_depth = 1L:2L, min_count = 3L,
                      smoothing = c("floor", "kneser_ney"),
                      prune = FALSE, folds = 3L)
  expect_true(any(grepl("^floor", tg$smoothing)))
  expect_true(any(grepl("^kneser_ney", tg$smoothing)))
})

test_that("tune_tree on i.i.d. data prefers smallest max_depth", {
  m <- mk_iid_data(seed = 7)
  tg <- tune_tree(m, max_depth = 1L:3L, min_count = 3L,
                      smoothing = "floor", prune = FALSE, folds = 4L)
  best <- attr(tg, "best")
  expect_equal(best$max_depth, min(tg$max_depth))
})

test_that("tune_tree print method runs without error", {
  m <- mk_iid_data()
  tg <- tune_tree(m, max_depth = 1L, min_count = 3L,
                      smoothing = "floor", prune = FALSE, folds = 3L)
  expect_output(print(tg), "<transitrees_tune>")
  expect_output(print(tg), "best")
})

test_that("tune_tree is reproducible across calls with same seed", {
  m <- mk_iid_data()
  a <- tune_tree(m, max_depth = 1L:2L, min_count = 3L,
                     smoothing = "floor", prune = FALSE, folds = 3L,
                     seed = 42L)
  b <- tune_tree(m, max_depth = 1L:2L, min_count = 3L,
                     smoothing = "floor", prune = FALSE, folds = 3L,
                     seed = 42L)
  expect_equal(a$perplexity, b$perplexity)
})

test_that("plot.transitrees_tune returns a ggplot", {
  m <- mk_iid_data()
  tg <- tune_tree(m, max_depth = 1L:2L, min_count = c(3L, 5L),
                      smoothing = "floor",
                      prune = c(FALSE, TRUE), folds = 3L)
  expect_s3_class(plot(tg), "ggplot")
})

test_that("tune_tree errors when k > n_sequences", {
  m <- matrix(sample(c("A","B"), 4 * 5, replace = TRUE), 4, 5)
  expect_error(tune_tree(m, folds = 10L), "Not enough sequences")
})

test_that("tune_tree reshapes long-format data (action =), not row-by-row", {
  set.seed(3)
  long <- data.frame(
    id = rep(letters[1:12], each = 6),
    o  = rep(1:6, times = 12),
    s  = sample(c("A","B","C"), 72, replace = TRUE),
    stringsAsFactors = FALSE)
  t1 <- tune_tree(long, actor = "id", order = "o", action = "s",
                      max_depth = 1L:2L, min_count = 2L, folds = 3L)
  ## equals tuning the explicitly reshaped wide frame
  wide <- prepare_input(long, actor = "id", order = "o", action = "s")
  t2 <- tune_tree(wide, max_depth = 1L:2L, min_count = 2L, folds = 3L)
  expect_equal(t1$perplexity, t2$perplexity)
  ## actor without action is rejected (not silently read as wide)
  expect_error(tune_tree(long, actor = "id"), "action =")
})
