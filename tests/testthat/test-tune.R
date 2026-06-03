# Tests for tune_pathtree

mk_iid_data <- function(n = 60, len = 14, seed = 1,
                        states = c("A","B","C","D")) {
  set.seed(seed)
  matrix(sample(states, n * len, replace = TRUE), n, len)
}

test_that("tune_pathtree returns a sorted pathtree_tune data.frame", {
  m <- mk_iid_data()
  tg <- tune_pathtree(m, max_depth = 1L:2L, min_count = c(3L, 5L),
                      smoothing = "floor",
                      prune = c(FALSE, TRUE), folds = 3L)
  expect_s3_class(tg, "pathtree_tune")
  expect_s3_class(tg, "data.frame")
  expect_named(tg, c("max_depth", "nmin", "smoothing", "prune",
                     "logLik", "n_scored", "perplexity", "n_nodes_avg"))
  expect_equal(nrow(tg), 2L * 2L * 1L * 2L)
  expect_equal(order(tg$perplexity, na.last = TRUE), seq_len(nrow(tg)))
  expect_false(is.null(attr(tg, "best")))
})

test_that("tune_pathtree compares smoothing methods cleanly", {
  m <- mk_iid_data()
  tg <- tune_pathtree(m, max_depth = 1L:2L, min_count = 3L,
                      smoothing = c("floor", "kneser_ney"),
                      prune = FALSE, folds = 3L)
  expect_true(any(grepl("^floor", tg$smoothing)))
  expect_true(any(grepl("^kneser_ney", tg$smoothing)))
})

test_that("tune_pathtree on i.i.d. data prefers smallest max_depth", {
  m <- mk_iid_data(seed = 7)
  tg <- tune_pathtree(m, max_depth = 1L:3L, min_count = 3L,
                      smoothing = "floor", prune = FALSE, folds = 4L)
  best <- attr(tg, "best")
  expect_equal(best$max_depth, min(tg$max_depth))
})

test_that("tune_pathtree print method runs without error", {
  m <- mk_iid_data()
  tg <- tune_pathtree(m, max_depth = 1L, min_count = 3L,
                      smoothing = "floor", prune = FALSE, folds = 3L)
  expect_output(print(tg), "<pathtree_tune>")
  expect_output(print(tg), "best")
})

test_that("tune_pathtree is reproducible across calls with same seed", {
  m <- mk_iid_data()
  a <- tune_pathtree(m, max_depth = 1L:2L, min_count = 3L,
                     smoothing = "floor", prune = FALSE, folds = 3L,
                     seed = 42L)
  b <- tune_pathtree(m, max_depth = 1L:2L, min_count = 3L,
                     smoothing = "floor", prune = FALSE, folds = 3L,
                     seed = 42L)
  expect_equal(a$perplexity, b$perplexity)
})

test_that("plot.pathtree_tune returns a ggplot", {
  m <- mk_iid_data()
  tg <- tune_pathtree(m, max_depth = 1L:2L, min_count = c(3L, 5L),
                      smoothing = "floor",
                      prune = c(FALSE, TRUE), folds = 3L)
  expect_s3_class(plot(tg), "ggplot")
})

test_that("tune_pathtree errors when k > n_sequences", {
  m <- matrix(sample(c("A","B"), 4 * 5, replace = TRUE), 4, 5)
  expect_error(tune_pathtree(m, folds = 10L), "Not enough sequences")
})
