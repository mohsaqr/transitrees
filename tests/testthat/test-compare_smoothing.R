# ---- Tests for compare_smoothing() ----

.cs_seqs <- function() {
  set.seed(1)
  replicate(60, sample(c("A", "B", "C"), 12, replace = TRUE),
            simplify = FALSE)
}

test_that("compare_smoothing returns a tidy 3-column data.frame", {
  res <- compare_smoothing(.cs_seqs(), max_depth = 3L, min_count = 5L)
  expect_s3_class(res, "data.frame")
  expect_named(res, c("smoothing", "n_nodes", "perplexity"))
  expect_equal(nrow(res), 5L)               # all five schemes by default
  expect_type(res$n_nodes, "integer")
  expect_type(res$perplexity, "double")
})

test_that("topology is invariant across smoothers (n_nodes all equal)", {
  res <- compare_smoothing(.cs_seqs(), max_depth = 3L, min_count = 5L)
  expect_equal(length(unique(res$n_nodes)), 1L)
})

test_that("rows follow the requested smoothing order and subset", {
  res <- compare_smoothing(.cs_seqs(),
                           smoothing = c("kneser_ney", "floor"),
                           max_depth = 2L, min_count = 5L)
  expect_equal(res$smoothing, c("kneser_ney", "floor"))
  expect_equal(nrow(res), 2L)
})

test_that("perplexity matches a hand-fitted tree", {
  seqs <- .cs_seqs()
  res  <- compare_smoothing(seqs, smoothing = "floor",
                            max_depth = 2L, min_count = 5L)
  ref  <- perplexity(context_tree(seqs, max_depth = 2L, min_count = 5L,
                                  smoothing = "floor"))
  expect_equal(res$perplexity, ref)
})

test_that("a non-character smoothing argument errors clearly", {
  expect_error(compare_smoothing(.cs_seqs(), smoothing = 1:3),
               "non-empty character vector")
  expect_error(compare_smoothing(.cs_seqs(), smoothing = character(0)),
               "non-empty character vector")
})

test_that("compare_smoothing on a fitted tree re-smooths (frozen topology)", {
  tr <- context_tree(.cs_seqs(), max_depth = 3L, min_count = 5L)
  res <- compare_smoothing(tr)
  expect_named(res, c("smoothing", "n_nodes", "perplexity"))
  ## topology frozen: every n_nodes equals the fitted tree's
  expect_true(all(res$n_nodes == n_nodes(tr)))
  ## perplexity matches a hand re-smooth
  ref <- perplexity(smooth_pathtree(tr, "kneser_ney"))
  expect_equal(res$perplexity[res$smoothing == "kneser_ney"], ref)
})
