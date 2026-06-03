# ---- Tests for compare_pruning() ----

.cp_tree <- function() {
  set.seed(1)
  seqs <- replicate(80, sample(c("A", "B", "C"), 14, replace = TRUE),
                    simplify = FALSE)
  context_tree(seqs, max_depth = 4L, min_count = 3L)
}

test_that("compare_pruning returns a tidy 3-column data.frame", {
  res <- compare_pruning(.cp_tree())
  expect_s3_class(res, "data.frame")
  expect_named(res, c("criterion", "n_nodes", "reduction_pct"))
  expect_equal(nrow(res), 4L)                    # all four criteria by default
  expect_type(res$n_nodes, "integer")
  expect_equal(res$criterion, c("G2", "KL", "AIC", "BIC"))
})

test_that("n_nodes / reduction_pct match a hand-pruned tree", {
  tree <- .cp_tree()
  res  <- compare_pruning(tree, criterion = "G2", alpha = 0.05)
  ref  <- length(prune_pathtree(tree, criterion = "G2",
                                alpha = 0.05)$nodes)
  expect_equal(res$n_nodes, ref)
  expect_equal(res$reduction_pct,
               round(100 * (1 - ref / length(tree$nodes)), 1))
})

test_that("criterion order and subset are honoured", {
  res <- compare_pruning(.cp_tree(), criterion = c("BIC", "G2"))
  expect_equal(res$criterion, c("BIC", "G2"))
  expect_equal(nrow(res), 2L)
})

test_that("reduction_pct is between 0 and 100 and never negative", {
  res <- compare_pruning(.cp_tree())
  expect_true(all(res$reduction_pct >= 0 & res$reduction_pct <= 100))
})

test_that("a non-pathtree or empty criterion errors clearly", {
  expect_error(compare_pruning(list()), "pathtree")
  expect_error(compare_pruning(.cp_tree(), criterion = character(0)),
               "non-empty character vector")
})
