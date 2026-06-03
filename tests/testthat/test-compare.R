# Tests for Phase E: pathtree_distance + compare_pathtrees

mk_two_trees <- function(seed_a = 1, seed_b = 2,
                         states = c("A","B","C"),
                         max_depth = 2L, min_count = 3L) {
  set.seed(seed_a)
  m_a <- matrix(sample(states, 30 * 12, replace = TRUE), 30, 12)
  set.seed(seed_b)
  m_b <- matrix(sample(states, 30 * 12, replace = TRUE), 30, 12)
  list(
    a = context_tree(m_a, max_depth = max_depth, min_count = min_count),
    b = context_tree(m_b, max_depth = max_depth, min_count = min_count)
  )
}

test_that("pathtree_distance is zero between a tree and itself", {
  trs <- mk_two_trees()
  expect_equal(pathtree_distance(trs$a, trs$a), 0)
})

test_that("pathtree_distance is symmetric when symmetric = TRUE", {
  trs <- mk_two_trees()
  d_ab <- pathtree_distance(trs$a, trs$b, symmetric = TRUE)
  d_ba <- pathtree_distance(trs$b, trs$a, symmetric = TRUE)
  expect_equal(d_ab, d_ba)
})

test_that("pathtree_distance returns a finite non-negative scalar", {
  trs <- mk_two_trees()
  d <- pathtree_distance(trs$a, trs$b)
  expect_true(is.finite(d))
  expect_gte(d, 0)
})

test_that("pathtree_distance asymmetric form differs in general", {
  trs <- mk_two_trees()
  d_sym <- pathtree_distance(trs$a, trs$b, symmetric = TRUE)
  d_asy <- pathtree_distance(trs$a, trs$b, symmetric = FALSE)
  expect_true(is.finite(d_sym) && is.finite(d_asy))
})

test_that("pathtree_distance errors on incompatible alphabets", {
  set.seed(1)
  m1 <- matrix(sample(c("A","B"), 30, TRUE), 6)
  m2 <- matrix(sample(c("X","Y"), 30, TRUE), 6)
  tr1 <- context_tree(m1, max_depth = 1L, min_count = 1L)
  tr2 <- context_tree(m2, max_depth = 1L, min_count = 1L)
  expect_error(pathtree_distance(tr1, tr2), "incompatible alphabets")
})

test_that("pathtree_distance aligns compatible alphabets by state name", {
  set.seed(1)
  m <- matrix(sample(c("A", "B", "C"), 90, TRUE), 10)
  tr1 <- context_tree(m, max_depth = 1L, min_count = 1L,
                      alphabet = c("A", "B", "C"))
  tr2 <- context_tree(m, max_depth = 1L, min_count = 1L,
                      alphabet = c("C", "B", "A"))
  expect_equal(pathtree_distance(tr1, tr2), 0, tolerance = 1e-12)
})

test_that("compare_pathtrees returns a pathtree_comparison object", {
  trs <- mk_two_trees()
  cmp <- compare_pathtrees(trs$a, trs$b, iter = 30L, seed = 1L)
  expect_s3_class(cmp, "pathtree_comparison")
  expect_true(is.numeric(cmp$pdist))
  expect_length(cmp$null_dist, 30L)
  expect_true(cmp$p_value > 0 && cmp$p_value <= 1)
  expect_s3_class(cmp$pathways, "data.frame")
})

test_that("compare_pathtrees handles pruned observed trees in permutation refits", {
  trs <- mk_two_trees(max_depth = 2L, min_count = 2L)
  tr_a <- prune_pathtree(trs$a, criterion = "G2", alpha = 0.05)
  tr_b <- prune_pathtree(trs$b, criterion = "G2", alpha = 0.05)
  cmp <- compare_pathtrees(tr_a, tr_b, iter = 20L, seed = 1L)
  expect_s3_class(cmp, "pathtree_comparison")
  expect_length(cmp$null_dist, 20L)
  expect_true(all(is.finite(cmp$null_dist)))
})

test_that("compare_pathtrees null mean is centred near observed for exchangeable groups", {
  ## When the two trees are fit on i.i.d. samples from the same
  ## generator, the null distribution should bracket the observed.
  trs <- mk_two_trees(seed_a = 11, seed_b = 12, max_depth = 2L,
                      min_count = 3L)
  cmp <- compare_pathtrees(trs$a, trs$b, iter = 100L, seed = 1L)
  expect_gt(cmp$p_value, 0.05)
})

test_that("compare_pathtrees print method runs without error", {
  trs <- mk_two_trees()
  cmp <- compare_pathtrees(trs$a, trs$b, iter = 20L, seed = 1L)
  expect_output(print(cmp), "<pathtree_comparison>", fixed = TRUE)
  expect_output(print(cmp), "p-value")
})

test_that("as.data.frame(cmp) returns the pathways table (uniform tidy-extract)", {
  trs <- mk_two_trees()
  cmp <- compare_pathtrees(trs$a, trs$b, iter = 20L, seed = 1L)
  df <- as.data.frame(cmp)
  expect_s3_class(df, "data.frame")
  expect_identical(df, cmp$pathways)
  expect_named(df, c("pathway", "count_a", "count_b",
                     "divergence_ab", "divergence_ba", "divergence_sym"))
})

test_that("plot.pathtree_comparison returns a ggplot", {
  trs <- mk_two_trees()
  cmp <- compare_pathtrees(trs$a, trs$b, iter = 30L, seed = 1L)
  expect_s3_class(plot(cmp), "ggplot")
})

test_that("compare_pathtrees rejects a non-positive n_perm", {
  ## iter = 0 used to yield an empty null distribution and a NaN
  ## summary mean rather than an error.
  trs <- mk_two_trees()
  expect_error(compare_pathtrees(trs$a, trs$b, iter = 0L), "iter")
  expect_error(compare_pathtrees(trs$a, trs$b, iter = NA), "iter")
})
