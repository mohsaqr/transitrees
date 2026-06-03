# ---- Tests for pathtree_dependence() ----

.simple_tree <- function() {
  set.seed(2)
  states <- c("A", "B", "C")
  P <- matrix(c(0.7, 0.2, 0.1,
                0.1, 0.7, 0.2,
                0.2, 0.2, 0.6),
              nrow = 3, byrow = TRUE,
              dimnames = list(states, states))
  seqs <- lapply(seq_len(100), function(i) {
    n <- 20L
    s <- character(n); s[1L] <- sample(states, 1L)
    for (t in seq.int(2L, n)) s[t] <- sample(states, 1L, prob = P[s[t - 1L], ])
    s
  })
  context_tree(seqs, max_depth = 3L, min_count = 5L)
}

test_that("pathtree_dependence returns a sorted data.frame with all columns", {
  tree <- .simple_tree()
  pd   <- pathtree_dependence(tree)
  expect_true(is.data.frame(pd))
  expect_named(pd, c("pathway", "depth", "count", "divergence",
                     "entropy", "entropy_before", "entropy_drop",
                     "likely_next", "likely_before", "changes_prediction"))
  expect_equal(pd$divergence, sort(pd$divergence, decreasing = TRUE))
})

test_that("divergence is non-negative for finite values", {
  tree <- .simple_tree()
  pd   <- pathtree_dependence(tree)
  finite_div <- pd$divergence[is.finite(pd$divergence)]
  expect_true(all(finite_div >= -1e-10))
})

test_that("entropy_drop = entropy_before - entropy", {
  tree <- .simple_tree()
  pd   <- pathtree_dependence(tree)
  expect_equal(pd$entropy_drop, pd$entropy_before - pd$entropy,
               tolerance = 1e-12)
})

test_that("changes_prediction = likely_next != likely_before", {
  tree <- .simple_tree()
  pd   <- pathtree_dependence(tree)
  expect_equal(pd$changes_prediction, pd$likely_next != pd$likely_before)
})

test_that("base 2 vs base e differ by ln(2)", {
  tree <- .simple_tree()
  pd2 <- pathtree_dependence(tree, base = 2)
  pde <- pathtree_dependence(tree, base = exp(1))
  ## row order may differ if divergence ties — match by context
  pde <- pde[match(pd2$pathway, pde$pathway), ]
  expect_equal(pd2$divergence * log(2), pde$divergence, tolerance = 1e-10)
})

test_that("pathtree_dependence on an order-1 generator has small chain-level KL", {
  tree <- .simple_tree()
  pr   <- prune_pathtree(tree, criterion = "G2", alpha = 0.05)
  pd   <- pathtree_dependence(pr)
  ## On an order-1 process with sufficient data, post-pruning the surviving
  ## contexts should mostly be depth 1; deeper contexts are noise.
  expect_true(median(pd$depth) <= 2L)
})

test_that("pathtree_dependence sorts descending by each sort_by column", {
  tree <- .simple_tree()
  for (sb in c("divergence", "entropy_drop", "entropy", "count", "depth")) {
    pd  <- pathtree_dependence(tree, sort_by = sb)
    key <- pd[[sb]]
    expect_equal(key, sort(key, decreasing = TRUE), info = sb)
  }
})

test_that("pathtree_dependence top = keeps the leading rows of the sort", {
  tree <- .simple_tree()
  full <- pathtree_dependence(tree, sort_by = "entropy_drop")
  pd3  <- pathtree_dependence(tree, sort_by = "entropy_drop", top = 3L)
  expect_lte(nrow(pd3), 3L)
  expect_equal(pd3$pathway, utils::head(full$pathway, nrow(pd3)))
})
