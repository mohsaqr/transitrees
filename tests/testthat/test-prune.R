# ---- Tests for prune_tree() ----

.medium_seqs <- function(seed = 1) {
  set.seed(seed)
  ## 80 sequences from a noisy order-1 process — pruning should be aggressive
  states <- c("A", "B", "C")
  P <- matrix(c(0.7, 0.2, 0.1,
                0.1, 0.7, 0.2,
                0.2, 0.2, 0.6),
              nrow = 3, byrow = TRUE,
              dimnames = list(states, states))
  lapply(seq_len(80), function(i) {
    seq_len_i <- 20L
    s <- character(seq_len_i)
    s[1L] <- sample(states, 1L)
    for (t in seq.int(2L, seq_len_i)) {
      s[t] <- sample(states, 1L, prob = P[s[t - 1L], ])
    }
    s
  })
}

test_that("prune_tree marks the tree as pruned", {
  tree <- context_tree(.medium_seqs(), max_depth = 3L, min_count = 5L)
  pr   <- prune_tree(tree, criterion = "G2", alpha = 0.05)
  expect_true(isTRUE(pr$pruned))
  expect_equal(pr$pruning$criterion, "G2")
})

test_that("pruning never exceeds the unpruned node count", {
  tree <- context_tree(.medium_seqs(), max_depth = 3L, min_count = 5L)
  for (crit in c("G2", "KL", "AIC", "BIC")) {
    pr <- prune_tree(tree, criterion = crit)
    expect_lte(length(pr$nodes), length(tree$nodes), label = crit)
  }
})

test_that("pruning preserves the root", {
  tree <- context_tree(.medium_seqs(), max_depth = 3L, min_count = 5L)
  pr   <- prune_tree(tree, criterion = "G2", alpha = 0.001)
  expect_true("<root>" %in% names(pr$nodes))
})

test_that("aggressive pruning of an order-1 process collapses to depth 1", {
  tree <- context_tree(.medium_seqs(), max_depth = 3L, min_count = 5L)
  ## alpha = 0.001 -> critical value = qchisq(0.999, 2) ~ 13.8, very strict
  pr   <- prune_tree(tree, criterion = "G2", alpha = 0.001)
  depths <- vapply(pr$nodes, function(x) x$depth, integer(1))
  ## On a pure order-1 generator, depth >= 2 contexts carry no genuine
  ## signal; strict pruning should leave the tree dominated by depths 0-1.
  expect_true(mean(depths <= 1L) > 0.5,
              info = sprintf("depths kept: %s",
                              paste(table(depths), collapse = ",")))
})

test_that("invalid criterion errors", {
  tree <- context_tree(.medium_seqs(), max_depth = 2L, min_count = 5L)
  expect_error(prune_tree(tree, criterion = "invalid"),
               regexp = "should be one of")
})

test_that("lowercase criterion is rejected (case-sensitive match.arg)", {
  tree <- context_tree(.medium_seqs(), max_depth = 2L, min_count = 5L)
  expect_error(prune_tree(tree, criterion = "g2"),
               regexp = "should be one of")
})
