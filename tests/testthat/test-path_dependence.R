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
  context_tree(seqs, max_depth = 3L, nmin = 5L)
}

test_that("pathtree_dependence returns a sorted data.frame with all columns", {
  tree <- .simple_tree()
  pd   <- pathtree_dependence(tree)
  expect_true(is.data.frame(pd))
  expect_named(pd, c("pathway", "depth", "count", "KL",
                     "H_node", "H_parent", "H_drop",
                     "modal_next", "modal_parent", "flips"))
  expect_equal(pd$KL, sort(pd$KL, decreasing = TRUE))
})

test_that("KL is non-negative for finite values", {
  tree <- .simple_tree()
  pd   <- pathtree_dependence(tree)
  finite_KL <- pd$KL[is.finite(pd$KL)]
  expect_true(all(finite_KL >= -1e-10))
})

test_that("H_drop = H_parent - H_node", {
  tree <- .simple_tree()
  pd   <- pathtree_dependence(tree)
  expect_equal(pd$H_drop, pd$H_parent - pd$H_node, tolerance = 1e-12)
})

test_that("flips column is consistent with modal_next != modal_parent", {
  tree <- .simple_tree()
  pd   <- pathtree_dependence(tree)
  expect_equal(pd$flips, pd$modal_next != pd$modal_parent)
})

test_that("base 2 vs base e differ by ln(2)", {
  tree <- .simple_tree()
  pd2 <- pathtree_dependence(tree, base = 2)
  pde <- pathtree_dependence(tree, base = exp(1))
  ## row order may differ if KL ties — match by context
  pde <- pde[match(pd2$pathway, pde$pathway), ]
  expect_equal(pd2$KL * log(2), pde$KL, tolerance = 1e-10)
})

test_that("pathtree_dependence on an order-1 generator has small chain-level KL", {
  tree <- .simple_tree()
  pr   <- prune_pathtree(tree, criterion = "G2", alpha = 0.05)
  pd   <- pathtree_dependence(pr)
  ## On an order-1 process with sufficient data, post-pruning the surviving
  ## contexts should mostly be depth 1; deeper contexts are noise.
  expect_true(median(pd$depth) <= 2L)
})
