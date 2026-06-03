library(testthat)
library(pathtree)

test_that("compare_pathtrees detects a real difference between different generators", {
  ## Create two different Markov chains
  set.seed(42)
  states <- c("A", "B")
  
  # Group A: A -> A (0.8), B -> B (0.8)
  trajs_a <- replicate(50, {
    s <- character(20)
    s[1] <- sample(states, 1)
    for (t in 2:20) {
      p <- if (s[t-1] == "A") c(0.8, 0.2) else c(0.2, 0.8)
      s[t] <- sample(states, 1, prob = p)
    }
    s
  }, simplify = FALSE)
  
  # Group B: A -> B (0.8), B -> A (0.8)
  trajs_b <- replicate(50, {
    s <- character(20)
    s[1] <- sample(states, 1)
    for (t in 2:20) {
      p <- if (s[t-1] == "A") c(0.2, 0.8) else c(0.8, 0.2)
      s[t] <- sample(states, 1, prob = p)
    }
    s
  }, simplify = FALSE)
  
  tr_a <- context_tree(trajs_a, max_depth = 1, min_count = 2)
  tr_b <- context_tree(trajs_b, max_depth = 1, min_count = 2)
  
  # Compare them
  cmp <- compare_pathtrees(tr_a, tr_b, iter = 100, seed = 1)
  
  # Expect significant p-value
  expect_lt(cmp$p_value, 0.05)
  expect_gt(cmp$pdist, 0.3) # Should be substantial distance
})

test_that("bootstrap_pathways correctly identifies non-informative pathways in noisy data", {
  ## Create data where A -> B and A -> C are equally likely (random)
  ## But by chance in a small sample, one might look slightly more likely.
  ## The bootstrap should show it's NOT informative.
  set.seed(123)
  trajs <- replicate(30, {
    s <- character(10)
    s[1] <- "A"
    for (t in 2:10) {
      s[t] <- sample(c("A", "B", "C"), 1)
    }
    s
  }, simplify = FALSE)
  
  tr <- context_tree(trajs, max_depth = 1, min_count = 2)
  b  <- bootstrap_pathways(tr, iter = 100, seed = 1)
  
  # Pathways should mostly NOT be informative
  # (Root is usually informative if alphabet is biased, but here it's uniform)
  # Check non-root pathways
  s_sub <- b$summary[b$summary$pathway != "(start)", ]
  if (nrow(s_sub) > 0) {
    expect_true(all(s_sub$informative_rate < 0.95))
    expect_true(all(!s_sub$informative))
  }
})

test_that("bootstrap_pathways handles very short sequences gracefully", {
  # Sequences shorter than max_depth
  trajs <- list(c("A"), c("B"), c("A", "B"))
  tr <- context_tree(trajs, max_depth = 2, min_count = 1)
  
  # This should run without error
  expect_error(b <- bootstrap_pathways(tr, iter = 20, seed = 1), NA)
  expect_s3_class(b, "pathtree_bootstrap")
})
