# CRAN gap-filling coverage: S3 methods and error paths not exercised by
# the existing suite. Base R + testthat only; small deterministic inputs.

.cov_tree <- function(seed = 1L, n = 40L, len = 10L,
                      states = c("A", "B", "C"),
                      max_depth = 2L, min_count = 3L) {
  set.seed(seed)
  m <- matrix(sample(states, n * len, replace = TRUE), n, len)
  context_tree(m, max_depth = max_depth, min_count = min_count)
}

# ---- bootstrap S3: print + summary (previously untested) ----------------

test_that("print.transitiontrees_bootstrap shows the header and returns invisibly", {
  b <- bootstrap_pathways(.cov_tree(), iter = 20L, seed = 1L)
  expect_output(print(b), "<transitiontrees_bootstrap>", fixed = TRUE)
  expect_output(print(b), "resamples")
  expect_output(print(b), "pathways")
  expect_invisible(print(b))
})

test_that("print.transitiontrees_bootstrap honours n and notes truncation", {
  b <- bootstrap_pathways(.cov_tree(), iter = 20L, seed = 1L)
  ## With n = 1 there should be a "more pathways" footer (the tree has
  ## many contexts).
  expect_output(print(b, n = 1L), "more pathways")
})

test_that("summary.transitiontrees_bootstrap returns the summary table verbatim", {
  b <- bootstrap_pathways(.cov_tree(), iter = 20L, seed = 1L)
  s <- summary(b)
  expect_s3_class(s, "data.frame")
  expect_identical(s, b$summary)
  ## summary() and as.data.frame() are interchangeable extractors.
  expect_identical(s, as.data.frame(b))
})

# ---- mine_sequences error / edge paths ----------------------------------

test_that("mine_sequences rejects an unknown 'which'", {
  tree <- .cov_tree()
  set.seed(2)
  new <- replicate(5, sample(c("A", "B", "C"), 8, replace = TRUE),
                   simplify = FALSE)
  expect_error(mine_sequences(tree, new, which = "weird"),
               "should be one of")
})

test_that("mine_sequences rejects a non-transitiontrees first argument", {
  expect_error(mine_sequences(list(), list(c("A", "B"))),
               "inherits")
})

test_that("mine_sequences on empty newdata returns the empty score schema", {
  tree <- .cov_tree()
  ## All tokens out of vocabulary -> nothing scorable -> empty tidy frame.
  out <- mine_sequences(tree, list(c("Z", "Z", "Z")))
  expect_s3_class(out, "data.frame")
  expect_named(out, c("sequence_id", "n_scored", "log_lik", "perplexity"))
})

# ---- score_sequences / score_positions error + edge paths ---------------

test_that("score_sequences errors on input with no usable sequences", {
  tree <- .cov_tree()
  expect_error(score_sequences(tree, list()),
               "No usable held-out sequences")
})

test_that("score_sequences rejects a non-transitiontrees first argument", {
  expect_error(score_sequences(list(), list(c("A", "B"))), "inherits")
})

test_that("score_positions(worst=) returns the least-probable positions in order", {
  tree <- .cov_tree()
  set.seed(3)
  new <- replicate(8, sample(c("A", "B", "C"), 8, replace = TRUE),
                   simplify = FALSE)
  full <- score_positions(tree, new)
  w3   <- score_positions(tree, new, worst = 3L)
  expect_lte(nrow(w3), 3L)
  ## ascending predicted_prob (most surprising first)
  expect_equal(w3$predicted_prob, sort(w3$predicted_prob))
  ## the worst rows are the global minima of predicted_prob
  expect_equal(w3$predicted_prob,
               utils::head(sort(full$predicted_prob), nrow(w3)))
})

# ---- query_pathway error path -------------------------------------------

test_that("query_pathway errors when next_state is outside the alphabet", {
  tree <- .cov_tree()
  expect_error(query_pathway(tree, "A", next_state = "Z"),
               "not in the tree's alphabet")
})

# ---- unknown smoothing method strings ------------------------------------

test_that("context_tree rejects an unknown smoothing method name", {
  set.seed(4)
  m <- matrix(sample(c("A", "B", "C"), 30 * 8, replace = TRUE), 30, 8)
  expect_error(context_tree(m, max_depth = 1L, smoothing = "bogus"),
               "Unknown smoothing method|one of")
})

test_that("smooth_tree rejects an unknown smoothing method name", {
  tree <- .cov_tree()
  expect_error(smooth_tree(tree, "bogus"),
               "Unknown smoothing method|one of")
})

# ---- subtree / n_nodes / tree_distance argument guards -------------------

test_that("subtree rejects a non-transitiontrees object", {
  expect_error(subtree(list(), "A"))
})

test_that("tree_distance is self-consistent for a pruned tree", {
  tree <- .cov_tree(max_depth = 3L, min_count = 2L)
  pr   <- prune_tree(tree, criterion = "G2", alpha = 0.05)
  expect_equal(tree_distance(pr, pr), 0)
})
