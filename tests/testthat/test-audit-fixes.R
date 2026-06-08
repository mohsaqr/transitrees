# Regression tests for the read-only audit findings.

mk_seqs <- function(seed = 1L, n = 40L) {
  set.seed(seed)
  replicate(n, sample(c("A", "B", "C"), 8, replace = TRUE), simplify = FALSE)
}

# ---- Weights: stored, consistent, and not silently dropped ----

test_that("weighted fits store their weights; uniform weights normalise away", {
  seqs <- mk_seqs()
  trw <- context_tree(seqs, max_depth = 2L, min_count = 2L,
                      weights = c(rep(5, 20), rep(1, 20)))
  expect_false(is.null(trw$weights))
  # uniform weights are equivalent to unweighted -> stored as NULL
  expect_null(context_tree(seqs, max_depth = 2L, weights = rep(2, 40))$weights)
})

test_that("nobs() agrees with logLik()'s nobs under weights (audit #2)", {
  seqs <- mk_seqs()
  trw <- context_tree(seqs, max_depth = 2L, min_count = 2L,
                      weights = c(rep(10, 20), rep(1, 20)))
  expect_identical(as.integer(nobs(trw)), attr(logLik(trw), "nobs"))
  expect_identical(model_fit(trw)$nobs, as.integer(nobs(trw)))
})

test_that("resampling/permutation surfaces reject weighted trees (audit #1)", {
  seqs <- mk_seqs()
  trw  <- context_tree(seqs, max_depth = 2L, min_count = 2L,
                       weights = c(rep(5, 20), rep(1, 20)))
  trw2 <- context_tree(mk_seqs(2), max_depth = 2L, min_count = 2L,
                       weights = c(rep(3, 20), rep(1, 20)))
  expect_error(bootstrap_pathways(trw, iter = 10L), "weighted")
  expect_error(compare_trees(trw, trw2, iter = 5L), "weighted")
  grp <- structure(list(x = trw, y = trw2),
                   class = c("transitiontrees_group", "list"),
                   group = "g")
  expect_error(compare_groups(grp, iter = 10L), "weighted")
})

test_that("all-zero weights are rejected (audit #5)", {
  expect_error(context_tree(mk_seqs(), weights = rep(0, 40)), "positive")
})

# ---- Argument validation (audit #3, #6) ----

test_that("prune_tree validates alpha (audit #3)", {
  tree <- context_tree(mk_seqs(), max_depth = 2L)
  expect_error(prune_tree(tree, alpha = 2),   "alpha")
  expect_error(prune_tree(tree, alpha = 0),   "alpha")
  expect_error(prune_tree(tree, alpha = -0.1), "alpha")
})

test_that(".suffix_chain validates alpha (audit #3)", {
  tree <- context_tree(mk_seqs(), max_depth = 3L, min_count = 3L)
  expect_error(transitiontrees:::.suffix_chain(tree, "A -> B", alpha = 2),
               "alpha")
})

test_that("score_positions validates worst (audit #6)", {
  tree <- context_tree(mk_seqs(), max_depth = 2L)
  expect_error(score_positions(tree, mk_seqs(2), worst = -1), "worst")
  expect_error(score_positions(tree, mk_seqs(2), worst = 0),  "worst")
})

# ---- tune_tree never reports a failed config as best (audit #4) ----

test_that("tune_tree excludes failed configs from best and warns", {
  set.seed(3)
  seqs <- replicate(40, sample(c("A", "B"), 6, replace = TRUE),
                    simplify = FALSE)
  # ymin = 0.6 on a binary alphabet makes the floor smoother invalid, so
  # every fold of that config errors; the valid 'laplace' config must win.
  expect_warning(
    tt <- tune_tree(seqs, max_depth = 1:2, min_count = 2L,
                    smoothing = list(list("laplace"),
                                     list("floor", ymin = 0.6)),
                    folds = 3L),
    "fold")
  best <- attr(tt, "best")
  expect_true("folds_failed" %in% names(tt))
  expect_false(is.null(best))                  # a valid config exists
  expect_identical(best$folds_failed, 0L)      # best never has failed folds
  # the invalid floor config is present but not selected
  expect_true(any(tt$folds_failed > 0L))
})

test_that("tune_tree returns best = NULL when every config fails", {
  set.seed(4)
  seqs <- replicate(30, sample(c("A", "B"), 6, replace = TRUE),
                    simplify = FALSE)
  ws <- testthat::capture_warnings(
    tt <- tune_tree(seqs, max_depth = 1L, min_count = 2L,
                    smoothing = list(list("floor", ymin = 0.6)), folds = 3L))
  expect_true(any(grepl("fold", ws)))
  expect_true(any(grepl("best.*NULL|NULL", ws)))
  expect_null(attr(tt, "best"))
})
