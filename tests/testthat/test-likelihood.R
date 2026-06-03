# Tests for Phase A: logLik / perplexity / score_*

mk_tree <- function(seed = 1L, n = 40L, len = 12L,
                    states = c("A","B","C"),
                    max_depth = 3L, min_count = 3L) {
  set.seed(seed)
  m <- matrix(sample(states, n * len, replace = TRUE), n, len)
  context_tree(m, max_depth = max_depth, min_count = min_count)
}

test_that("logLik.transitrees returns a logLik object with df and nobs", {
  tr <- mk_tree()
  ll <- logLik(tr)
  expect_s3_class(ll, "logLik")
  expect_true(is.finite(as.numeric(ll)))
  expect_true(!is.null(attr(ll, "nobs")))
  expect_true(!is.null(attr(ll, "df")))
  expect_identical(attr(ll, "df"),
                   length(tr$nodes) * (length(tr$alphabet) - 1L))
})

test_that("nobs.transitrees returns total tokens", {
  tr <- mk_tree()
  expect_identical(nobs(tr), as.integer(tr$n_obs))
})

test_that("AIC and BIC compute via logLik attributes", {
  tr <- mk_tree()
  ll <- logLik(tr)
  expect_identical(AIC(tr), -2 * as.numeric(ll) + 2 * attr(ll, "df"))
  expect_identical(
    BIC(tr),
    -2 * as.numeric(ll) + log(attr(ll, "nobs")) * attr(ll, "df"))
})

test_that("in-sample logLik attributes nobs == n_obs (every token scored)", {
  tr <- mk_tree()
  expect_identical(attr(logLik(tr), "nobs"), as.integer(tr$n_obs))
})

test_that("in-sample logLik <= 0 (probabilities <= 1 ⇒ log <= 0)", {
  tr <- mk_tree()
  expect_lte(as.numeric(logLik(tr)), 0)
})

test_that("perplexity equals exp(-ll/n)", {
  tr <- mk_tree()
  ll <- logLik(tr); n <- attr(ll, "nobs")
  expect_equal(perplexity(tr), exp(-as.numeric(ll) / n))
})

test_that("perplexity of unsmoothed tree on i.i.d. data is near alphabet size", {
  set.seed(42)
  states <- c("A","B","C","D")
  m <- matrix(sample(states, 60 * 20, replace = TRUE), 60, 20)
  ## max_depth 0 -> only root (marginal); ymin = 0 -> exact MLE.
  tr <- context_tree(m, max_depth = 0L, min_count = 1L,
                     smoothing = list("floor", ymin = 0))
  pp <- perplexity(tr)
  expect_gt(pp, 3.5); expect_lt(pp, 4.5)
})

test_that("held-out perplexity >= in-sample perplexity on same data", {
  ## Tree fits its training data better than a fresh resample of the
  ## same generating distribution.
  set.seed(7)
  states <- c("A","B","C")
  train <- matrix(sample(states, 50 * 12, replace = TRUE), 50, 12)
  test  <- matrix(sample(states, 50 * 12, replace = TRUE), 50, 12)
  tr <- context_tree(train, max_depth = 2L, min_count = 3L)
  pp_in  <- perplexity(tr)
  pp_out <- perplexity(tr, newdata = test)
  expect_gte(pp_out, pp_in - 1e-8)
})

test_that("held-out logLik attributes match score_positions row count", {
  tr <- mk_tree()
  set.seed(123)
  test <- matrix(sample(tr$alphabet, 30 * 10, replace = TRUE), 30, 10)
  ll <- logLik(tr, newdata = test)
  pos <- score_positions(tr, newdata = test)
  expect_identical(attr(ll, "nobs"), as.integer(nrow(pos)))
  expect_equal(as.numeric(ll), sum(pos$log_lik))
})

test_that("score_sequences logLik sums match overall held-out logLik", {
  tr <- mk_tree()
  set.seed(11)
  test <- matrix(sample(tr$alphabet, 20 * 10, replace = TRUE), 20, 10)
  seq_tab <- score_sequences(tr, test)
  ll <- as.numeric(logLik(tr, newdata = test))
  expect_equal(sum(seq_tab$log_lik), ll)
  expect_true(all(seq_tab$perplexity > 0))
})

test_that("score_positions returns expected columns", {
  tr <- mk_tree()
  set.seed(99)
  test <- matrix(sample(tr$alphabet, 4 * 5, replace = TRUE), 4, 5)
  pos <- score_positions(tr, test)
  expect_named(pos, c("sequence_id", "position", "matched_context",
                      "observed", "predicted_prob", "log_lik"))
  expect_true(all(pos$predicted_prob > 0 & pos$predicted_prob <= 1))
})

test_that("logLik on out-of-vocab held-out data skips unscorable tokens", {
  tr <- mk_tree()
  test <- matrix(c("A","B","C","Z","C","A"), 2, 3)  # "Z" not in alphabet
  expect_silent(ll <- logLik(tr, newdata = test))
  expect_lt(attr(ll, "nobs"), length(test))
})
