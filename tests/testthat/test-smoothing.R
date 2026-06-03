# Tests for Phase B: smoothing variants and re-smoothing

mk_seq_data <- function(n = 40, len = 12, seed = 1,
                        states = c("A","B","C","D")) {
  set.seed(seed)
  matrix(sample(states, n * len, replace = TRUE), n, len)
}

test_that("Laplace smoothing returns a probability vector summing to 1", {
  counts <- c(3, 0, 1, 0)
  p <- pathtree:::.smooth_laplace(counts, alpha = 1)
  expect_equal(sum(p), 1)
  expect_true(all(p > 0))
})

test_that("Laplace with alpha = 0 reduces to MLE", {
  counts <- c(3, 1, 4)
  expect_equal(pathtree:::.smooth_laplace(counts, alpha = 0),
               counts / sum(counts))
})

test_that("Kneser-Ney returns a probability vector summing to 1", {
  counts <- c(5, 0, 2, 0)
  parent <- c(0.25, 0.25, 0.25, 0.25)
  p <- pathtree:::.smooth_kneser_ney(counts, parent, discount = 0.75)
  expect_equal(sum(p), 1)
  expect_true(all(p >= 0))
})

test_that("Kneser-Ney with discount 0 reduces to MLE", {
  counts <- c(5, 1, 2, 0)
  parent <- c(0.25, 0.25, 0.25, 0.25)
  p <- pathtree:::.smooth_kneser_ney(counts, parent, discount = 0)
  expect_equal(p, counts / sum(counts))
})

test_that("Witten-Bell sums to 1 and reduces to parent when no obs", {
  parent <- c(0.5, 0.3, 0.2)
  expect_equal(pathtree:::.smooth_witten_bell(c(0,0,0), parent), parent)
  p <- pathtree:::.smooth_witten_bell(c(2, 1, 0), parent)
  expect_equal(sum(p), 1)
})

test_that("Jelinek-Mercer mixes MLE and parent", {
  parent <- c(0.5, 0.3, 0.2)
  counts <- c(3, 0, 1)
  p <- pathtree:::.smooth_jelinek_mercer(counts, parent, lambda = 0.5)
  mle <- counts / sum(counts)
  expect_equal(p, 0.5 * mle + 0.5 * parent)
  expect_equal(sum(p), 1)
})

test_that("context_tree() accepts all smoothing methods and sums to 1", {
  m <- mk_seq_data()
  for (s in c("floor", "laplace", "kneser_ney", "witten_bell",
              "jelinek_mercer")) {
    tr <- context_tree(m, max_depth = 2L, min_count = 3L, smoothing = s)
    expect_s3_class(tr, "pathtree")
    sums <- vapply(tr$nodes, function(n) sum(n$prob), numeric(1))
    expect_true(all(abs(sums - 1) < 1e-9),
                info = paste("smoothing =", s))
    expect_identical(tr$smoothing$method, s)
  }
})

test_that("invalid smoothing hyperparameters error early", {
  m <- mk_seq_data()
  expect_error(context_tree(m, smoothing = list("floor", ymin = -0.1)),
               "ymin")
  expect_error(context_tree(m, smoothing = list("laplace", alpha = -1)),
               "alpha")
  expect_error(context_tree(m, smoothing = list("kneser_ney",
                                                discount = 1.5)),
               "discount")
  expect_error(context_tree(m, smoothing = list("jelinek_mercer",
                                                lambda = 1.5)),
               "lambda")
  expect_error(smooth_pathtree(context_tree(m, max_depth = 1L),
                               list("jelinek_mercer", lambda = -0.1)),
               "lambda")
})

test_that("valid boundary smoothing hyperparameters are accepted", {
  m <- mk_seq_data()
  expect_s3_class(context_tree(m, smoothing = list("floor", ymin = 0)),
                  "pathtree")
  expect_s3_class(context_tree(m, smoothing = list("laplace", alpha = 0)),
                  "pathtree")
  expect_s3_class(context_tree(m, smoothing = list("kneser_ney",
                                                   discount = 1)),
                  "pathtree")
  expect_s3_class(context_tree(m, smoothing = list("jelinek_mercer",
                                                   lambda = 0)),
                  "pathtree")
  expect_s3_class(context_tree(m, smoothing = list("jelinek_mercer",
                                                   lambda = 1)),
                  "pathtree")
})

test_that("smoothing = 'floor' default reproduces v0.0 behaviour", {
  ## Equivalence to ymin floor: tree built with default args matches
  ## tree built with explicit smoothing = 'floor'.
  m <- mk_seq_data()
  a <- context_tree(m, max_depth = 2L, min_count = 3L)
  b <- context_tree(m, max_depth = 2L, min_count = 3L, smoothing = "floor")
  for (ctx in names(a$nodes))
    expect_equal(a$nodes[[ctx]]$prob, b$nodes[[ctx]]$prob)
})

test_that("smooth_pathtree() preserves topology, swaps probs", {
  m <- mk_seq_data()
  tr  <- context_tree(m, max_depth = 2L, min_count = 3L)
  tr2 <- smooth_pathtree(tr, "laplace")
  expect_identical(names(tr$nodes), names(tr2$nodes))
  expect_identical(tr$edges, tr2$edges)
  expect_identical(tr2$smoothing$method, "laplace")
  for (ctx in names(tr$nodes))
    expect_identical(tr$nodes[[ctx]]$counts, tr2$nodes[[ctx]]$counts)
})

test_that("smooth_pathtree() with floor reproduces tree's original probs", {
  m <- mk_seq_data()
  tr  <- context_tree(m, max_depth = 2L, min_count = 3L,
                      smoothing = list("floor", ymin = 0.001))
  tr2 <- smooth_pathtree(tr, list("floor", ymin = 0.001))
  for (ctx in names(tr$nodes))
    expect_equal(tr$nodes[[ctx]]$prob, tr2$nodes[[ctx]]$prob)
})

test_that("Kneser-Ney never produces an infinite held-out perplexity", {
  set.seed(2025)
  states <- c("A","B","C","D","E")
  train  <- matrix(sample(states, 30 * 8,  replace = TRUE), 30, 8)
  test   <- matrix(sample(states, 30 * 8,  replace = TRUE), 30, 8)
  tr_kn  <- context_tree(train, max_depth = 3L, min_count = 1L,
                         smoothing = "kneser_ney")
  pp_kn <- perplexity(tr_kn, newdata = test)
  expect_true(is.finite(pp_kn))
  expect_gt(pp_kn, 1)
})

test_that("smooth_pathtree() walks top-down (parent updated before child)", {
  m <- mk_seq_data(seed = 3)
  tr <- context_tree(m, max_depth = 2L, min_count = 3L)
  tr2 <- smooth_pathtree(tr, list("jelinek_mercer", lambda = 1))
  ## With lambda = 1, every node's prob should equal its parent's
  ## smoothed prob. The root has no parent (uniform fallback).
  k <- length(tr$alphabet)
  for (ctx in names(tr2$nodes)) {
    if (identical(ctx, pathtree:::.ROOT)) {
      expect_equal(tr2$nodes[[ctx]]$prob, rep(1 / k, k))
    } else {
      par <- pathtree:::.pt_parent_ctx(ctx)
      expect_equal(tr2$nodes[[ctx]]$prob,
                   tr2$nodes[[par]]$prob)
    }
  }
})
