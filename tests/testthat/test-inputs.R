# Tests for Phase F: stslist input + per-sequence weights

mk_mat <- function(n = 30, len = 12, seed = 1,
                   states = c("A","B","C")) {
  set.seed(seed)
  matrix(sample(states, n * len, replace = TRUE), n, len)
}

test_that("weights = rep(1, n) gives byte-identical fit to unweighted", {
  m <- mk_mat()
  tr1 <- context_tree(m, max_depth = 2L, nmin = 3L)
  tr2 <- context_tree(m, max_depth = 2L, nmin = 3L,
                      weights = rep(1, nrow(m)))
  ## Unweighted is integer; weighted is numeric — compare values.
  expect_equal(names(tr1$nodes), names(tr2$nodes))
  for (ctx in names(tr1$nodes)) {
    expect_equal(as.numeric(tr1$nodes[[ctx]]$counts),
                 as.numeric(tr2$nodes[[ctx]]$counts))
    expect_equal(tr1$nodes[[ctx]]$prob,
                 tr2$nodes[[ctx]]$prob)
  }
})

test_that("weighted fit equals unweighted fit on row-duplicated data", {
  ## Duplicating a row twice should be equivalent to weighting that row 3.
  set.seed(11)
  m <- mk_mat()
  dup_idx <- c(seq_len(nrow(m)), 1L, 1L)   # row 1 appears 3x in total
  m_dup <- m[dup_idx, , drop = FALSE]
  tr_dup <- context_tree(m_dup, max_depth = 2L, nmin = 1L)

  w <- rep(1, nrow(m)); w[1L] <- 3
  tr_w   <- context_tree(m, max_depth = 2L, nmin = 1L, weights = w)

  expect_equal(sort(names(tr_dup$nodes)), sort(names(tr_w$nodes)))
  for (ctx in names(tr_dup$nodes)) {
    expect_equal(as.numeric(tr_dup$nodes[[ctx]]$counts),
                 as.numeric(tr_w$nodes[[ctx]]$counts))
  }
})

test_that("weights validation errors on mismatched length", {
  m <- mk_mat()
  expect_error(context_tree(m, max_depth = 1L, nmin = 1L,
                            weights = c(1, 2, 3)),
               "length equal to number of input")
})

test_that("weights validation errors on negative values", {
  m <- mk_mat()
  w <- rep(1, nrow(m)); w[1L] <- -1
  expect_error(context_tree(m, max_depth = 1L, nmin = 1L, weights = w),
               "non-negative")
})

test_that("stslist input matches as.matrix() input on the same data", {
  skip_if_not_installed("TraMineR")
  m <- mk_mat()
  ## TraMineR seqdef warns on default; suppress for the test.
  sts <- suppressMessages(suppressWarnings(
    TraMineR::seqdef(m)
  ))
  tr_mat <- context_tree(m,   max_depth = 2L, nmin = 3L)
  tr_sts <- context_tree(sts, max_depth = 2L, nmin = 3L)
  expect_equal(sort(tr_mat$alphabet), sort(tr_sts$alphabet))
  expect_equal(sort(names(tr_mat$nodes)), sort(names(tr_sts$nodes)))
  for (ctx in names(tr_mat$nodes)) {
    expect_equal(as.numeric(tr_mat$nodes[[ctx]]$counts),
                 as.numeric(tr_sts$nodes[[ctx]]$counts))
  }
})

test_that("stslist with weights attribute is auto-detected", {
  skip_if_not_installed("TraMineR")
  m <- mk_mat()
  w <- runif(nrow(m), 0.5, 2)
  sts <- suppressMessages(suppressWarnings(
    TraMineR::seqdef(m, weights = w)
  ))
  tr_auto <- context_tree(sts, max_depth = 2L, nmin = 1L)
  tr_explicit <- context_tree(m, max_depth = 2L, nmin = 1L, weights = w)
  for (ctx in names(tr_auto$nodes)) {
    expect_equal(as.numeric(tr_auto$nodes[[ctx]]$counts),
                 as.numeric(tr_explicit$nodes[[ctx]]$counts))
  }
})
