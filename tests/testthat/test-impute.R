# Tests for impute_sequences()

mk_imp_tree <- function(seed = 1L) {
  set.seed(seed)
  seqs <- replicate(80, sample(c("A", "B", "C"), 10, replace = TRUE),
                    simplify = FALSE)
  context_tree(seqs, max_depth = 2L)
}

test_that("impute_sequences fills internal gaps and preserves shape (list)", {
  tree <- mk_imp_tree()
  gappy <- list(c("A", NA, "C"), c("B", "B", NA, "A"))
  out <- impute_sequences(tree, gappy)
  expect_type(out, "list")
  expect_length(out, 2L)
  expect_false(anyNA(out[[1]]))
  expect_false(anyNA(out[[2]]))
  expect_true(all(unlist(out) %in% tree$alphabet))
})

test_that("impute_sequences leaves trailing padding untouched", {
  tree <- mk_imp_tree()
  # internal gap at pos 2, trailing padding at pos 4-5
  m <- matrix(c("A", NA, "B", NA, NA), nrow = 1)
  out <- impute_sequences(tree, m)
  expect_false(is.na(out[1, 2]))      # internal gap filled
  expect_true(is.na(out[1, 4]))       # trailing padding left as NA
  expect_true(is.na(out[1, 5]))
  expect_true(is.matrix(out))
  expect_equal(dim(out), dim(m))
})

test_that("impute_sequences method = 'modal' is deterministic", {
  tree <- mk_imp_tree()
  gappy <- list(c("A", NA, "C"))
  a <- impute_sequences(tree, gappy, method = "modal")
  b <- impute_sequences(tree, gappy, method = "modal")
  expect_identical(a, b)
})

test_that("impute_sequences method = 'prob' is reproducible with a seed", {
  tree <- mk_imp_tree()
  gappy <- list(c("A", NA, NA, "C"), c("B", NA, "A"))
  a <- impute_sequences(tree, gappy, method = "prob", seed = 42L)
  b <- impute_sequences(tree, gappy, method = "prob", seed = 42L)
  expect_identical(a, b)
})

test_that("impute_sequences leaves fully observed sequences unchanged", {
  tree <- mk_imp_tree()
  full <- list(c("A", "B", "C"))
  expect_identical(impute_sequences(tree, full), list(c("A", "B", "C")))
})

test_that("impute_sequences handles a data.frame and a bare vector", {
  tree <- mk_imp_tree()
  df <- as.data.frame(matrix(c("A", NA, "C", "B", "A", NA), nrow = 2,
                             byrow = TRUE), stringsAsFactors = FALSE)
  out_df <- impute_sequences(tree, df)
  expect_s3_class(out_df, "data.frame")
  expect_false(is.na(out_df[1, 2]))   # internal gap row 1
  expect_true(is.na(out_df[2, 3]))    # trailing padding row 2

  out_vec <- impute_sequences(tree, c("A", NA, "B"))
  expect_false(anyNA(out_vec))
  expect_length(out_vec, 3L)
})

test_that("impute_sequences returns an all-missing sequence unchanged", {
  tree <- mk_imp_tree()
  out <- impute_sequences(tree, list(c(NA, NA, NA)))
  expect_true(all(is.na(out[[1]])))
})

test_that("impute_sequences(method='prob', seed=) restores the caller's RNG (#6)", {
  tree <- mk_imp_tree()
  set.seed(99); before <- runif(1)
  set.seed(99)
  invisible(impute_sequences(tree, list(c("A", NA, "C")),
                             method = "prob", seed = 7L))
  after <- runif(1)               # RNG stream must be untouched by the seed= call
  expect_equal(before, after)
})
