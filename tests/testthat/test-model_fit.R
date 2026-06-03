# ---- Tests for model_fit() ----

.mf_tree <- function() {
  set.seed(1)
  context_tree(replicate(60, sample(c("A","B","C"), 12, replace = TRUE),
                         simplify = FALSE),
               max_depth = 2L, min_count = 3L)
}

test_that("model_fit returns one tidy row with the six scalars", {
  res <- model_fit(.mf_tree())
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1L)
  expect_named(res, c("logLik", "df", "nobs", "AIC", "BIC", "perplexity"))
})

test_that("model_fit reproduces the standard generics exactly", {
  tree <- .mf_tree()
  res  <- model_fit(tree)
  expect_equal(res$logLik,     as.numeric(logLik(tree)))
  expect_equal(res$nobs,       as.integer(nobs(tree)))
  expect_equal(res$AIC,        AIC(tree))
  expect_equal(res$BIC,        BIC(tree))
  expect_equal(res$perplexity, perplexity(tree))
})

test_that("model_fit(newdata) is out-of-sample", {
  tree <- .mf_tree()
  test <- replicate(10, sample(c("A","B","C"), 12, replace = TRUE),
                    simplify = FALSE)
  res  <- model_fit(tree, newdata = test)
  expect_equal(res$perplexity, perplexity(tree, newdata = test))
  expect_equal(res$logLik, as.numeric(logLik(tree, newdata = test)))
})

test_that("model_fit on a group returns one row per group, tagged", {
  wide <- data.frame(
    T1 = c("A","B","A","B","A","B","A","B"),
    T2 = c("B","A","B","A","B","A","B","A"),
    T3 = c("A","A","B","B","A","B","B","A"),
    stringsAsFactors = FALSE)
  g   <- context_tree(wide, max_depth = 1L, min_count = 1L,
                      group = rep(c("x","y"), each = 4L))
  res <- model_fit(g)
  expect_true("group" %in% names(res))
  expect_setequal(res$group, c("x", "y"))
  expect_equal(nrow(res), 2L)
})

# ---- n_nodes() ----

test_that("n_nodes counts contexts and matches length(tree$nodes)", {
  tr <- .mf_tree()
  expect_equal(n_nodes(tr), length(tr$nodes))
  expect_type(n_nodes(tr), "integer")
})

test_that("n_nodes on a group returns one named count per group", {
  wide <- data.frame(T1 = c("A","B","A","B"), T2 = c("B","A","B","A"),
                     stringsAsFactors = FALSE)
  g <- context_tree(wide, max_depth = 1L, min_count = 1L,
                    group = c("x","x","y","y"))
  nn <- n_nodes(g)
  expect_named(nn, c("x","y"))
  expect_equal(unname(nn), c(length(g$x$nodes), length(g$y$nodes)))
})
