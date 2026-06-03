# ---- Tests for predict() and generate_sequences() ----

.three_state_seqs <- function(seed = 1) {
  set.seed(seed)
  states <- c("A", "B", "C")
  P <- matrix(c(0.7, 0.2, 0.1,
                0.1, 0.7, 0.2,
                0.2, 0.2, 0.6),
              nrow = 3, byrow = TRUE,
              dimnames = list(states, states))
  lapply(seq_len(60), function(i) {
    n <- 18L
    s <- character(n)
    s[1L] <- sample(states, 1L)
    for (t in seq.int(2L, n)) {
      s[t] <- sample(states, 1L, prob = P[s[t - 1L], ])
    }
    s
  })
}

test_that("predict returns a probability matrix summing to 1 per row", {
  tree <- context_tree(.three_state_seqs(), max_depth = 2L, min_count = 5L)
  preds <- predict(tree, newdata = list(c("A","B"), c("B","C"), c("C","A")))
  expect_true(is.matrix(preds))
  expect_equal(ncol(preds), 3L)
  expect_equal(rowSums(preds), rep(1, 3), tolerance = 1e-10)
  expect_true(all(preds >= 0))
})

test_that("predict(type = 'class') returns modal predictions", {
  tree <- context_tree(.three_state_seqs(), max_depth = 2L, min_count = 5L)
  classes <- predict(tree,
                     newdata = list(c("A","A"), c("B","B"), c("C","C")),
                     type = "class")
  expect_type(classes, "character")
  expect_length(classes, 3L)
  expect_true(all(classes %in% tree$alphabet))
})

test_that("predict on a list-wrapped single history returns a 1xk matrix", {
  ## Schema-stable: container input (list/df/matrix) always yields a matrix.
  tree <- context_tree(.three_state_seqs(), max_depth = 2L, min_count = 5L)
  m <- predict(tree, newdata = list(c("A","B")))
  expect_true(is.matrix(m))
  expect_equal(dim(m), c(1L, 3L))
  expect_equal(colnames(m), tree$alphabet)
  expect_equal(sum(m), 1, tolerance = 1e-10)
})

test_that("predict on a bare character vector returns a named vector", {
  ## Interactive shortcut: bare vector input collapses to a named vector.
  tree <- context_tree(.three_state_seqs(), max_depth = 2L, min_count = 5L)
  v <- predict(tree, newdata = c("A","B"))
  expect_true(is.numeric(v))
  expect_false(is.matrix(v))
  expect_equal(length(v), 3L)
  expect_named(v, tree$alphabet)
})

test_that("predict falls back to root for empty histories", {
  tree <- context_tree(.three_state_seqs(), max_depth = 2L, min_count = 5L)
  m <- predict(tree, newdata = list(character(0)))
  expect_equal(unname(m[1L, ]), tree$nodes[["<root>"]]$prob,
               tolerance = 1e-12)
})

test_that("predict accepts a data.frame and a matrix newdata (-> matrix)", {
  tree <- context_tree(.three_state_seqs(), max_depth = 2L, min_count = 5L)
  df <- data.frame(t1 = c("A", "B"), t2 = c("B", "C"),
                   stringsAsFactors = FALSE)
  pm <- predict(tree, newdata = df)
  expect_true(is.matrix(pm))
  expect_equal(dim(pm), c(2L, 3L))
  expect_equal(unname(rowSums(pm)), c(1, 1), tolerance = 1e-10)
  ## matrix input gives the same probabilities as the equivalent data.frame
  mm <- predict(tree, newdata = as.matrix(df))
  expect_equal(unname(pm), unname(mm))
  ## type = "class" on a container returns one modal label per row
  cl <- predict(tree, newdata = df, type = "class")
  expect_length(cl, 2L)
  expect_true(all(cl %in% tree$alphabet))
})

test_that("predict rejects an unsupported newdata type", {
  tree <- context_tree(.three_state_seqs(), max_depth = 2L, min_count = 5L)
  expect_error(predict(tree, newdata = 1:3),
               "must be a list, data.frame, matrix, or character vector")
})

test_that("generate_sequences returns the right shape", {
  tree <- context_tree(.three_state_seqs(), max_depth = 2L, min_count = 5L)
  set.seed(42)
  m <- generate_sequences(tree, n = 5L, length = 10L)
  expect_equal(dim(m), c(5L, 10L))
  expect_true(all(m %in% tree$alphabet))
})

test_that("generate_sequences honours start", {
  tree <- context_tree(.three_state_seqs(), max_depth = 2L, min_count = 5L)
  set.seed(42)
  m <- generate_sequences(tree, n = 4L, length = 5L,
                          start = c("A","A","B","C"))
  expect_equal(m[, 1L], c("A","A","B","C"))
})

test_that("simulate.pathtree is a generate_sequences alias with nsim", {
  tree <- context_tree(.three_state_seqs(), max_depth = 2L, min_count = 5L)
  m <- simulate(tree, nsim = 4L, seed = 42L, length = 6L,
                start = c("A", "B", "C", "A"))
  expect_equal(dim(m), c(4L, 6L))
  expect_equal(m[, 1L], c("A", "B", "C", "A"))
  set.seed(42)
  g <- generate_sequences(tree, n = 4L, length = 6L,
                          start = c("A", "B", "C", "A"))
  expect_equal(m, g)
})

test_that("generate_sequences validates n and length", {
  tree <- context_tree(.three_state_seqs(), max_depth = 2L, min_count = 5L)
  expect_error(generate_sequences(tree, n = 2L, length = 0L), "'length'")
  expect_error(generate_sequences(tree, n = 0L, length = 5L), "'n'")
  expect_error(generate_sequences(tree, n = 2L, length = NA), "'length'")
  ## length == 1 returns just the start column (no descending loop)
  m1 <- generate_sequences(tree, n = 3L, length = 1L)
  expect_equal(dim(m1), c(3L, 1L))
  expect_false(anyNA(m1))
})
