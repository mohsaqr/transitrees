# ---- Tests for context_tree() ----

.simple_seqs <- function() {
  list(
    c("A", "B", "A", "B", "A"),
    c("B", "A", "B", "A", "B"),
    c("A", "A", "B", "B", "A"),
    c("B", "B", "A", "A", "B"),
    c("A", "B", "B", "A", "A"),
    c("B", "A", "A", "B", "B")
  )
}

# ---- structure ----

test_that("context_tree returns ctxtree with correct fields", {
  tree <- context_tree(.simple_seqs(), max_depth = 2L, nmin = 1L)
  expect_s3_class(tree, "pathtree")
  expect_named(tree, c("nodes", "edges", "alphabet",
                       "max_depth", "nmin",
                       "n_seq", "n_obs", "smoothing",
                       "pruned", "pruning", "data"))
  expect_setequal(tree$alphabet, c("A", "B"))
  expect_equal(tree$n_seq, 6L)
  expect_false(isTRUE(tree$pruned))
})

test_that("root node carries the marginal next-state distribution", {
  tree <- context_tree(.simple_seqs(), max_depth = 2L, nmin = 1L)
  root <- tree$nodes[["<root>"]]
  expect_false(is.null(root))
  expect_equal(sum(root$prob), 1, tolerance = 1e-12)
  expect_equal(length(root$prob), 2L)
  expect_equal(root$depth, 0L)
})

test_that("root node is retained even when nmin exceeds observations", {
  tree <- context_tree(.simple_seqs(), max_depth = 2L, nmin = 1000L)
  expect_true("<root>" %in% names(tree$nodes))
  expect_equal(length(tree$nodes), 1L)
  expect_equal(tree$max_depth, 0L)
})

test_that("node probabilities sum to 1 (after smoothing)", {
  tree <- context_tree(.simple_seqs(), max_depth = 3L, nmin = 1L)
  for (info in tree$nodes) {
    expect_equal(sum(info$prob), 1, tolerance = 1e-10)
  }
})

test_that("nmin filters out rare contexts", {
  tree_loose  <- context_tree(.simple_seqs(), max_depth = 2L, nmin = 1L)
  tree_strict <- context_tree(.simple_seqs(), max_depth = 2L, nmin = 5L)
  expect_lte(length(tree_strict$nodes), length(tree_loose$nodes))
})

test_that("depth respects max_depth", {
  tree <- context_tree(.simple_seqs(), max_depth = 2L, nmin = 1L)
  depths <- vapply(tree$nodes, function(x) x$depth, integer(1))
  expect_true(all(depths <= 2L))
})

# ---- input dispatch ----

test_that("context_tree accepts a character matrix", {
  m <- do.call(rbind, lapply(.simple_seqs(),
                              function(x) c(x, rep(NA, 6L - length(x)))))
  tree <- context_tree(m, max_depth = 2L, nmin = 1L)
  expect_s3_class(tree, "pathtree")
  expect_setequal(tree$alphabet, c("A", "B"))
})

test_that("context_tree accepts a wide data.frame", {
  df <- as.data.frame(do.call(rbind,
        lapply(.simple_seqs(),
               function(x) c(x, rep(NA, 6L - length(x))))),
        stringsAsFactors = FALSE)
  tree <- context_tree(df, max_depth = 2L, nmin = 1L)
  expect_s3_class(tree, "pathtree")
})

test_that("context_tree errors on unusable input", {
  expect_error(context_tree(matrix(1:6, 2, 3)),
               regexp = "data.*must be|wide data\\.frame")
  expect_error(context_tree(123L), regexp = "data.*must be|wide data\\.frame")
})

# ---- print/summary ----

test_that("print and summary dispatch correctly", {
  tree <- context_tree(.simple_seqs(), max_depth = 2L, nmin = 1L)
  expect_output(print(tree), "<pathtree>")
  s <- summary(tree)
  expect_s3_class(s, "summary.pathtree")
  expect_true(is.data.frame(s$table))
  expect_true(any(s$table$pathway == "(root)"))
  expect_output(print(s), "pathtree summary")
})

test_that("summary(tree)$table carries the canonical 7-column schema", {
  tree <- context_tree(.simple_seqs(), max_depth = 2L, nmin = 1L)
  s <- summary(tree)
  expect_named(s$table, c("pathway", "depth", "count",
                          "modal_next", "prob_next", "KL", "flips"))
  ## And it sorts by (depth, -count) — structural tree order.
  expect_equal(s$table$depth, sort(s$table$depth))
})

test_that("as.data.frame.pathtree is the canonical tidy view", {
  tree <- context_tree(.simple_seqs(), max_depth = 2L, nmin = 1L)
  df <- as.data.frame(tree)
  expect_s3_class(df, "data.frame")
  expect_named(df, c("pathway", "depth", "count",
                     "modal_next", "prob_next", "KL", "flips"))
  expect_equal(nrow(df), length(tree$nodes))
  expect_true(all(df$prob_next > 0 & df$prob_next <= 1))
  ## Identical to pathtree_pathways()
  expect_identical(df, pathtree_pathways(tree))
})
