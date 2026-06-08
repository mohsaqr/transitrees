# ---- Tests for tree_pathways() and friends ----

.three_state_tree <- function(seed = 1) {
  set.seed(seed)
  states <- c("A", "B", "C")
  P <- matrix(c(0.7, 0.2, 0.1,
                0.1, 0.7, 0.2,
                0.2, 0.2, 0.6),
              nrow = 3, byrow = TRUE,
              dimnames = list(states, states))
  seqs <- lapply(seq_len(80), function(i) {
    n <- 18L
    s <- character(n); s[1L] <- sample(states, 1L)
    for (t in seq.int(2L, n)) {
      s[t] <- sample(states, 1L, prob = P[s[t - 1L], ])
    }
    s
  })
  context_tree(seqs, max_depth = 3L, min_count = 5L)
}

# ---- structure ----

test_that("tree_pathways() returns the expected columns", {
  tree <- .three_state_tree()
  pw   <- tree_pathways(tree)
  expect_true(is.data.frame(pw))
  expect_named(pw, c("pathway", "depth", "count",
                     "likely_next", "next_probability", "divergence",
                     "changes_prediction"))
})

test_that("tree_pathways() includes (root) as a pathway of length 0", {
  tree <- .three_state_tree()
  pw   <- tree_pathways(tree)
  expect_true("(start)" %in% pw$pathway)
  root_row <- pw[pw$pathway == "(start)", ]
  expect_equal(root_row$depth, 0L)
  expect_true(is.na(root_row$divergence))
})

test_that("tree_pathways() in arrow notation, never with > separator", {
  tree <- .three_state_tree()
  pw   <- tree_pathways(tree)
  non_root <- pw$pathway[pw$pathway != "(start)"]
  ## Arrow notation present in any multi-state pathway
  multi <- non_root[grepl(" ", non_root)]
  expect_true(all(grepl(" -> ", multi)))
  expect_true(all(!grepl(" > ", multi)))
})

test_that("tree_pathways() default sort is by count, descending", {
  tree <- .three_state_tree()
  pw   <- tree_pathways(tree)
  expect_equal(pw$count, sort(pw$count, decreasing = TRUE))
})

test_that("sort_by 'divergence' sorts non-NA descending", {
  tree <- .three_state_tree()
  pw   <- tree_pathways(tree, sort_by = "divergence")
  finite_div <- pw$divergence[!is.na(pw$divergence)]
  expect_equal(finite_div, sort(finite_div, decreasing = TRUE))
  ## Any NA (only the root) should be the last row
  na_pos <- which(is.na(pw$divergence))
  if (length(na_pos))
    expect_true(all(na_pos >= nrow(pw) - length(na_pos) + 1L))
})

test_that("min_count filters out rare pathways", {
  tree <- .three_state_tree()
  all_pw <- tree_pathways(tree, min_count = 1L)
  big    <- tree_pathways(tree, min_count = 30L)
  expect_lte(nrow(big), nrow(all_pw))
  expect_true(all(big$count >= 30L))
})

test_that("empty tree_pathways() result is schema-stable", {
  tree <- .three_state_tree()
  empty <- tree_pathways(tree, min_count = 10000L)
  expect_true(is.data.frame(empty))
  expect_equal(nrow(empty), 0L)
  expect_named(empty, c("pathway", "depth", "count",
                        "likely_next", "next_probability", "divergence",
                        "changes_prediction"))
})

# ---- common_pathways ----

test_that("common_pathways returns top n by count", {
  tree <- .three_state_tree()
  cp <- common_pathways(tree, top = 5L)
  expect_equal(nrow(cp), 5L)
  expect_equal(cp$count, sort(cp$count, decreasing = TRUE))
})

test_that("common_pathways(depth = k) restricts to depth k", {
  tree <- .three_state_tree()
  cp <- common_pathways(tree, top = 99L, depth = 2L)
  expect_true(all(cp$depth == 2L))
})

# ---- divergent_pathways ----

test_that("divergent_pathways returns top n by divergence", {
  tree <- .three_state_tree()
  dp <- divergent_pathways(tree, top = 4L, min_count = 5L)
  expect_lte(nrow(dp), 4L)
  expect_true(all(!is.na(dp$divergence)))
  expect_equal(dp$divergence, sort(dp$divergence, decreasing = TRUE))
})

test_that("divergent_pathways(flips_only = TRUE) returns only flipping rows", {
  tree <- .three_state_tree()
  dp <- divergent_pathways(tree, flips_only = TRUE,
                           min_count = 5L, top = 99L)
  if (nrow(dp) > 0L)
    expect_true(all(dp$changes_prediction, na.rm = TRUE))
})

# ---- sharp_pathways ----

test_that("sharp_pathways returns top n by next_probability", {
  tree <- .three_state_tree()
  sp <- sharp_pathways(tree, top = 5L, min_count = 5L)
  expect_lte(nrow(sp), 5L)
  expect_equal(sp$next_probability,
               sort(sp$next_probability, decreasing = TRUE))
})

# ---- tree_pathways is a plain function (no longer an S3 generic) ----

test_that("tree_pathways is a plain exported function", {
  expect_true(is.function(tree_pathways))
  ## The generic was retired to avoid collision with other packages' pathways().
  expect_false(isTRUE(any(grepl("UseMethod", deparse(body(tree_pathways))))))
})
