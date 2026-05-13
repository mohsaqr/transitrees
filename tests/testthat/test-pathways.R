# ---- Tests for pathtree_pathways() and friends ----

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
  context_tree(seqs, max_depth = 3L, nmin = 5L)
}

# ---- structure ----

test_that("pathtree_pathways() returns the expected columns", {
  tree <- .three_state_tree()
  pw   <- pathtree_pathways(tree)
  expect_true(is.data.frame(pw))
  expect_named(pw, c("pathway", "depth", "count",
                     "modal_next", "prob_next", "KL", "flips"))
})

test_that("pathtree_pathways() includes (root) as a pathway of length 0", {
  tree <- .three_state_tree()
  pw   <- pathtree_pathways(tree)
  expect_true("(root)" %in% pw$pathway)
  root_row <- pw[pw$pathway == "(root)", ]
  expect_equal(root_row$depth, 0L)
  expect_true(is.na(root_row$KL))
})

test_that("pathtree_pathways() in arrow notation, never with > separator", {
  tree <- .three_state_tree()
  pw   <- pathtree_pathways(tree)
  non_root <- pw$pathway[pw$pathway != "(root)"]
  ## Arrow notation present in any multi-state pathway
  multi <- non_root[grepl(" ", non_root)]
  expect_true(all(grepl(" -> ", multi)))
  expect_true(all(!grepl(" > ", multi)))
})

test_that("pathtree_pathways() default sort is by count, descending", {
  tree <- .three_state_tree()
  pw   <- pathtree_pathways(tree)
  expect_equal(pw$count, sort(pw$count, decreasing = TRUE))
})

test_that("pathtree_pathways(sort_by = 'KL') sorts non-NA KL descending", {
  tree <- .three_state_tree()
  pw   <- pathtree_pathways(tree, sort_by = "KL")
  finite_kl <- pw$KL[!is.na(pw$KL)]
  expect_equal(finite_kl, sort(finite_kl, decreasing = TRUE))
  ## Any NA (only the root) should be the last row
  na_pos <- which(is.na(pw$KL))
  if (length(na_pos))
    expect_true(all(na_pos >= nrow(pw) - length(na_pos) + 1L))
})

test_that("min_count filters out rare pathways", {
  tree <- .three_state_tree()
  all_pw <- pathtree_pathways(tree, min_count = 1L)
  big    <- pathtree_pathways(tree, min_count = 30L)
  expect_lte(nrow(big), nrow(all_pw))
  expect_true(all(big$count >= 30L))
})

test_that("empty pathtree_pathways() result is schema-stable", {
  tree <- .three_state_tree()
  empty <- pathtree_pathways(tree, min_count = 10000L)
  expect_true(is.data.frame(empty))
  expect_equal(nrow(empty), 0L)
  expect_named(empty, c("pathway", "depth", "count",
                        "modal_next", "prob_next", "KL", "flips"))
})

# ---- common_pathways ----

test_that("common_pathways returns top n by count", {
  tree <- .three_state_tree()
  cp <- common_pathways(tree, n = 5L)
  expect_equal(nrow(cp), 5L)
  expect_equal(cp$count, sort(cp$count, decreasing = TRUE))
})

test_that("common_pathways(depth = k) restricts to depth k", {
  tree <- .three_state_tree()
  cp <- common_pathways(tree, n = 99L, depth = 2L)
  expect_true(all(cp$depth == 2L))
})

# ---- divergent_pathways ----

test_that("divergent_pathways returns top n by KL", {
  tree <- .three_state_tree()
  dp <- divergent_pathways(tree, n = 4L, min_count = 5L)
  expect_lte(nrow(dp), 4L)
  expect_true(all(!is.na(dp$KL)))
  expect_equal(dp$KL, sort(dp$KL, decreasing = TRUE))
})

test_that("divergent_pathways(flips_only = TRUE) returns only flipping rows", {
  tree <- .three_state_tree()
  dp <- divergent_pathways(tree, flips_only = TRUE,
                           min_count = 5L, n = 99L)
  if (nrow(dp) > 0L) expect_true(all(dp$flips, na.rm = TRUE))
})

# ---- sharp_pathways ----

test_that("sharp_pathways returns top n by prob_next", {
  tree <- .three_state_tree()
  sp <- sharp_pathways(tree, n = 5L, min_count = 5L)
  expect_lte(nrow(sp), 5L)
  expect_equal(sp$prob_next, sort(sp$prob_next, decreasing = TRUE))
})

# ---- pathtree_pathways is a plain function (no longer an S3 generic) ----

test_that("pathtree_pathways is a plain exported function", {
  expect_true(is.function(pathtree_pathways))
  ## The generic was retired to avoid collision with Nestimate::pathways.
  expect_false(isTRUE(any(grepl("UseMethod", deparse(body(pathtree_pathways))))))
})
