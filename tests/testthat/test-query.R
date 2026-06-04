# Tests for Phase D: query_pathway / subtree / pathway_exists

mk_simple <- function(seed = 1L) {
  set.seed(seed)
  m <- matrix(sample(c("A","B","C"), 30 * 12, replace = TRUE), 30, 12)
  context_tree(m, max_depth = 2L, min_count = 3L)
}

test_that("query_pathway returns probability vector for existing pathway", {
  tr <- mk_simple()
  ctx <- names(tr$nodes)[2L]      # any non-root context
  p <- query_pathway(tr, ctx)
  expect_equal(sum(p), 1)
  expect_named(p, tr$alphabet)
})

test_that("query_pathway accepts arrow string and character vector equally", {
  tr <- mk_simple()
  p_str <- query_pathway(tr, "A")
  p_vec <- query_pathway(tr, c("A"))
  expect_equal(p_str, p_vec)
})

test_that("query_pathway with next_state returns scalar probability", {
  tr <- mk_simple()
  v <- query_pathway(tr, "A", next_state = "B")
  expect_length(v, 1L)
  expect_equal(v, query_pathway(tr, "A")[["B"]])
})

test_that("query_pathway exact = TRUE returns NA on missing pathway", {
  tr <- mk_simple()
  fake <- "Z -> Z -> Z"
  expect_true(all(is.na(query_pathway(tr, fake, exact = TRUE))))
})

test_that("query_pathway exact = FALSE falls back to suffix match", {
  tr <- mk_simple()
  ## The fall-back uses .ct_match_context: an arbitrary history is
  ## resolved to the deepest matching suffix or the root.
  p <- query_pathway(tr, c("Z","A"), exact = FALSE)
  expect_equal(sum(p), 1)
})

test_that("pathway_exists returns TRUE/FALSE correctly", {
  tr <- mk_simple()
  expect_true(pathway_exists(tr, "(start)"))   # new display label
  expect_true(pathway_exists(tr, "(root)"))    # legacy label still accepted
  expect_true(pathway_exists(tr, names(tr$nodes)[2L]))
  expect_false(pathway_exists(tr, "Z -> Z -> Z"))
})

test_that("subtree() returns a strict descendant of the original tree", {
  tr <- mk_simple()
  ## Find a depth-1 node to root the subtree at
  d1 <- vapply(tr$nodes, function(n) n$depth, integer(1)) == 1L
  pathway <- names(tr$nodes)[d1][1L]
  sub <- subtree(tr, pathway)
  expect_s3_class(sub, "transitiontrees")
  expect_true(all(names(sub$nodes) %in% names(tr$nodes)))
  expect_true(pathway %in% names(sub$nodes))
  expect_identical(attr(sub, "local_root"), pathway)
})

test_that("subtree() of a leaf has just one node and no edges", {
  tr <- mk_simple()
  leaves <- setdiff(names(tr$nodes), tr$edges$parent)
  if (length(leaves) > 0L) {
    sub <- subtree(tr, leaves[1L])
    expect_length(sub$nodes, 1L)
    expect_equal(nrow(sub$edges), 0L)
  }
})

test_that("subtree() errors on missing pathway", {
  tr <- mk_simple()
  expect_error(subtree(tr, "ZZZ"), "not a node")
})
