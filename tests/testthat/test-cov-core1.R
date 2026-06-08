# Targeted coverage tests for currently-uncovered branches.
# Base R + testthat only; deterministic, small.

# ---- shared fixtures ----------------------------------------------------

make_tree <- function() {
  set.seed(1)
  m <- matrix(sample(c("A", "B", "C"), 400, replace = TRUE), nrow = 40)
  context_tree(m, max_depth = 3L, min_count = 2L)
}

# ---- path_dependence.R --------------------------------------------------

test_that("tree_dependence on a root-only tree returns the empty schema", {
  # max_depth = 0 yields a tree with only the root node, so there are no
  # non-root contexts and tree_dependence() takes the empty-data.frame path
  # (.pt_empty_dependence_df()).
  set.seed(2)
  m <- matrix(sample(c("A", "B", "C"), 300, replace = TRUE), nrow = 30)
  tr0 <- context_tree(m, max_depth = 0L, min_count = 1L)
  expect_length(tr0$nodes, 1L)

  dep <- tree_dependence(tr0)
  expect_s3_class(dep, "data.frame")
  expect_equal(nrow(dep), 0L)
  expect_identical(
    names(dep),
    c("pathway", "depth", "count", "divergence", "entropy",
      "entropy_before", "entropy_drop", "likely_next", "likely_before",
      "changes_prediction")
  )
  expect_type(dep$changes_prediction, "logical")
})

test_that("tree_dependence skips a context whose suffix-parent is absent", {
  # A subtree drops the root, so its local-root node (here "A") has a
  # suffix-parent (.ROOT) that is no longer present; that row is skipped
  # (path_dependence.R parent-missing guard).
  tr <- make_tree()
  st <- subtree(tr, "A")
  expect_false(transitiontrees:::.ROOT %in% names(st$nodes))

  dep <- tree_dependence(st)
  expect_s3_class(dep, "data.frame")
  # The local root "A" is excluded because its parent is missing.
  expect_false("A" %in% dep$pathway)
  expect_true(nrow(dep) >= 1L)
})

# ---- pathways.R ---------------------------------------------------------

test_that("tree_pathways yields NA KL/flips when the parent is missing", {
  # Same subtree situation: node "A" has no parent in the subtree, so KL and
  # flips fall to NA (pathways.R par_info-NULL branch).
  tr <- make_tree()
  st <- subtree(tr, "A")
  pw <- tree_pathways(st, min_count = 1L)
  row_a <- pw[pw$pathway == "A", , drop = FALSE]
  expect_equal(nrow(row_a), 1L)
  expect_true(is.na(row_a$divergence))
  expect_true(is.na(row_a$changes_prediction))
})

# ---- prune.R ------------------------------------------------------------

test_that("prune_tree validates threshold", {
  tr <- make_tree()
  expect_error(prune_tree(tr, criterion = "KL", threshold = -1),
               "threshold")
})

test_that("prune_tree on a subtree protects an orphaned local root", {
  # Aggressive KL pruning collapses the whole subtree; when the local root
  # "A" becomes a leaf its suffix-parent (.ROOT) is absent, so the parent
  # guard skips it instead of removing it (prune.R parent-missing branch).
  tr <- make_tree()
  st <- subtree(tr, "A")
  pr <- prune_tree(st, criterion = "KL", threshold = 100)
  expect_s3_class(pr, "transitiontrees")
  expect_identical(names(pr$nodes), "A")
})

# ---- predict.R ----------------------------------------------------------

test_that("predict returns a uniform row when the matched context is absent", {
  # On a subtree the root is gone; a history matching no subtree node falls
  # back to .ROOT, which is missing, triggering the uniform-distribution
  # guard in predict.transitiontrees().
  tr <- make_tree()
  st <- subtree(tr, "A")
  p <- predict(st, list(c("B")), type = "prob")
  expect_true(all(abs(p[1, ] - 1 / length(st$alphabet)) < 1e-9))
})

test_that("generate_sequences rejects a mis-sized start vector", {
  tr <- make_tree()
  expect_error(
    generate_sequences(tr, n = 3L, length = 4L, start = c("A", "B")),
    "length"
  )
})

# ---- query.R ------------------------------------------------------------

test_that("query_pathway normalises an empty pathway to the root", {
  tr <- make_tree()
  root_prob <- query_pathway(tr, character(0))
  expect_named(root_prob, tr$alphabet)
  expect_equal(sum(root_prob), 1, tolerance = 1e-8)
  expect_equal(query_pathway(tr, NULL), root_prob)
})

# ---- prepare_input.R (.pt_parse_time) -----------------------------------

test_that(".pt_parse_time handles Date, ISO, formatted and bad input", {
  pt <- transitiontrees:::.pt_parse_time
  # Date branch
  d <- pt(as.Date(c("2020-01-01", "2020-01-02")))
  expect_s3_class(d, "POSIXct")
  # character ISO (no format)
  iso <- pt(c("2020-01-01 09:00:00", "2020-01-01 10:00:00"))
  expect_s3_class(iso, "POSIXct")
  # character with explicit format
  fmt <- pt(c("01/02/2020", "03/04/2020"), format = "%m/%d/%Y")
  expect_s3_class(fmt, "POSIXct")
  # unparseable character -> errors during parsing (base as.POSIXct throws
  # for the no-format ISO path before the package's own message is reached)
  expect_error(pt(c("not a date", "still not")))
})

test_that("prepare_input parses a Date time column end to end", {
  long <- data.frame(
    user  = c("a", "a", "b", "b"),
    day   = as.Date(c("2020-01-01", "2020-01-02", "2020-01-01", "2020-01-03")),
    state = c("X", "Y", "Y", "Z"),
    stringsAsFactors = FALSE
  )
  wide <- prepare_input(long, actor = "user", time = "day", action = "state")
  expect_s3_class(wide, "data.frame")
  expect_true(nrow(wide) >= 1L)
})

# ---- smoothing.R --------------------------------------------------------

test_that("smoothing kernels return the back-off on empty counts", {
  z <- c(0, 0, 0)
  # floor: NULL parent -> uniform; supplied parent -> parent unchanged
  expect_equal(transitiontrees:::.smooth_floor(z, NULL, ymin = 0.001),
               rep(1 / 3, 3))
  expect_equal(transitiontrees:::.smooth_floor(z, c(.2, .3, .5)),
               c(.2, .3, .5))
  expect_equal(transitiontrees:::.smooth_kneser_ney(z, c(.2, .3, .5)),
               c(.2, .3, .5))
  expect_equal(transitiontrees:::.smooth_jelinek_mercer(z, c(.2, .3, .5)),
               c(.2, .3, .5))
})

test_that(".ct_smooth_dispatch errors on an unknown resolved method", {
  expect_error(
    transitiontrees:::.ct_smooth_dispatch(list(method = "bogus"), c(1, 2, 3)),
    "Unknown smoothing method"
  )
})

test_that(".pt_resolve_smoothing rejects bad specs", {
  # list with an unknown method
  expect_error(transitiontrees:::.pt_resolve_smoothing(list("bogus")),
               "Smoothing method must be one of")
  # neither character nor list
  expect_error(transitiontrees:::.pt_resolve_smoothing(42),
               "must be a method name")
})

test_that(".pt_validate_smoothing checks the floor rule and unknown methods", {
  expect_error(
    transitiontrees:::.pt_validate_smoothing(
      list(method = "floor", ymin = 0.01, rule = "nonsense")),
    "rule"
  )
  expect_error(
    transitiontrees:::.pt_validate_smoothing(list(method = "bogus")),
    "Unknown smoothing method"
  )
})

test_that("suffix-walk helpers handle the root sentinel", {
  root <- transitiontrees:::.ROOT
  # parent of the root is NA
  expect_true(is.na(transitiontrees:::.pt_parent_ctx(root)))
  # walking up from a depth-1 context whose root is absent yields NULL
  nodes <- list(A = list(prob = c(0.5, 0.3, 0.2)))
  expect_null(transitiontrees:::.pt_parent_prob(nodes, "A"))
  # the root itself has no parent prob
  expect_null(transitiontrees:::.pt_parent_prob(list(), root))
})
