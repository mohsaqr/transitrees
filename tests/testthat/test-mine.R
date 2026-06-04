# Tests for mine_contexts() and mine_sequences()

mk_mine_tree <- function(seed = 1L) {
  set.seed(seed)
  seqs <- replicate(120, sample(c("A", "B", "C"), 10, replace = TRUE),
                    simplify = FALSE)
  context_tree(seqs, max_depth = 2L)
}

test_that("mine_contexts returns the canonical schema sorted by prob", {
  tree <- mk_mine_tree()
  out <- mine_contexts(tree, state = "A")
  expect_s3_class(out, "data.frame")
  expect_identical(names(out),
                   c("pathway", "depth", "count", "state", "prob", "is_modal"))
  expect_true(all(out$state == "A"))
  expect_true(all(diff(out$prob) <= 1e-12))   # descending
  expect_true(all(out$prob >= 0 & out$prob <= 1))
})

test_that("mine_contexts filters on min_prob / max_prob", {
  tree <- mk_mine_tree()
  hi <- mine_contexts(tree, state = "A", min_prob = 0.4)
  expect_true(all(hi$prob >= 0.4))
  lo <- mine_contexts(tree, state = "A", max_prob = 0.3)
  expect_true(all(lo$prob <= 0.3))
  band <- mine_contexts(tree, state = "A", min_prob = 0.2, max_prob = 0.5)
  expect_true(all(band$prob >= 0.2 & band$prob <= 0.5))
})

test_that("mine_contexts is_modal flags the context's modal state", {
  tree <- mk_mine_tree()
  out <- mine_contexts(tree, state = "A")
  # a row is_modal iff "A" is that context's argmax next state
  chk <- vapply(out$pathway, function(p) {
    ctx <- if (p == "(start)") "<root>" else p
    tree$alphabet[which.max(tree$nodes[[ctx]]$prob)] == "A"
  }, logical(1))
  expect_equal(out$is_modal, unname(chk))
})

test_that("mine_contexts errors on a state outside the alphabet", {
  tree <- mk_mine_tree()
  expect_error(mine_contexts(tree, state = "Z"), "alphabet")
  expect_error(mine_contexts(tree), "alphabet")
})

test_that("mine_contexts min_count returns an empty typed frame when nothing qualifies", {
  tree <- mk_mine_tree()
  out <- mine_contexts(tree, state = "A", min_count = 1e6)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 0L)
  expect_identical(names(out),
                   c("pathway", "depth", "count", "state", "prob", "is_modal"))
})

test_that("mine_sequences ranks by perplexity in both directions", {
  tree <- mk_mine_tree()
  set.seed(2)
  new <- replicate(30, sample(c("A", "B", "C"), 10, replace = TRUE),
                   simplify = FALSE)
  surprising <- mine_sequences(tree, new, n = 5, which = "surprising")
  expected   <- mine_sequences(tree, new, n = 5, which = "expected")
  expect_equal(nrow(surprising), 5L)
  expect_identical(names(surprising),
                   c("sequence_id", "n_scored", "log_lik", "perplexity"))
  expect_true(all(diff(surprising$perplexity) <= 1e-9))   # high -> low
  expect_true(all(diff(expected$perplexity)   >= -1e-9))  # low -> high
  expect_gte(surprising$perplexity[1], expected$perplexity[1])
})
