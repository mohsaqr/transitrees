# Tests for plot_trajectories()

mk_traj_tree <- function(seed = 1L) {
  set.seed(seed)
  seqs <- replicate(120, sample(c("A", "B", "C"), 8, replace = TRUE),
                    simplify = FALSE)
  context_tree(seqs, max_depth = 3L, min_count = 3L)
}

test_that("plot_trajectories returns a ggplot for both measures", {
  skip_if_not_installed("ggforce")
  tree <- mk_traj_tree()
  expect_s3_class(plot_trajectories(tree, "frequency"), "ggplot")
  expect_s3_class(plot_trajectories(tree, "predictability"), "ggplot")
})

test_that("plot_trajectories renders without error", {
  skip_if_not_installed("ggforce")
  tree <- mk_traj_tree()
  p <- plot_trajectories(tree, "predictability")
  pdf(NULL); on.exit(dev.off()); expect_silent(print(p))
})

test_that("plot_trajectories errors when no prefix clears min_count", {
  skip_if_not_installed("ggforce")
  tree <- mk_traj_tree()
  expect_error(plot_trajectories(tree, min_count = 1e6), "min_count")
})

test_that("plot_trajectories accepts a pruned tree", {
  skip_if_not_installed("ggforce")
  tree <- mk_traj_tree()
  pruned <- prune_tree(tree)
  expect_s3_class(plot_trajectories(pruned, "predictability"), "ggplot")
})

test_that("plot_trajectories validates its input", {
  expect_error(plot_trajectories(list(1, 2, 3)), "transitiontrees")
})

# ---- data-correctness of the underlying trajectory table (audit) ----

test_that(".trajectory_data has correct schema, counts, connectivity, probs", {
  tree <- mk_traj_tree()
  seqs <- transitiontrees:::.ct_traj(tree$data)
  D <- transitiontrees:::.trajectory_data(tree, min_count = 4L)

  expect_named(D, c("node", "parent", "depth", "count", "last",
                    "eprob"))
  ## root row: depth 0, NA parent/last/eprob
  root <- D[D$node == "(start)", ]
  expect_identical(nrow(root), 1L)
  expect_identical(root$depth, 0L)
  expect_true(is.na(root$eprob))

  ## prefix counts match a manual recount of the sequences
  manual <- table(unlist(lapply(seqs, function(s)
    vapply(seq_along(s), function(k) paste(s[seq_len(k)], collapse = " -> "),
           character(1)))))
  nb <- D[D$node != "(start)", ]
  expect_equal(nb$count, as.integer(manual[nb$node]))
  expect_true(all(nb$count >= 4L))             # min_count honoured

  ## connectivity: every node's parent is itself a node in the table
  expect_true(all(nb$parent %in% D$node))

  ## predictability values are valid probabilities and match query_pathway
  expect_true(all(nb$eprob >= 0 & nb$eprob <= 1))
  deep <- nb[lengths(strsplit(nb$node, " -> ", fixed = TRUE)) >= 2L, ]
  if (nrow(deep)) {
    r <- deep[1, ]
    toks   <- strsplit(r$node, " -> ", fixed = TRUE)[[1]]
    hist   <- paste(utils::head(toks, -1L), collapse = " -> ")
    expect_equal(r$eprob,
                 as.numeric(query_pathway(tree, hist)[[r$last]]))
  }
})

test_that(".trajectory_data drops prefixes whose parent was filtered out", {
  tree <- mk_traj_tree()
  D <- transitiontrees:::.trajectory_data(tree, min_count = 4L)
  ## tree is connected: no non-root node points at a missing parent
  nb <- D[D$node != "(start)", ]
  expect_length(setdiff(nb$parent, D$node), 0L)
})
