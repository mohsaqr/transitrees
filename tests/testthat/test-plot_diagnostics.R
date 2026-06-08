# Tests for plot_distributions() and plot_predictive()

mk_diag_tree <- function(seed = 1L) {
  set.seed(seed)
  seqs <- replicate(100, sample(c("A", "B", "C"), 10, replace = TRUE),
                    simplify = FALSE)
  context_tree(seqs, max_depth = 2L)
}

test_that("plot_distributions returns a ggplot", {
  tree <- mk_diag_tree()
  p <- plot_distributions(tree, top = 6)
  expect_s3_class(p, "ggplot")
})

test_that("plot_distributions accepts an explicit context set", {
  tree <- mk_diag_tree()
  p <- plot_distributions(tree, contexts = c("(start)", "A"))
  expect_s3_class(p, "ggplot")
})

test_that("plot_distributions errors when no contexts qualify", {
  tree <- mk_diag_tree()
  expect_error(plot_distributions(tree, min_count = 1e6), "No contexts")
})

test_that("plot_predictive returns a ggplot for all types", {
  tree <- mk_diag_tree()
  set.seed(3)
  new <- replicate(12, sample(c("A", "B", "C"), 10, replace = TRUE),
                   simplify = FALSE)
  expect_s3_class(plot_predictive(tree, new, type = "position"), "ggplot")
  expect_s3_class(plot_predictive(tree, new, type = "ecdf"), "ggplot")
  expect_s3_class(plot_predictive(tree, new, type = "logloss"), "ggplot")
})

test_that("plot_predictive(logloss) equals -log2(predicted_prob)", {
  tree <- mk_diag_tree()
  set.seed(3)
  new <- replicate(12, sample(c("A", "B", "C"), 10, replace = TRUE),
                   simplify = FALSE)
  sp <- score_positions(tree, new)
  expect_equal(-sp$log_lik / log(2), -log2(sp$predicted_prob))
})

test_that("plot_pruning returns a ggplot and renders", {
  set.seed(5)
  seqs <- replicate(120, sample(c("A", "B", "C"), 12, replace = TRUE),
                    simplify = FALSE)
  tree <- context_tree(seqs, max_depth = 3L, min_count = 3L)
  p <- plot_pruning(tree, "A -> B -> C")
  expect_s3_class(p, "ggplot")
  pdf(NULL); on.exit(dev.off()); expect_silent(print(p))
})

test_that(".suffix_chain walks full chain to root with correct schema", {
  set.seed(5)
  seqs <- replicate(120, sample(c("A", "B", "C"), 12, replace = TRUE),
                    simplify = FALSE)
  tree <- context_tree(seqs, max_depth = 3L, min_count = 3L)
  ch <- transitiontrees:::.suffix_chain(tree, "A -> B -> C")
  expect_named(ch, c("L", "context", "label", "state", "prob", "count",
                     "g2", "diverges", "pruned", "retained", "status"))
  ## one row per (context, state); chain reaches the root
  expect_true(0L %in% ch$L)
  expect_setequal(ch$status[ch$L == 0L], "root")
  ## probabilities at each context sum to 1
  agg <- tapply(ch$prob, ch$context, sum)
  expect_true(all(abs(agg - 1) < 1e-8))
  ## an informative context must diverge; pruned/retained must not
  nb <- ch[!duplicated(ch$context) & ch$L > 0L, ]
  expect_true(all(nb$diverges[nb$status == "informative"]))
  expect_true(all(!nb$diverges[nb$status == "pruned"]))
  expect_true(all(!nb$diverges[nb$status == "retained"]))
})

test_that(".suffix_chain separates informative from retained ancestors", {
  ## A non-diverging ancestor that survives only because a DEEPER context
  ## diverges must be labelled 'retained', not 'informative'. (Codex review.)
  set.seed(1)
  seqs <- replicate(150, sample(c("A", "B", "C"), sample(6:12, 1),
                                replace = TRUE), simplify = FALSE)
  tree <- context_tree(seqs, max_depth = 3L, min_count = 3L)
  keys <- setdiff(names(tree$nodes), "<root>")
  deep <- keys[lengths(strsplit(keys, " -> ", fixed = TRUE)) == 3L]
  nb <- NULL
  for (p in deep) {
    ch <- transitiontrees:::.suffix_chain(tree, p)
    cand <- ch[!duplicated(ch$context) & ch$L > 0L, ]
    if (any(cand$status == "retained") && any(cand$status == "informative")) {
      nb <- cand; break
    }
  }
  expect_false(is.null(nb))                       # fixture must exhibit both
  expect_true(all(nb$diverges[nb$status == "informative"]))
  ## the retained ancestors do NOT themselves diverge — the bug was calling
  ## them 'kept (adds memory)'.
  expect_true(all(!nb$diverges[nb$status == "retained"]))
  expect_true(any(nb$status == "retained" & !nb$diverges))
})

test_that("plot_pruning errors when the pathway has no fitted context", {
  tree <- mk_diag_tree()  # max_depth 2
  expect_error(transitiontrees:::.suffix_chain(tree, "Z"),
               "not in the tree")
})

test_that("plot_pruning errors when the requested context is absent but a suffix exists", {
  ## A depth-1 tree cannot hold the 2-move context 'A -> B', though 'B'
  ## exists. The chain must NOT silently start at 'B'. (Codex review.)
  set.seed(3)
  seqs <- replicate(80, sample(c("A", "B", "C"), 8, replace = TRUE),
                    simplify = FALSE)
  tree <- context_tree(seqs, max_depth = 1L, min_count = 3L)
  expect_error(plot_pruning(tree, "A -> B"), "not in the tree")
  expect_error(transitiontrees:::.suffix_chain(tree, "A -> B"),
               "deepest fitted suffix")
})

test_that("plot_predictive errors when nothing scores", {
  tree <- mk_diag_tree()
  # all-out-of-alphabet held-out data -> zero scored positions
  expect_error(plot_predictive(tree, list(c("Z", "Z", "Z"))),
               "No held-out positions")
})

test_that("plot_difference returns a ggplot for a two-group object", {
  set.seed(1)
  gx <- replicate(40, sample(c("A","B","C"), 8, replace = TRUE,
                             prob = c(.2,.6,.2)), simplify = FALSE)
  gy <- replicate(40, sample(c("A","B","C"), 8, replace = TRUE,
                             prob = c(.2,.2,.6)), simplify = FALSE)
  grp <- context_tree(c(gx, gy), group = rep(c("x","y"), each = 40),
                      max_depth = 2L, min_count = 3L)
  expect_s3_class(plot_difference(grp), "ggplot")
  expect_s3_class(plot_difference(grp, depth = 1), "ggplot")      # order-1
  # significance overlay from a comparison
  cmp <- compare_groups(grp, iter = 99L)
  expect_s3_class(plot_difference(grp, comparison = cmp), "ggplot")
})

test_that("plot_difference requires naming two groups when there are >2", {
  set.seed(2)
  g <- context_tree(
    replicate(90, sample(c("A","B","C"), 8, replace = TRUE), simplify = FALSE),
    group = rep(c("p","q","r"), each = 30), max_depth = 1L, min_count = 3L)
  expect_error(plot_difference(g), "name the two")
  expect_s3_class(plot_difference(g, groups = c("p", "r")), "ggplot")
})

test_that("plot_difference layout = tree colours the phylogram by group difference", {
  set.seed(1)
  gx <- replicate(50, sample(c("A","B","C"), 8, replace = TRUE,
                             prob = c(.2,.6,.2)), simplify = FALSE)
  gy <- replicate(50, sample(c("A","B","C"), 8, replace = TRUE,
                             prob = c(.2,.2,.6)), simplify = FALSE)
  grp <- context_tree(c(gx, gy), group = rep(c("x","y"), each = 50),
                      max_depth = 2L, min_count = 3L)
  p <- plot_difference(grp, layout = "tree")
  expect_s3_class(p, "ggplot")
  pdf(NULL); on.exit(dev.off()); expect_silent(print(p))
})

test_that("plot_difference measure = residual/probability both render", {
  set.seed(1)
  gx <- replicate(40, sample(c("A","B","C"), 8, replace = TRUE,
                             prob = c(.2,.6,.2)), simplify = FALSE)
  gy <- replicate(40, sample(c("A","B","C"), 8, replace = TRUE,
                             prob = c(.2,.2,.6)), simplify = FALSE)
  grp <- context_tree(c(gx, gy), group = rep(c("x","y"), each = 40),
                      max_depth = 2L, min_count = 3L)
  expect_s3_class(plot_difference(grp, measure = "residual"), "ggplot")
  expect_s3_class(plot_difference(grp, measure = "probability"), "ggplot")
})
