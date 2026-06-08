# Targeted coverage tests (round 2) for currently-uncovered branches in
# context_tree.R, bootstrap.R, compare.R, compare_groups.R, likelihood.R,
# impute.R, tune.R. Base R + testthat (+ ggplot2 class checks) only;
# deterministic, small, fast.

ROOT <- "<root>"

# ---- context_tree.R: network-object grouping ----------------------------

# A bare list with a 2-D $data slot is duck-typed as a netobject.
make_netobj <- function(t1, t2, with_nodes = TRUE) {
  o <- list(data = data.frame(t1 = t1, t2 = t2, stringsAsFactors = FALSE))
  if (with_nodes)
    o$nodes <- data.frame(id = 1:3, label = c("A", "B", "C"),
                          stringsAsFactors = FALSE)
  o
}

test_that("a named list of netobjects fits as a transitiontrees_group", {
  # .ct_is_group() duck-typing arm (return TRUE for a named list whose
  # every element is a netobject), grouped dispatch, and the union-alphabet
  # branch of .ct_group_alphabet().
  no1 <- make_netobj(c("A", "B", "A"), c("B", "C", "B"))
  no2 <- make_netobj(c("C", "A"), c("A", "C"))
  expect_true(transitiontrees:::.ct_is_group(list(g1 = no1, g2 = no2)))

  g <- context_tree(list(g1 = no1, g2 = no2), max_depth = 1L, min_count = 1L)
  expect_s3_class(g, "transitiontrees_group")
  expect_named(g, c("g1", "g2"))
  # shared alphabet is the union of each element's declared node set
  expect_identical(g[["g1"]]$alphabet, c("A", "B", "C"))
})

test_that("an unnamed netobject_group gets default group names", {
  # .ct_group_elements() supplies group1/group2 when names are absent.
  no1 <- make_netobj(c("A", "B"), c("B", "C"))
  no2 <- make_netobj(c("C", "A"), c("A", "C"))
  gg  <- structure(list(no1, no2), class = c("netobject_group", "list"))
  g   <- context_tree(gg, max_depth = 1L, min_count = 1L)
  expect_named(g, c("group1", "group2"))
})

test_that("a grouped object with no declared alphabet derives one", {
  # .ct_group_alphabet() else-branch: elements carry no node/label set, so
  # the union is empty and the alphabet is derived from the sequences.
  no1 <- make_netobj(c("A", "B"), c("B", "A"), with_nodes = FALSE)
  no2 <- make_netobj(c("B", "A"), c("A", "B"), with_nodes = FALSE)
  gg  <- structure(list(no1, no2), class = c("netobject_group", "list"))
  g   <- context_tree(gg, max_depth = 1L, min_count = 1L)
  expect_setequal(g[["group1"]]$alphabet, c("A", "B"))
})

test_that("an explicit alphabet is honoured for a grouped object", {
  # .ct_group_alphabet() early return when 'alphabet' is supplied.
  no1 <- make_netobj(c("A", "B"), c("B", "C"))
  no2 <- make_netobj(c("C", "A"), c("A", "C"))
  g <- context_tree(list(g1 = no1, g2 = no2), alphabet = c("A", "B", "C", "D"),
                    max_depth = 1L, min_count = 1L)
  expect_identical(g[["g1"]]$alphabet, c("A", "B", "C", "D"))
})

# ---- context_tree.R: .ct_group_split_by + .ct_split_weights -------------

make_meta_netobj <- function() {
  list(
    data = data.frame(t1 = c("A", "B", "A", "C"),
                      t2 = c("B", "A", "C", "A"),
                      stringsAsFactors = FALSE),
    metadata = data.frame(grp = c("x", "x", "y", "y"),
                          stringsAsFactors = FALSE))
}

test_that("group= names a metadata column of a network object", {
  no <- make_meta_netobj()
  g  <- context_tree(no, group = "grp", max_depth = 1L, min_count = 1L)
  expect_named(g, c("x", "y"))
})

test_that("group= as a per-sequence vector splits a network object", {
  no <- make_meta_netobj()
  g  <- context_tree(no, group = c("x", "x", "y", "y"),
                     max_depth = 1L, min_count = 1L)
  expect_length(g, 2L)
})

test_that("group= naming an absent metadata column errors", {
  no <- make_meta_netobj()
  expect_error(context_tree(no, group = "nope", max_depth = 1L),
               "not a column of the object's")
})

test_that(".ct_split_weights validates weights in the group path", {
  seqs <- list(c("A", "B", "C"), c("B", "A", "C"),
               c("C", "A", "B"), c("A", "C", "B"))
  grp  <- c("a", "a", "b", "b")
  expect_error(
    context_tree(seqs, group = grp, weights = c("p", "q", "r", "s"),
                 max_depth = 1L, min_count = 1L),
    "must be a numeric vector")
  expect_error(
    context_tree(seqs, group = grp, weights = c(1, 2, NA, 4),
                 max_depth = 1L, min_count = 1L),
    "must not contain NA")
  expect_error(
    context_tree(seqs, group = grp, weights = c(1, 2, -1, 4),
                 max_depth = 1L, min_count = 1L),
    "non-negative")
})

# ---- context_tree.R: small helpers --------------------------------------

test_that(".ct_coerce unwraps a network object to its sequence frame", {
  no <- list(data = data.frame(t1 = c("A", "B"), t2 = c("B", "C"),
                               stringsAsFactors = FALSE))
  out <- transitiontrees:::.ct_coerce(no)
  expect_s3_class(out, "data.frame")
  expect_equal(dim(out), c(2L, 2L))
})

test_that(".ct_kl is Inf when p>0 where q=0, and .ct_g2 is 0 on empty", {
  expect_identical(transitiontrees:::.ct_kl(c(1, 0), c(0, 1)), Inf)
  expect_identical(transitiontrees:::.ct_g2(c(0, 0, 0), c(.5, .3, .2)), 0)
})

# ---- context_tree.R: long-format group/block column checks --------------

test_that("a long-mode group naming an absent column errors", {
  long <- data.frame(id = c("a", "a", "b", "b"),
                     act = c("X", "Y", "Y", "Z"),
                     stringsAsFactors = FALSE)
  expect_error(
    context_tree(long, action = "act", actor = "id", group = "nocol",
                 max_depth = 1L),
    "is not a column of 'data'")
})

# ---- context_tree.R: block length validation in the group path ----------

test_that("block length must match the number of sequences", {
  seqs <- list(c("A", "B", "C"), c("B", "A", "C"),
               c("C", "A", "B"), c("A", "C", "B"))
  expect_error(
    context_tree(seqs, group = c("a", "a", "b", "b"), block = c(1, 2),
                 max_depth = 1L, min_count = 1L),
    "one entry per input sequence")
})

# ---- context_tree.R: coercion / weight validation (single fit) ----------

test_that("an all-empty input errors after coercion", {
  expect_error(
    context_tree(list(c(NA, NA), c("", "")), max_depth = 1L),
    "No usable sequences after coercion")
})

test_that("non-numeric weights are rejected in the single-fit path", {
  expect_error(
    context_tree(list(c("A", "B"), c("B", "C"), c("C", "A")),
                 weights = c("p", "q", "r"), max_depth = 1L, min_count = 1L),
    "must be a numeric vector")
})

# ---- context_tree.R: print / summary methods ----------------------------

test_that("print shows the pruning banner for a pruned tree", {
  set.seed(11)
  m  <- matrix(sample(c("A", "B", "C"), 400, replace = TRUE), nrow = 40)
  tr <- context_tree(m, max_depth = 3L, min_count = 2L)
  pr <- prune_tree(tr, criterion = "G2", alpha = 0.05)
  expect_true(isTRUE(pr$pruned))
  out <- utils::capture.output(print(pr))
  expect_true(any(grepl("pruned by", out)))
})

test_that("print shows the subtree banner for a subtree", {
  set.seed(12)
  m  <- matrix(sample(c("A", "B", "C"), 400, replace = TRUE), nrow = 40)
  tr <- context_tree(m, max_depth = 3L, min_count = 2L)
  st <- subtree(tr, "A")
  out <- utils::capture.output(print(st))
  expect_true(any(grepl("subtree of", out)))

  # local_root == .ROOT prints the root label (.ROOT branch of the banner)
  st2 <- st
  attr(st2, "local_root") <- ROOT
  out2 <- utils::capture.output(print(st2))
  expect_true(any(grepl("subtree of: \\(start\\)", out2)))
})

test_that("print truncates a tree with more lines than max_lines", {
  set.seed(13)
  m  <- matrix(sample(c("A", "B", "C"), 600, replace = TRUE), nrow = 50)
  tr <- context_tree(m, max_depth = 4L, min_count = 1L)
  out <- utils::capture.output(print(tr, max_lines = 2L))
  expect_true(any(grepl("more nodes", out)))
})

test_that("print.summary.transitiontrees truncates extra rows", {
  set.seed(14)
  m  <- matrix(sample(c("A", "B", "C"), 400, replace = TRUE), nrow = 40)
  tr <- context_tree(m, max_depth = 3L, min_count = 2L)
  out <- utils::capture.output(print(summary(tr), n = 2L))
  expect_true(any(grepl("more rows", out)))
})

test_that("witten_bell smoothing prints its bare method label", {
  # .pt_smoothing_label() returns the method name alone when there are no
  # extra hyperparameters.
  set.seed(15)
  m  <- matrix(sample(c("A", "B", "C"), 200, replace = TRUE), nrow = 20)
  tr <- context_tree(m, max_depth = 1L, min_count = 2L,
                     smoothing = "witten_bell")
  out <- utils::capture.output(print(tr))
  expect_true(any(grepl("witten_bell", out)))
})

# ---- context_tree.R: as.data.frame.transitiontrees_group empties --------

test_that("as.data.frame on a group with all-empty members returns schema", {
  seqs <- list(c("A", "B", "C"), c("B", "A", "C"),
               c("C", "A", "B"), c("A", "C", "B"))
  g <- context_tree(seqs, group = c("a", "a", "b", "b"),
                    max_depth = 1L, min_count = 1L)
  # forcing a huge min_count makes every member's pathway table empty
  out <- as.data.frame(g, min_count = 1e6)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 0L)
  expect_true("group" %in% names(out))
})

# ---- bootstrap.R helpers ------------------------------------------------

test_that(".pt_pathway_depth maps labels to depths", {
  d <- transitiontrees:::.pt_pathway_depth(c("(start)", "A -> B", "A"))
  expect_identical(d, c(0L, 2L, 1L))
})

test_that(".pt_resample_count_matrices handles an empty-context depth", {
  # short sequences leave depth 1 with no contexts -> the 0-row branch
  pc <- transitiontrees:::.pt_precompute_for_boot(
    list(c("A"), c("B")), max_depth = 2L, alphabet = c("A", "B"))
  mats <- transitiontrees:::.pt_resample_count_matrices(pc, c(1L, 2L))
  expect_equal(nrow(mats[[2L]]), 0L)
})

test_that(".pt_pathway_stats_from_counts handles ptot=0 and KL=Inf", {
  layout <- list(pw_depth = c(1L, 1L), pw_ctx_idx = c(1L, 2L),
                 parent_depth = c(0L, 0L), parent_ctx_idx = c(1L, 2L))
  cbd <- list(matrix(c(0, 0, 0, 5), nrow = 2, byrow = TRUE),
              matrix(c(2, 0, 3, 0), nrow = 2, byrow = TRUE))
  s <- transitiontrees:::.pt_pathway_stats_from_counts(cbd, layout)
  expect_true(is.na(s$KL[1]))   # parent total 0 -> skipped
  expect_identical(s$KL[2], Inf) # p>0 where parent prob 0 -> Inf
})

test_that("bootstrap_pathways errors when the tree has no pathways", {
  fake <- structure(
    list(nodes = list(), data = list(c("A", "B"), c("B", "A")),
         alphabet = c("A", "B"), max_depth = 0L, weights = NULL),
    class = "transitiontrees")
  expect_error(bootstrap_pathways(fake, iter = 2L),
               "no pathways")
})

test_that("bootstrap_pathways runs with a progress bar", {
  seqs <- replicate(20, sample(c("A", "B", "C"), 8, replace = TRUE),
                    simplify = FALSE)
  tr <- context_tree(seqs, max_depth = 1L, min_count = 2L)
  out <- utils::capture.output(
    boot <- bootstrap_pathways(tr, iter = 3L, progress = TRUE))
  expect_s3_class(boot, "transitiontrees_bootstrap")
})

test_that("plot.transitiontrees_bootstrap uses the default min_stability", {
  seqs <- replicate(40, sample(c("A", "B", "C"), 10, replace = TRUE),
                    simplify = FALSE)
  tr   <- context_tree(seqs, max_depth = 1L, min_count = 2L)
  boot <- bootstrap_pathways(tr, iter = 40L, stability_threshold = 0.01)
  p <- plot(boot)            # min_stability = NULL -> x$stability_threshold
  expect_s3_class(p, "ggplot")
})

test_that("plot.transitiontrees_bootstrap errors when nothing passes", {
  seqs <- replicate(40, sample(c("A", "B", "C"), 10, replace = TRUE),
                    simplify = FALSE)
  tr   <- context_tree(seqs, max_depth = 1L, min_count = 2L)
  boot <- bootstrap_pathways(tr, iter = 40L)
  expect_error(plot(boot, min_stability = 2),
               "No pathways meet")
})

# ---- compare.R ----------------------------------------------------------

test_that("tree_distance returns 0 when total weight is zero", {
  # a shared context whose count is 0 in both trees -> total weight 0
  zero_node <- list("<root>" = list(n = 0, prob = c(0.5, 0.5)))
  t1 <- structure(list(nodes = zero_node, alphabet = c("A", "B")),
                  class = "transitiontrees")
  t2 <- structure(list(nodes = zero_node, alphabet = c("A", "B")),
                  class = "transitiontrees")
  expect_identical(tree_distance(t1, t2), 0)
})

test_that("compare_trees errors when a pruned tree lacks pruning metadata", {
  set.seed(21)
  m1 <- matrix(sample(c("A", "B", "C"), 200, replace = TRUE), nrow = 20)
  m2 <- matrix(sample(c("A", "B", "C"), 200, replace = TRUE), nrow = 20)
  ta <- context_tree(m1, max_depth = 2L, min_count = 2L)
  tb <- context_tree(m2, max_depth = 2L, min_count = 2L)
  ta$pruned <- TRUE; ta$pruning <- NULL   # corrupt the metadata
  expect_error(compare_trees(ta, tb, iter = 3L),
               "missing pruning metadata")
})

# ---- compare_groups.R helpers + branches --------------------------------

test_that(".cg_entropy_bits and .cg_usage_g2 handle empty/zero input", {
  expect_identical(transitiontrees:::.cg_entropy_bits(c(0, 0)), 0)
  expect_identical(transitiontrees:::.cg_usage_g2(c(0, 0), c(0, 0)), 0)
})

test_that(".cg_context_counts substitutes an empty trajectory list", {
  cc <- transitiontrees:::.cg_context_counts(
    list(g1 = list(c("A", "B")), g2 = list()),
    grp_names = c("g1", "g2"), contexts = ROOT, depths = 0L,
    alphabet = c("A", "B"))
  expect_length(cc, 1L)
  expect_equal(unname(rowSums(cc[[1L]])["g2"]), 0)
})

test_that("compare_groups adds the root when member nodes lack it", {
  set.seed(22)
  seqs <- replicate(40, sample(c("A", "B", "C"), 8, replace = TRUE),
                    simplify = FALSE)
  g <- context_tree(seqs, group = rep(c("x", "y"), each = 20),
                    max_depth = 1L, min_count = 1L)
  # strip the root node from each member so the context union omits it
  g[["x"]]$nodes[[ROOT]] <- NULL
  g[["y"]]$nodes[[ROOT]] <- NULL
  cmp <- compare_groups(g, iter = 49L)
  expect_true("(start)" %in% cmp$pathways$pathway)
})

test_that("compare_groups handles a single-context (root-only) union", {
  # huge min_count leaves only the root node in each member tree
  set.seed(23)
  seqs <- replicate(40, sample(c("A", "B", "C"), 6, replace = TRUE),
                    simplify = FALSE)
  g <- context_tree(seqs, group = rep(c("x", "y"), each = 20),
                    max_depth = 1L, min_count = 1000L)
  expect_identical(names(g[["x"]]$nodes), ROOT)
  cmp <- compare_groups(g, iter = 19L)
  expect_s3_class(cmp, "transitiontrees_group_comparison")
  expect_equal(nrow(cmp$pathways), 1L)
})

test_that("print.transitiontrees_group_comparison truncates extra pathways", {
  set.seed(24)
  gx <- replicate(40, sample(c("A", "B", "C"), 8, replace = TRUE,
                             prob = c(.2, .6, .2)), simplify = FALSE)
  gy <- replicate(40, sample(c("A", "B", "C"), 8, replace = TRUE,
                             prob = c(.2, .2, .6)), simplify = FALSE)
  g <- context_tree(c(gx, gy), group = rep(c("x", "y"), each = 40),
                    max_depth = 1L, min_count = 1L)
  cmp <- compare_groups(g, iter = 49L)
  out <- utils::capture.output(print(cmp, n = 1L))
  expect_true(any(grepl("more", out)))
})

test_that("plot divergence errors when all behavioral divergence is zero", {
  # identical groups -> JSD is 0 everywhere -> nothing to plot
  set.seed(25)
  base <- replicate(20, sample(c("A", "B", "C"), 8, replace = TRUE),
                    simplify = FALSE)
  g <- context_tree(c(base, base), group = rep(c("x", "y"), each = 20),
                    max_depth = 1L, min_count = 1L)
  cmp <- compare_groups(g, iter = 19L)
  expect_error(plot(cmp, style = "divergence"),
               "positive behavioral divergence")
})

# ---- likelihood.R -------------------------------------------------------

test_that("in-sample logLik handles a node with zero unique counts", {
  # node "A" is fully absorbed by its child "Z -> A" -> unique counts 0
  tr <- context_tree(list(c("Z", "A", "B"), c("Z", "A", "B")),
                     max_depth = 2L, min_count = 1L)
  expect_equal(sum(transitiontrees:::.pt_unique_counts(tr)[["A"]]), 0)
  ll <- logLik(tr)
  expect_true(is.finite(as.numeric(ll)))
})

test_that("perplexity and model_fit are NA when nothing is scored", {
  tr <- context_tree(list(c("A", "B", "C"), c("B", "C", "A")),
                     max_depth = 1L, min_count = 1L)
  oov <- list(c("Z", "Z"), c("Y", "Y"))   # all out-of-vocabulary
  expect_true(is.na(perplexity(tr, newdata = oov)))
  mf <- model_fit(tr, newdata = oov)
  expect_equal(mf$nobs, 0L)
  expect_true(is.na(mf$logLik))
  expect_true(is.na(mf$AIC))
  expect_true(is.na(mf$perplexity))
})

# ---- impute.R -----------------------------------------------------------

test_that("impute_sequences rejects an unsupported newdata type", {
  tr <- context_tree(list(c("A", "B", "C"), c("B", "C", "A")),
                     max_depth = 1L, min_count = 1L)
  expect_error(impute_sequences(tr, 42),
               "must be a list, matrix, data.frame, or character vector")
})

test_that("impute_sequences restores a missing global RNG seed", {
  tr <- context_tree(list(c("A", "B", "C"), c("B", "C", "A")),
                     max_depth = 1L, min_count = 1L)
  if (exists(".Random.seed", envir = globalenv()))
    rm(".Random.seed", envir = globalenv())
  out <- impute_sequences(tr, list(c("A", NA, "C")), method = "prob",
                          seed = 1L)
  expect_length(out, 1L)
  expect_false(anyNA(out[[1L]]))
})

# ---- tune.R -------------------------------------------------------------

test_that(".pt_smoothing_grid accepts a single explicit spec list", {
  set.seed(31)
  m <- matrix(sample(c("A", "B", "C"), 30 * 8, replace = TRUE), 30, 8)
  res <- tune_tree(m, max_depth = 1L, min_count = 2L,
                   smoothing = list("floor", ymin = 0.001),
                   prune = FALSE, folds = 3L)
  expect_s3_class(res, "transitiontrees_tune")
  expect_true(all(grepl("floor", res$smoothing)))
})

test_that("tune_tree rejects a non-list, non-character smoothing grid", {
  set.seed(32)
  m <- matrix(sample(c("A", "B", "C"), 30 * 8, replace = TRUE), 30, 8)
  expect_error(
    tune_tree(m, max_depth = 1L, smoothing = 5L, prune = FALSE, folds = 3L),
    "must be a character vector or a list")
})
