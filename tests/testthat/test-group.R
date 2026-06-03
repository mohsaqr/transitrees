# ---- Tests for grouped fits: context_tree() -> transitrees_group ----
#
# Objects are fabricated (class-stamped) so the suite needs neither
# Nestimate nor tna installed, exactly like the single-object tna tests
# in test-context_tree.R.

.grp_wide <- function() {
  data.frame(
    T1 = c("A","B","A","B","A","B","A","B"),
    T2 = c("B","A","B","A","B","A","B","A"),
    T3 = c("A","A","B","B","A","B","B","A"),
    stringsAsFactors = FALSE)
}
.grp_labels <- c("x","x","x","x","y","y","y","y")

.mk_netobj <- function(wide, meta = NULL) {
  structure(list(data = wide,
                 nodes = data.frame(id = 1:2, label = c("A","B"),
                                    stringsAsFactors = FALSE),
                 metadata = meta),
            class = c("netobject", "cograph_network"))
}

# ---- a grouped object passed directly -----------------------------------

test_that("a netobject_group fits to a transitrees_group", {
  w <- .grp_wide()
  nog <- structure(
    list(x = .mk_netobj(w[1:4, ]), y = .mk_netobj(w[5:8, ])),
    class = "netobject_group")
  g <- context_tree(nog, max_depth = 2L, min_count = 1L)
  expect_s3_class(g, "transitrees_group")
  expect_named(g, c("x", "y"))
  expect_s3_class(g$x, "transitrees")
  expect_s3_class(g$y, "transitrees")
})

test_that("a fabricated group_tna (integer-coded) decodes and fits", {
  mk_tna <- function(m) structure(list(labels = c("A","B"), data = m),
                                  class = "tna")
  gtna <- structure(
    list(g1 = mk_tna(matrix(c(1L,2L,1L, 2L,1L,2L, 1L,1L,2L),
                            nrow = 3, byrow = TRUE)),
         g2 = mk_tna(matrix(c(2L,1L,2L, 1L,2L,1L, 2L,2L,1L),
                            nrow = 3, byrow = TRUE))),
    class = "group_tna")
  g <- context_tree(gtna, max_depth = 1L, min_count = 1L)
  expect_s3_class(g, "transitrees_group")
  expect_named(g, c("g1", "g2"))
  expect_setequal(g$g1$alphabet, c("A","B"))   # decoded, not "1"/"2"
})

# ---- group = on a single object -----------------------------------------

test_that("group = a metadata column splits a netobject", {
  no <- .mk_netobj(.grp_wide(),
                   meta = data.frame(grp = .grp_labels,
                                     stringsAsFactors = FALSE))
  g <- context_tree(no, max_depth = 2L, min_count = 1L, group = "grp")
  expect_s3_class(g, "transitrees_group")
  expect_named(g, c("x", "y"))
  expect_identical(attr(g, "group"), "grp")
})

test_that("group = a vector splits a wide frame, shared alphabet", {
  g <- context_tree(.grp_wide(), max_depth = 1L, min_count = 1L,
                    group = .grp_labels)
  expect_s3_class(g, "transitrees_group")
  expect_setequal(g$x$alphabet, g$y$alphabet)   # union alphabet shared
  expect_setequal(g$x$alphabet, c("A","B"))
})

# ---- error handling ------------------------------------------------------

test_that("grouped object plus group = is rejected", {
  nog <- structure(list(x = .mk_netobj(.grp_wide()[1:4, ]),
                        y = .mk_netobj(.grp_wide()[5:8, ])),
                   class = "netobject_group")
  expect_error(context_tree(nog, group = "grp"),
               "either a grouped object or 'group ='")
})

test_that("a length-mismatched group vector errors", {
  expect_error(context_tree(.grp_wide(), group = c("x","y")),
               "one entry per input sequence")
})

test_that("a string group on a plain frame errors with guidance", {
  expect_error(context_tree(.grp_wide(), group = "grp"),
               "only valid for a network object")
})

test_that("a plain ragged list is NOT treated as a group", {
  lst <- list(c("A","B","A","B"), c("B","A","B","A"))
  out <- context_tree(lst, max_depth = 1L, min_count = 1L)
  expect_s3_class(out, "transitrees")          # single tree, not a group
  expect_false(inherits(out, "transitrees_group"))
})

# ---- transitrees_group methods ---------------------------------------------

test_that("as.data.frame.transitrees_group tags rows with group", {
  g <- context_tree(.grp_wide(), max_depth = 1L, min_count = 1L,
                    group = .grp_labels)
  df <- as.data.frame(g)
  expect_s3_class(df, "data.frame")
  expect_true("group" %in% names(df))
  expect_setequal(unique(df$group), c("x", "y"))
})

test_that("print.transitrees_group returns invisibly", {
  g <- context_tree(.grp_wide(), max_depth = 1L, min_count = 1L,
                    group = .grp_labels)
  expect_invisible(print(g))
})

# ---- compare_trees() accepts a 2-group transitrees_group ---------------

test_that("compare_trees(group) compares the pair", {
  g   <- context_tree(.grp_wide(), max_depth = 1L, min_count = 1L,
                      group = .grp_labels)
  cmp <- compare_trees(g, iter = 20L, seed = 1L)
  expect_s3_class(cmp, "transitrees_comparison")
  ## identical to the explicit two-tree call
  cmp2 <- compare_trees(g$x, g$y, iter = 20L, seed = 1L)
  expect_equal(cmp$pdist, cmp2$pdist)
})

test_that("compare_trees rejects a non-pairwise group", {
  one <- structure(list(only = context_tree(.grp_wide(), max_depth = 1L,
                                            min_count = 1L)),
                   class = c("transitrees_group", "list"))
  expect_error(compare_trees(one), "exactly 2 groups")
})

test_that("compare_trees rejects a group plus tree_b", {
  g <- context_tree(.grp_wide(), max_depth = 1L, min_count = 1L,
                    group = .grp_labels)
  expect_error(compare_trees(g, g$x), "not a group plus")
})

test_that("compare_trees still requires tree_b for plain trees", {
  tr <- context_tree(.grp_wide(), max_depth = 1L, min_count = 1L)
  expect_error(compare_trees(tr), "'tree_b' is required")
})

# ---- prune / smooth dispatch over a group -------------------------------

test_that("prune_tree on a group prunes each, keeps the wrapper", {
  g  <- context_tree(.grp_wide(), max_depth = 2L, min_count = 1L,
                     group = .grp_labels)
  pg <- prune_tree(g, criterion = "G2", alpha = 0.05)
  expect_s3_class(pg, "transitrees_group")
  expect_named(pg, names(g))
  expect_true(all(vapply(pg, function(t) isTRUE(t$pruned), logical(1))))
  ## identical to pruning each member by hand
  expect_equal(pg$x$nodes,
               prune_tree(g$x, criterion = "G2", alpha = 0.05)$nodes)
})

test_that("smooth_tree on a group re-smooths each, keeps the wrapper", {
  g  <- context_tree(.grp_wide(), max_depth = 2L, min_count = 1L,
                     group = .grp_labels)
  sg <- smooth_tree(g, "kneser_ney")
  expect_s3_class(sg, "transitrees_group")
  expect_named(sg, names(g))
  expect_equal(sg$x$smoothing$method, "kneser_ney")
})

# ---- regression: grouped fits preserve alphabet + weights ---------------

test_that("group = preserves a netobject's declared alphabet (absent state)", {
  ## $nodes declares a state ("C") that never appears in the observed
  ## rows; a single fit keeps it, so a grouped split must too.
  no <- structure(list(
    data     = .grp_wide(),
    nodes    = data.frame(id = 1:3, label = c("A","B","C"),
                          stringsAsFactors = FALSE),
    metadata = data.frame(grp = .grp_labels, stringsAsFactors = FALSE)),
    class = c("netobject", "cograph_network"))
  single  <- context_tree(no, max_depth = 1L, min_count = 1L)
  grouped <- context_tree(no, max_depth = 1L, min_count = 1L, group = "grp")
  expect_true("C" %in% single$alphabet)
  expect_true("C" %in% grouped$x$alphabet)
  expect_true("C" %in% grouped$y$alphabet)
  expect_setequal(grouped$x$alphabet, single$alphabet)
})

test_that("group = splits and forwards weights to each per-group fit", {
  w   <- .grp_wide()
  grp <- .grp_labels
  wt  <- rep(1, length(grp)); wt[1] <- 4L          # up-weight one x-row
  g_un <- context_tree(w, max_depth = 1L, min_count = 1L, group = grp)
  g_w  <- context_tree(w, max_depth = 1L, min_count = 1L, group = grp,
                       weights = wt)
  ## per-group result equals a direct weighted fit of that group's rows
  direct_x <- context_tree(w[grp == "x", ], max_depth = 1L, min_count = 1L,
                           weights = wt[grp == "x"])
  expect_equal(g_w$x$nodes, direct_x$nodes)
  ## weighting actually changed the x fit, and left y untouched
  expect_false(identical(g_w$x$nodes[["<root>"]]$n,
                         g_un$x$nodes[["<root>"]]$n))
  expect_equal(g_w$y$nodes, g_un$y$nodes)
})

test_that("group = validates the weights length against total sequences", {
  expect_error(
    context_tree(.grp_wide(), max_depth = 1L, min_count = 1L,
                 group = .grp_labels, weights = c(1, 2, 3)),
    "length equal to number of input")
})

test_that("weights with a grouped object is rejected, not silently dropped", {
  nog <- structure(list(x = .mk_netobj(.grp_wide()[1:4, ]),
                        y = .mk_netobj(.grp_wide()[5:8, ])),
                   class = "netobject_group")
  expect_error(context_tree(nog, weights = rep(1, 8)),
               "not supported with a grouped object")
})
