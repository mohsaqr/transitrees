# Tests for compare_groups()

mk_diff_group <- function(seed = 1L) {
  set.seed(seed)
  gx <- replicate(60, sample(c("A", "B", "C"), 8, replace = TRUE,
                             prob = c(.34, .5, .16)), simplify = FALSE)
  gy <- replicate(60, sample(c("A", "B", "C"), 8, replace = TRUE,
                             prob = c(.34, .16, .5)), simplify = FALSE)
  context_tree(c(gx, gy), group = rep(c("x", "y"), each = 60),
               max_depth = 1L, min_count = 3L)
}

mk_same_group <- function(seed = 2L) {
  set.seed(seed)
  ga <- replicate(60, sample(c("A", "B", "C"), 8, replace = TRUE),
                  simplify = FALSE)
  gb <- replicate(60, sample(c("A", "B", "C"), 8, replace = TRUE),
                  simplify = FALSE)
  context_tree(c(ga, gb), group = rep(c("a", "b"), each = 60),
               max_depth = 1L, min_count = 3L)
}

test_that("compare_groups returns the expected object and schema", {
  cmp <- compare_groups(mk_diff_group(), iter = 99L)
  expect_s3_class(cmp, "transitiontrees_group_comparison")
  expect_true(all(c("pathways", "omnibus", "distance_matrix", "groups",
                    "iter", "seed") %in% names(cmp)))
  expect_s3_class(cmp$pathways, "data.frame")
  # dynamic per-group columns present
  expect_true(all(c("count_x", "count_y", "modal_x", "modal_y") %in%
                    names(cmp$pathways)))
  expect_true(all(c("jsd_bits", "jsd_p", "jsd_padj",
                    "usage_g2", "usage_p", "usage_padj", "flips") %in%
                    names(cmp$pathways)))
})

test_that("compare_groups flags a real difference and clears the negative control", {
  diff <- compare_groups(mk_diff_group(), iter = 499L)
  same <- compare_groups(mk_same_group(), iter = 499L)
  # behavioral omnibus: significant for the real difference, not for identical groups
  expect_lt(diff$omnibus$p_value[diff$omnibus$axis == "behavioral"], 0.05)
  expect_gt(same$omnibus$p_value[same$omnibus$axis == "behavioral"], 0.05)
  # the "A" context (designed to differ) has a flip and a small adjusted p
  row_A <- diff$pathways[diff$pathways$pathway == "A", ]
  expect_true(row_A$flips)
  expect_lt(row_A$jsd_padj, 0.05)
})

test_that("compare_groups jsd_bits is non-negative and sorted descending", {
  cmp <- compare_groups(mk_diff_group(), iter = 99L)
  expect_true(all(cmp$pathways$jsd_bits >= 0))
  expect_true(all(diff(cmp$pathways$jsd_bits) <= 1e-9))
})

test_that("compare_groups distance matrix is symmetric with zero diagonal", {
  cmp <- compare_groups(mk_diff_group(), iter = 49L)
  dm <- cmp$distance_matrix
  expect_equal(dim(dm), c(2L, 2L))
  expect_equal(diag(dm), c(x = 0, y = 0))
  expect_equal(dm[1, 2], dm[2, 1])
  expect_gt(dm[1, 2], 0)
})

test_that("compare_groups is reproducible with a fixed seed", {
  a <- compare_groups(mk_diff_group(), iter = 99L, seed = 7L)
  b <- compare_groups(mk_diff_group(), iter = 99L, seed = 7L)
  expect_equal(a$pathways$jsd_p, b$pathways$jsd_p)
  expect_equal(a$omnibus$p_value, b$omnibus$p_value)
})

test_that("compare_groups handles three groups", {
  set.seed(3)
  g <- context_tree(
    replicate(90, sample(c("A", "B", "C"), 8, replace = TRUE),
              simplify = FALSE),
    group = rep(c("p", "q", "r"), each = 30), max_depth = 1L, min_count = 3L)
  cmp <- compare_groups(g, iter = 99L)
  expect_length(cmp$groups, 3L)
  expect_equal(dim(cmp$distance_matrix), c(3L, 3L))
  expect_true(all(c("count_p", "count_q", "count_r") %in%
                    names(cmp$pathways)))
})

test_that("compare_groups errors on a single group", {
  set.seed(4)
  g1 <- context_tree(
    replicate(20, sample(c("A", "B"), 6, replace = TRUE), simplify = FALSE),
    group = rep("only", 20), max_depth = 1L, min_count = 1L)
  expect_error(compare_groups(g1), "at least two groups")
})

test_that("compare_groups S3 methods work", {
  cmp <- compare_groups(mk_diff_group(), iter = 99L)
  expect_output(print(cmp), "transitiontrees_group_comparison")
  expect_output(summary(cmp), "Group comparison")
  expect_s3_class(as.data.frame(cmp), "data.frame")
  expect_s3_class(plot(cmp), "ggplot")
  expect_s3_class(plot(cmp, style = "matrix"), "ggplot")
})

test_that("compare_groups validates block length", {
  cmp_obj <- mk_diff_group()
  expect_error(compare_groups(cmp_obj, iter = 9L, block = c("a", "b")),
               "one entry per sequence")
})

test_that("stratified permutation keeps within-block effects but controls between-block confounds", {
  mk <- function(early_seq, late_seq, n_stu, per, between = FALSE) {
    seqs <- list(); student <- character(); phase <- character()
    for (s in seq_len(n_stu)) {
      sid <- paste0("s", s)
      if (between) {
        ph <- if (s <= n_stu / 2) "early" else "late"
        sq <- if (ph == "early") early_seq else late_seq
        for (k in seq_len(per)) {
          seqs <- c(seqs, list(sq)); student <- c(student, sid)
          phase <- c(phase, ph)
        }
      } else {
        for (k in seq_len(per)) {
          seqs <- c(seqs, list(early_seq)); student <- c(student, sid)
          phase <- c(phase, "early")
        }
        for (k in seq_len(per)) {
          seqs <- c(seqs, list(late_seq)); student <- c(student, sid)
          phase <- c(phase, "late")
        }
      }
    }
    list(seqs = seqs, student = student, phase = phase)
  }
  pval <- function(d, use_block) {
    grp <- context_tree(d$seqs, group = d$phase, max_depth = 1L,
                        min_count = 2L)
    blk <- if (use_block)
      c(d$student[d$phase == "early"], d$student[d$phase == "late"]) else NULL
    cmp <- compare_groups(grp, iter = 499L, block = blk)
    cmp$omnibus$p_value[cmp$omnibus$axis == "behavioral"]
  }
  set.seed(1)
  within <- mk(c("A","B","A","B"), c("A","C","A","C"), 24, 2, between = FALSE)
  betw   <- mk(c("A","B","A","B"), c("A","C","A","C"), 24, 3, between = TRUE)

  # a real within-student effect survives the stratified null
  expect_lt(pval(within, TRUE), 0.05)
  # a pure between-student confound is NOT called significant by the
  # stratified test (it warns that it is degenerate), though the naive
  # (unstratified) test falsely fires
  expect_gt(suppressWarnings(pval(betw, TRUE)),  0.50)
  expect_lt(pval(betw, FALSE), 0.05)
})

test_that("stratified comparison prints its scheme", {
  d_seqs <- replicate(40, sample(c("A","B","C"), 6, replace = TRUE),
                      simplify = FALSE)
  lab <- rep(c("x","y"), each = 20)
  grp <- context_tree(d_seqs, group = lab, max_depth = 1L, min_count = 2L)
  blk <- rep(seq_len(20), times = 2)   # 20 blocks, one x + one y each
  cmp <- compare_groups(grp, iter = 49L, block = blk)
  expect_true(isTRUE(cmp$stratified))
  expect_output(print(cmp), "stratified")
})

test_that("all-singleton blocks warn that stratification is degenerate", {
  set.seed(11)
  seqs <- replicate(40, sample(c("A","B","C"), 6, replace = TRUE),
                    simplify = FALSE)
  grp  <- context_tree(seqs, group = rep(c("x","y"), each = 20),
                       max_depth = 1L, min_count = 2L)
  blk  <- as.character(seq_along(seqs))     # one sequence per block
  expect_warning(compare_groups(grp, iter = 49L, block = blk), "degenerate")
})

test_that("singleton blocks beside mixed blocks do not recycle (sample length-1 fix)", {
  set.seed(11)
  seqs <- replicate(6, sample(c("A","B","C"), 6, replace = TRUE),
                    simplify = FALSE)
  grp  <- context_tree(seqs, group = rep(c("x","y"), each = 3),
                       max_depth = 1L, min_count = 1L)
  # pooled order x,x,x,y,y,y; blocks 1,2 are mixed, 3 and 4 are singletons
  blk  <- c(1, 2, 3, 1, 2, 4)
  expect_silent(compare_groups(grp, iter = 49L, block = blk))
})

test_that("context_tree carries block; compare_groups uses it automatically", {
  set.seed(12)
  long <- data.frame(
    id    = rep(1:40, each = 3),
    step  = rep(1:3, times = 40),
    move  = sample(c("A","B","C"), 120, replace = TRUE),
    grp   = rep(ifelse(1:40 <= 20, "x", "y"), each = 3),
    subj  = rep(rep(1:10, times = 4), each = 3),   # each subj spans both groups
    stringsAsFactors = FALSE)
  g <- context_tree(long, actor = "id", action = "move", order = "step",
                    group = "grp", block = "subj", max_depth = 1L,
                    min_count = 2L)
  expect_s3_class(g, "transitiontrees_group")
  expect_false(is.null(attr(g[[1L]], "block")))   # block lives on member trees
  cmp <- compare_groups(g, iter = 99L)
  expect_true(isTRUE(cmp$stratified))             # picked up the stored block
})

test_that("group=/block= as columns align under time/session sessionization (#1)", {
  long <- data.frame(
    id   = rep(c("u1","u2"), each = 4),
    t    = as.POSIXct("2020-01-01") + c(0,10,5000,5010, 0,10,5000,5010),
    move = c("A","B","A","C","B","A","C","A"),
    grp  = rep(c("x","y"), each = 4),
    subj = rep(c("s1","s2"), each = 4),
    stringsAsFactors = FALSE)
  g <- context_tree(long, actor = "id", time = "t", action = "move",
                    group = "grp", block = "subj",
                    max_depth = 1L, min_count = 1L, time_threshold = 900)
  expect_setequal(names(g), c("x", "y"))         # was empty before the fix
  expect_false(is.null(attr(g[[1L]], "block")))
})

test_that("stored block stays aligned when an empty sequence is dropped (#2)", {
  m   <- rbind(c("A","B","C"), c(NA,NA,NA), c("B","A","B"), c("C","C","A"))
  g   <- context_tree(m, group = c("x","x","y","y"),
                      block = c("s1","s2","s1","s2"),
                      max_depth = 1L, min_count = 1L)
  expect_error(compare_groups(g, iter = 19L), NA)   # no length-mismatch error
})

test_that("usage_g2 is NA for the root and absent from its FDR family (#5)", {
  cmp <- compare_groups(mk_diff_group(), iter = 99L)
  root <- cmp$pathways[cmp$pathways$pathway == "(start)", ]
  expect_true(is.na(root$usage_g2))
  expect_true(is.na(root$usage_padj))
})
