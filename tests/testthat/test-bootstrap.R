# Tests for bootstrap_pathways(): correctness, stability properties,
# reproducibility, and validation against simulated structure.

mk_boot_tree <- function(n = 40L, L = 12L, alpha = c("A","B","C"),
                         seed = 1L, max_depth = 2L, min_count = 3L) {
  set.seed(seed)
  m <- matrix(sample(alpha, n * L, replace = TRUE), n, L)
  context_tree(m, max_depth = max_depth, min_count = min_count)
}

test_that("bootstrap_pathways returns a transitrees_bootstrap with the expected slots", {
  tr <- mk_boot_tree()
  b  <- bootstrap_pathways(tr, iter = 20L, seed = 42L)
  expect_s3_class(b, "transitrees_bootstrap")
  expect_true(all(c("summary", "pathways_orig", "M_count",
                    "M_next_probability", "M_divergence",
                    "M_changes_prediction",
                    "iter", "stat", "consistency_range",
                    "stability_threshold") %in% names(b)))
  expect_equal(b$iter, 20L)
  ## Lock the loosened defaults (wide band, high stability rate,
  ## relaxed informativeness gate).
  expect_equal(b$consistency_range, c(0.5, 1.5))
  expect_equal(b$stability_threshold, 0.95)
  expect_equal(b$informative_threshold, 0.80)
  expect_equal(nrow(b$M_count), 20L)
  expect_equal(nrow(b$M_next_probability), 20L)
  expect_equal(nrow(b$M_divergence), 20L)
  expect_equal(nrow(b$M_changes_prediction), 20L)
  expect_equal(ncol(b$M_count), nrow(b$summary))
})

test_that("summary table has the canonical column vocabulary", {
  tr <- mk_boot_tree()
  b  <- bootstrap_pathways(tr, iter = 20L, seed = 42L)
  expect_true(all(c("pathway", "depth", "count", "next_probability",
                    "divergence", "G2",
                    "p_stability", "stability_rate", "stable",
                    "informative_rate", "informative",
                    "flip_consistency",
                    "mean_count", "ci_count_lo", "ci_count_hi",
                    "mean_divergence", "ci_divergence_lo", "ci_divergence_hi",
                    "mean_G2", "ci_G2_lo", "ci_G2_hi") %in%
                  names(b$summary)))
  ## Generic p_value, sig, appearance_rate were removed.
  expect_false("p_value"         %in% names(b$summary))
  expect_false("sig"             %in% names(b$summary))
  expect_false("appearance_rate" %in% names(b$summary))
})

test_that("summary table has symmetric sd_* columns (one per stat)", {
  tr <- mk_boot_tree()
  b  <- bootstrap_pathways(tr, iter = 30L, seed = 42L)
  expect_true(all(c("sd_count", "sd_next_probability", "sd_divergence",
                    "sd_G2") %in% names(b$summary)))
  ## SDs are non-negative or NA
  for (col in c("sd_count", "sd_next_probability", "sd_divergence",
                "sd_G2")) {
    v <- b$summary[[col]]
    v <- v[!is.na(v)]
    expect_true(all(v >= 0), info = col)
  }
})

test_that("keep_resamples = FALSE drops the M_* matrices but keeps the summary", {
  tr <- mk_boot_tree()
  b  <- bootstrap_pathways(tr, iter = 30L, seed = 1L,
                           keep_resamples = FALSE)
  expect_true(is.null(b$M_count))
  expect_true(is.null(b$M_next_probability))
  expect_true(is.null(b$M_divergence))
  expect_true(is.null(b$M_G2))
  expect_true(is.null(b$M_changes_prediction))
  expect_s3_class(b$summary, "data.frame")
  expect_true(nrow(b$summary) > 0L)
  ## plot_pathway_resamples errors helpfully when resamples were dropped
  expect_error(plot_pathway_resamples(b),
               "keep_resamples = TRUE")
})

test_that("summary + pathways_orig both lead with the canonical pathway schema", {
  tr <- mk_boot_tree()
  b  <- bootstrap_pathways(tr, iter = 25L, seed = 1L)
  canonical <- c("pathway", "depth", "count",
                 "likely_next", "next_probability", "divergence",
                 "changes_prediction")
  ## Both tables expose the same canonical leading columns.
  expect_equal(names(b$summary)[seq_along(canonical)], canonical)
  expect_equal(names(b$pathways_orig)[seq_along(canonical)], canonical)
  ## pathways_orig carries G2 right after the canonical block.
  expect_equal(names(b$pathways_orig)[length(canonical) + 1L], "G2")
})

test_that("keep_resamples = TRUE (default) retains M_* matrices", {
  tr <- mk_boot_tree()
  b  <- bootstrap_pathways(tr, iter = 30L, seed = 1L)
  expect_true(is.matrix(b$M_count))
  expect_true(is.matrix(b$M_next_probability))
  expect_true(is.matrix(b$M_divergence))
  expect_true(is.matrix(b$M_G2))
  expect_true(is.matrix(b$M_changes_prediction))
})

test_that("p_stability and stability_rate are in [0, 1]", {
  tr <- mk_boot_tree()
  b  <- bootstrap_pathways(tr, iter = 50L, seed = 1L)
  sr <- b$summary$stability_rate
  sr <- sr[!is.na(sr)]
  expect_true(all(sr >= 0))
  expect_true(all(sr <= 1))
  ps <- b$summary$p_stability
  ps <- ps[!is.na(ps)]
  expect_true(all(ps >= 0))
  expect_true(all(ps <= 1))
})

test_that("p_stability matches the corrected outside-band count", {
  tr <- mk_boot_tree()
  b  <- bootstrap_pathways(tr, iter = 50L, seed = 1L,
                           stat = "count")
  s <- b$summary
  for (j in seq_len(nrow(s))) {
    pathway <- s$pathway[j]
    observed <- s$count[j]
    lo <- min(observed * b$consistency_range)
    hi <- max(observed * b$consistency_range)
    boot_values <- b$M_count[, pathway]
    n_out <- sum(boot_values < lo | boot_values > hi, na.rm = TRUE)
    expect_equal(s$p_stability[j], (n_out + 1) / (b$iter + 1),
                 tolerance = 1e-12)
  }
})

test_that("stricter stability threshold cannot increase stable pathways", {
  tr <- mk_boot_tree()
  loose <- bootstrap_pathways(tr, iter = 50L, seed = 1L,
                              stability_threshold = 0.60)
  strict <- bootstrap_pathways(tr, iter = 50L, seed = 1L,
                               stability_threshold = 0.90)
  expect_lte(sum(strict$summary$stable, na.rm = TRUE),
             sum(loose$summary$stable, na.rm = TRUE))
})

test_that("seed makes the bootstrap reproducible", {
  tr <- mk_boot_tree()
  b1 <- bootstrap_pathways(tr, iter = 30L, seed = 7L)
  b2 <- bootstrap_pathways(tr, iter = 30L, seed = 7L)
  expect_equal(b1$summary, b2$summary)
  expect_equal(b1$M_count, b2$M_count)
})

test_that("different seeds give different resamples", {
  tr <- mk_boot_tree()
  b1 <- bootstrap_pathways(tr, iter = 30L, seed = 1L)
  b2 <- bootstrap_pathways(tr, iter = 30L, seed = 2L)
  expect_false(identical(b1$M_count, b2$M_count))
})

test_that("alphabet is locked across resamples", {
  ## Even when a small bootstrap sample misses a state, the resampled
  ## tree's alphabet must match the input tree's, so all resamples
  ## align to the same union of pathways.
  tr <- mk_boot_tree(n = 30L, L = 8L, alpha = c("A","B","C","D"),
                     min_count = 2L)
  b  <- bootstrap_pathways(tr, iter = 30L, seed = 1L)
  expect_true(all(c("A","B","C","D") %in% tr$alphabet))
  ## All resamples should be able to produce pathways over the locked
  ## alphabet; we verify by checking every column's name parses to
  ## states all in the original alphabet.
  parse_states <- function(p) {
    if (p == "(start)") return(character(0))
    strsplit(p, " -> ", fixed = TRUE)[[1]]
  }
  used <- unique(unlist(lapply(colnames(b$M_count), parse_states)))
  expect_true(all(used %in% tr$alphabet))
})

test_that("missing pathway in a resample contributes count = 0 / NA elsewhere", {
  tr <- mk_boot_tree(n = 25L, L = 10L, max_depth = 3L, min_count = 4L)
  b  <- bootstrap_pathways(tr, iter = 40L, seed = 1L)
  ## Find a pathway that did NOT appear in every resample.
  rates <- colMeans(b$M_count > 0)
  partial <- which(rates < 1 & rates > 0)
  if (length(partial) > 0L) {
    j <- partial[[1L]]
    miss_rows <- which(b$M_count[, j] == 0)
    expect_true(all(is.na(b$M_next_probability[miss_rows, j])))
    expect_true(all(is.na(b$M_divergence[miss_rows, j])))
    expect_true(all(is.na(b$M_changes_prediction[miss_rows, j])))
  } else {
    succeed()
  }
})

test_that("validation: simulated deterministic next-state is informative", {
  ## Construct sequences where context "A" deterministically predicts
  ## "B": stability_rate near 1 AND informative_rate near 1.
  set.seed(42)
  n <- 80; L <- 14
  trajs <- replicate(n, {
    out <- character(L)
    out[1] <- sample(c("A","B","C"), 1)
    for (t in 2:L) {
      out[t] <- if (out[t - 1] == "A") "B" else
                  sample(c("A","B","C"), 1)
    }
    out
  }, simplify = FALSE)
  tr <- context_tree(trajs, max_depth = 2L, min_count = 3L,
                     smoothing = list("floor", ymin = 0))
  b  <- bootstrap_pathways(tr, iter = 80L, seed = 1L)
  row_A <- b$summary[b$summary$pathway == "A", ]
  expect_equal(nrow(row_A), 1L)
  expect_lt(row_A$p_stability,    0.25)
  expect_gt(row_A$stability_rate,   0.75)
  expect_gt(row_A$informative_rate, 0.95)
  expect_true(row_A$informative)
  ## CI on next_probability should sit near 1
  expect_gt(row_A$ci_next_probability_lo, 0.7)
})

test_that("validation: random data leaves few pathways informative", {
  set.seed(7)
  n <- 80; L <- 14
  trajs <- replicate(n,
    sample(c("A","B","C"), L, replace = TRUE),
    simplify = FALSE)
  tr <- context_tree(trajs, max_depth = 1L, min_count = 3L)
  b  <- bootstrap_pathways(tr, iter = 100L, seed = 1L)
  ## Random data should mostly be stable on count (counts reproduce)
  ## but very few pathways should be informative.
  s <- b$summary[b$summary$pathway != "(start)", ]
  expect_true(mean(s$informative, na.rm = TRUE) <= 0.5)
})

test_that("stable frequent pathways in iid data are not called informative", {
  set.seed(2)
  trajs <- replicate(200L,
    sample(c("A", "B", "C"), 30L, replace = TRUE),
    simplify = FALSE)
  tr <- context_tree(trajs, max_depth = 1L, min_count = 1L,
                     smoothing = list("floor", ymin = 0))
  b <- bootstrap_pathways(tr, iter = 100L, seed = 1L,
                          stat = "count")
  s <- b$summary[b$summary$pathway != "(start)", ]
  expect_true(all(s$stable))
  expect_true(all(!s$informative))
  expect_true(all(s$informative_rate < b$informative_threshold))
})

test_that("validation: count CI for high-frequency pathway brackets the original count", {
  tr <- mk_boot_tree(n = 80L, L = 14L, max_depth = 1L, min_count = 5L)
  b  <- bootstrap_pathways(tr, iter = 100L, seed = 1L)
  ## Look at the highest-count non-root pathway
  s <- b$summary[b$summary$pathway != "(start)", ]
  s <- s[order(-s$count), ]
  top <- s[1L, ]
  ## CI should bracket the original count (95% level → coverage ~95%)
  expect_lte(top$ci_count_lo, top$count)
  expect_gte(top$ci_count_hi, top$count)
})

test_that("stable flag follows p_stability threshold", {
  tr <- mk_boot_tree()
  b  <- bootstrap_pathways(tr, iter = 50L, seed = 1L,
                            stability_threshold = 0.6)
  s <- b$summary
  expect_true(all(s$p_stability[!is.na(s$p_stability)] >= 0))
  expect_true(all(s$p_stability[!is.na(s$p_stability)] <= 1))
  expect_true(all(s$stability_rate[!is.na(s$stability_rate)] >= 0))
  expect_true(all(s$stability_rate[!is.na(s$stability_rate)] <= 1))
  ## stable flag must match the documented p_stability criterion.
  expected_stable <- !is.na(s$p_stability) & s$p_stability < 0.4
  expect_equal(s$stable, expected_stable)
})

test_that("simulation guards against stability false negatives and false positives", {
  ## Most sequences repeatedly carry A -> B, so the pathway count is
  ## almost unchanged under sequence bootstrap. A few high-leverage
  ## carrier sequences repeatedly carry X -> Y, making X genuinely
  ## informative but unstable at the sample level.
  stable <- replicate(75L, rep(c("A", "B"), 20L), simplify = FALSE)
  carriers <- replicate(5L, rep(c("X", "Y"), 20L), simplify = FALSE)
  tr <- context_tree(c(stable, carriers), max_depth = 1L, min_count = 1L,
                     smoothing = list("floor", ymin = 0))
  ## Pin the tolerance band so this mechanism test (stable vs
  ## carrier-driven) is decoupled from the package default band.
  b <- bootstrap_pathways(tr, iter = 200L, seed = 1L, stat = "count",
                          consistency_range = c(0.75, 1.25))

  row_A <- b$summary[b$summary$pathway == "A", ]
  row_X <- b$summary[b$summary$pathway == "X", ]
  expect_equal(nrow(row_A), 1L)
  expect_equal(nrow(row_X), 1L)

  ## False negative guard: a genuinely stable, common pathway should
  ## not be missed.
  expect_lt(row_A$p_stability, 0.05)
  expect_true(row_A$stable)
  expect_gt(row_A$informative_rate, 0.95)

  ## False positive guard: a deterministic but carrier-driven pathway
  ## should not be called stable just because its next state is sharp.
  expect_gt(row_X$p_stability, 0.25)
  expect_false(row_X$stable)
  expect_gt(row_X$informative_rate, 0.95)
  expect_true(row_X$informative)
})

test_that("stat = 'next_probability' flags probability instability when counts are stable", {
  ## Every sequence contributes the same number of A contexts, so count
  ## is stable. But whole-sequence resampling changes the mix of
  ## A -> B carrier sequences and A -> C carrier sequences, so the
  ## modal next-state probability is unstable.
  b_carriers <- replicate(5L, rep(c("A", "B"), 20L), simplify = FALSE)
  c_carriers <- replicate(5L, rep(c("A", "C"), 20L), simplify = FALSE)
  tr <- context_tree(c(b_carriers, c_carriers),
                     max_depth = 1L, min_count = 1L,
                     smoothing = list("floor", ymin = 0))
  ## Pin the tolerance band so the count-stable / prob-unstable
  ## contrast is decoupled from the package default band.
  b_count <- bootstrap_pathways(tr, iter = 200L, seed = 1L,
                                stat = "count",
                                consistency_range = c(0.75, 1.25))
  b_prob <- bootstrap_pathways(tr, iter = 200L, seed = 1L,
                               stat = "next_probability",
                               consistency_range = c(0.75, 1.25))

  row_count <- b_count$summary[b_count$summary$pathway == "A", ]
  row_prob <- b_prob$summary[b_prob$summary$pathway == "A", ]
  expect_true(row_count$stable)
  expect_lt(row_count$p_stability, 0.05)
  expect_false(row_prob$stable)
  expect_gt(row_prob$p_stability, 0.25)
})

test_that("stat = 'divergence' makes near-zero divergence instability explicit", {
  set.seed(2)
  trajs <- replicate(200L,
    sample(c("A", "B", "C"), 30L, replace = TRUE),
    simplify = FALSE)
  tr <- context_tree(trajs, max_depth = 1L, min_count = 1L,
                     smoothing = list("floor", ymin = 0))
  b <- bootstrap_pathways(tr, iter = 100L, seed = 1L, stat = "divergence")
  s <- b$summary[b$summary$pathway != "(start)", ]
  expect_true(all(is.finite(s$divergence)))
  expect_true(all(s$divergence < 0.01))
  expect_true(all(s$p_stability > 0.25))
  expect_true(all(!s$stable))
})

test_that("changing stat changes which statistic stability is computed on", {
  tr <- mk_boot_tree()
  b_count <- bootstrap_pathways(tr, iter = 30L, seed = 1L,
                                 stat = "count")
  b_prob  <- bootstrap_pathways(tr, iter = 30L, seed = 1L,
                                 stat = "next_probability")
  expect_equal(b_count$stat, "count")
  expect_equal(b_prob$stat,  "next_probability")
  ## Different stats produce different stability_rates in general.
  ## (Same per-pathway count/prob distributions but different
  ## tolerance bands → different rates.)
  expect_false(identical(b_count$summary$stability_rate,
                         b_prob$summary$stability_rate))
})

test_that("as.data.frame(boot) returns the summary table (uniform tidy-extract)", {
  tr <- mk_boot_tree()
  b  <- bootstrap_pathways(tr, iter = 20L, seed = 1L)
  df <- as.data.frame(b)
  expect_s3_class(df, "data.frame")
  expect_identical(df, b$summary)
})

test_that("plot.transitrees_bootstrap returns a ggplot", {
  tr <- mk_boot_tree()
  b  <- bootstrap_pathways(tr, iter = 30L, seed = 1L)
  p  <- plot(b, min_stability = 0)
  expect_s3_class(p, "ggplot")
})

test_that("input tree without $data raises an informative error", {
  tr <- mk_boot_tree()
  tr$data <- NULL
  expect_error(bootstrap_pathways(tr, iter = 5L),
               "tree\\$data is missing")
})

test_that("validation errors on bad arguments", {
  tr <- mk_boot_tree()
  expect_error(bootstrap_pathways(tr, iter = 1L), "iter")
  expect_error(bootstrap_pathways(tr, iter = 5L,
                                   stability_threshold = 0))
  expect_error(bootstrap_pathways(tr, iter = 5L,
                                   consistency_range = c(1.5, 0.5)))
})

test_that("plot_pathway_resamples renders a ggplot for each stat", {
  tr <- mk_boot_tree()
  b  <- bootstrap_pathways(tr, iter = 30L, seed = 1L)  # keep_resamples = TRUE
  for (st in c("count", "next_probability", "divergence", "G2"))
    expect_s3_class(plot_pathway_resamples(b, stat = st, top = 3L), "ggplot")
})

test_that("plot_pathway_resamples errors on an unknown pathway", {
  tr <- mk_boot_tree()
  b  <- bootstrap_pathways(tr, iter = 20L, seed = 1L)
  expect_error(plot_pathway_resamples(b, pathways = "Z -> Z -> Z"),
               "Unknown pathway")
})
