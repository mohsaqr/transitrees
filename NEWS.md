# pathtree 0.1.0

Initial CRAN release.

## Core fitting

* `context_tree()` fits a variable-depth pathway tree (prediction
  suffix tree; Ron, Singer & Tishby 1996) from a wide character
  matrix, data.frame, list of character vectors, or TraMineR
  `stslist`.
* Five smoothing schemes are available via the unified `smoothing`
  argument: `"floor"`, `"laplace"`, `"kneser_ney"`, `"witten_bell"`,
  and `"jelinek_mercer"`. Hyperparameters can be passed as
  `list(method, ...)`.
* `prune()` supports four pruning criteria: likelihood-ratio `G2`,
  Kullback-Leibler, `AIC`, and `BIC`.
* `smooth_pathtree()` re-smooths a fitted tree without refitting
  the count tensor.
* `compare_smoothing()` fits the tree under several smoothing schemes
  (all five by default) with `max_depth`/`nmin` held fixed and returns
  a tidy one-row-per-scheme table of `n_nodes` and in-sample
  `perplexity` — a one-call replacement for the manual `lapply()` loop.
* `compare_pruning()` prunes a fitted tree under several criteria (all
  four by default) with `alpha`/`threshold` held fixed and returns a
  tidy one-row-per-criterion table of `n_nodes` and `reduction_pct` —
  the pruning analogue of `compare_smoothing()`.
* `n_nodes()` is a small accessor for the number of contexts in a tree
  (an intuitive `length(tree$nodes)`); returns one count per group for a
  `pathtree_group`.
* `compare_smoothing()` now also accepts an already-fitted `pathtree`:
  it re-smooths the tree under each scheme (topology frozen, no
  re-count) instead of refitting from data — handy for sweeping
  smoothers on a pruned model.
* `model_fit()` bundles the standard fit scalars (`logLik`, `df`,
  `nobs`, `AIC`, `BIC`, `perplexity`) into one tidy row — a one-call
  replacement for `logLik(); nobs(); AIC(); BIC(); perplexity()`. Takes
  optional `newdata` for held-out evaluation and returns one row per
  group for a `pathtree_group`.
* **Plain-English output columns.** The pathway tables, the dependence
  table, the comparison breakdown, and the bootstrap summary now use
  readable column names: `modal_next` -> `likely_next`, `prob_next` ->
  `next_probability`, `KL` -> `divergence`, `flips` ->
  `changes_prediction`; `pathtree_dependence()`'s `H_node`/`H_parent`/
  `H_drop`/`modal_parent` -> `entropy`/`entropy_before`/`entropy_drop`/
  `likely_before`; comparison `KL_ab`/`KL_ba`/`sym_KL` ->
  `divergence_ab`/`divergence_ba`/`divergence_sym`; the bootstrap
  `mean_*`/`ci_*`/`M_*` columns follow suit. Argument *values* keep
  backward-compatible aliases: `sort_by = "KL"` and
  `bootstrap_pathways(stat = "prob_next" / "KL")` still work.
* `context_tree()` gains grouped fits. Pass a grouped family object
  (Nestimate `netobject_group`, tna `group_tna`, or any named list of
  family objects) or a single dataset plus `group =` (a `$metadata`
  column name, or a per-sequence vector) and it fits one tree per
  group over a shared alphabet, returning a new `pathtree_group` (named
  list of trees, with `print` and `as.data.frame` methods).
  `compare_pathtrees()` accepts a 2-group `pathtree_group` directly:
  `compare_pathtrees(context_tree(net_group))`. Previously a grouped
  object was silently mis-fitted into a meaningless tree.
* The root context is now retained even when `nmin` exceeds all
  observed counts, producing a well-formed root-only tree instead
  of an empty object.

## Pathway-centric API

* `pathways()`, `common_pathways()`, `divergent_pathways()`, and
  `sharp_pathways()` rank pathways by frequency, KL from the
  suffix-parent, or modal-flip status.
* `path_dependence()` exposes the per-context KL diagnostic table.
* `query_pathway()`, `pathway_exists()`, and `subtree()` provide
  tree introspection.

## Prediction and scoring

* `predict.pathtree()` returns the next-state probability or top-k
  classification for new partial sequences.
* `simulate.pathtree()` is the standard S3 simulation method; it
  delegates to `generate_sequences()`.
* `logLik()`, `nobs()`, `AIC()`, `BIC()`, `perplexity()`,
  `score_sequences()`, and `score_positions()` provide a complete
  predictive-evaluation toolchain.

## Resampling and comparison

* `tune_pathtree()` performs k-fold cross-validated hyperparameter
  tuning across `max_depth`, `nmin`, smoothing scheme, and prune
  on/off.
* `bootstrap_pathways()` reports `p_stability` (probability the
  pathway statistic falls outside the consistency band) alongside
  `stability_rate`, plus an `informative_rate` based on the per-
  resample `G2` against the chi-square reference cutoff.
* `compare_pathtrees()` runs a permutation test for two-tree
  divergence; permutation refits respect each observed tree's
  hyperparameters and pruning state.
* `pathtree_distance()` computes pairwise KL between two trees and
  aligns probability vectors by alphabet name (so trees with the
  same alphabet in different orders compare correctly).

## Visualisation

* `plot.pathtree()` renders a pure-ggplot2 dendrogram coloured by
  KL divergence with modal-flip ring markers.
* `plot.pathtree()` with `style = "icicle"` renders a static
  `ggraph` partition diagram (suggested dependency).
* `plot_pathways()`, `plot_divergence()`, `plot_pathway_resamples()`,
  and the `plot.pathtree_bootstrap()` /
  `plot.pathtree_comparison()` / `plot.pathtree_tune()` methods
  cover bootstrap, comparison, and tuning diagnostics.

## Validation

* Equivalence-tested at machine precision against the archived
  `PST` package (Gabadinho & Ritschard 2013) across 30 random
  configurations: counts exact, probabilities within 1.11e-16.
