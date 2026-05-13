# pathtree — Design Conventions

This file is the authoritative reference for `pathtree`'s public API
shape. All v0.1 work follows these rules. Departures must be justified
in this file before being shipped.

## 1. The unifying object

A **pathtree** is a fitted variable-depth context tree:

- a list of **nodes** keyed by pathway string
- a flat **edges** data.frame for traversal
- the alphabet, the fit hyperparameters, and a `pruned` flag
- class: `c("pathtree")`

Everything else in the package is a function that consumes a `pathtree`
and returns either: another `pathtree` (e.g. `prune()`), a tidy
data.frame (the pathway-centric API), a scalar / numeric vector
(scoring), or a graphics object (plotting).

## 2. Naming conventions

### Verbs (functions that act)

`verb_pathtree()` for any non-generic verb that acts on a tree, so
the public surface does not collide with sibling packages
(`tna::prune`, `rpart::prune`, `Nestimate::path_dependence`,
`Nestimate::pathways`):

```r
context_tree()         # fit
prune_pathtree()       # was prune()
smooth_pathtree()      # was smooth()
tune_pathtree()        # CV-based hyperparameter selection
query_pathway()        # look up a specific pathway
compare_pathtrees()    # KL / permutation between two trees
bootstrap_pathways()   # bootstrap pathway uncertainty
```

### Nouns (functions that retrieve / summarise)

Plural for collections, singular for one thing. Nouns that collide
with sibling packages get a `pathtree_` prefix:

```r
pathtree_pathways(tree)     # all pathways, tidy (was pathways())
common_pathways(tree, n)    # top n by frequency
divergent_pathways(tree, n) # top n by KL
sharp_pathways(tree, n)     # top n by predictive sharpness
pathtree_dependence(tree)   # per-context KL diagnostic (was path_dependence())
subtree(tree, pathway)      # extract subtree rooted at a pathway
```

### S3 generics (use the standard ones)

If R has a generic, use it. Never invent a parallel name.

```r
predict(tree, newdata)
simulate(tree, nsim, seed)        # alias of generate_sequences()
logLik(tree, newdata = NULL)      # in-sample if NULL, held-out otherwise
nobs(tree)
AIC(tree)  / BIC(tree)            # via logLik + nobs
print(tree)
summary(tree)
plot(tree, style = ...)
```

`generate_sequences()` is kept as the user-facing alias because
"simulate" suggests stochastic-process language some readers won't
parse the same way; they should be byte-identical.

## 3. Output shapes

### Pathway tables

Always a `data.frame` (not tibble — keep base-R only). Column order:

```
pathway | depth | count | modal_next | prob_next | KL | flips
```

Sorted by the leading numeric column (`count` for
`pathtree_pathways()`, `KL` for `divergent_pathways`, `prob_next` for
`sharp_pathways`). Sort direction is descending unless the column is
a length / cost. The empty case returns a schema-stable 0-row
data.frame with the same columns.

`pathtree_dependence()` extends this with `H_node`, `H_parent`,
`H_drop`, `modal_parent` to support information-theoretic diagnostics.

### Tuning grid

`tune_pathtree()` returns a `data.frame` with one row per grid point:

```
max_depth | nmin | smoothing | prune | logLik | n_scored | perplexity | n_nodes_avg
```

`smoothing` is rendered as a label (`"floor"`, `"floor(ymin=0.001)"`,
`"kneser_ney(discount=0.5)"`, …) so multi-method or
hyperparameter-sweep grids are unambiguous in one column.
`n_nodes_avg` is the mean tree size across folds (the post-prune size
when `prune = TRUE`). Sorted by `perplexity` ascending; the minimum
row is exposed via `attr(grid, "best")`.

### Scoring

Scalar: `logLik(tree, newdata)`, `perplexity(tree, newdata)`.

Per-sequence: `score_sequences(tree, newdata)` returns a data.frame
with `sequence_id, n_scored, log_lik, perplexity` — `n_scored` rather
than `n_obs` because positions whose observed state is outside the
training alphabet are dropped from the score, not counted as obs.

Per-position: `score_positions(tree, newdata)` returns a data.frame
with `sequence_id, position, matched_context, observed,
predicted_prob, log_lik`.

### Comparisons

`compare_pathtrees(tree_a, tree_b, n_perm = 200)` returns a
`pathtree_comparison` object:

```
$pdist        # symmetric (or asymmetric) KL, scalar
$null_dist    # numeric vector, length n_perm
$p_value      # one-sided p
$pathways     # data.frame with per-pathway divergence
              # (columns: pathway, count_a, count_b, KL_ab, KL_ba, sym_KL)
```

## 4. Smoothing

`context_tree(..., smoothing = c("floor", "laplace", "kneser_ney", "witten_bell", "jelinek_mercer"))`

- `"floor"` (default) — current behaviour, `ymin` floor + renormalise
- `"laplace"` — additive-α; `alpha` argument, default 1
- `"kneser_ney"` — back-off along the suffix path with discount `D` and
  continuation probabilities; canonical for PSTs (Begleiter, El-Yaniv
  & Yona 2004 §3); `discount` argument, default 0.75
- `"witten_bell"` — interpolation with novelty-based weight
- `"jelinek_mercer"` — fixed-λ interpolation between order-k and
  order-(k-1); `lambda` argument

Smoothing is applied *once* during `context_tree()`. Re-smoothing a
fitted tree is `smooth_pathtree(tree, smoothing)` — returns a new tree.

## 5. Inputs

Three accepted shapes:

1. **Wide character matrix or data.frame** (rows = sequences, columns
   = time-steps). NAs allowed; trailing NAs are treated as
   end-of-sequence.
2. **List of character vectors** (each element is one ragged
   sequence). NAs allowed inside.
3. **`stslist`** (TraMineR) — extracted via `as.matrix(stslist)`;
   alphabet honoured; weights honoured if present.

Numeric matrices are *not* accepted as input — they collide with
"square numeric matrix = transition matrix" in the network-toolkit
ecosystem. Numeric inputs must be cast to character explicitly.

## 6. S3 method conventions

Every S3 class has at minimum `print`, `summary`, and a `plot` (when
visualisation is meaningful). `summary` returns its own S3 class
named `summary.<class>` so `print(summary(x))` dispatches correctly.

Public S3 classes:

- `pathtree` — fitted tree
- `summary.pathtree` — tidy summary
- `pathtree_tune` — output of `tune_pathtree()`
- `pathtree_comparison` — output of `compare_pathtrees()`
- `pathtree_bootstrap` — output of `bootstrap_pathways()`

## 7. Plotting

`plot.pathtree(tree, style = c("dendrogram", "icicle"))`:

- **`"dendrogram"`** (default) — pure-ggplot2 radial dendrogram. Root
  at the centre, leaves on the outer ring; node size = count,
  fill = state, hairline elbow edges. No optional-package dependency.
- **`"icicle"`** — circular partition (sunburst) via `ggraph` (in
  Suggests). Errors with an install hint if `ggraph`/`tidygraph` are
  missing.

`plot_pathways()` and `plot_divergence()` are pathway-centric views
that complement (don't replace) the tree plot.

The previous `"interactive"` (`collapsibleTree`) style was retired in
favour of the radial dendrogram so the default render works in any
R context (script, knitr, headless) without an htmlwidget dependency.

## 8. Hard rules

- Pure base R + ggplot2 in core. `collapsibleTree` and `ggraph` /
  `tidygraph` only loaded inside the relevant plot branches via
  `requireNamespace()`; they are Suggests, not Imports.
- No Rcpp. No data.table. No tidyverse anywhere.
- All vectorised operations use `vapply` / `tapply` / matrix indexing.
- No for loops in hot paths (training loop in EM is the one exception
  in Nestimate; pathtree has no equivalent).
- Roxygen documentation on every exported function. `@return` clause
  mandatory.
- Tests for every public function in `tests/testthat/`.
- Equivalence tests for every estimator in
  `local_testing_and_equivalence/`, gated by
  `PATHTREE_EQUIV_TESTS=true`, target machine epsilon when an external
  reference exists.

## 9. What `pathtree` is not

- Not a hidden-Markov-model package — that's `seqHMM`, `depmixS4`.
- Not a sequence-distance / clustering package — that's `TraMineR`.
- Not an HONEM / higher-order-network embedding package — that's the
  `Nestimate::build_honem()` / `Saebi 2020` family.
- Not a fixed-order Markov chain package — that's `markovchain`.

`pathtree` is the variable-order PST toolkit. Anything that doesn't
serve that focus belongs in another package.
