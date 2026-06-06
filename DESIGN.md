# transitiontrees — Design Conventions

This file is the authoritative reference for `transitiontrees`'s public API
shape. All v0.1 work follows these rules. Departures must be justified
in this file before being shipped.

## 1. The unifying object

A **transitiontrees** is a fitted variable-depth context tree:

- a list of **nodes** keyed by pathway string
- a flat **edges** data.frame for traversal
- the alphabet, the fit hyperparameters, and a `pruned` flag
- class: `c("transitiontrees")`

Everything else in the package is a function that consumes a `transitiontrees`
and returns either: another `transitiontrees` (e.g. `prune()`), a tidy
data.frame (the pathway-centric API), a scalar / numeric vector
(scoring), or a graphics object (plotting).

## 2. Naming conventions

### Verbs (functions that act)

`verb_transitiontrees()` for any non-generic verb that acts on a tree, so
the public surface does not collide with sibling packages
(`tna::prune`, `rpart::prune`, `Nestimate::path_dependence`,
`Nestimate::pathways`):

```r
context_tree()         # fit
prune_tree()       # was prune()
compare_pruning()      # prune under N criteria, tidy size/reduction table
smooth_tree()      # was smooth()
compare_smoothing()    # fit under N smoothers, tidy size/perplexity table
tune_tree()        # CV-based hyperparameter selection
query_pathway()        # look up a specific pathway
tree_distance()    # scalar (a)symmetric KL between two trees
compare_trees()    # KL / permutation between two trees
bootstrap_pathways()   # bootstrap pathway uncertainty
```

### Nouns (functions that retrieve / summarise)

Plural for collections, singular for one thing. Nouns that collide
with sibling packages get a `transitiontrees_` prefix:

```r
tree_pathways(tree)       # all pathways, tidy (was pathways())
common_pathways(tree, top)    # top-`top` by frequency
divergent_pathways(tree, top) # top-`top` by KL
sharp_pathways(tree, top)     # top-`top` by predictive sharpness
tree_dependence(tree)     # per-context KL diagnostic (was path_dependence())
subtree(tree, pathway)        # extract subtree rooted at a pathway
n_nodes(tree)                 # context (node) count; accessor for length(tree$nodes)
```

### Argument vocabulary (one concept = one name)

Argument names are consistent across the package and aligned with the
sibling package **Nestimate** wherever an equivalent exists. The same
concept always gets the same argument name:

| Concept | Argument | Notes / Nestimate match |
|---|---|---|
| The input sequence data | `data` | matches `build_network(data=)` |
| A fitted tree (non-generic verbs) | `tree` | the transitiontrees object |
| A fitted tree (S3 generics) | `x` / `object` | required by `predict`/`logLik`/`plot`/`print` generics |
| Held-out data to score | `newdata` | standard R scoring convention |
| Minimum observation count | `min_count` | matches Nestimate `min_count` |
| How many ranked rows to return | `top` | on `common_/divergent_/sharp_pathways()` |
| Resamples / permutations | `iter` | matches Nestimate `iter` |
| CV folds | `folds` | on `tune_tree()` |
| Significance level | `alpha` | matches Nestimate `alpha` |
| CI tail probability | `ci_level` | matches Nestimate `ci_level` |
| RNG seed | `seed` | matches Nestimate `seed` |
| Grouping | `group` | matches Nestimate `group` |
| Tree depth / context order | `max_depth` | the depth cap |

`generate_sequences(tree, n =, length =)` keeps `n` for the **number of
sequences to generate** — a generation count, not a "top-`top`" of a
ranked table — and `simulate()` keeps the stats-generic `nsim`.

These names are the only accepted spelling — there are no deprecated
aliases. The stored object **fields** use the internal names (e.g.
`tree$nmin`, `bootstrap$alpha_g2`, the tuning-grid `nmin` column); the
object contract is independent of the user-facing argument vocabulary.

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
pathway | depth | count | likely_next | next_probability | divergence | changes_prediction
```

Plain-English column names (renamed from the original
`modal_next`/`prob_next`/`KL`/`flips` for readability). `divergence` is
the Kullback-Leibler divergence from the parent context's prediction;
`changes_prediction` is the modal-flip flag. Sorted by the leading
numeric column (`count` for `tree_pathways()`, `divergence` for
`divergent_pathways`, `next_probability` for `sharp_pathways`). Sort
direction is descending unless the column is a length / cost. The empty
case returns a schema-stable 0-row data.frame with the same columns.

Argument *values* use the canonical names only: `sort_by` takes
`c("count", "divergence", "depth")`; `bootstrap_pathways(stat=)` takes
`c("count", "next_probability", "divergence")` (plus `"G2"` for
`plot_pathway_resamples()`).

`tree_dependence()` extends this with `entropy` (Shannon entropy of
the next-state distribution), `entropy_before` (the parent's entropy),
`entropy_drop` (`entropy_before - entropy`), and `likely_before` (the
parent's most likely next state) to support information-theoretic
diagnostics. `compare_trees()`'s breakdown uses `count_a`,
`count_b`, `divergence_ab`, `divergence_ba`, `divergence_sym`;
`bootstrap_pathways()` uses the same leading schema plus
`mean_divergence` / `ci_divergence_*` / `M_divergence` etc.

### Tuning grid

`tune_tree()` returns a `data.frame` with one row per grid point:

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

One-call bundle: `model_fit(tree, newdata = NULL)` returns a tidy
one-row `data.frame` (`logLik | df | nobs | AIC | BIC | perplexity`),
out-of-sample when `newdata` is given, one row per group for a
`transitiontrees_group`.

Per-sequence: `score_sequences(tree, newdata)` returns a data.frame
with `sequence_id, n_scored, log_lik, perplexity` — `n_scored` rather
than `n_obs` because positions whose observed state is outside the
training alphabet are dropped from the score, not counted as obs.

Per-position: `score_positions(tree, newdata)` returns a data.frame
with `sequence_id, position, matched_context, observed,
predicted_prob, log_lik`.

### Comparisons

`compare_trees(tree_a, tree_b, iter = 200)` returns a
`transitiontrees_comparison` object:

```
$pdist        # symmetric (or asymmetric) KL, scalar
$null_dist    # numeric vector, length iter
$p_value      # one-sided p
$pathways     # data.frame with per-pathway divergence
              # (columns: pathway, count_a, count_b,
              #  divergence_ab, divergence_ba, divergence_sym)
```

## 4. Smoothing

`context_tree(..., smoothing = c("floor", "laplace", "kneser_ney", "witten_bell", "jelinek_mercer"))`

- `"floor"` (default) — `ymin` floor, with a `rule` argument:
  - `"interpolate"` (default) — PST-compatible: a distribution with a
    zero-count state is shifted toward uniform,
    `p_i = (1 - k*ymin)*p_i + ymin`, so each zero lands at exactly
    `ymin`; fully observed distributions stay raw MLE. This is what
    gives full-depth machine-precision parity with `PST` (see
    `PARITY.md`).
  - `"cap"` — the original transitiontrees rule: clamp every probability up to
    `ymin` and renormalise. Opt-in via
    `list("floor", ymin = .., rule = "cap")`.
- `"laplace"` — additive-α; `alpha` argument, default 1
- `"kneser_ney"` — back-off along the suffix path with discount `D` and
  continuation probabilities; canonical for PSTs (Begleiter, El-Yaniv
  & Yona 2004 §3); `discount` argument, default 0.75
- `"witten_bell"` — interpolation with novelty-based weight
- `"jelinek_mercer"` — fixed-λ interpolation between order-k and
  order-(k-1); `lambda` argument

Smoothing is applied *once* during `context_tree()`. Re-smoothing a
fitted tree is `smooth_tree(tree, smoothing)` — returns a new tree.

## 5. Inputs

Four accepted shapes:

1. **Wide character matrix or data.frame** (rows = sequences, columns
   = time-steps). NAs allowed; trailing NAs are treated as
   end-of-sequence.
2. **List of character vectors** (each element is one ragged
   sequence). NAs allowed inside.

   For both 1 and 2, missing cells (`NA` or `""`) are **dropped** from
   a sequence — whether trailing or internal — leaving the observed
   states in order; the gap is closed, never coerced into a literal
   `"NA"` state. Both shapes share one cleaning rule (`.ct_clean_seq`)
   so identical data yields an identical tree regardless of container.
   A sequence is kept as long as it has ≥ 1 observed state (it then
   contributes to the root marginal even if it admits no transition).
3. **`stslist`** (TraMineR) — extracted via `as.matrix(stslist)`;
   alphabet honoured; weights honoured if present.
4. **Dynalytics / `mohsaqr`-family model object**, taken directly —
   a Nestimate `netobject` (`c("netobject", "cograph_network")`), a
   cograph network object (`"cograph_network"`), or a `tna` object
   (`"tna"`; also what `codyna::to_tna()` returns). Detection is by
   class **or structurally** (any list-like object exposing a 2-D
   `$data`/`$sequences`/`$seqdata` slot), so a future sibling that
   follows the convention needs no code change. The state-label map
   is resolved across family conventions — `$nodes` id/label table,
   positional `$labels`, or a `labels`/`alphabet` attribute — and an
   integer-coded frame (tna) is decoded through it. So any object in
   the
   family is accepted *as is*: the sequence frame is **extracted**
   from wherever the upstream kept it — `$data` first (the documented
   handoff), then `$sequences` / `$seqdata` / a `"data"`/`"sequences"`
   attribute / an embedded `netobject` — so the caller need not know
   the storage slot. The frame is a non-empty wide trailing-NA-padded
   frame (one row per session), either **character** (Nestimate
   netobject) **or integer-coded** (tna sequences, surfaced by
   `cograph::as_cograph(<tna>)`): an integer-coded frame is *decoded*
   through the `$nodes` id/label table (`label[match(code, id)]`), NA
   preserved — not rejected. Positional `$labels` are read as
   **1-based** (code k = state k); a `0` or otherwise out-of-range
   code, or an id-table code absent from `$nodes`, is an explicit
   error, not a silent drop/recycle. The "reject numeric matrices" rule
   disambiguates a *top-level* square transition matrix and does not
   apply inside a contractual sequence slot (the object's transition
   matrix lives in `$weights`, not `$data`). `$nodes$label` (the
   network's canonical node set) becomes the default `alphabet` unless
   the caller passes one explicitly. The sequence unit is whatever the
   upstream builder's `session=` defined; it is encoded in the object,
   not re-chosen here. This is the contract boundary with the sibling
   network packages: the network object *is* the handoff format.
   Unwrapping happens before row-count / weight detection so those see
   the frame, not the raw slot list.

   A network object that carries **no `$data` sequences** (a pure
   graph: nodes/edges/weights only — i.e. an aggregated transition
   network) is **rejected with an explicit, instructive error**, not
   silently coerced. The original sequences cannot be recovered from
   edge weights — this is the *same* invariant as the numeric-matrix
   rejection below, applied to the object form.

**Grouped inputs.** A grouped family object — Nestimate
`netobject_group`, tna `group_tna`, or any *named* list whose every
element is a single family object — fits one tree per element and
returns a `transitiontrees_group`. Equivalently, a single dataset plus
`group =` (a `$metadata` column name for a network object, or a vector
with one entry per sequence) is split and batch-fitted. All groups
share one alphabet (the union, or an explicit `alphabet =`) so the
trees stay comparable. A grouped object passed where a single object is
expected is detected by class or structurally — it is **not** silently
mistaken for a ragged list.

Numeric matrices are *not* accepted as input — they collide with
"square numeric matrix = transition matrix" in the network-toolkit
ecosystem. Numeric inputs must be cast to character explicitly. The
guiding rule for both this and the pure-graph case: **transitiontrees fits on
sequences, never on aggregated transitions.**

## 6. S3 method conventions

Every S3 class has at minimum `print`, `summary`, and a `plot` (when
visualisation is meaningful). `summary` returns its own S3 class
named `summary.<class>` so `print(summary(x))` dispatches correctly.

Public S3 classes:

- `transitiontrees` — fitted tree
- `transitiontrees_group` — named list of `transitiontrees`s, one per group, from a
  grouped fit (`context_tree(..., group =)` or a grouped family object).
  Follows the Nestimate `*_group` convention: same-key same-order named
  list, plain `lapply` (a single-group failure aborts the batch), outer
  class `c("transitiontrees_group", "list")`. Has `print` and `as.data.frame`
  (row-binds per-group pathway tables with a leading `group` column);
  `compare_trees()` accepts a 2-group one directly.
- `summary.transitiontrees` — tidy summary
- `transitiontrees_tune` — output of `tune_tree()`
- `transitiontrees_comparison` — output of `compare_trees()`
- `transitiontrees_bootstrap` — output of `bootstrap_pathways()`

## 7. Plotting

`plot.transitiontrees(tree, style = c("horizontal", "dendrogram", "icicle", "interactive"))`:

- **`"horizontal"`** (default) — pure-ggplot2 left-to-right phylogram.
  Root at the left, depth increasing rightward, smooth curved branches,
  each leaf labelled with its full arrow-form pathway. Internal-node
  names sit below their nodes. No optional-package dependency; the
  layout to use when citing specific pathways inline.
- **`"dendrogram"`** — pure-ggplot2 radial dendrogram. Root at the
  centre, leaves on the outer ring; node size = count, fill = state.
  No optional-package dependency.
- **`"icicle"`** — circular partition (sunburst) via `ggraph` (in
  Suggests). Errors with an install hint if `ggraph`/`tidygraph` are
  missing.
- **`"interactive"`** — draggable/zoomable hierarchical htmlwidget via
  `visNetwork` (in Suggests). Same encoding as the static styles —
  node size = count, edge width = child's count (flow), fill = last
  state — plus hover tooltips. Errors with an install hint if
  `visNetwork` is missing. `point_size_range` / `edge_size_range`
  control the pixel sizing (defaults `c(10, 45)` / `c(1, 10)`).

The two pure-ggplot2 styles (`dendrogram`, `horizontal`) carry no
optional-package dependency, so the default render works in any R
context (script, knitr, headless). `icicle` and `interactive` each
gate on their Suggests package and error with an install hint when it
is absent.

`plot_pathways()` and `plot_divergence()` are pathway-centric views
that complement (don't replace) the tree plot.

## 8. Hard rules

- Pure base R + ggplot2 in core. `visNetwork` and `ggraph` /
  `tidygraph` only loaded inside the relevant plot branches via
  `requireNamespace()`; they are Suggests, not Imports.
- No Rcpp. No data.table. No tidyverse anywhere.
- All vectorised operations use `vapply` / `tapply` / matrix indexing.
- No for loops in hot paths (training loop in EM is the one exception
  in Nestimate; transitiontrees has no equivalent).
- Roxygen documentation on every exported function. `@return` clause
  mandatory.
- Tests for every public function in `tests/testthat/`.
- Equivalence tests for every estimator in
  `local_testing_and_equivalence/`, gated by
  `TRANSITREES_EQUIV_TESTS=true`, target machine epsilon when an external
  reference exists.

## 9. What `transitiontrees` is not

- Not a hidden-Markov-model package — that's `seqHMM`, `depmixS4`.
- Not a sequence-distance / clustering package — that's `TraMineR`.
- Not an HONEM / higher-order-network embedding package — that's the
  `Nestimate::build_honem()` / `Saebi 2020` family.
- Not a fixed-order Markov chain package — that's `markovchain`.

`transitiontrees` is the variable-order PST toolkit. Anything that doesn't
serve that focus belongs in another package.
