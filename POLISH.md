# pathtree v0.1.1 — Polish Pass

A function-by-function audit of the public surface for *dead-simple
calling, tidy outputs, and consistent conventions*. All changes are
backward-compatible at the **call** level except where noted.

## Inventory and verdicts

| Function | Audit | Change |
|---|---|---|
| `context_tree()` | First-call OK (`context_tree(data)` works); 11 args feels heavy but each has a sensible default and most users will never touch the smoothing kwargs. | none |
| `prune()` | `prune(tree)` works (default G²). | none |
| `predict.pathtree()` | Standard S3 `predict` shape; output matrix/named vector. | none |
| `generate_sequences()` | `n = 1` is a useless default — the function returns a 1×L matrix. | **default bumped to `n = 5`** |
| `pathways()` | Tidy data.frame; column names mostly OK but `length` was confusable (vs base R `length()`) and the underlying concept is **depth**. | **column renamed `length` → `depth`** |
| `common_pathways()` | Took a `length =` filter arg; same naming issue. | **arg renamed `length` → `depth`** |
| `divergent_pathways()` / `sharp_pathways()` | Wrap `pathways()`; inherit the rename automatically. | none |
| `path_dependence()` | Tidy data.frame, but column names diverged from `pathways()`: `context` vs `pathway`, `n` vs `count`. | **renamed `context` → `pathway`, `n` → `count`** |
| `logLik.pathtree()` / `nobs.pathtree()` | Standard `stats` generics; AIC / BIC work for free. | none |
| `perplexity()` | Took a `base =` arg that was *mathematically a no-op* (perplexity is base-invariant). Confusing. | **`base` arg removed** |
| `score_sequences()` / `score_positions()` | Tidy data.frames. | none (already polished by `simplify` pass) |
| `smooth_pathtree()` | Many smoothing-method kwargs but each scheme genuinely needs its own. | none |
| `tune_pathtree()` | Returns `pathtree_tune` data.frame with `print` method, but no `plot()` — users couldn't visualise their grid. | **added `plot.pathtree_tune()`** |
| `query_pathway()` / `subtree()` / `pathway_exists()` | Already simple; subtree adds `local_root` attribute. | none |
| `pathtree_distance()` | Bare scalar. | none |
| `compare_pathtrees()` | Returns `pathtree_comparison` with `print` method, no `plot()`. | **added `plot.pathtree_comparison()`** |
| `plot.pathtree()` | Two ggplot styles after the dependency cleanup. | none |
| `plot_pathways()` / `plot_divergence()` | Tidy lollipop / heatmap. | none |
| `summary.pathtree()` | Old summary table used `context` / `n` — same drift as `path_dependence`. | **table columns renamed to `pathway, depth, count, modal_next, prob_next`** |
| `print.pathtree()` | Old `(ctx) [n=23 p=0.5/0.3/0.2]` rendering was hard to scan; the per-node prob-vector slash list took the eye away from the modal next state — the thing users actually want at a glance. | **rewritten** — aligned columns, shows count + modal next + its probability, lists the smoothing scheme used |
| `print.summary.pathtree()` | Same header re-styled. | banner format aligned with `print.pathtree` |
| **No `as.data.frame` method** | A user with a fitted tree shouldn't have to call `pathways()` to get a flat node table — that's `as.data.frame`'s job. | **added `as.data.frame.pathtree()`** returning `pathway, depth, count, modal_next, prob_next` |

## What "tidy output" means in pathtree

Every tidy data.frame returned by the package now uses the same column
vocabulary:

```
pathway     character   arrow-notation pathway, e.g. "A -> B"; (root) for the root
depth       integer     pathway length (depth in the tree)
count       numeric     observed count of the pathway in training
modal_next  character   the alphabet symbol with highest predicted probability
prob_next   numeric     P(modal_next | pathway)
KL          numeric     KL divergence vs. (k-1)-suffix; NA at root
flips       logical     TRUE iff modal_next changes between this pathway and parent
```

These are the **canonical column names**. They appear in
`pathways()`, `path_dependence()`, `as.data.frame.pathtree()`, and the
`summary.pathtree()` table. They are stable.

## What "dead simple to call" means

Every public function has at most one required argument. All other
parameters have defaults that produce a useful result on a typical
dataset:

```r
context_tree(data)                  # fits with sensible defaults
prune(tree)                         # G² with alpha = 0.05
pathways(tree)                      # all pathways, sorted by count
common_pathways(tree)               # top 10
divergent_pathways(tree)            # top 10 by KL
sharp_pathways(tree)                # top 10 by sharpness
path_dependence(tree)               # full diagnostic table
logLik(tree); AIC(tree); BIC(tree)  # standard model-comparison
perplexity(tree)                    # in-sample
generate_sequences(tree)            # 5 sequences of length 10
plot(tree)                          # dendrogram (ggplot)
as.data.frame(tree)                 # one-row-per-node tidy table
```

## What "every result has a plot()" means

| Result class | `print()` | `plot()` |
|---|---|---|
| `pathtree` | tree skeleton | dendrogram (default) / icicle (Suggests) |
| `summary.pathtree` | banner + node table | — |
| `pathtree_tune` | head + chosen line | **NEW**: perplexity surface, faceted by smoothing × prune, star = best |
| `pathtree_comparison` | observed + p + top divergent pathways | **NEW**: null-distribution histogram + observed line + p-value annotation |

## Backward compatibility

- **`pathways()$length` → `$depth`**: any code reading the `length`
  column directly will need to update. `length()` (the base function)
  on the data.frame still works.
- **`common_pathways(length = k)` → `common_pathways(depth = k)`**:
  same rename of a kwarg.
- **`path_dependence()$context` / `$n` → `$pathway` / `$count`**.
- **`summary(tree)$table` columns**: `context, n, modal` →
  `pathway, count, modal_next` (plus new `prob_next`).
- **`perplexity(..., base = ...)`**: argument removed. Was a no-op.
- **`generate_sequences(n = 1)` → `n = 5`** by default.
- **`print.pathtree()` output** has changed visual format; existing
  scripts that grep its output may need updating.

## Tests

Tests for the renamed columns updated. Three new tests added
(`as.data.frame.pathtree`, `plot.pathtree_tune`, `plot.pathtree_comparison`).
**Suite: 273 tests, all pass.** **PST 0.94.1 equivalence**: still
machine precision (1.11e-16). **`R CMD check`: 0 errors, 0 warnings,
0 notes.**
