# Numerical parity of `transitiontrees`

`transitiontrees` is an independent implementation of the variable-order
prediction-suffix-tree model. Its fit and predictive surface are
cross-validated, at machine precision, against independent external
reference implementations by other authors:

- a **canonical probabilistic-suffix-tree reference** — used to check
  counts, node probabilities at all depths, context queries, per-position
  prediction, `logLik`, and topology;
- an **independent first-order Markov reference** — used as an order-1
  cross-check on the transition probabilities;
- a **standard long → sequence reshaper** — the reference for
  `transitiontrees`'s own `prepare_input()`.

The equivalence suites live **outside** the package (so `R CMD check`
never needs the references) and are gated by an environment variable:

```bash
TRANSITREES_EQUIV_TESTS=true Rscript -e \
  'testthat::test_dir("local_testing_and_equivalence")'
```

## What is at parity (and the tolerance)

| transitiontrees | reference | result | tolerance |
|---|---|---|---|
| `context_tree()` counts | suffix-tree fit | exact, **all depths** | integer-exact (0) |
| `context_tree()` node probabilities | suffix-tree fit | exact, **all depths** | measured 0; asserted ≤1e-10 |
| `query_pathway()` | reference query | exact | asserted ≤1e-10 |
| per-position predictive prob (`score_positions`) | reference per-position prediction | exact | measured 0; asserted ≤1e-9 |
| `nobs()` | reference observation count | exact | integer-exact |
| node set / topology | reference node names | exact set | exact |
| `logLik()` | reference log-likelihood | exact | **1e-10 rel.** (machine) |
| order-1 transition probs | first-order Markov reference | exact | 1e-12 |
| `prepare_input()` (long->sequence) | standard reshaper | identical (see note) | exact |
| `plot_pruning()` chain distributions | reference per-suffix query | exact | **0** |
| `plot_pruning()` chain counts | reference node counts | exact | integer-exact |
| `plot_predictive(type="logloss")` | reference log-loss view | exact | **~1e-15** |

**`prepare_input()` note.** Driven over **420 argument combinations**
(every combination of actor [none / single / two columns], time [none /
timestamp], order, action, and `time_threshold`) on three long-format
example datasets:

- **Time-based cases** — the timestamp/session logic (a new session
  starts when the gap exceeds `time_threshold`): **byte-identical** to
  the reference reshaper, including row order.
- **No-time cases** — pure actor/order grouping: the **set of sequences
  is identical**, but row order can differ. The reference sorts the
  session id by its native column type (integer numerically, multi-column
  by interaction-factor level) while transitiontrees sorts it as a string.
  Row order does not affect `context_tree()` (it fits on the sequence
  set).

The fit is validated against the suffix-tree reference across 30
randomised configurations (alphabet size 2–4, varied sequence
count/length, `max_depth`, `min_count`, `ymin`): **count error 0,
probability error 0** over all shared contexts (6 tests / 209 assertions,
0 failures).

**Published-example check.** The published worked `s1` example
(reference Table 1) reproduces exactly: `query` at contexts `e`, `a`,
`a-b`, `a-b-a`, `a-b-a-a` all match the reference to 0 (e.g.
`a-b-a` → 0.600 / 0.400).

**logLik note.** Both the *per-position* predictive probabilities **and**
the aggregate scalar `logLik()` are machine-exact against the reference
(max relative error ~1e-14 over the 30 configs); there is no longer any
rounding gap, so the scalar is asserted at 1e-10. (The reference exposes
`nobs`/`predict`/`logLik` as S4 methods; when `transitiontrees` is
attached its S3 methods shadow those generics, so the harness calls the
reference's methods via `getMethod()`.)

## The floor-smoothing alignment (2026-06-03)

Reaching full-depth probability parity required matching the reference
`floor` smoother exactly. It shifts a distribution that contains a
zero-count state toward uniform:

```
p_i = (1 - k * ymin) * p_i_raw + ymin      (k = alphabet size)
```

applied **only** to distributions with at least one zero entry; fully
observed distributions stay raw MLE. transitiontrees's `floor` now uses this
rule by default (`rule = "interpolate"`). The original transitiontrees rule —
clamp each probability up to `ymin` and renormalise — remains available
as an opt-in:

```r
context_tree(data, smoothing = list("floor", ymin = 0.001, rule = "cap"))
```

The two rules differ only at sparse nodes (a state with zero count) and
only by O(`ymin`); on fully observed distributions they are identical.
The earlier "machine-precision" claim was, in fact, only validated on
depth ≤1 (a harness bug compared deep contexts under mismatched key
formats); the suite now compares **every** depth.

## G² pruning, measured (2026-06-07)

`prune_tree(criterion = "G2", alpha)` and the reference `prune(gain =
"G2", C)` use the same likelihood-ratio gain but differ in two documented
ways, both now quantified over the gated configs:

1. **Degrees of freedom.** transitiontrees uses the textbook df = k−1 for a
   k-category multinomial LRT — it keeps a child when
   `2·N·KL(child‖parent) > qchisq(1−alpha, k−1)`. The reference compares
   `N·KL` (= G²/2) against a raw cutoff `C`, and its worked example sets
   `C95 = qchisq(0.95, 1)/2` — i.e. the reference **hardcodes df = 1**. So
   `prune_tree(alpha = 0.05)` equals the reference with
   `C = qchisq(.95, k−1)/2`, **not** the worked example's
   `C = qchisq(.95, 1)/2`, except at k = 2 where both are df 1.
2. **Empirical vs smoothed + recursion.** transitiontrees's G² uses the
   *empirical* child distribution (counts/n); the reference uses the
   *smoothed* node prob — an O(`ymin`) difference that can flip a node
   sitting exactly on the cutoff. The reference also prunes top-down
   level-by-level vs transitiontrees's bottom-up collapse, and cannot prune
   an unsmoothed (`ymin = 0`) tree at all.

Measured over 30 randomised configs (alphabet 2–4):

| Cutoff fed to the reference | mean Jaccard | exact-match configs |
|---|---|---|
| worked-example default `qchisq(.95, 1)/2` (df 1) | 0.63 | 10 / 30 |
| matched `qchisq(.95, k−1)/2` (df k−1) | **0.98** | **28 / 30** |

Under the df = 1 cutoff transitiontrees's kept set is **always a subset** of
the reference's (it never keeps a node the reference drops), and at k = 2
the two agree **exactly**. Once the df is matched, agreement is exact on
the vast majority; the residual is the O(`ymin`) smoothed-vs-empirical
tie-break. Conclusion: the two pruners are equivalent up to (a) the df
convention (transitiontrees = statistically-correct k−1; the reference =
hardcoded 1) and (b) an O(`ymin`) borderline tie-break — not an
algorithmic divergence.

## Diagnostic plots (2026-06-07)

The two ggplot diagnostics that mirror the reference's base-graphics
plots are verified on the **numbers they draw**, not pixels:

- **`plot_pruning()`** walks a pathway's suffix chain. Its per-context
  distributions are bit-identical to the reference query (error 0) and
  its counts equal the reference node counts exactly. Its keep/prune
  colouring uses `prune_tree`'s G² rule; the suite asserts the divergence
  flags agree with the reference G² (matched df) on at least 90% of chain
  nodes (the most recent run measured 22/22), the only possible
  disagreement being the same O(`ymin`) smoothed-vs-empirical tie-break
  documented above.
- **`plot_predictive(type = "logloss")`** is
  \eqn{-\log_2 P(\mathrm{observed})}; since the per-position predictive
  probabilities are machine-exact vs the reference, the log-loss is too
  (max abs err ~1e-15).

## What is deliberately NOT at parity

These are transitiontrees's own designs with no equivalent reference
algorithm, so they are validated by transitiontrees's own unit tests, not
against an external reference:

- **`prune_tree()`** — familywise-α G²/KL/AIC/BIC with bottom-up
  collapse. Not a blind "not at parity": the G² gain is **measured**
  against the reference pruner and the relationship is now characterised
  exactly (see "G² pruning, measured" above). KL/AIC/BIC, and the
  bottom-up recursion, remain transitiontrees's own design.
- **`tree_distance()`** — count-weighted symmetric KL; the reference
  distance is a different quantity with a required level argument.
- **`tune_tree()`** — k-fold CV; different split semantics from the
  reference tuner.
- **`generate_sequences()` / `simulate()`** — stochastic; different RNG
  consumption, so comparable only in distribution, not draw-for-draw.
- The other four smoothers (`laplace`, `kneser_ney`, `witten_bell`,
  `jelinek_mercer`) are transitiontrees additions with no reference
  counterpart.

## Files

The gated equivalence suites live in `local_testing_and_equivalence/`:
a suffix-tree parity suite, an order-1 first-order-Markov cross-check, and
a `prepare_input()` long→sequence reshaping check. They run only with
`TRANSITREES_EQUIV_TESTS=true` and are excluded from the package build.
