# Numerical parity of `transitiontrees`

`transitiontrees` is an independent implementation of the variable-order
prediction-suffix-tree model. Its fit and predictive surface are
cross-validated, at machine precision, against three external reference
packages by other authors:

- **`PST`** (Gabadinho & Ritschard, CRAN-archived) — the canonical
  Probabilistic Suffix Tree package.
- **`markovchain`** (Spedicato) — a first-order Markov chain estimator,
  used as an independent order-1 check.
- **`tna`** (Saqr et al.) — its `prepare_data()` long->sequence reshaper
  is the reference for transitiontrees's own `prepare_input()`.

The equivalence suites live **outside** the package (so `R CMD check`
never needs the references) and are gated by an environment variable:

```bash
TRANSITREES_EQUIV_TESTS=true Rscript -e \
  'testthat::test_dir("local_testing_and_equivalence")'
```

## What is at parity (and the tolerance)

| transitiontrees | reference | result | tolerance |
|---|---|---|---|
| `context_tree()` counts | `PST::pstree` | exact, **all depths** | integer-exact (0) |
| `context_tree()` node probabilities | `PST::pstree` | exact, **all depths** | **0** (≤1e-16) |
| `query_pathway()` | `PST::query` | exact | 1e-10 |
| per-position predictive prob (`score_positions`) | `PST::predict(decomp=TRUE)` | exact | **0** (≤1e-9) |
| `nobs()` | `PST::nobs` | exact | integer-exact |
| node set / topology | `PST::nodenames` | exact set | exact |
| `logLik()` | `PST::logLik` | matches | 1e-4 rel. (see note) |
| order-1 transition probs | `markovchain::markovchainFit` | exact | 1e-12 |
| `prepare_input()` (long->sequence) | `tna::prepare_data` | identical (see note) | exact |

**`prepare_input()` note.** Driven over **420 argument combinations**
(every combination of actor [none / single / two columns], time [none /
timestamp], order, action, and `time_threshold`) on three long-format
Nestimate datasets (`group_regulation_long`, `ai_long`, `human_long`):

- **Time-based cases** — the timestamp/session logic (a new session
  starts when the gap exceeds `time_threshold`): **byte-identical** to
  tna, including row order.
- **No-time cases** — pure actor/order grouping: the **set of sequences
  is identical**, but row order can differ. tna sorts the session id by
  its native column type (integer numerically, multi-column by
  interaction-factor level) while transitiontrees sorts it as a string. Row
  order does not affect `context_tree()` (it fits on the sequence set).

The fit is validated across 30 randomised configurations (alphabet size
2–4, varied sequence count/length, `max_depth`, `min_count`, `ymin`):
**count error 0, probability error 0** over all 1222 shared contexts.

**logLik note.** The *per-position* predictive probabilities are
machine-exact (the `score_positions` row above), so the sequence
likelihoods are identical. PST's aggregate `logLik()` scalar differs by
~1e-5 only because of internal rounding in that method, not because the
model scores anything differently — hence the looser tolerance on the
scalar and the exact tolerance on the per-position probabilities.

## The floor-smoothing alignment (2026-06-03)

Reaching full-depth probability parity required matching PST's `floor`
smoother exactly. PST shifts a distribution that contains a zero-count
state toward uniform:

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

## What is deliberately NOT at parity

These are transitiontrees's own designs with no equivalent algorithm in PST, so
they are validated by transitiontrees's own unit tests, not against PST:

- **`prune_tree()`** — familywise-α G²/KL/AIC/BIC with bottom-up
  collapse. PST's `prune(gain, C)` uses a different recursion; the
  surviving node sets diverge on some configs, and `PST::prune()` cannot
  prune an unsmoothed (`ymin = 0`) tree at all.
- **`tree_distance()`** — count-weighted symmetric KL; PST's `pdist`
  is a different distance with a required level argument.
- **`tune_tree()`** — k-fold CV; different split semantics from
  PST's `tune`.
- **`generate_sequences()` / `simulate()`** — stochastic; different RNG
  consumption, so comparable only in distribution, not draw-for-draw.
- The other four smoothers (`laplace`, `kneser_ney`, `witten_bell`,
  `jelinek_mercer`) are transitiontrees additions with no PST counterpart.

## Files

- `local_testing_and_equivalence/test-equiv-PST.R` — PST parity suite.
- `local_testing_and_equivalence/test-equiv-markovchain.R` — order-1
  cross-check against `markovchain`.
- `local_testing_and_equivalence/test-equiv-tna-prepare.R` —
  `prepare_input()` vs `tna::prepare_data()`.
