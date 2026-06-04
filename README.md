# transitiontrees

Predictive pathway discovery in categorical sequence data.

`transitiontrees` fits a variable-depth pathway tree (prediction suffix tree;
Ron, Singer & Tishby 1996) from sequences and exposes a tidy,
pathway-centric API: most common pathways, most predictively divergent
pathways, modal-flip diagnostics, and next-state predictions.

## Why this package

Sequence-analysis tooling tends to either count k-grams (descriptive)
or fit fixed-order Markov chains (one-size-fits-all memory). Real
trajectories — student engagement weeks, patient pathways, clickstreams,
play sequences — usually have memory that varies *by context*: some
states predict the next state on their own, others only when paired
with a longer history. `transitiontrees` makes that variable-depth structure
the central object and reports it as a ranked list of pathways the
data actually supports.

An earlier R implementation, `PST` 0.94.1 (Gabadinho & Ritschard
2013, *Journal of Statistical Software* **53**(3)), was archived from
CRAN on 2025-11-27. `transitiontrees` is an **independent implementation** of
the same model — not a replacement for, or fork of, any package:

- pure base R + ggplot2, no Rcpp, no data.table, no tidyverse;
- tidy `data.frame` outputs by default, in one canonical schema so every
  table joins cleanly;
- validated at machine precision against three external references by
  other authors — `PST` (counts, node probabilities at all depths,
  `query`, per-position `predict`, `logLik`, topology), `markovchain`
  (order-1 transitions), and `tna::prepare_data` (the long → sequence
  reshaper behind `prepare_input()`); see `PARITY.md`;
- a pathway-centric API (`tree_pathways()`, `common_pathways()`,
  `divergent_pathways()`, `sharp_pathways()`) that ranks trajectories
  by frequency, predictive divergence, or modal-flip — the structure
  domain experts actually want to read off the model;
- naming chosen so `library(transitiontrees)` does not collide with sibling
  packages in the Dynalytics ecosystem: the exports that would have
  shadowed `tna::prune` / `Nestimate::pathways` / `Nestimate::path_dependence`
  are namespaced as `prune_tree()` / `tree_pathways()` /
  `tree_dependence()`;
- a full predictive-modelling toolchain: `logLik` / `AIC` / `BIC` /
  `perplexity` / `score_sequences` / `score_positions` / `model_fit`,
  five smoothing schemes (`floor`, `laplace`, `kneser_ney`,
  `witten_bell`, `jelinek_mercer`), k-fold cross-validated `tune_tree()`,
  and permutation-tested two-tree `compare_trees()`;
- sequence tools: `impute_sequences()` fills internal gaps from the
  fitted tree; `mine_contexts()` / `mine_sequences()` scan for contexts
  where a state is (un)usually likely and for the best/worst-fit
  held-out sequences;
- flexible ingestion — a wide character matrix/data.frame, a list of
  character vectors, a long event log (via `actor` / `time` / `action` /
  `order` arguments, or the standalone `prepare_input()`), a TraMineR
  `stslist`, or a sibling-package network object — plus per-sequence
  weights for social-science workflows;
- four plot styles — `horizontal` (default), `dendrogram`, `icicle`,
  and an `interactive` (visNetwork) view — with count-sized nodes and
  flow-sized edges.

## Installation

```r
# GitHub (development version)
remotes::install_github("mohsaqr/transitiontrees")
```

## Usage

```r
library(transitiontrees)

# --- Fit ------------------------------------------------------------
# From a wide character matrix / data.frame, or a list of character
# vectors. `min_count` is the minimum occurrences a context needs to
# get its own node.
tree <- context_tree(seqs, max_depth = 4, min_count = 5, smoothing = "floor")

# Smoothing as an explicit list (hyperparameters), or a different scheme:
# context_tree(seqs, max_depth = 4, smoothing = list("floor", ymin = 0.001))
# context_tree(seqs, max_depth = 4, smoothing = "kneser_ney")

# From a long event log directly (no manual reshaping):
# context_tree(events, actor = "id", action = "move", order = "step",
#              max_depth = 3)
# or build the wide frame yourself:
# wide <- prepare_input(events, actor = "id", time = "timestamp",
#                       action = "move")

# --- Inspect --------------------------------------------------------
print(tree)
summary(tree)
model_fit(tree)                          # logLik, df, nobs, AIC, BIC, perplexity

# --- Prune ----------------------------------------------------------
# Collapse contexts whose extra depth is not a significant gain over
# their parent (likelihood-ratio G^2 test at familywise alpha = 0.05).
pruned <- prune_tree(tree, criterion = "G2", alpha = 0.05)

# --- The pathway-centric API ---------------------------------------
tree_pathways(pruned)                       # all pathways, sorted by count
common_pathways(pruned,    top = 8)         # top by frequency
divergent_pathways(pruned, top = 6)         # top by KL from shorter history
divergent_pathways(pruned, flips_only = TRUE)
sharp_pathways(pruned,     top = 5)         # most deterministic continuations
tree_dependence(pruned)                     # per-context entropy/KL diagnostic

# --- Predict & evaluate --------------------------------------------
predict(pruned, newdata = list(c("A","B","B"), c("A","A","C")))
predict(pruned, c("A","B"), type = "class")   # modal next state
logLik(pruned); AIC(pruned); BIC(pruned)
perplexity(pruned, newdata = test_seqs)
score_sequences(pruned, newdata = test_seqs)

# Generate sequences from the fitted tree
generate_sequences(pruned, n = 100, length = 20)
simulate(pruned, nsim = 100, length = 20)     # S3 alias

# --- Cross-validated tuning ----------------------------------------
tg <- tune_tree(seqs, max_depth = 1:4, min_count = c(3, 5),
                smoothing = c("floor", "kneser_ney"),
                prune = c(FALSE, TRUE), folds = 5)
attr(tg, "best")

# --- Two-tree comparison with a permutation test -------------------
tree_a <- context_tree(group_a)
tree_b <- context_tree(group_b)
compare_trees(tree_a, tree_b, iter = 200)

# --- Bootstrap pathway reliability ---------------------------------
boot <- bootstrap_pathways(pruned, iter = 1000, stat = "count")
summary(boot)[, c("pathway", "p_stability", "stable",
                  "informative_rate", "informative")]
plot(boot)

# --- Tree introspection --------------------------------------------
query_pathway(pruned, "A -> B")             # full predicted distribution
pathway_exists(pruned, c("A","B","C"))
sub <- subtree(pruned, "A")                 # restrict to descendants of "A"

# --- Impute & mine --------------------------------------------------
impute_sequences(pruned, gappy_seqs)                  # fill internal NA gaps
mine_contexts(pruned, state = "B", min_prob = 0.6)    # contexts where B is likely
mine_sequences(pruned, test_seqs, which = "surprising")  # worst-fit sequences

# --- Visualisation --------------------------------------------------
plot(pruned)                                # horizontal tree (default)
plot(pruned, style = "dendrogram")          # pure-ggplot dendrogram
plot(pruned, style = "icicle")              # ggraph partition diagram
plot(pruned, style = "interactive")         # visNetwork (HTML widget)
plot_pathways(pruned)                       # next-move probability heatmap
plot_divergence(pruned)                     # per-context KL lollipop
plot_distributions(pruned)                  # per-context next-state bars
plot_predictive(pruned, test_seqs)          # held-out confidence diagnostics
```

## Bundled datasets

Four datasets ship with the package for examples and tests:

| Dataset | Shape | Notes |
|---|---|---|
| `trajectories` | 138 × 15 wide character matrix | engagement states over time |
| `group_regulation_long` | long event log, POSIXct time | regulation-of-learning events |
| `ai_long` | long event log, Unix time + session id | AI-prompting moves |
| `engagement` | TraMineR `stslist` | weekly engagement sequences |

## Pathway notation

Pathways are reported in arrow notation (`A -> B -> C`), matching the
convention used elsewhere in the Dynalytics ecosystem (`Nestimate`,
`tna`, `cograph`). The leftmost state is the *oldest*; the next-state
prediction is conditional on the trajectory ending at the rightmost
state. The root context (the marginal next-state distribution) is shown
as `(start)`.

## Bootstrap interpretation

`bootstrap_pathways()` reports stability and informativeness separately.
`p_stability` is the bootstrap-estimated probability that a pathway
statistic falls outside the chosen consistency band under sequence-level
resampling. Small values mean the pathway rarely fails the reproducibility
criterion. `informative_rate` is the fraction of resamples where the
pathway's empirical G² statistic exceeds the chi-square reference cutoff
against its suffix-parent context.

Read the two flags together:

| stable | informative | Interpretation |
|---|---|---|
| TRUE | TRUE | reproducible and predictively distinctive pathway |
| TRUE | FALSE | reproducible pathway, but not predictively distinctive |
| FALSE | TRUE | sharp/divergent pathway carried by an unstable subset of sequences |
| FALSE | FALSE | weak or sample-fragile pathway |

## References

Begleiter, R., El-Yaniv, R., Yona, G. (2004). On prediction using
variable-order Markov models. *Journal of Artificial Intelligence
Research*, **22**, 385–421.

Gabadinho, A., Ritschard, G. (2013). Searching for typical life
trajectories applied to childbirth histories. In: *Gendered Life
Courses Between Standardization and Individualization*, 287–312.

Ron, D., Singer, Y., Tishby, N. (1996). The power of amnesia: learning
probabilistic automata with variable memory length. *Machine Learning*,
**25**, 117–149.

Willems, F.M.J., Shtarkov, Y.M., Tjalkens, T.J. (1995). The context-tree
weighting method: basic properties. *IEEE Transactions on Information
Theory*, **41**, 653–664.

## License

GPL-3.
