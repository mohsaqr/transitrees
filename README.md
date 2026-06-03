# transitrees

Predictive pathway discovery in categorical sequence data.

`transitrees` fits a variable-depth pathway tree (prediction suffix tree;
Ron, Singer & Tishby 1996) from sequences and exposes a tidy,
pathway-centric API: most common pathways, most predictively divergent
pathways, modal-flip diagnostics, and next-state predictions.

## Why this package

Sequence-analysis tooling tends to either count k-grams (descriptive)
or fit fixed-order Markov chains (one-size-fits-all memory). Real
trajectories — student engagement weeks, patient pathways, clickstreams,
play sequences — usually have memory that varies *by context*: some
states predict the next state on their own, others only when paired
with a longer history. `transitrees` makes that variable-depth structure
the central object and reports it as a ranked list of pathways the
data actually supports.

An earlier R implementation, `PST` 0.94.1 (Gabadinho & Ritschard
2013, *Journal of Statistical Software* **53**(3)), was archived from
CRAN on 2025-11-27. `transitrees` is an **independent implementation** of
the same model — not a replacement for, or fork of, any package:

- pure base R + ggplot2, no Rcpp;
- tidy data.frame outputs by default;
- validated at machine precision against two external references —
  `PST` (fit, `query`, `predict`, `logLik`, topology) and `markovchain`
  (order-1); see `PARITY.md`;
- a pathway-centric API (`tree_pathways()`, `common_pathways()`,
  `divergent_pathways()`, `sharp_pathways()`) that ranks trajectories
  by frequency, predictive divergence, or modal-flip — the structure
  domain experts actually want to read off the model;
- naming chosen so `library(transitrees)` does not collide with sibling
  packages in the `mohsaqr` family: the three exports that would have
  shadowed `tna::prune` / `Nestimate::pathways` / `Nestimate::path_dependence`
  are namespaced as `prune_tree()` / `tree_pathways()` /
  `tree_dependence()`;
- a full predictive-modelling toolchain: `logLik` / `AIC` / `BIC` /
  `perplexity` / `score_*`, five smoothing schemes
  (`floor`, `laplace`, `kneser_ney`, `witten_bell`, `jelinek_mercer`),
  k-fold cross-validated `tune_tree()`, and permutation-tested
  two-tree `compare_trees()`;
- TraMineR `stslist` ingestion + per-sequence weights for social-
  science workflows;
- two ggplot-based plot styles — a pure-ggplot dendrogram (default;
  KL fill, modal-flip ring, count-sized nodes) and a static `ggraph`
  icicle for publication-ready figures.

## Installation

```r
# Github (development version)
remotes::install_github("mohsaqr/transitrees")
```

## Usage

```r
library(transitrees)

# Fit a pathway tree from a wide character matrix or data.frame
tree <- context_tree(seqs, max_depth = 4, nmin = 5, smoothing = "floor")
# explicit hyperparameter:
# context_tree(seqs, max_depth = 4, smoothing = list("floor", ymin = 0.001))
# different scheme:
# context_tree(seqs, max_depth = 4, smoothing = "kneser_ney")

# Inspect
print(tree)
summary(tree)
plot(tree)

# Prune by likelihood-ratio G^2 with familywise alpha = 0.05
pruned <- prune_tree(tree, criterion = "G2", alpha = 0.05)

# The pathway-centric API
tree_pathways(pruned)                   # all pathways, sorted by count
common_pathways(pruned,    n = 8)           # top by frequency
divergent_pathways(pruned, n = 6)           # top by KL from shorter history
divergent_pathways(pruned, flips_only = TRUE)
sharp_pathways(pruned,     n = 5)           # most deterministic continuations

# Per-context KL diagnostic — full table
tree_dependence(pruned)

# Predict next state for new partial sequences
predict(pruned, newdata = list(c("A","B","B"), c("A","A","C")))

# Generate sequences from the fitted tree
generate_sequences(pruned, n = 100, length = 20)
simulate(pruned, nsim = 100, length = 20)  # S3 alias

# Predictive evaluation
logLik(pruned)
AIC(pruned); BIC(pruned)
perplexity(pruned, newdata = test_seqs)
score_sequences(pruned, newdata = test_seqs)

# Cross-validated hyperparameter tuning
tg <- tune_tree(seqs, max_depth = 1:4, nmin = c(3, 5),
                    smoothing = c("floor", "kneser_ney"),
                    prune = c(FALSE, TRUE), k = 5)
attr(tg, "best")

# Two-tree comparison with permutation test
tree_a <- context_tree(group_a)
tree_b <- context_tree(group_b)
compare_trees(tree_a, tree_b, n_perm = 200)

# Bootstrap pathway reliability
boot <- bootstrap_pathways(pruned, iter = 1000, stat = "count")
summary(boot)[, c("pathway", "p_stability", "stable",
                  "informative_rate", "informative")]
plot(boot)

# Tree introspection
query_pathway(pruned, "A -> B")          # full predicted distribution
pathway_exists(pruned, c("A","B","C"))
sub <- subtree(pruned, "A")              # restrict to descendants of "A"

# Visualisation
plot(pruned)                             # pure-ggplot dendrogram (default)
plot(pruned, style = "icicle")           # ggraph partition diagram
```

## Pathway notation

Pathways are reported in arrow notation (`A -> B -> C`), matching the
convention used elsewhere in the `mohsaqr` package family (`Nestimate`,
`tna`, `cograph`). The leftmost state is the *oldest*; the next-state
prediction is conditional on the trajectory ending at the rightmost
state.

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
