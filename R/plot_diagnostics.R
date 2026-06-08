# ---- Distribution and predictive-diagnostic plots ----

#' Per-Context Next-State Distributions
#'
#' @description
#' Small-multiples bar chart of the full next-state distribution
#' \code{P(next | context)}, one panel per context. A per-context
#' probability display: where \code{\link{plot_pathways}()} renders the same
#' numbers as a heatmap, this shows each context's distribution as its
#' own panel, with the modal bar highlighted.
#'
#' @param tree A \code{transitiontrees}.
#' @param contexts Character vector of pathway strings to show, or
#'   \code{NULL} (default) to take the \code{top} most frequent
#'   contexts.
#' @param top Integer. Number of contexts to show when \code{contexts}
#'   is \code{NULL}. Default 12.
#' @param min_count Integer. Drop contexts below this count. Default 1.
#'
#' @return A ggplot object.
#'
#' @examples
#' \donttest{
#' seqs <- replicate(60, sample(c("A", "B", "C"), 10, replace = TRUE),
#'                   simplify = FALSE)
#' tree <- context_tree(seqs, max_depth = 2L)
#' plot_distributions(tree, top = 9)
#' }
#' @export
plot_distributions <- function(tree, contexts = NULL, top = 12L,
                               min_count = 1L) {
  stopifnot(inherits(tree, "transitiontrees"))
  alpha <- tree$alphabet
  pw <- tree_pathways(tree, min_count = min_count, sort_by = "count")
  pw <- if (is.null(contexts))
    utils::head(pw, as.integer(top)) else
    pw[pw$pathway %in% as.character(contexts), , drop = FALSE]
  if (nrow(pw) == 0L)
    stop("No contexts to plot.", call. = FALSE)

  long <- do.call(rbind, lapply(pw$pathway, function(p) {
    ctx  <- if (identical(p, .ROOT_LABEL)) .ROOT else p
    info <- tree$nodes[[ctx]]
    data.frame(pathway  = p,
               state    = factor(alpha, levels = alpha),
               prob     = info$prob,
               is_modal = seq_along(alpha) == which.max(info$prob),
               stringsAsFactors = FALSE)
  }))
  long$pathway <- factor(long$pathway, levels = pw$pathway)

  ggplot2::ggplot(long, ggplot2::aes(x = .data$state, y = .data$prob,
                                     fill = .data$is_modal)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::facet_wrap(~ pathway) +
    ggplot2::scale_fill_manual(
      values = c(`TRUE` = "#08519c", `FALSE` = "#9ecae1"),
      guide = "none") +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(
      x = "next state", y = "P(next | context)",
      title = "Per-context next-state distributions",
      subtitle = sprintf("%d contexts; modal next state highlighted",
                         nrow(pw))) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40", size = 9),
      axis.text.x   = ggplot2::element_text(angle = 45, hjust = 1),
      strip.text    = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank())
}

#' Predictive Diagnostics for Held-Out Scoring
#'
#' @description
#' Visual diagnostics for how a fitted tree scores held-out sequences,
#' built on \code{\link{score_positions}()}:
#' \describe{
#'   \item{\code{type = "position"}}{predicted probability of the
#'     observed next state against position in the sequence — where the
#'     model is confident vs. surprised as a sequence unfolds.}
#'   \item{\code{type = "ecdf"}}{the empirical cumulative distribution of
#'     those predicted probabilities — a calibration-style view of how
#'     often the model assigns high vs. low probability to what actually
#'     happened.}
#'   \item{\code{type = "logloss"}}{the per-position log-loss
#'     \eqn{-\log_2 P(\mathrm{observed})} against position — a
#'     per-position log-loss view. Lower is better
#'     (0 = certain and correct); the dashed line is the mean log-loss
#'     over all scored positions.}
#' }
#'
#' @param tree A \code{transitiontrees}.
#' @param newdata Held-out sequence data in any format accepted by
#'   \code{context_tree()}.
#' @param type One of \code{"position"} (default), \code{"ecdf"}, or
#'   \code{"logloss"}.
#'
#' @return A ggplot object.
#'
#' @examples
#' \donttest{
#' fit  <- replicate(60, sample(c("A", "B", "C"), 10, replace = TRUE),
#'                   simplify = FALSE)
#' tree <- context_tree(fit, max_depth = 2L)
#' new  <- replicate(15, sample(c("A", "B", "C"), 10, replace = TRUE),
#'                   simplify = FALSE)
#' plot_predictive(tree, new, type = "position")
#' plot_predictive(tree, new, type = "ecdf")
#' plot_predictive(tree, new, type = "logloss")
#' }
#' @export
plot_predictive <- function(tree, newdata,
                            type = c("position", "ecdf", "logloss")) {
  stopifnot(inherits(tree, "transitiontrees"))
  type <- match.arg(type)
  sp <- score_positions(tree, newdata)
  if (nrow(sp) == 0L)
    stop("No held-out positions to score.", call. = FALSE)

  base_theme <- ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40", size = 9),
      panel.grid.minor = ggplot2::element_blank())

  if (type == "position") {
    sp$sequence_id <- factor(sp$sequence_id)
    return(
      ggplot2::ggplot(sp, ggplot2::aes(x = .data$position,
                                       y = .data$predicted_prob)) +
        ggplot2::geom_line(ggplot2::aes(group = .data$sequence_id),
                           colour = "grey70", alpha = 0.5) +
        ggplot2::geom_point(ggplot2::aes(colour = .data$predicted_prob),
                            size = 2) +
        ggplot2::scale_colour_gradient(low = "#D55E00", high = "#0072B2",
                                       limits = c(0, 1), name = "P(obs)") +
        ggplot2::ylim(0, 1) +
        ggplot2::labs(
          x = "position in sequence",
          y = "predicted probability of observed move",
          title = "Predictive confidence by position",
          subtitle = "low (orange) = the model was surprised here") +
        base_theme)
  }

  if (type == "logloss") {
    ## -log2 P(observed): per-position log-loss.
    sp$logloss <- -sp$log_lik / log(2)
    mean_ll <- mean(sp$logloss)
    ## per-position mean (how surprise evolves as sequences unfold) — a
    ## real trend, unlike connecting points across unrelated sequences.
    pos_mean <- stats::aggregate(logloss ~ position, data = sp, FUN = mean)
    return(
      ggplot2::ggplot(sp, ggplot2::aes(x = .data$position,
                                       y = .data$logloss)) +
        ggplot2::geom_jitter(ggplot2::aes(colour = .data$logloss),
                             width = 0.18, height = 0, size = 1.6,
                             alpha = 0.55) +
        ggplot2::geom_line(data = pos_mean,
                           ggplot2::aes(x = .data$position,
                                        y = .data$logloss),
                           colour = "grey20", linewidth = 0.9) +
        ggplot2::geom_point(data = pos_mean,
                            ggplot2::aes(x = .data$position,
                                         y = .data$logloss),
                            colour = "grey20", size = 1.8) +
        ggplot2::geom_hline(yintercept = mean_ll, linetype = "dashed",
                            colour = "#D55E00") +
        ggplot2::annotate("text", x = max(sp$position), y = mean_ll,
                          label = sprintf(" mean %.2f bits", mean_ll),
                          hjust = 1, vjust = -0.5, size = 3,
                          colour = "#D55E00") +
        ggplot2::scale_colour_gradient(low = "#0072B2", high = "#D55E00",
                                       name = "log-loss\n(bits)") +
        ggplot2::labs(
          x = "position in sequence",
          y = expression(-log[2] ~ P(observed ~ move)),
          title = "Predictive log-loss by position",
          subtitle = paste0("lower = better; points = positions, ",
                            "black line = per-position mean, ",
                            "dashed = overall mean")) +
        base_theme)
  }

  ggplot2::ggplot(sp, ggplot2::aes(x = .data$predicted_prob)) +
    ggplot2::stat_ecdf(geom = "step", colour = "#0072B2", linewidth = 1) +
    ggplot2::xlim(0, 1) +
    ggplot2::labs(
      x = "predicted probability of observed move",
      y = "cumulative fraction of positions",
      title = "Predictive probability ECDF",
      subtitle = sprintf("%d held-out positions scored", nrow(sp))) +
    base_theme
}

# ---- Pruning illustration (suffix-chain view) -----------------------

#' Suffix-chain pruning data for one pathway (internal)
#'
#' Walks the suffix chain of \code{pathway} — the full context, then the
#' context with its oldest (leftmost) move dropped, and so on down to the
#' root — and, at each link, evaluates the \code{G2} gain of a context
#' against its one-shorter parent under the same rule as
#' \code{\link{prune_tree}} (criterion \code{"G2"}). The requested pathway
#' must itself be a fitted context; the function errors rather than
#' silently starting from a shorter suffix. Each context is classified into
#' one of three states: \code{"informative"} (its own \code{G2} clears the
#' cutoff, so it adds memory over its parent), \code{"pruned"} (it and every
#' deeper context fail to diverge — the redundant tail), or \code{"retained"}
#' (its own \code{G2} is below the cutoff but a deeper context diverges, so
#' \code{prune_tree()} keeps it only as a structural bridge).
#'
#' @return A long \code{data.frame}, one row per (context, state):
#'   \code{L} (memory length), \code{context}, \code{label} (the
#'   most-recent move), \code{state}, \code{prob}, \code{count},
#'   \code{g2}, \code{diverges} (own G2 clears the cutoff), \code{pruned}
#'   (redundant tail), \code{retained} (kept only for a deeper context),
#'   \code{status} (\code{"informative"} / \code{"retained"} /
#'   \code{"pruned"} / \code{"root"}).
#' @noRd
.suffix_chain <- function(tree, pathway, alpha = 0.05) {
  stopifnot(inherits(tree, "transitiontrees"),
            is.character(pathway), length(pathway) == 1L)
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha >= 1)
    stop("'alpha' must be a single number in (0, 1).", call. = FALSE)
  states <- tree$alphabet
  k      <- length(states)
  g2_cv  <- stats::qchisq(1 - alpha, df = k - 1L)

  toks <- if (pathway %in% c(.ROOT_LABEL, .ROOT)) character(0) else
    strsplit(pathway, " -> ", fixed = TRUE)[[1L]]
  n_tok <- length(toks)

  ## Contexts by memory length L = n_tok .. 1 (deepest first), then root.
  keys <- character(0); Ls <- integer(0)
  if (n_tok >= 1L) {
    Ls   <- seq.int(n_tok, 1L)
    keys <- vapply(Ls, function(d)
      paste(toks[(n_tok - d + 1L):n_tok], collapse = " -> "), character(1))
  }
  present <- keys %in% names(tree$nodes)
  ## The requested pathway must itself be a fitted context. If only a shorter
  ## suffix exists, the chain would silently start below the requested depth
  ## and the plot would answer a different memory question than the one asked,
  ## so refuse rather than quietly shorten.
  if (n_tok >= 1L && !present[1L]) {
    deepest <- if (any(present)) keys[which(present)[1L]] else .ROOT_LABEL
    stop(sprintf(
      "Context \"%s\" is not in the tree (max_depth = %d%s). The deepest fitted suffix is \"%s\"; call plot_pruning() with a pathway whose full context exists.",
      pathway, tree$max_depth,
      if (isTRUE(tree$pruned)) ", pruned" else "", deepest),
      call. = FALSE)
  }
  ## Suffix property: if the deepest context exists, so do all its suffixes.
  ## Guard against any gap regardless, so the chain is always contiguous from
  ## the requested pathway down to the root.
  if (n_tok >= 1L && !all(present)) {
    missing <- keys[!present]
    stop(sprintf(
      "Suffix chain of \"%s\" is not contiguous: context(s) %s are missing from the tree.",
      pathway, paste(sprintf('"%s"', missing), collapse = ", ")),
      call. = FALSE)
  }
  keys <- c(keys, .ROOT)
  Ls   <- c(Ls, 0L)

  nodes <- lapply(keys, function(kk) tree$nodes[[kk]])
  m     <- length(keys)

  ## G2 of each context vs its parent (the next, one-shorter context).
  g2 <- c(vapply(seq_len(m - 1L), function(i)
    .ct_g2(nodes[[i]]$counts, nodes[[i + 1L]]$prob), numeric(1)),
    NA_real_)
  ## A context is INFORMATIVE if its own G2 clears the cutoff (it predicts
  ## differently from its one-shorter parent). It is PRUNED if it, and every
  ## deeper context in the chain, fail to diverge (the redundant tail). The
  ## remainder are RETAINED: their own G2 is below the cutoff, but a deeper
  ## context diverges, so prune_tree() keeps them structurally as a bridge to
  ## that descendant. Only INFORMATIVE contexts genuinely add memory.
  diverges    <- g2 > g2_cv                          # NA for the root
  pruned_tail <- c(cumprod(!diverges[-m]) == 1, NA)  # contiguous tail from leaf
  status <- ifelse(is.na(diverges), "root",
             ifelse(diverges, "informative",
              ifelse(pruned_tail, "pruned", "retained")))
  pruned   <- status == "pruned"                     # back-compat logical
  retained <- status == "retained"
  pruned[is.na(diverges)]   <- NA
  retained[is.na(diverges)] <- NA
  labels <- ifelse(keys == .ROOT, .ROOT_LABEL,
                   vapply(keys, function(kk)
                     utils::tail(strsplit(kk, " -> ", fixed = TRUE)[[1L]],
                                 1L), character(1)))
  display <- ifelse(keys == .ROOT, .ROOT_LABEL, keys)

  rows <- lapply(seq_len(m), function(i)
    data.frame(L        = Ls[i],
               context  = display[i],
               label    = labels[i],
               state    = factor(states, levels = states),
               prob     = as.numeric(nodes[[i]]$prob),
               count    = sum(nodes[[i]]$counts),
               g2       = g2[i],
               diverges = diverges[i],
               pruned   = pruned[i],
               retained = retained[i],
               status   = status[i],
               stringsAsFactors = FALSE))
  do.call(rbind, rows)
}

#' Illustrate Pruning Along a Pathway's Suffix Chain
#'
#' @description
#' A suffix-chain pruning view: take one pathway and show, side by
#' side, the next-state distribution at every context along its suffix
#' chain — the full context, then the context with its oldest move
#' dropped, and so on down to the root — marking which contexts
#' \code{\link{prune_tree}} (criterion \code{"G2"}) keeps versus prunes.
#' It answers "how much memory does this pathway actually need?": each
#' context is drawn as its own panel (deepest memory on the left, root on
#' the right) and classified into three states by opacity. \strong{Solid}
#' contexts are \emph{informative} — their own \eqn{G^2} clears the cutoff,
#' so they add predictive information over their one-shorter parent.
#' \strong{Mid-opacity} contexts are \emph{retained}: their own \eqn{G^2} is
#' below the cutoff, but a deeper context diverges, so \code{prune_tree()}
#' keeps them only as a structural bridge — they do not themselves add
#' memory. \strong{Faded} contexts are \emph{pruned} (the redundant tail).
#' The panel title carries the full context, the decision, and the
#' \eqn{G^2}. The requested pathway must be a fitted context; the function
#' errors rather than silently plotting a shorter suffix.
#'
#' @param tree A \code{transitiontrees} (typically \emph{unpruned}, so the
#'   full chain is visible).
#' @param pathway A single pathway string in arrow form
#'   (\code{"A -> B -> C"}, oldest on the left).
#' @param alpha Significance level for the \code{G2} keep/prune decision.
#'   Default 0.05.
#'
#' @return A ggplot object.
#'
#' @details
#' The keep/prune decision is exactly \code{prune_tree}'s G2 rule
#' (\eqn{2N \cdot \mathrm{KL} > \chi^2_{1-\alpha, k-1}}); the cumulative
#' \code{pruned} flag follows the leaf-up amnesia rule. The distributions
#' and counts shown are the same node values reported throughout the
#' pathway API (see \code{PARITY.md}).
#'
#' @examples
#' \donttest{
#' seqs <- replicate(80, sample(c("A","B","C"), 12, replace = TRUE),
#'                   simplify = FALSE)
#' tree <- context_tree(seqs, max_depth = 3L, min_count = 3L)
#' plot_pruning(tree, "A -> B -> C")
#' }
#' @seealso \code{\link{prune_tree}}, \code{\link{plot_distributions}}.
#' @export
plot_pruning <- function(tree, pathway, alpha = 0.05) {
  stopifnot(inherits(tree, "transitiontrees"))
  long <- .suffix_chain(tree, pathway, alpha = alpha)
  states <- tree$alphabet
  pal <- .okabe_ito(length(states))

  ## One panel per context along the chain, deepest (full memory) on the
  ## left, root on the right. The decision is encoded by opacity so it never
  ## competes with the next-state fill colours; the full context, the
  ## decision, and G2 live in the panel title. Three states are kept
  ## distinct: informative contexts add memory over their parent; retained
  ## contexts do not, but survive because a deeper context does; pruned
  ## contexts are the redundant tail.
  node <- long[!duplicated(long$context), , drop = FALSE]
  ord  <- node$context[order(-node$L)]
  tag  <- c(informative = "adds memory",
            retained    = "kept for a deeper context",
            pruned      = "pruned (redundant)")
  mk_lab <- function(ctx, L, status, g2)
    if (status == "root")
      sprintf("L=%d   %s\n(root  -  no parent)", L, ctx)
    else
      sprintf("L=%d   %s\n%s   -   G\u00b2 = %.1f", L, ctx, tag[[status]], g2)
  labmap <- stats::setNames(
    mapply(mk_lab, node$context, node$L, node$status, node$g2), node$context)
  long$panel <- factor(labmap[long$context], levels = labmap[ord])

  ## Highlight the modal next state in each context. Opacity encodes the
  ## decision: informative solid, retained mid, pruned faded.
  long$modal <- stats::ave(long$prob, long$context,
                           FUN = function(p) seq_along(p) == which.max(p)) == 1
  op_map <- c(informative = 1, retained = 0.7, pruned = 0.45, root = 1)
  long$opacity <- unname(op_map[long$status])

  ggplot2::ggplot(long, ggplot2::aes(x = .data$state, y = .data$prob,
                                     fill = .data$state)) +
    ggplot2::geom_col(ggplot2::aes(alpha = .data$opacity,
                                   colour = .data$modal),
                      width = 0.85, linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$prob)),
                       vjust = -0.4, size = 2.7, colour = "grey25") +
    ggplot2::facet_wrap(~ .data$panel, nrow = 1) +
    ggplot2::scale_fill_manual(values = pal, name = "next state") +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "grey20",
                                            `FALSE` = NA), guide = "none") +
    ggplot2::scale_alpha_identity() +
    ggplot2::scale_y_continuous(limits = c(0, 1.08),
                                breaks = seq(0, 1, 0.25)) +
    ggplot2::labs(
      x = NULL, y = "P(next move | context)",
      title = sprintf("How much memory does \u201c%s\u201d need?", pathway),
      subtitle = paste0("each panel = one context along the suffix chain ",
                        "(deepest left \u2192 root right); solid = adds memory, ",
                        "mid = kept for a deeper context, faded = pruned")) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40", size = 9),
      strip.text    = ggplot2::element_text(size = 8.5, lineheight = 1.1),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      legend.position = "right")
}

#' Okabe-Ito colorblind-safe qualitative palette (internal)
#' @noRd
.okabe_ito <- function(n) {
  pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
           "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000")
  if (n <= length(pal)) pal[seq_len(n)] else
    grDevices::hcl.colors(n, "Dark 3")
}

#' Difference (Subtraction) Map Between Two Groups
#'
#' @description
#' A per-context map of how two groups differ in their next-state
#' predictions: one row per shared context, one column per next state.
#' By default each cell is the Pearson standardized residual of the first
#' group against the no-difference null (red = more than expected, blue =
#' less; \code{|r| > 2} notable), which is support-aware and decomposes
#' the per-context \eqn{G^2}; \code{measure = "probability"} shows the raw
#' \code{P(group1) - P(group2)} instead.
#'
#' @param group A \code{transitiontrees_group}.
#' @param groups Optional length-2 character vector naming the two groups
#'   to subtract (\code{group1 - group2}). Defaults to the first two; it
#'   is required when the object has more than two groups.
#' @param depth Integer or \code{NULL}. Restrict to contexts of this
#'   depth (\code{depth = 1} gives the order-1 transition-matrix
#'   difference). \code{NULL} (default) uses all shared contexts.
#' @param min_count Integer. Drop contexts whose count in \emph{either}
#'   group is below this. Default 1.
#' @param comparison Optional \code{transitiontrees_group_comparison} (from
#'   \code{\link{compare_groups}()}). When supplied, contexts whose
#'   behavioral difference is significant (\code{jsd_padj < alpha}) are
#'   starred, so the map shows which differences survived the
#'   permutation test.
#' @param alpha Numeric. FDR threshold for the significance stars when
#'   \code{comparison} is given. Default 0.05.
#' @param annotate Logical. Print the signed difference in each cell
#'   (\code{layout = "tile"} only). Default \code{TRUE}.
#' @param layout One of \code{"tile"} (default; per-context heatmap over
#'   all shared contexts) or \code{"tree"} (the horizontal context-tree
#'   phylogram drawn on a pooled backbone, with each node and branch
#'   coloured by which group reaches that context more — red for the first
#'   group, blue for the second; node size = pooled count).
#' @param measure For \code{layout = "tile"}, what each cell encodes:
#'   \code{"residual"} (default) is the Pearson standardized residual of
#'   the first group against the no-group null (observed vs expected from
#'   the context's margins) — support-aware and decomposing the per-context
#'   \eqn{G^2}; \code{"probability"} is the raw next-state probability
#'   difference \code{P(group1) - P(group2)}. Ignored by \code{"tree"}.
#'
#' @return A ggplot object.
#'
#' @examples
#' \donttest{
#' gx <- replicate(60, sample(c("A","B","C"), 8, replace = TRUE,
#'                            prob = c(.2,.6,.2)), simplify = FALSE)
#' gy <- replicate(60, sample(c("A","B","C"), 8, replace = TRUE,
#'                            prob = c(.2,.2,.6)), simplify = FALSE)
#' grp <- context_tree(c(gx, gy), group = rep(c("x","y"), each = 60),
#'                     max_depth = 2L)
#' plot_difference(grp)                       # residual heatmap
#' plot_difference(grp, layout = "tree")      # difference on the tree map
#' }
#' @seealso \code{\link{compare_groups}} for the significance test.
#' @export
plot_difference <- function(group, groups = NULL, depth = NULL,
                            min_count = 1L, comparison = NULL,
                            alpha = 0.05, annotate = TRUE,
                            layout = c("tile", "tree"),
                            measure = c("residual", "probability")) {
  stopifnot(inherits(group, "transitiontrees_group"))
  layout  <- match.arg(layout)
  measure <- match.arg(measure)
  gn <- names(group)
  if (is.null(groups)) {
    if (length(gn) != 2L)
      stop("The object has ", length(gn), " groups; name the two to ",
           "subtract via 'groups = c(\"a\", \"b\")'.", call. = FALSE)
    groups <- gn[1:2]
  }
  if (length(groups) != 2L || !all(groups %in% gn))
    stop("'groups' must name two groups of the object.", call. = FALSE)
  a <- group[[groups[1L]]]; b <- group[[groups[2L]]]
  alphabet <- a$alphabet

  if (layout == "tree")
    return(.diff_tree(a, b, groups, alphabet))

  ctx <- intersect(names(a$nodes), names(b$nodes))
  if (!is.null(depth))
    ctx <- ctx[vapply(ctx, .cg_depth, integer(1)) == as.integer(depth)]
  ctx <- ctx[vapply(ctx, function(cc)
    a$nodes[[cc]]$n >= min_count && b$nodes[[cc]]$n >= min_count,
    logical(1))]
  if (length(ctx) == 0L)
    stop("No shared contexts to compare (try a lower 'min_count' or ",
         "another 'depth').", call. = FALSE)

  ## significance stars from the permutation comparison, if supplied
  sig_label <- function(p) p
  if (!is.null(comparison)) {
    stopifnot(inherits(comparison, "transitiontrees_group_comparison"))
    cp <- comparison$pathways
    sig <- cp$pathway[cp$jsd_padj < alpha]
    sig_label <- function(p) ifelse(p %in% sig, paste0(p, " *"), p)
  }

  pretty <- function(cc) {
    lab <- if (identical(cc, .ROOT)) .ROOT_LABEL else cc
    sig_label(lab)
  }
  ## per-cell value: Pearson residual of group 1 against the no-group
  ## null (observed vs expected from the context's margins), or the raw
  ## probability difference.
  cell_vals <- function(cc) {
    if (measure == "probability")
      return(a$nodes[[cc]]$prob - b$nodes[[cc]]$prob)
    o <- rbind(a$nodes[[cc]]$counts, b$nodes[[cc]]$counts)
    e <- .cg_expected(o)
    r <- (o[1L, ] - e[1L, ]) / sqrt(e[1L, ])
    r[!is.finite(r)] <- 0
    r
  }
  long <- do.call(rbind, lapply(ctx, function(cc) {
    v <- cell_vals(cc)
    data.frame(pathway    = pretty(cc),
               next_state = factor(alphabet, levels = alphabet),
               val        = v, abs_tot = sum(abs(v)),
               stringsAsFactors = FALSE)
  }))
  ## most-different contexts at the top
  row_order <- unique(long$pathway[order(long$abs_tot)])
  long$pathway <- factor(long$pathway, levels = row_order)

  fmt   <- if (measure == "probability") "%+.2f" else "%+.1f"
  thresh <- if (measure == "probability") 0.01 else 0.5
  fill_name <- if (measure == "probability")
    sprintf("P(%s) - P(%s)", groups[1L], groups[2L]) else
    sprintf("Pearson\nresidual\n(%s)", groups[1L])
  M <- max(abs(long$val)); if (M == 0) M <- 1

  p <- ggplot2::ggplot(long, ggplot2::aes(.data$next_state, .data$pathway,
                                          fill = .data$val)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4)
  if (isTRUE(annotate))
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = ifelse(abs(.data$val) >= thresh,
                                  sprintf(fmt, .data$val), "")),
      size = 3)
  sub <- if (measure == "probability")
    sprintf("red = more likely in %s, blue = more likely in %s",
            groups[1L], groups[2L]) else
    sprintf(paste0("observed vs expected (no-group null): red = %s does ",
                   "the move MORE than expected, blue = less; |r| > 2 notable"),
            groups[1L])
  p +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      limits = c(-M, M), name = fill_name) +
    ggplot2::labs(
      x = "next state", y = NULL,
      title = sprintf("Difference map: %s vs %s", groups[1L], groups[2L]),
      subtitle = sprintf("%s%s", sub,
        if (is.null(comparison)) "" else
          sprintf("  |  * = FDR < %.2f", alpha))) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40", size = 9),
      axis.text.y   = ggplot2::element_text(family = "mono"),
      panel.grid    = ggplot2::element_blank())
}

#' Difference rendered on the horizontal context-tree layout: a pooled
#' backbone whose nodes and edges are coloured by which group reaches a
#' context more (prevalence early - late).
#' @noRd
.diff_tree <- function(a, b, groups, alphabet) {
  trajs_a <- .ct_traj(a$data); trajs_b <- .ct_traj(b$data)
  pooled <- context_tree(c(trajs_a, trajs_b),
                         max_depth = a$max_depth, min_count = a$nmin)
  if (length(pooled$nodes) == 0L)
    stop("Pooled tree is empty; nothing to plot.", call. = FALSE)

  contexts <- names(pooled$nodes)
  depths   <- vapply(contexts, .cg_depth, integer(1))
  cc <- .cg_context_counts(stats::setNames(list(trajs_a, trajs_b), groups),
                           groups, contexts, depths, alphabet)
  Ng <- vapply(groups, function(g) sum(cc[[.ROOT]][g, ]), numeric(1))
  delta <- vapply(contexts, function(ctx)
    sum(cc[[ctx]][groups[1L], ]) / Ng[1L] -
    sum(cc[[ctx]][groups[2L], ]) / Ng[2L], numeric(1))
  names(delta) <- contexts

  layout  <- .ct_horizontal_layout(pooled)
  e_paths <- .ct_horizontal_edge_paths(layout, pooled$edges)
  layout$delta  <- delta[layout$context]
  e_paths$delta <- delta[pooled$edges$child][e_paths$edge_id]

  is_root <- layout$context == .ROOT
  body_df <- layout[!is_root, , drop = FALSE]
  root_df <- layout[ is_root, , drop = FALSE]
  is_leaf <- !(layout$context %in% unique(pooled$edges$parent))
  leaf_df <- layout[is_leaf & !is_root, , drop = FALSE]
  leaf_df$label <- leaf_df$context
  internal_df <- layout[!is_leaf & !is_root, , drop = FALSE]
  internal_df$abbr <- vapply(internal_df$context, .pt_last_state,
                             character(1))

  M <- max(abs(layout$delta), na.rm = TRUE); if (!is.finite(M) || M == 0) M <- 1
  div_fill <- ggplot2::scale_fill_gradient2(
    low = "#2166AC", mid = "grey92", high = "#B2182B", midpoint = 0,
    limits = c(-M, M),
    name = sprintf("prevalence\n%s - %s", groups[1L], groups[2L]))
  div_col <- ggplot2::scale_colour_gradient2(
    low = "#2166AC", mid = "grey75", high = "#B2182B", midpoint = 0,
    limits = c(-M, M), guide = "none")

  ggplot2::ggplot() +
    ggplot2::geom_path(
      data = e_paths,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$edge_id,
                   linewidth = .data$edge_weight, colour = .data$delta),
      lineend = "round", linejoin = "round") +
    ggplot2::geom_point(
      data = body_df,
      ggplot2::aes(x = .data$x, y = .data$y, size = .data$n,
                   fill = .data$delta),
      shape = 21, colour = "grey25", stroke = 0.2, na.rm = TRUE) +
    ggplot2::geom_text(
      data = internal_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$abbr),
      vjust = 1, nudge_y = -0.05, size = 2.7, colour = "grey20",
      fontface = "bold", na.rm = TRUE) +
    ggplot2::geom_text(
      data = leaf_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      hjust = 0, nudge_x = 0.10, size = 3, colour = "grey20",
      family = "mono") +
    ggplot2::geom_point(
      data = root_df, ggplot2::aes(x = .data$x, y = .data$y),
      shape = 16, size = 12, colour = "black") +
    ggplot2::geom_text(
      data = root_df, ggplot2::aes(x = .data$x, y = .data$y),
      label = .ROOT_LABEL, hjust = 1, nudge_x = -0.14, size = 3,
      colour = "grey25", fontface = "bold") +
    div_fill + div_col +
    ggplot2::scale_size_continuous(range = c(4, 14), name = "context count") +
    ggplot2::scale_linewidth_continuous(range = c(0.3, 2.5),
                                        name = "edge flow") +
    ggplot2::scale_x_continuous(
      breaks = seq.int(0L, pooled$max_depth, by = 1L),
      expand = ggplot2::expansion(mult = c(0.17, 0.32))) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40", size = 9),
      legend.position = "right",
      axis.title.y = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid   = ggplot2::element_blank()) +
    ggplot2::labs(
      x = "depth",
      title = sprintf("Difference tree: %s vs %s", groups[1L], groups[2L]),
      subtitle = sprintf(
        "node & edge colour = which phase reaches the context more (red = %s, blue = %s); size = pooled count",
        groups[1L], groups[2L]))
}
