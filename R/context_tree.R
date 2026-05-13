# ---- Core: context_tree() and friends ----

#' Internal sentinel for the root context.
#' @noRd
.ROOT <- "<root>"

#' @noRd
.ct_coerce <- function(data) {
  if (inherits(data, "stslist")) {
    ## TraMineR sequence object: extract the wide character matrix.
    m <- as.matrix(data); storage.mode(m) <- "character"
    return(as.data.frame(m, stringsAsFactors = FALSE))
  }
  if (is.matrix(data) && !is.numeric(data)) {
    return(as.data.frame(data, stringsAsFactors = FALSE))
  }
  if (is.list(data) && !is.data.frame(data)) {
    ## list of character vectors → ragged trajectories; keep as list
    return(lapply(data, as.character))
  }
  if (is.data.frame(data)) return(data)
  stop("'data' must be a wide data.frame, character/logical matrix, ",
       "list of character vectors, or 'stslist'.", call. = FALSE)
}

#' @noRd
.ct_data_weights <- function(data) {
  ## Pull TraMineR-style weights when present; NULL otherwise.
  if (inherits(data, "stslist")) {
    w <- attr(data, "weights")
    if (!is.null(w)) return(as.numeric(w))
  }
  NULL
}

#' @noRd
.ct_data_n_rows <- function(data) {
  if (inherits(data, "stslist") || is.matrix(data) || is.data.frame(data))
    return(nrow(data))
  if (is.list(data)) return(length(data))
  NA_integer_
}

#' @noRd
.ct_traj <- function(data) {
  ## Attaches an "idx" attribute mapping surviving sequences back to
  ## their original row index — needed to realign per-sequence weights
  ## after empty/short rows are dropped.
  if (is.list(data) && !is.data.frame(data)) {
    out <- lapply(data, function(x) {
      x <- as.character(x)
      x[!is.na(x) & nzchar(x)]
    })
    keep <- which(lengths(out) > 0L)
    res <- out[keep]
    attr(res, "idx") <- keep
    return(res)
  }
  m <- as.matrix(data)
  storage.mode(m) <- "character"
  trajs <- lapply(seq_len(nrow(m)), function(i) {
    r <- m[i, ]
    keep <- !is.na(r) & nzchar(r)
    if (!any(keep)) return(character(0))
    last <- max(which(keep))
    r[seq_len(last)]
  })
  keep <- which(lengths(trajs) >= 2L)
  res <- trajs[keep]
  attr(res, "idx") <- keep
  res
}

#' @noRd
.ct_alphabet <- function(trajs) {
  sort(unique(unlist(trajs, use.names = FALSE)))
}

#' @noRd
.ct_count_table <- function(trajs, depth, alphabet, weights = NULL) {
  ## Vectorised k-gram counter. Returns a named list keyed by context
  ## (string with " -> " separator); the .ROOT key denotes the
  ## marginal distribution. Values are numeric vectors indexed by
  ## alphabet position (integer when weights = NULL).
  k <- depth
  na_size <- length(alphabet)
  weighted <- !is.null(weights)

  if (k == 0L) {
    pooled <- unlist(trajs, use.names = FALSE)
    if (!weighted) {
      counts <- tabulate(match(pooled, alphabet), nbins = na_size)
    } else {
      pooled_w <- unlist(mapply(function(t, w) rep(w, length(t)),
                                 trajs, weights, SIMPLIFY = FALSE),
                          use.names = FALSE)
      idx <- match(pooled, alphabet)
      counts <- vapply(seq_len(na_size),
                       function(j) sum(pooled_w[idx == j], na.rm = TRUE),
                       numeric(1))
    }
    return(setNames(list(counts), .ROOT))
  }

  per_traj <- lapply(seq_along(trajs), function(i) {
    traj <- trajs[[i]]
    L <- length(traj)
    if (L < k + 1L) return(NULL)
    starts <- seq_len(L - k)
    ctx <- vapply(starts, function(t)
      paste(traj[t:(t + k - 1L)], collapse = " -> "),
      character(1))
    list(ctx = ctx,
         nxt = traj[starts + k],
         w   = if (weighted) rep(weights[i], length(ctx)) else NULL)
  })
  ctx_vec <- unlist(lapply(per_traj, `[[`, "ctx"), use.names = FALSE)
  nxt_vec <- unlist(lapply(per_traj, `[[`, "nxt"), use.names = FALSE)
  if (weighted)
    w_vec <- unlist(lapply(per_traj, `[[`, "w"), use.names = FALSE)
  if (length(ctx_vec) == 0L) return(list())

  nxt_idx <- match(nxt_vec, alphabet)
  keep    <- !is.na(nxt_idx)
  ctx_vec <- ctx_vec[keep]
  nxt_idx <- nxt_idx[keep]
  if (weighted) w_vec <- w_vec[keep]

  unique_ctx <- unique(ctx_vec)
  out <- lapply(unique_ctx, function(c) {
    msk <- ctx_vec == c
    if (!weighted) {
      tabulate(nxt_idx[msk], nbins = na_size)
    } else {
      vapply(seq_len(na_size),
             function(j) sum(w_vec[msk & nxt_idx == j]),
             numeric(1))
    }
  })
  setNames(out, unique_ctx)
}

#' @noRd
.ct_kl <- function(p, q) {
  ## D_KL(p || q) in nats; defines 0 log 0 := 0 and returns Inf when
  ## p_i > 0 and q_i == 0.
  msk <- p > 0
  if (any(p[msk] > 0 & q[msk] == 0)) return(Inf)
  sum(p[msk] * log(p[msk] / q[msk]))
}

#' @noRd
.ct_g2 <- function(N_child, p_parent) {
  ## Likelihood-ratio G^2 statistic for the child's count vector tested
  ## against the parent's distribution. Df = (alphabet - 1) for a leaf.
  n <- sum(N_child)
  if (n == 0) return(0)
  expected <- n * p_parent
  msk <- N_child > 0
  2 * sum(N_child[msk] * log(N_child[msk] / pmax(expected[msk],
                                                  .Machine$double.eps)))
}

#' Fit a Prediction Suffix Tree from Categorical Sequence Data
#'
#' @description
#' Estimates a variable-depth context tree (PST; Ron, Singer & Tishby
#' 1996) from a collection of sequences. Each internal node represents
#' a context (string of recent states); each leaf carries a smoothed
#' conditional distribution over the next state. The tree is grown to
#' \code{max_depth}, then optionally pruned via
#' \code{\link{prune_pathtree}()}.
#'
#' @param data Sequence data: wide data.frame / character matrix
#'   (rows = trajectories, columns = time-steps), list of character
#'   vectors, or TraMineR \code{stslist}. Numeric matrices are
#'   rejected; cast to character explicitly.
#' @param max_depth Integer. Maximum context length the tree may
#'   represent. Default 5.
#' @param nmin Integer. Minimum number of times a context must occur
#'   to receive its own node. Default 5. Contexts seen fewer than
#'   \code{nmin} times are absorbed into their parent.
#' @param smoothing Smoothing specification: a method name as a string
#'   (uses defaults for that method's hyperparameters) or a list of
#'   the form \code{list(method, ...kwargs)} for explicit hyperparameters.
#'   Methods: \code{"floor"} (default; \code{ymin = 0.001}),
#'   \code{"laplace"} (\code{alpha = 1}), \code{"kneser_ney"}
#'   (\code{discount = 0.75}), \code{"witten_bell"},
#'   \code{"jelinek_mercer"} (\code{lambda = 0.5}).
#' @param alphabet Character vector. Optional. Override the data-derived
#'   alphabet (useful when the test set may include states unseen in
#'   training).
#' @param weights Numeric vector of per-sequence weights, length equal
#'   to the number of input rows / list elements. If \code{NULL}
#'   (default) and \code{data} is a TraMineR \code{stslist} carrying
#'   weights, those are auto-detected.
#'
#' @return A \code{pathtree} object: a list with components
#' \describe{
#'   \item{nodes}{Named list of node descriptors. Names are context
#'     strings (e.g. \code{"A -> B"}); the root is keyed by the literal
#'     sentinel \code{"<root>"}. Each entry has \code{depth} (integer),
#'     \code{counts} (numeric vector indexed by alphabet),
#'     \code{prob} (smoothed probability vector), and \code{n} (sum of
#'     counts).}
#'   \item{edges}{data.frame with columns \code{parent}, \code{child},
#'     \code{symbol} for fast tree traversal.}
#'   \item{alphabet}{character vector of states.}
#'   \item{max_depth}{Integer. The fitted depth (may be less than the
#'     requested \code{max_depth} if data are short).}
#'   \item{nmin}{the chosen min-count threshold.}
#'   \item{smoothing}{resolved smoothing list (\code{method} +
#'     hyperparameters).}
#'   \item{n_seq, n_obs}{number of sequences and observations.}
#'   \item{pruned}{Logical. \code{TRUE} after
#'     \code{\link{prune_pathtree}()} has been applied.}
#'   \item{pruning}{When \code{pruned} is \code{TRUE}, a list capturing
#'     the criterion / alpha / threshold used; otherwise \code{NULL}.}
#'   \item{data}{The cleaned trajectories (a list of character vectors)
#'     retained for downstream bootstrap and permutation routines.}
#' }
#'
#' @details
#' Construction is bottom-up via k-gram counting. The root holds the
#' marginal next-state distribution; depth-k nodes hold the
#' next-state distribution conditional on the most recent k states.
#' Nodes whose total count falls below \code{nmin} are not created.
#' All nodes start "live" and unpruned; \code{\link{prune_pathtree}()}
#' decides which to retain.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' seqs <- replicate(50, sample(c("A","B","C"), 12, replace = TRUE),
#'                   simplify = FALSE)
#' tree <- context_tree(seqs, max_depth = 3)
#' tree
#' summary(tree)
#' }
#'
#' @references
#' Ron, D., Singer, Y., Tishby, N. (1996). The power of amnesia:
#' learning probabilistic automata with variable memory length.
#' \emph{Machine Learning}, 25, 117-149.
#'
#' @export
context_tree <- function(data,
                         max_depth = 5L,
                         nmin = 5L,
                         smoothing = "floor",
                         alphabet  = NULL,
                         weights   = NULL) {
  n_rows_in <- .ct_data_n_rows(data)
  if (is.null(weights)) weights <- .ct_data_weights(data)

  data  <- .ct_coerce(data)
  trajs <- .ct_traj(data)
  if (length(trajs) == 0L)
    stop("No usable sequences after coercion.", call. = FALSE)

  if (!is.null(weights)) {
    if (length(weights) != n_rows_in)
      stop("'weights' must have length equal to number of input ",
           "sequences (got ", length(weights), ", expected ",
           n_rows_in, ").", call. = FALSE)
    if (any(weights < 0))
      stop("'weights' must be non-negative.", call. = FALSE)
    weights <- weights[attr(trajs, "idx")]
  }

  alphabet  <- if (is.null(alphabet)) .ct_alphabet(trajs) else alphabet
  max_depth <- as.integer(max_depth)
  nmin      <- as.integer(nmin)
  sm        <- .pt_resolve_smoothing(smoothing)

  ## Smooth top-down so each node's parent prob is available when the
  ## node itself is computed (kneser_ney/witten_bell/jelinek_mercer
  ## interpolate with the parent).
  nodes <- list()
  for (d in seq.int(0L, max_depth)) {
    counts_d <- .ct_count_table(trajs, depth = d, alphabet = alphabet,
                                  weights = weights)
    keep <- vapply(counts_d, sum, numeric(1)) >= nmin
    if (d == 0L && .ROOT %in% names(counts_d)) keep[.ROOT] <- TRUE
    counts_d <- counts_d[keep]
    if (length(counts_d) == 0L) {
      max_depth <- d - 1L
      break
    }
    for (ctx in names(counts_d)) {
      counts <- counts_d[[ctx]]
      parent_prob <- if (d == 0L) NULL else .pt_parent_prob(nodes, ctx)
      nodes[[ctx]] <- list(
        depth  = d,
        counts = counts,
        prob   = .ct_smooth_dispatch(sm, counts, parent_prob),
        n      = sum(counts)
      )
    }
  }

  ## Build parent/child relationships:
  ## a node "x_1 > ... > x_k" has parent "x_2 > ... > x_k" (drop the leftmost).
  ## The root carries the .ROOT sentinel as its name.
  edges <- do.call(rbind, lapply(names(nodes), function(ctx) {
    if (identical(ctx, .ROOT)) return(NULL)
    parts <- strsplit(ctx, " -> ", fixed = TRUE)[[1L]]
    parent <- if (length(parts) == 1L) .ROOT else
      paste(parts[-1L], collapse = " -> ")
    if (!parent %in% names(nodes)) return(NULL)
    data.frame(parent = parent, child = ctx, symbol = parts[1L],
               stringsAsFactors = FALSE)
  }))
  if (is.null(edges))
    edges <- data.frame(parent = character(0), child = character(0),
                        symbol = character(0), stringsAsFactors = FALSE)

  structure(
    list(
      nodes      = nodes,
      edges      = edges,
      alphabet   = alphabet,
      max_depth  = max_depth,
      nmin       = nmin,
      n_seq      = length(trajs),
      n_obs      = sum(lengths(trajs)),
      smoothing  = sm,
      pruned     = FALSE,
      pruning    = NULL,
      data       = trajs
    ),
    class = "pathtree"
  )
}

#' Coerce a Pathtree to a Tidy Data Frame
#'
#' @description
#' Returns the canonical tidy node table — identical to
#' \code{pathtree_pathways(tree)}. Lets users do
#' \code{as.data.frame(tree)} and immediately filter, sort, or export
#' with base-R idioms.
#'
#' @param x A \code{pathtree}.
#' @param row.names,optional Ignored.
#' @param ... Forwarded to \code{\link{pathtree_pathways}()}.
#'
#' @return A data.frame with columns \code{pathway}, \code{depth},
#'   \code{count}, \code{modal_next}, \code{prob_next}, \code{KL},
#'   \code{flips}. See \code{\link{pathtree_pathways}}.
#' @export
as.data.frame.pathtree <- function(x, row.names = NULL,
                                    optional = FALSE, ...) {
  pathtree_pathways(x, ...)
}

#' Print a Context Tree
#'
#' @param x A \code{pathtree}.
#' @param max_lines Integer. Maximum tree-rendering lines. Default 25.
#' @param digits Integer. Probability digits. Default 2.
#' @param ... Ignored.
#' @return \code{x} invisibly.
#' @export
print.pathtree <- function(x, max_lines = 25L, digits = 2L, ...) {
  alpha <- x$alphabet
  sm_label <- .pt_smoothing_label(x$smoothing)
  cat(sprintf(
    "<pathtree>  %d nodes, depth <= %d, %d states  [%s]\n",
    length(x$nodes), x$max_depth, length(alpha),
    if (isTRUE(x$pruned)) "pruned" else "unpruned"
  ))
  cat(sprintf("  alphabet : %s\n", paste(alpha, collapse = ", ")))
  cat(sprintf("  fit on   : %d sequences, %d observations\n",
              x$n_seq, x$n_obs))
  cat(sprintf("  smoothing: %s   nmin = %d\n", sm_label, x$nmin))
  if (length(x$nodes) == 0L) return(invisible(x))

  render <- function(parent, prefix) {
    children <- sort(x$edges$child[x$edges$parent == parent])
    for (i in seq_along(children)) {
      child  <- children[i]; last <- i == length(children)
      branch <- if (last) "`-- " else "|-- "
      info   <- x$nodes[[child]]
      modal  <- alpha[which.max(info$prob)]
      label  <- if (identical(child, .ROOT)) "(root)" else
        strsplit(child, " -> ", fixed = TRUE)[[1L]][[1L]]
      lab <- sprintf("%-8s  n=%-5s  -> %s (%.*f)",
                     label, format(as.integer(info$n)), modal, digits,
                     max(info$prob))
      cat(prefix, branch, lab, "\n", sep = "")
      render(child, paste0(prefix, if (last) "    " else "|   "))
    }
  }
  root <- x$nodes[[.ROOT]]
  if (!is.null(root)) {
    modal <- alpha[which.max(root$prob)]
    cat(sprintf("(root)    n=%-5s  -> %s (%.*f)\n",
                format(as.integer(root$n)), modal, digits,
                max(root$prob)))
  }
  out <- utils::capture.output(render(.ROOT, ""))
  if (length(out) > max_lines) {
    cat(paste(out[seq_len(max_lines)], collapse = "\n"), "\n")
    cat(sprintf("... %d more nodes (use as.data.frame(x) or summary(x))\n",
                length(out) - max_lines))
  } else if (length(out) > 0L) {
    cat(paste(out, collapse = "\n"), "\n")
  }
  invisible(x)
}

#' Summary of a Context Tree
#'
#' @param object A \code{pathtree}.
#' @param ... Ignored.
#' @return A \code{summary.pathtree} object. The \code{$table} slot is
#'   the canonical pathway data.frame from
#'   \code{\link{pathtree_pathways}} (columns \code{pathway},
#'   \code{depth}, \code{count}, \code{modal_next}, \code{prob_next},
#'   \code{KL}, \code{flips}), re-sorted by \code{(depth, -count)} so
#'   the structural tree order is read top-to-bottom.
#' @export
summary.pathtree <- function(object, ...) {
  tbl <- pathtree_pathways(object, min_count = 1L)
  if (nrow(tbl) > 0L) {
    tbl <- tbl[order(tbl$depth, -tbl$count), , drop = FALSE]
    rownames(tbl) <- NULL
  }
  structure(
    list(table = tbl, n_states = length(object$alphabet),
         max_depth = object$max_depth, n_nodes = length(object$nodes),
         pruned = isTRUE(object$pruned)),
    class = "summary.pathtree"
  )
}

#' @export
print.summary.pathtree <- function(x, n = 10L, ...) {
  cat(sprintf("<pathtree summary>  %d nodes, depth <= %d, %d states  [%s]\n\n",
              x$n_nodes, x$max_depth, x$n_states,
              if (x$pruned) "pruned" else "unpruned"))
  tbl <- x$table
  print(utils::head(tbl, n), row.names = FALSE)
  if (nrow(tbl) > n)
    cat(sprintf("# ... %d more rows (use as.data.frame(tree) for the full table)\n",
                nrow(tbl) - n))
  invisible(x)
}
