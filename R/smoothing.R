# ---- Smoothing schemes for transitreess ----
#
# Implements:
#   .smooth_floor          - floor MLE at ymin (PST-compatible "interpolate"
#                            rule by default; "cap" = original renormalise)
#   .smooth_laplace        - additive-alpha
#   .smooth_kneser_ney     - back-off along suffix path (parent-prob
#                            approximation; see notes)
#   .smooth_witten_bell    - novelty-weighted interpolation
#   .smooth_jelinek_mercer - fixed-lambda interpolation
#   .ct_smooth_dispatch    - switch dispatcher
#   smooth_tree()      - re-smooth a fitted tree without refitting

#' Floor smoothing of a count vector to a probability distribution.
#'
#' Two rules:
#' \itemize{
#'   \item \code{"interpolate"} (default, PST-compatible): a distribution
#'     that has at least one zero-count state is shifted toward uniform,
#'     \eqn{p_i = (1 - k\,y_{min}) p_i + y_{min}} (k = alphabet size), so
#'     each zero state lands at exactly \code{ymin}. Distributions with
#'     every state observed are left as the raw MLE. This is the floor
#'     used by the archived \pkg{PST} package.
#'   \item \code{"cap"}: clamp every probability up to \code{ymin} and
#'     renormalise (\code{pmax(p, ymin) / sum(...)}). transitrees's original
#'     rule; kept available for back-compatibility.
#' }
#' @noRd
.smooth_floor <- function(counts, parent_prob = NULL, ...,
                          ymin = 0.001, rule = c("interpolate", "cap")) {
  rule <- match.arg(rule)
  k <- length(counts); n <- sum(counts)
  if (n == 0) return(if (is.null(parent_prob)) rep(1 / k, k) else parent_prob)
  p <- counts / n
  if (ymin <= 0) return(p)
  if (rule == "cap") {
    q <- pmax(p, ymin)
    return(q / sum(q))
  }
  ## "interpolate": shift a zero-containing distribution toward uniform.
  ## The coefficient (1 - k*ymin) must stay positive, else probabilities
  ## go negative; require ymin < 1/k.
  if (k * ymin >= 1)
    stop("floor smoothing with rule = \"interpolate\" needs ymin < 1/k ",
         "(k = alphabet size = ", k, "); got ymin = ", ymin,
         ". Lower ymin, or use rule = \"cap\".", call. = FALSE)
  if (any(p == 0)) p <- (1 - k * ymin) * p + ymin
  p
}

#' @noRd
.smooth_laplace <- function(counts, parent_prob = NULL, ...,
                            alpha = 1) {
  k <- length(counts)
  (counts + alpha) / (sum(counts) + alpha * k)
}

#' @noRd
.smooth_kneser_ney <- function(counts, parent_prob = NULL, ...,
                               discount = 0.75) {
  k <- length(counts); n <- sum(counts)
  if (is.null(parent_prob)) parent_prob <- rep(1 / k, k)
  if (n == 0) return(parent_prob)
  n_pos  <- sum(counts > 0)
  high   <- pmax(counts - discount, 0) / n
  back_w <- (discount * n_pos) / n
  high + back_w * parent_prob
}

#' @noRd
.smooth_witten_bell <- function(counts, parent_prob = NULL, ...) {
  k <- length(counts); n <- sum(counts)
  if (is.null(parent_prob)) parent_prob <- rep(1 / k, k)
  if (n == 0) return(parent_prob)
  n_pos  <- sum(counts > 0)
  lambda <- n_pos / (n_pos + n)
  mle    <- counts / n
  (1 - lambda) * mle + lambda * parent_prob
}

#' @noRd
.smooth_jelinek_mercer <- function(counts, parent_prob = NULL, ...,
                                   lambda = 0.5) {
  k <- length(counts); n <- sum(counts)
  if (is.null(parent_prob)) parent_prob <- rep(1 / k, k)
  if (n == 0) return(parent_prob)
  mle <- counts / n
  (1 - lambda) * mle + lambda * parent_prob
}

#' @noRd
.ct_smooth_dispatch <- function(smoothing, counts, parent_prob = NULL) {
  ## `smoothing` is a resolved list from .pt_resolve_smoothing().
  switch(smoothing$method,
    floor          = .smooth_floor(counts, parent_prob,
                                   ymin = smoothing$ymin,
                                   rule = smoothing$rule),
    laplace        = .smooth_laplace(counts, parent_prob,
                                     alpha = smoothing$alpha),
    kneser_ney     = .smooth_kneser_ney(counts, parent_prob,
                                        discount = smoothing$discount),
    witten_bell    = .smooth_witten_bell(counts, parent_prob),
    jelinek_mercer = .smooth_jelinek_mercer(counts, parent_prob,
                                            lambda = smoothing$lambda),
    stop("Unknown smoothing method: '", smoothing$method, "'.",
         call. = FALSE))
}

#' @noRd
.pt_smoothing_defaults <- list(
  floor          = list(method = "floor",          ymin = 0.001,
                        rule = "interpolate"),
  laplace        = list(method = "laplace",        alpha = 1),
  kneser_ney     = list(method = "kneser_ney",     discount = 0.75),
  witten_bell    = list(method = "witten_bell"),
  jelinek_mercer = list(method = "jelinek_mercer", lambda = 0.5)
)

#' @noRd
.pt_resolve_smoothing <- function(spec) {
  ## Accepts a method name (character scalar) or a list whose first
  ## element / "method" entry names the scheme; remaining named entries
  ## override the method's defaults.
  if (is.character(spec) && length(spec) == 1L) {
    if (!spec %in% names(.pt_smoothing_defaults))
      stop("Unknown smoothing method: '", spec, "'. Choose one of: ",
           paste(names(.pt_smoothing_defaults), collapse = ", "),
           call. = FALSE)
    return(.pt_validate_smoothing(.pt_smoothing_defaults[[spec]]))
  }
  if (is.list(spec)) {
    method <- if (!is.null(spec$method)) spec$method else spec[[1L]]
    if (!is.character(method) ||
        !method %in% names(.pt_smoothing_defaults))
      stop("Smoothing method must be one of: ",
           paste(names(.pt_smoothing_defaults), collapse = ", "),
           call. = FALSE)
    nm <- names(spec)
    overrides <- if (is.null(nm)) list() else
      spec[nzchar(nm) & nm != "method"]
    return(.pt_validate_smoothing(
      modifyList(.pt_smoothing_defaults[[method]], overrides)))
  }
  stop("'smoothing' must be a method name (character) or a list ",
       "(method, ...kwargs).", call. = FALSE)
}

#' @noRd
.pt_validate_smoothing <- function(sm) {
  scalar_num <- function(x) is.numeric(x) && length(x) == 1L &&
    !is.na(x) && is.finite(x)
  bad <- function(arg, rule) {
    stop("Invalid smoothing parameter '", arg, "': expected ",
         rule, ".", call. = FALSE)
  }
  switch(sm$method,
    floor = {
      if (!scalar_num(sm$ymin) || sm$ymin < 0)
        bad("ymin", "a finite number >= 0")
      if (!is.null(sm$rule) && !(length(sm$rule) == 1L &&
          sm$rule %in% c("interpolate", "cap")))
        bad("rule", 'either "interpolate" or "cap"')
    },
    laplace = {
      if (!scalar_num(sm$alpha) || sm$alpha < 0)
        bad("alpha", "a finite number >= 0")
    },
    kneser_ney = {
      if (!scalar_num(sm$discount) || sm$discount < 0 || sm$discount > 1)
        bad("discount", "a finite number in [0, 1]")
    },
    witten_bell = NULL,
    jelinek_mercer = {
      if (!scalar_num(sm$lambda) || sm$lambda < 0 || sm$lambda > 1)
        bad("lambda", "a finite number in [0, 1]")
    },
    stop("Unknown smoothing method: '", sm$method, "'.", call. = FALSE))
  sm
}

#' @noRd
.pt_last_state <- function(ctx) {
  if (identical(ctx, .ROOT)) return(.ROOT_LABEL)
  strsplit(ctx, " -> ", fixed = TRUE)[[1L]][[1L]]
}

#' @noRd
.pt_children_of <- function(tree) {
  if (nrow(tree$edges) == 0L) return(list())
  split(tree$edges$child, tree$edges$parent)
}

#' @noRd
.pt_parent_ctx <- function(ctx) {
  if (identical(ctx, .ROOT)) return(NA_character_)
  parts <- strsplit(ctx, " -> ", fixed = TRUE)[[1L]]
  if (length(parts) == 1L) return(.ROOT)
  paste(parts[-1L], collapse = " -> ")
}

#' @noRd
.pt_parent_prob <- function(nodes, ctx) {
  ## Walks up the suffix path to the nearest ancestor present in
  ## `nodes`; needed because nmin can drop intermediate ancestors.
  if (identical(ctx, .ROOT)) return(NULL)
  cur <- ctx
  while (TRUE) {
    cur <- .pt_parent_ctx(cur)
    if (is.na(cur)) return(NULL)
    if (cur %in% names(nodes)) return(nodes[[cur]]$prob)
    if (identical(cur, .ROOT)) return(NULL)
  }
}

#' Re-Smooth a Fitted Pathtree
#'
#' @description
#' Replaces every node's probability vector with a new smoothing
#' scheme without refitting the tree. Walks nodes top-down by depth so
#' each node's parent is re-smoothed before its children read it.
#'
#' @param tree A \code{transitrees}.
#' @param smoothing Smoothing specification: either a method name as a
#'   string (uses defaults for that method's hyperparameters) or a list
#'   of the form \code{list(method, ...kwargs)} for explicit
#'   hyperparameters. Available methods: \code{"floor"} (\code{ymin =
#'   0.001}), \code{"laplace"} (\code{alpha = 1}), \code{"kneser_ney"}
#'   (\code{discount = 0.75}), \code{"witten_bell"},
#'   \code{"jelinek_mercer"} (\code{lambda = 0.5}).
#'
#' @return A new \code{transitrees} with re-smoothed probabilities. Counts
#'   and topology are unchanged.
#'
#' @details
#' For \code{"kneser_ney"} the canonical continuation-distribution
#' formulation requires per-state \emph{type counts}. transitrees does
#' not track these; the implementation uses the parent's smoothed
#' probability as the back-off distribution, an approximation
#' discussed in Begleiter, El-Yaniv & Yona (2004), \emph{JAIR} 22, §3.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' m  <- matrix(sample(c("A","B","C"), 200, TRUE), 20)
#' tr <- context_tree(m, max_depth = 2L, min_count = 3L)
#' smooth_tree(tr, "kneser_ney")
#' smooth_tree(tr, list("kneser_ney", discount = 0.5))
#' }
#' @export
smooth_tree <- function(tree, smoothing = "floor") {
  ## A transitrees_group re-smooths each member, preserving the wrapper.
  if (inherits(tree, "transitrees_group")) {
    out <- lapply(tree, smooth_tree, smoothing = smoothing)
    return(structure(out, class = class(tree),
                     group = attr(tree, "group")))
  }
  stopifnot(inherits(tree, "transitrees"))
  sm <- .pt_resolve_smoothing(smoothing)

  ord <- order(vapply(tree$nodes, function(x) x$depth, integer(1)))
  new_nodes <- tree$nodes
  for (ctx in names(tree$nodes)[ord]) {
    parent_prob <- .pt_parent_prob(new_nodes, ctx)
    new_nodes[[ctx]]$prob <- .ct_smooth_dispatch(
      smoothing   = sm,
      counts      = new_nodes[[ctx]]$counts,
      parent_prob = parent_prob
    )
  }
  tree$nodes     <- new_nodes
  tree$smoothing <- sm
  tree
}

#' Compare Smoothing Schemes on One Dataset
#'
#' @description
#' Fits a context tree under several smoothing schemes — holding
#' \code{max_depth}, \code{nmin} and every other argument fixed — and
#' returns a tidy one-row-per-scheme comparison of tree size and
#' in-sample perplexity. A convenience wrapper over repeated
#' \code{\link{context_tree}} calls that collapses the usual five-line
#' \code{lapply()} loop into a single call.
#'
#' @details
#' The perplexity reported is \strong{in-sample} (computed on the
#' fitting data), so it rewards memorisation and must \emph{not} be used
#' to pick a smoother — use \code{\link{tune_tree}()} for
#' out-of-sample selection. The point of this table is the side-by-side
#' view and the invariance of \code{n_nodes} across schemes: smoothing
#' changes the \emph{probabilities} inside the tree, never \emph{which}
#' contexts exist (topology is set by \code{nmin}, not by the smoother).
#'
#' @param data Either sequence data in any form accepted by
#'   \code{\link{context_tree}} (wide matrix / data.frame, list of
#'   character vectors, TraMineR \code{stslist}, or a
#'   \code{mohsaqr}-family network object) — fitted afresh under each
#'   scheme — \strong{or} an already-fitted \code{transitrees}, which is
#'   \emph{re-smoothed} under each scheme (topology frozen, no
#'   re-count; e.g. to sweep smoothers on a pruned tree).
#' @param smoothing Character vector of smoothing-method names to
#'   compare. Defaults to all five: \code{"floor"}, \code{"laplace"},
#'   \code{"kneser_ney"}, \code{"witten_bell"}, \code{"jelinek_mercer"}.
#' @param ... Further arguments passed to \code{\link{context_tree}}
#'   (e.g. \code{max_depth}, \code{nmin}, \code{alphabet}), held fixed
#'   across every scheme. Ignored when \code{data} is a fitted tree.
#'
#' @return A \code{data.frame} with one row per scheme (in the order
#'   given by \code{smoothing}) and columns \code{smoothing} (method
#'   name), \code{n_nodes} (tree size) and \code{perplexity} (in-sample).
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' seqs <- replicate(50, sample(c("A", "B", "C"), 12, replace = TRUE),
#'                   simplify = FALSE)
#' compare_smoothing(seqs, max_depth = 3L, min_count = 5L)
#' compare_smoothing(seqs, smoothing = c("floor", "kneser_ney"),
#'                   max_depth = 2L)
#' }
#'
#' @seealso \code{\link{smooth_tree}} to re-smooth a fitted tree
#'   without re-counting; \code{\link{tune_tree}} for
#'   cross-validated selection.
#' @export
compare_smoothing <- function(data,
                              smoothing = c("floor", "laplace", "kneser_ney",
                                            "witten_bell", "jelinek_mercer"),
                              ...) {
  if (!is.character(smoothing) || length(smoothing) < 1L)
    stop("'smoothing' must be a non-empty character vector of method names.",
         call. = FALSE)
  ## A fitted tree is re-smoothed (topology frozen, no re-count); raw
  ## data is fitted afresh under each scheme.
  fits <- if (inherits(data, "transitrees"))
    lapply(smoothing, function(s) smooth_tree(data, s))
  else
    lapply(smoothing, function(s) context_tree(data, smoothing = s, ...))
  data.frame(
    smoothing  = smoothing,
    n_nodes    = vapply(fits, n_nodes, integer(1)),
    perplexity = vapply(fits, perplexity, numeric(1)),
    stringsAsFactors = FALSE,
    row.names  = NULL)
}
