# ---- Core: context_tree() and friends ----

#' Internal sentinel for the root context. This value is the pinned
#' object key (parity contract) and never changes.
#' @noRd
.ROOT <- "<root>"

#' User-facing display label for the root context, shown in pathway
#' tables, plots, print methods, and accepted by \code{query_pathway()}.
#' Distinct from the internal sentinel \code{.ROOT}.
#' @noRd
.ROOT_LABEL <- "(start)"

#' Is this a Dynalytics / \code{mohsaqr}-family model object?
#'
#' pathtree is a sibling of \code{Nestimate}, \code{cograph},
#' \code{tna}, \code{codyna}, \code{temporal}, \code{Saqrlab},
#' \code{Snakeplot}, so it takes their model objects directly. Detected
#' by known class — \code{netobject} (Nestimate),
#' \code{cograph_network} (cograph), \code{tna} (tna; also what
#' \code{codyna::to_tna()} produces) — \emph{or} structurally: any
#' list-like object carrying a 2-D sequence frame in a
#' \code{$data}/\code{$sequences}/\code{$seqdata} slot. The structural
#' arm means a new sibling that follows the same convention works with
#' no code change here. Whether it can actually be fitted still depends
#' on a usable sequence frame being present — see
#' \code{.ct_unwrap_netobject()}.
#' @noRd
.ct_is_netobject <- function(x) {
  if (inherits(x, c("netobject", "cograph_network", "tna")))
    return(TRUE)
  ## Structural duck-typing for any other family model: a list (but
  ## not a bare data.frame) exposing a rectangular sequence slot.
  if (is.list(x) && !is.data.frame(x)) {
    for (nm in c("data", "sequences", "seqdata")) {
      v <- x[[nm]]
      if (!is.null(v) && length(dim(v)) == 2L &&
          NROW(v) > 0L && NCOL(v) > 0L)
        return(TRUE)
    }
  }
  FALSE
}

#' Unwrap a Nestimate / cograph network object into its sequence frame.
#'
#' Objects in the \code{mohsaqr} network family
#' (\code{c("netobject", "cograph_network")} from
#' \code{Nestimate::build_network()}, or \code{c("cograph_network",
#' "list")} from \code{cograph::cograph()}) optionally carry a
#' \code{$data} slot: the wide, trailing-NA-padded character data.frame
#' \code{context_tree()} consumes — one row per session, one column per
#' step. \code{$nodes$label} holds the canonical alphabet (the network's
#' node set), surfaced so the fitted tree's symbol set matches the
#' network even when a code never appears within a counted window.
#'
#' The sequence frame is extracted "as is" by scanning every place a
#' \code{mohsaqr} network object is known to keep it — the \code{$data}
#' slot first (the documented handoff), then other plausible
#' slots/attributes (\code{$sequences}, \code{$seqdata}, an embedded
#' netobject, or a \code{"data"}/\code{"sequences"} attribute) — so a
#' caller can hand over any family object without knowing where the
#' upstream stashed it.
#'
#' A sequence slot may hold an \strong{integer-coded} frame plus a
#' \code{$nodes} id/label table (this is exactly how \code{tna} stores
#' sequences, surfaced by \code{cograph::as_cograph(<tna>)}). The
#' "reject numeric matrices" rule is about disambiguating a
#' \emph{top-level} square transition matrix — it does \strong{not}
#' apply here, because a network object's sequence slot is
#' contractually sequences (its transition matrix lives in
#' \code{$weights}). So a numeric frame found in the slot is decoded,
#' not rejected: integer codes are mapped through \code{$nodes}
#' (\code{label[match(code, id)]}) when a label table exists, else
#' cast to character; \code{NA} is preserved as end-of-sequence.
#'
#' Only a \emph{pure graph} projection — nodes/edges/weights with the
#' sequences provably absent everywhere (e.g. the output of
#' \code{cograph::cograph(<netobject>)}, whose constructor nulls the
#' data slot) — cannot yield a variable-order tree: the original
#' sequences are unrecoverable from edge weights. Only then does this
#' error, with explicit guidance, rather than fabricating sequences.
#' @return list(data = <data.frame>, alphabet = <character or NULL>).
#' @noRd
.ct_seq_frame_ok <- function(d) {
  ## Any non-empty rectangular frame in a contractual sequence slot
  ## counts — numeric/integer included (tna codes states as integers).
  (is.data.frame(d) || is.matrix(d)) &&
    length(dim(d)) == 2L && nrow(d) > 0L && ncol(d) > 0L
}

#' Resolve a code -> label decoder across the family's conventions.
#'
#' Sibling packages disagree on where the state-label map lives:
#' cograph/Nestimate use a \code{$nodes} id/label table; \code{tna}
#' uses a positional \code{$labels} vector (code k = state k) and also
#' stamps \code{attr(data, "labels")} / \code{"alphabet")}. We try them
#' in that order and fall back to no decoding.
#' @return list(ids = <codes> or NULL, labels = <character> or NULL).
#'   \code{ids = NULL} with non-NULL \code{labels} means positional.
#' @noRd
.ct_resolve_labels <- function(x, hit) {
  nodes <- x$nodes
  if (is.data.frame(nodes) && !is.null(nodes$label) && !is.null(nodes$id))
    return(list(ids = nodes$id, labels = as.character(nodes$label)))
  for (src in list(x$labels, x$alphabet,
                   attr(hit, "labels"), attr(hit, "alphabet"),
                   attr(hit, "levels")))
    if (is.character(src) && length(src) > 0L)
      return(list(ids = NULL, labels = src))   # positional
  list(ids = NULL, labels = NULL)
}

#' Decode a (possibly integer-coded) sequence frame to a character
#' data.frame, remapping codes through whatever label map the family
#' object carries.
#' @noRd
.ct_decode_seq_frame <- function(x, hit) {
  m <- as.matrix(hit)
  if (is.numeric(m)) {
    lk  <- .ct_resolve_labels(x, hit)
    if (!is.null(lk$labels)) {
      codes <- as.vector(m)
      if (is.null(lk$ids)) {
        ## Positional: code k = state k, 1-based (the tna convention).
        ## A 0 (or any out-of-range code) would index labels[0] / past
        ## the end and silently drop or recycle, corrupting the frame.
        bad <- !is.na(codes) &
          (codes < 1L | codes > length(lk$labels) | codes != floor(codes))
        if (any(bad))
          stop("Integer-coded sequence frame has positional state ",
               "code(s) outside 1..", length(lk$labels), " (the ",
               "label set): ", paste(sort(unique(codes[bad])),
                                      collapse = ", "), ". pathtree ",
               "reads positional $labels as 1-based (code k = state ",
               "k); recode to 1-based or supply a $nodes id/label ",
               "table.", call. = FALSE)
        idx <- codes
      } else {
        idx <- match(codes, lk$ids)                  # id table
        bad <- !is.na(codes) & is.na(idx)
        if (any(bad))
          stop("Integer-coded sequence frame has code(s) absent from ",
               "the $nodes id table: ",
               paste(sort(unique(codes[bad])), collapse = ", "), ".",
               call. = FALSE)
      }
      out <- matrix(lk$labels[idx], nrow = nrow(m), ncol = ncol(m),
                    dimnames = dimnames(m))
    } else {
      out <- m
      storage.mode(out) <- "character"               # 1,2 -> "1","2"
    }
  } else {
    out <- m
    storage.mode(out) <- "character"                 # already char
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}

.ct_unwrap_netobject <- function(x) {
  cls <- paste(class(x), collapse = "/")

  ## Candidate sequence locations, in priority order. $data is the
  ## documented handoff; the rest are defensive so the family object
  ## is accepted "as is" wherever the upstream kept the sequences.
  cand <- list(x$data, x$sequences, x$seqdata,
               attr(x, "data"), attr(x, "sequences"))
  ## An embedded netobject (some wrappers nest the original).
  emb <- Filter(function(z) inherits(z, "netobject"), x)
  if (length(emb)) cand <- c(cand, lapply(emb, `[[`, "data"))

  hit <- Find(.ct_seq_frame_ok, cand)
  if (is.null(hit))
    stop("This ", cls, " object carries no sequence data anywhere ",
         "(no usable $data / $sequences / $seqdata / embedded ",
         "netobject). It looks like a pure graph - an aggregated ",
         "transition network - and pathtree fits on raw sequences, ",
         "not on aggregated transitions: the original sequences ",
         "cannot be recovered from edge weights (the same reason ",
         "numeric transition matrices are rejected). Pass the ",
         "sequence-bearing object instead, e.g. the netobject from ",
         "Nestimate::build_network(..., format = \"wide\") or ",
         "cograph::as_cograph(<netobject>), or the original wide ",
         "sequence data.frame.", call. = FALSE)

  d <- .ct_decode_seq_frame(x, hit)
  ## Default alphabet = the object's declared node/label set, so the
  ## fitted tree's symbol space matches the source model.
  alpha <- .ct_resolve_labels(x, hit)$labels
  list(data = d, alphabet = alpha)
}

#' Is this a grouped family object (a named list of fittable objects)?
#'
#' Recognises Nestimate's \code{netobject_group} and \pkg{tna}'s
#' \code{group_tna} by class, plus any \emph{named} list whose every
#' element is itself a single fittable family object (netobject / tna /
#' cograph network). A bare ragged list of character vectors is
#' \strong{not} a group (its elements are not family objects), and a
#' single network object is not a group (it is detected as a single).
#' @noRd
.ct_is_group <- function(x) {
  if (inherits(x, c("netobject_group", "group_tna", "pathtree_group")))
    return(TRUE)
  if (is.list(x) && !is.data.frame(x) && !.ct_is_netobject(x) &&
      length(x) > 0L) {
    nm <- names(x)
    if (!is.null(nm) && all(nzchar(nm)) &&
        all(vapply(x, .ct_is_netobject, logical(1))))
      return(TRUE)
  }
  FALSE
}

#' Normalise a grouped object to a plain named list of elements.
#' @noRd
.ct_group_elements <- function(x) {
  els <- unclass(x)
  attr(els, "group") <- NULL
  if (is.null(names(els)) || any(!nzchar(names(els))))
    names(els) <- paste0("group", seq_along(els))
  els
}

#' Shared alphabet across a group's elements: the explicit \code{alphabet}
#' if supplied, else the union of each element's declared symbol set.
#' @noRd
.ct_group_alphabet <- function(elements, alphabet) {
  if (!is.null(alphabet)) return(alphabet)
  als <- lapply(elements, function(e)
    if (.ct_is_netobject(e)) .ct_unwrap_netobject(e)$alphabet else NULL)
  als <- sort(unique(unlist(als, use.names = FALSE)))
  if (length(als)) als else NULL
}

#' Split one dataset into a named list of subsets by a grouping factor.
#'
#' \code{group} is either a character scalar naming a column of a network
#' object's \code{$metadata}, or a vector with one entry per input
#' sequence. Returns the subsets, the resolved grouping-variable name,
#' and the overall (shared) alphabet so every subtree spans the same
#' symbols.
#' @noRd
.ct_group_split_by <- function(data, group) {
  var <- NA_character_
  if (.ct_is_netobject(data)) {
    base <- .ct_unwrap_netobject(data)$data
    if (is.character(group) && length(group) == 1L) {
      meta <- data$metadata
      if (is.null(meta) || is.null(meta[[group]]))
        stop("group = \"", group, "\" is not a column of the object's ",
             "$metadata.", call. = FALSE)
      g <- meta[[group]]; var <- group
    } else g <- group
  } else {
    if (is.character(group) && length(group) == 1L)
      stop("A single-string 'group' names a $metadata column and is ",
           "only valid for a network object; for a matrix / data.frame ",
           "/ list, pass a grouping vector with one entry per sequence.",
           call. = FALSE)
    base <- data
    g    <- group
  }
  n <- .ct_data_n_rows(base)
  if (length(g) != n)
    stop("'group' must have one entry per input sequence (got ",
         length(g), ", expected ", n, ").", call. = FALSE)
  g   <- as.factor(g)
  idx <- split(seq_len(n), g)
  subsets <- lapply(idx, function(ix)
    if (is.list(base) && !is.data.frame(base)) base[ix]
    else base[ix, , drop = FALSE])
  ## Shared alphabet: honour a network object's declared node set (so a
  ## valid state absent from the observed rows is not lost), exactly as a
  ## single-object netobject fit does; otherwise derive it from the
  ## observed sequences.
  declared <- if (.ct_is_netobject(data))
                .ct_unwrap_netobject(data)$alphabet else NULL
  overall  <- if (!is.null(declared)) declared
              else .ct_alphabet(.ct_traj(.ct_coerce(base)))
  list(subsets = subsets, var = var, alphabet = overall, idx = idx, n = n)
}

#' Validate a top-level \code{weights} vector and split it across the
#' groups produced by \code{\link{.ct_group_split_by}} (same row index
#' as \code{$subsets}). \code{NULL} weights yield a same-length list of
#' \code{NULL}s so each per-group fit is unweighted.
#' @noRd
.ct_split_weights <- function(weights, parts) {
  if (is.null(weights)) {
    out <- vector("list", length(parts$subsets))
    names(out) <- names(parts$subsets)
    return(out)
  }
  if (!is.numeric(weights))
    stop("'weights' must be a numeric vector.", call. = FALSE)
  if (length(weights) != parts$n)
    stop("'weights' must have length equal to number of input ",
         "sequences (got ", length(weights), ", expected ",
         parts$n, ").", call. = FALSE)
  if (anyNA(weights))
    stop("'weights' must not contain NA.", call. = FALSE)
  if (any(weights < 0))
    stop("'weights' must be non-negative.", call. = FALSE)
  lapply(parts$idx, function(ix) weights[ix])
}

#' Assemble a list of fitted trees into a \code{pathtree_group}.
#' @noRd
.ct_as_group <- function(trees, group_var) {
  structure(trees, class = c("pathtree_group", "list"),
            group = if (is.null(group_var)) NA_character_ else group_var)
}

#' @noRd
.ct_coerce <- function(data) {
  if (.ct_is_netobject(data)) {
    ## Nestimate netobject / cograph network object: take its $data
    ## sequence frame (errors with guidance if it is a pure graph with
    ## no sequences). context_tree() normally unwraps upstream so
    ## weights/row-count see the frame; this branch keeps the
    ## documented dispatch list self-consistent.
    return(.ct_unwrap_netobject(data)$data)
  }
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
       "list of character vectors, 'stslist', or a sequence-bearing ",
       "Nestimate 'netobject' / cograph network object.",
       call. = FALSE)
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

#' Drop missing / empty cells from one sequence, leaving observed
#' states in order. NA and "" are treated identically whether they fall
#' at the end of a row (end-of-sequence padding) or internally (a gap in
#' the recording): the gap is closed, never turned into a literal "NA"
#' state. This is the single cleaning rule shared by the list and the
#' matrix/data.frame paths so identical data yields identical trees
#' regardless of input shape.
#' @noRd
.ct_clean_seq <- function(x) {
  x <- as.character(x)
  x[!is.na(x) & nzchar(x)]
}

#' @noRd
.ct_traj <- function(data) {
  ## Attaches an "idx" attribute mapping surviving sequences back to
  ## their original row index — needed to realign per-sequence weights
  ## after empty/short rows are dropped.
  if (is.list(data) && !is.data.frame(data)) {
    out <- lapply(data, .ct_clean_seq)
  } else {
    m <- as.matrix(data)
    storage.mode(m) <- "character"
    out <- lapply(seq_len(nrow(m)), function(i) .ct_clean_seq(m[i, ]))
  }
  keep <- which(lengths(out) > 0L)
  res <- out[keep]
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
#' @param data Sequence data in any of these forms: a wide data.frame /
#'   character matrix (rows = trajectories, columns = time-steps), a
#'   list of character vectors, a TraMineR \code{stslist}, or a
#'   \strong{Dynalytics / \code{mohsaqr}-family model object}, taken
#'   directly: a Nestimate \code{netobject}
#'   (\code{Nestimate::build_network()}), a \code{cograph} network
#'   object (\code{cograph::cograph()} / \code{as_cograph()}), or a
#'   \code{tna} object (the \pkg{tna} package; also what
#'   \code{codyna::to_tna()} returns). Any other family object that
#'   follows the same convention (a \code{$data}/\code{$sequences}/
#'   \code{$seqdata} sequence slot) is detected structurally and works
#'   with no special-casing. For these objects the sequence frame is
#'   extracted from wherever the upstream stored it; an integer-coded
#'   frame (\code{tna}) is decoded through the object's label map
#'   (\code{$nodes} id/label table, positional \code{$labels}, or a
#'   \code{labels}/\code{alphabet} attribute), and that label set
#'   becomes the default alphabet so the tree shares the model's
#'   symbol space. A \emph{pure graph} object that carries no
#'   sequences anywhere (an aggregated transition network) is rejected
#'   with guidance — sequences cannot be recovered from edge weights,
#'   the same reason numeric transition matrices are rejected; route
#'   to the sequence-bearing object explicitly.
#' @param max_depth Integer. Maximum context length the tree may
#'   represent. Default 5.
#' @param min_count Integer. Minimum number of times a context must occur
#'   to receive its own node. Default 5. Contexts seen fewer than
#'   \code{min_count} times are absorbed into their parent.
#' @param smoothing Smoothing specification: a method name as a string
#'   (uses defaults for that method's hyperparameters) or a list of
#'   the form \code{list(method, ...kwargs)} for explicit hyperparameters.
#'   Methods: \code{"floor"} (default; \code{ymin = 0.001}),
#'   \code{"laplace"} (\code{alpha = 1}), \code{"kneser_ney"}
#'   (\code{discount = 0.75}), \code{"witten_bell"},
#'   \code{"jelinek_mercer"} (\code{lambda = 0.5}). The \code{"floor"}
#'   method also takes \code{rule}: \code{"interpolate"} (default, the
#'   PST-compatible floor — a distribution with a zero-count state is
#'   shifted toward uniform so each zero lands at exactly \code{ymin})
#'   or \code{"cap"} (clamp every probability up to \code{ymin} and
#'   renormalise), e.g. \code{list("floor", ymin = 0.001, rule = "cap")}.
#' @param alphabet Character vector. Optional. Override the data-derived
#'   alphabet (useful when the test set may include states unseen in
#'   training).
#' @param weights Numeric vector of per-sequence weights, length equal
#'   to the number of input rows / list elements. If \code{NULL}
#'   (default) and \code{data} is a TraMineR \code{stslist} carrying
#'   weights, those are auto-detected.
#' @param group Optional grouping for a \strong{batch fit}. Either a
#'   character scalar naming a column of a network object's
#'   \code{$metadata}, or a vector with one entry per input sequence.
#'   When supplied (or when \code{data} is itself a grouped family
#'   object such as a Nestimate \code{netobject_group} or a \pkg{tna}
#'   \code{group_tna}), \code{context_tree()} fits one tree per group
#'   over a shared alphabet and returns a \code{pathtree_group} (a named
#'   list of \code{pathtree}s). Default \code{NULL} (single tree).
#'
#' @return For a single fit, a \code{pathtree} object (described below).
#'   For a grouped fit (\code{group =} supplied, or a grouped family
#'   object passed in) a \code{pathtree_group}: a named list of
#'   \code{pathtree}s, one per group, in the group's key order, with its
#'   own \code{print} and \code{as.data.frame} methods. A single
#'   \code{pathtree} is a list with components
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
#' Nodes whose total count falls below \code{min_count} are not created.
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
                         min_count = 5L,
                         smoothing = "floor",
                         alphabet  = NULL,
                         weights   = NULL,
                         group     = NULL) {
  if (!is.numeric(max_depth) || length(max_depth) != 1L ||
      is.na(max_depth) || max_depth < 0)
    stop("'max_depth' must be a single non-negative integer.",
         call. = FALSE)
  if (!is.numeric(min_count) || length(min_count) != 1L ||
      is.na(min_count) || min_count < 1)
    stop("'min_count' must be a single integer >= 1.", call. = FALSE)

  ## ---- Grouped fit -> a pathtree_group (one tree per group) ----------
  ## A grouped family object (netobject_group / group_tna / a named list
  ## of family objects) dispatches over its elements; an explicit
  ## group = (a metadata column name or a per-sequence vector) splits a
  ## single dataset. Both share one alphabet so the trees are
  ## comparable, and return a pathtree_group.
  if (.ct_is_group(data)) {
    if (!is.null(group))
      stop("Pass either a grouped object or 'group =', not both.",
           call. = FALSE)
    if (!is.null(weights))
      stop("'weights' is not supported with a grouped object; each ",
           "element carries its own weights (e.g. an 'stslist'), or fit ",
           "each group separately.", call. = FALSE)
    elements <- .ct_group_elements(data)
    alpha    <- .ct_group_alphabet(elements, alphabet)
    trees <- lapply(elements, function(e)
      context_tree(e, max_depth = max_depth, min_count = min_count,
                   smoothing = smoothing, alphabet = alpha))
    gv <- attr(data, "group")
    return(.ct_as_group(trees,
      if (is.character(gv) && length(gv) == 1L) gv else NA_character_))
  }
  if (!is.null(group)) {
    parts      <- .ct_group_split_by(data, group)
    alpha      <- if (is.null(alphabet)) parts$alphabet else alphabet
    ## Split a supplied weights vector across groups by the same row
    ## index so each per-group fit stays weighted (was silently dropped).
    w_by_group <- .ct_split_weights(weights, parts)
    trees <- Map(function(d, w)
      context_tree(d, max_depth = max_depth, min_count = min_count,
                   smoothing = smoothing, alphabet = alpha, weights = w),
      parts$subsets, w_by_group)
    return(.ct_as_group(trees, parts$var))
  }

  if (.ct_is_netobject(data)) {
    ## Unwrap before row-count / weight detection so they see the
    ## sequence frame, not the 18-slot network list. The network's
    ## node set becomes the default alphabet unless the caller
    ## overrides it explicitly.
    no <- .ct_unwrap_netobject(data)
    data <- no$data
    if (is.null(alphabet)) alphabet <- no$alphabet
  }
  n_rows_in <- .ct_data_n_rows(data)
  if (is.null(weights)) weights <- .ct_data_weights(data)

  data  <- .ct_coerce(data)
  trajs <- .ct_traj(data)
  if (length(trajs) == 0L)
    stop("No usable sequences after coercion.", call. = FALSE)

  if (!is.null(weights)) {
    if (!is.numeric(weights))
      stop("'weights' must be a numeric vector.", call. = FALSE)
    if (length(weights) != n_rows_in)
      stop("'weights' must have length equal to number of input ",
           "sequences (got ", length(weights), ", expected ",
           n_rows_in, ").", call. = FALSE)
    if (anyNA(weights))
      stop("'weights' must not contain NA.", call. = FALSE)
    if (any(weights < 0))
      stop("'weights' must be non-negative.", call. = FALSE)
    weights <- weights[attr(trajs, "idx")]
  }

  alphabet  <- if (is.null(alphabet)) .ct_alphabet(trajs) else alphabet
  max_depth <- as.integer(max_depth)
  min_count <- as.integer(min_count)
  sm        <- .pt_resolve_smoothing(smoothing)

  ## Smooth top-down so each node's parent prob is available when the
  ## node itself is computed (kneser_ney/witten_bell/jelinek_mercer
  ## interpolate with the parent).
  nodes <- list()
  for (d in seq.int(0L, max_depth)) {
    counts_d <- .ct_count_table(trajs, depth = d, alphabet = alphabet,
                                  weights = weights)
    keep <- vapply(counts_d, sum, numeric(1)) >= min_count
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
      nmin       = min_count,
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
#'   \code{count}, \code{likely_next}, \code{next_probability},
#'   \code{divergence}, \code{changes_prediction}. See
#'   \code{\link{pathtree_pathways}}.
#' @export
as.data.frame.pathtree <- function(x, row.names = NULL,
                                    optional = FALSE, ...) {
  pathtree_pathways(x, ...)
}

#' Print a Group of Context Trees
#'
#' @param x A \code{pathtree_group} (named list of \code{pathtree}s).
#' @param ... Ignored.
#' @return \code{x} invisibly.
#' @export
print.pathtree_group <- function(x, ...) {
  gv <- attr(x, "group")
  by <- if (!is.null(gv) && !is.na(gv)) sprintf(" by '%s'", gv) else ""
  cat(sprintf("<pathtree_group>  %d groups%s\n", length(x), by))
  for (nm in names(x)) {
    t <- x[[nm]]
    cat(sprintf("  %-14s %3d nodes, depth <= %d, %d seq, %d obs%s\n",
                nm, length(t$nodes), t$max_depth, t$n_seq, t$n_obs,
                if (isTRUE(t$pruned)) "  [pruned]" else ""))
  }
  invisible(x)
}

#' Coerce a Group of Trees to One Tidy Data Frame
#'
#' @description
#' Row-binds each group's \code{\link{pathtree_pathways}} table, tagged
#' with a leading \code{group} column, so the whole batch is one tidy
#' frame ready to filter, facet, or join.
#'
#' @param x A \code{pathtree_group}.
#' @param row.names,optional Ignored.
#' @param ... Forwarded to \code{\link{pathtree_pathways}()}.
#'
#' @return A data.frame: the canonical pathway columns with a leading
#'   \code{group} column identifying the source tree.
#' @export
as.data.frame.pathtree_group <- function(x, row.names = NULL,
                                         optional = FALSE, ...) {
  parts <- lapply(names(x), function(nm) {
    d <- pathtree_pathways(x[[nm]], ...)
    if (nrow(d) == 0L) return(NULL)
    cbind(group = nm, d, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, parts)
  if (is.null(out))
    return(cbind(group = character(0), .pt_empty_pathways_df(),
                 stringsAsFactors = FALSE))
  rownames(out) <- NULL
  out
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
  cat(sprintf("  smoothing: %s   min_count = %d\n", sm_label, x$nmin))
  if (isTRUE(x$pruned) && !is.null(x$pruning))
    cat(sprintf("  pruned by: %s   alpha = %s\n",
                x$pruning$criterion, format(x$pruning$alpha)))
  local_root <- attr(x, "local_root")
  if (!is.null(local_root))
    cat(sprintf("  subtree of: %s\n",
                if (identical(local_root, .ROOT)) .ROOT_LABEL else local_root))
  if (length(x$nodes) == 0L) return(invisible(x))

  render <- function(parent, prefix) {
    children <- sort(x$edges$child[x$edges$parent == parent])
    for (i in seq_along(children)) {
      child  <- children[i]; last <- i == length(children)
      branch <- if (last) "`-- " else "|-- "
      info   <- x$nodes[[child]]
      modal  <- alpha[which.max(info$prob)]
      label  <- if (identical(child, .ROOT)) .ROOT_LABEL else
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
    cat(sprintf("%-8s  n=%-5s  -> %s (%.*f)\n",
                .ROOT_LABEL, format(as.integer(root$n)), modal, digits,
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
#'   \code{depth}, \code{count}, \code{likely_next},
#'   \code{next_probability}, \code{divergence},
#'   \code{changes_prediction}), re-sorted by \code{(depth, -count)} so
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
