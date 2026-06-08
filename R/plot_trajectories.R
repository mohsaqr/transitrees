# ---- Forward trajectory tree ----

#' Forward trajectory-tree node table (internal, testable)
#'
#' Builds the forward prefix tree: every prefix occurring at least
#' \code{min_count} times becomes a node, connected to its one-shorter
#' prefix (the root for first moves). \code{count} is the prefix
#' frequency; \code{eprob} is the model's \eqn{P(\text{last move} \mid
#' \text{history})} read via \code{query_pathway()} (the root distribution
#' for first moves). Prefixes whose parent did not survive the count
#' filter are dropped so the returned tree is always connected to the root.
#'
#' @return A \code{data.frame} with columns \code{node}, \code{parent},
#'   \code{depth}, \code{count}, \code{last}, \code{eprob}.
#' @noRd
.trajectory_data <- function(tree, min_count = 4L) {
  stopifnot(inherits(tree, "transitiontrees"))
  states <- tree$alphabet
  seqs   <- .ct_traj(tree$data)
  if (length(seqs) == 0L)
    stop("Tree carries no sequences to draw.", call. = FALSE)

  prefixes <- unlist(lapply(seqs, function(s)
    vapply(seq_along(s), function(k) paste(s[seq_len(k)], collapse = " -> "),
           character(1))), use.names = FALSE)
  pc   <- table(prefixes)
  keep <- names(pc)[pc >= as.integer(min_count)]
  if (length(keep) == 0L)
    stop(sprintf("No prefix occurs at least %d times; lower min_count.",
                 as.integer(min_count)), call. = FALSE)

  parent_of <- function(p) {
    m <- strsplit(p, " -> ", fixed = TRUE)[[1L]]
    if (length(m) == 1L) .ROOT_LABEL else paste(m[-length(m)], collapse = " -> ")
  }
  last_of <- function(p)
    utils::tail(strsplit(p, " -> ", fixed = TRUE)[[1L]], 1L)

  node  <- c(.ROOT_LABEL, keep)
  depth <- c(0L, lengths(strsplit(keep, " -> ", fixed = TRUE)))
  count <- c(length(seqs), as.integer(pc[keep]))
  par   <- c(NA_character_, vapply(keep, parent_of, character(1)))
  last  <- c(NA_character_, vapply(keep, last_of,   character(1)))
  names(depth) <- names(count) <- names(par) <- names(last) <- node
  ## drop any prefix whose parent was filtered out (keep the tree connected)
  ok   <- node == .ROOT_LABEL | par %in% node
  node <- node[ok]; depth <- depth[node]; count <- count[node]
  par  <- par[node]; last <- last[node]

  root_prob <- stats::setNames(as.numeric(tree$nodes[[.ROOT]]$prob), states)
  pred_of <- function(n) {
    if (n == .ROOT_LABEL) return(NA_real_)
    p <- par[[n]]
    if (p == .ROOT_LABEL) return(unname(root_prob[[last[[n]]]]))
    hist <- utils::tail(strsplit(p, " -> ", fixed = TRUE)[[1L]], tree$max_depth)
    pr <- tryCatch(query_pathway(tree, paste(hist, collapse = " -> ")),
                   error = function(e) NULL)
    if (is.null(pr)) return(NA_real_)
    val <- pr[[last[[n]]]]
    if (is.null(val)) NA_real_ else as.numeric(val)
  }
  eprob <- vapply(node, pred_of, numeric(1))

  data.frame(node = node, parent = unname(par[node]), depth = depth[node],
             count = count[node], last = last[node], eprob = eprob[node],
             stringsAsFactors = FALSE, row.names = NULL)
}

#' Forward Trajectory Tree (Prefix Tree)
#'
#' @description
#' Where \code{\link{plot.transitiontrees}} draws the fitted context tree
#' \emph{backwards} (each node is a suffix, the most-recent move), this
#' draws the same prompts \emph{forwards}: a prefix tree that starts at a
#' common root and follows each sequence move by move in time. The one
#' tree can be coloured two ways:
#' \describe{
#'   \item{\code{measure = "frequency"}}{node and edge colour and width
#'     encode how many sequences walk each path.}
#'   \item{\code{measure = "predictability"}}{colour encodes
#'     \eqn{P(\text{move} \mid \text{history})} from the model
#'     (\code{tree}); edge width still encodes flow.}
#' }
#' Higher values are drawn darker.
#'
#' @param tree A \code{transitiontrees} (pass a pruned tree to read
#'   predictability off the pruned model).
#' @param measure One of \code{"frequency"} (default) or
#'   \code{"predictability"}.
#' @param min_count Integer. Keep only prefixes occurring at least this
#'   many times. Default 4.
#'
#' @return A ggplot object.
#'
#' @details
#' The predictability of a node's last move is the model's conditional
#' probability of that move given the preceding history, truncated to the
#' tree's \code{max_depth} and read via \code{\link{query_pathway}()}
#' (the empty history for a first move uses the root distribution). This
#' is the forward-reading complement to the backward context tree; it is
#' a visualisation, so it depends on \code{ggforce} (in Suggests) for the
#' rounded node glyphs and errors with an install hint if it is missing.
#'
#' @examples
#' \donttest{
#' seqs <- replicate(120, sample(c("A", "B", "C"), 8, replace = TRUE),
#'                   simplify = FALSE)
#' tree   <- context_tree(seqs, max_depth = 3L, min_count = 3L)
#' pruned <- prune_tree(tree)
#' if (requireNamespace("ggforce", quietly = TRUE)) {
#'   plot_trajectories(tree,   measure = "frequency")
#'   plot_trajectories(pruned, measure = "predictability")
#' }
#' }
#' @seealso \code{\link{plot.transitiontrees}} for the backward context
#'   tree, \code{\link{plot_pruning}} for the suffix-chain view.
#' @export
plot_trajectories <- function(tree,
                              measure = c("frequency", "predictability"),
                              min_count = 4L) {
  stopifnot(inherits(tree, "transitiontrees"))
  measure <- match.arg(measure)
  if (!requireNamespace("ggforce", quietly = TRUE))
    stop("plot_trajectories() needs the 'ggforce' package; install it with ",
         "install.packages(\"ggforce\").", call. = FALSE)

  states <- tree$alphabet
  D     <- .trajectory_data(tree, min_count = min_count)
  node  <- D$node
  par   <- stats::setNames(D$parent, D$node)
  last  <- stats::setNames(D$last,   D$node)
  depth <- stats::setNames(D$depth,  D$node)
  count <- stats::setNames(D$count,  D$node)
  eprob <- stats::setNames(D$eprob,  D$node)

  ## --- vertical layout: leaves stack, parents centre on children ------
  kids <- split(node[node != .ROOT_LABEL], par[node != .ROOT_LABEL])
  yenv <- new.env(); leaf_y <- 0
  assign_y <- function(nd) {
    k <- kids[[nd]]
    y <- if (is.null(k)) { leaf_y <<- leaf_y + 1; leaf_y }
         else mean(vapply(k, assign_y, numeric(1)))
    assign(nd, y, envir = yenv); y
  }
  assign_y(.ROOT_LABEL)
  ypos <- vapply(node, function(n) get(n, envir = yenv), numeric(1))

  L <- data.frame(node = node, x = depth[node], y = ypos[node],
                  depth = depth[node], count = count[node],
                  last = last[node], eprob = eprob[node],
                  stringsAsFactors = FALSE)
  E <- L[L$node != .ROOT_LABEL, , drop = FALSE]
  E$px <- L$x[match(par[E$node], L$node)]
  E$py <- L$y[match(par[E$node], L$node)]

  ## cosine-smoothed curved edges
  tt <- seq(0, 1, length.out = 40); sm <- (1 - cos(pi * tt)) / 2
  elb <- do.call(rbind, lapply(seq_len(nrow(E)), function(i)
    data.frame(grp = i, x = E$px[i] + tt * (E$x[i] - E$px[i]),
               y = E$py[i] + sm * (E$y[i] - E$py[i]))))

  NB <- L[L$node != .ROOT_LABEL, , drop = FALSE]
  RT <- L[L$node == .ROOT_LABEL, , drop = FALSE]

  ## value driving colour: count or predictability
  nodeval <- if (measure == "frequency") NB$count else NB$eprob
  lim     <- if (measure == "frequency") NULL else c(0, 1)
  num     <- if (measure == "frequency") sprintf("n=%d", NB$count) else
    sprintf("%.0f%%", 100 * NB$eprob)
  leg     <- if (measure == "frequency") "sequences" else "P(move | history)"
  title   <- if (measure == "frequency")
    "Trajectory tree - by frequency" else
    "Trajectory tree - by predictability"
  sub     <- if (measure == "frequency")
    "move inside, count below; colour & width = sequences (higher = darker)" else
    "colour = P(move | history) from the model; width = sequences"

  ## fixed-size rounded node glyphs (variable size triggers a ggforce fill bug)
  HWX <- 0.17; HWY <- 0.26; RAD <- grid::unit(4, "mm")
  sqf <- function(df) do.call(rbind, lapply(seq_len(nrow(df)), function(i)
    data.frame(id = df$node[i], x = df$x[i] + c(-1, 1, 1, -1) * HWX,
               y = df$y[i] + c(-1, -1, 1, 1) * HWY,
               val = df$valcol[i], stringsAsFactors = FALSE)))

  d <- NB; d$valcol <- nodeval; poly <- sqf(d)
  rt <- RT; rt$valcol <- NA_real_
  vn  <- (nodeval - min(nodeval)) / (max(nodeval) - min(nodeval) + 1e-12)
  txt <- ifelse(vn > 0.45, "white", "grey15")
  ecol <- nodeval[match(E$node, NB$node)][elb$grp]
  ew   <- E$count[elb$grp]

  ggplot2::ggplot() +
    ggplot2::geom_path(
      data = elb, ggplot2::aes(x = .data$x, y = .data$y, group = .data$grp,
                               colour = ecol, linewidth = ew),
      lineend = "round") +
    ggplot2::scale_linewidth_continuous(range = c(0.4, 3.6),
                                        name = "flow (sequences)") +
    ggforce::geom_shape(
      data = poly, ggplot2::aes(x = .data$x, y = .data$y, group = .data$id,
                                fill = .data$val),
      radius = RAD, colour = "grey30", linewidth = 0.3) +
    ggforce::geom_shape(
      data = sqf(transform(rt, valcol = NA_real_)),
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$id),
      radius = RAD, fill = "grey15") +
    ggplot2::geom_text(data = NB, ggplot2::aes(x = .data$x, y = .data$y,
                                               label = .data$last),
                       colour = txt, size = 2.5, fontface = "bold") +
    ggplot2::geom_text(data = data.frame(NB, num = num),
                       ggplot2::aes(x = .data$x, y = .data$y, label = .data$num),
                       vjust = 1, nudge_y = -0.42, size = 2.5,
                       colour = "grey30") +
    ggplot2::scale_fill_viridis_c(option = "C", direction = -1, limits = lim,
                                  name = leg) +
    ggplot2::scale_colour_viridis_c(option = "C", direction = -1, limits = lim,
                                    guide = "none") +
    ggplot2::scale_x_continuous(breaks = 0:max(L$depth)) +
    ggplot2::labs(x = "move number (time ->)", y = NULL,
                  title = title, subtitle = sub) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40", size = 9))
}
