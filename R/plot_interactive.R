# ---- Interactive HTML-widget tree (visNetwork) ----

#' Linear-scale a numeric vector into a [lo, hi] pixel range.
#' @noRd
.pt_scale_range <- function(v, lo, hi) {
  v   <- as.numeric(v)
  rng <- range(v, na.rm = TRUE)
  if (!is.finite(rng[1L]) || rng[1L] == rng[2L])
    return(rep((lo + hi) / 2, length(v)))
  lo + (v - rng[1L]) / (rng[2L] - rng[1L]) * (hi - lo)
}

#' @noRd
.plot_interactive <- function(tree,
                              point_size_range = c(10, 45),
                              edge_size_range  = c(1, 10),
                              width = NULL, height = NULL, ...) {
  if (!requireNamespace("visNetwork", quietly = TRUE))
    stop("plot(style = 'interactive') requires the 'visNetwork' ",
         "package. Install it, or use style = 'dendrogram' / ",
         "'horizontal' / 'icicle'.", call. = FALSE)
  if (length(tree$nodes) == 0L)
    stop("Cannot plot an empty tree.", call. = FALSE)

  alphabet <- tree$alphabet
  pal      <- .pt_state_palette(alphabet)
  ids      <- names(tree$nodes)

  ## ---- Nodes: size = count, fill = last state, tooltip = distribution.
  node_rows <- lapply(ids, function(ctx) {
    info    <- tree$nodes[[ctx]]
    is_root <- identical(ctx, .ROOT)
    last    <- if (is_root) .ROOT_LABEL else .pt_last_state(ctx)

    modal_idx <- which.max(info$prob)
    modal     <- alphabet[modal_idx]
    modal_p   <- info$prob[modal_idx]
    dist_html <- paste(sprintf("&nbsp;&nbsp;%s: %.3f", alphabet, info$prob),
                       collapse = "<br/>")

    title <- sprintf(
      paste0("<b>Pathway</b>: %s<br/><b>Depth</b>: %d<br/>",
             "<b>Count</b>: %d<br/><b>Modal next</b>: %s (%.3f)<br/>",
             "<b>Next-state distribution</b>:<br/>%s"),
      if (is_root) .ROOT_LABEL else ctx, info$depth, as.integer(info$n),
      modal, modal_p, dist_html)

    data.frame(
      id    = ctx,
      label = last,
      count = as.numeric(info$n),
      level = info$depth,
      color = if (is_root) "#222222" else unname(pal[[last]]),
      title = title,
      stringsAsFactors = FALSE
    )
  })
  nodes <- do.call(rbind, node_rows)
  nodes$size <- .pt_scale_range(nodes$count,
                                point_size_range[1L], point_size_range[2L])

  ## ---- Edges: width = child's count ("flow" down the branch).
  e <- tree$edges
  if (is.null(e) || nrow(e) == 0L) {
    edges <- data.frame(from = character(0), to = character(0),
                        width = numeric(0), stringsAsFactors = FALSE)
  } else {
    child_count <- vapply(e$child,
                          function(c) as.numeric(tree$nodes[[c]]$n), numeric(1))
    edges <- data.frame(
      from  = e$parent,
      to    = e$child,
      width = .pt_scale_range(child_count,
                              edge_size_range[1L], edge_size_range[2L]),
      stringsAsFactors = FALSE
    )
  }

  visNetwork::visNetwork(nodes, edges, width = width, height = height) |>
    visNetwork::visHierarchicalLayout(direction = "LR",
                                      sortMethod = "directed",
                                      levelSeparation = 200,
                                      nodeSpacing = 120) |>
    visNetwork::visNodes(shape = "dot",
                         font = list(size = 16, face = "sans-serif"),
                         borderWidth = 1,
                         color = list(border = "#444444",
                                      highlight = list(border = "#000000"))) |>
    visNetwork::visEdges(color = list(color = "#cfcfcf", highlight = "#555555"),
                         smooth = list(enabled = TRUE, type = "cubicBezier",
                                       roundness = 0.5),
                         arrows = "to") |>
    visNetwork::visOptions(highlightNearest = list(enabled = TRUE,
                                                   degree = 1, hover = TRUE)) |>
    visNetwork::visInteraction(dragView = FALSE, zoomView = FALSE,
                               dragNodes = TRUE, hover = TRUE)
}
