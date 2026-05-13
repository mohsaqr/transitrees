# ---- Interactive HTML-widget tree (collapsibleTree) ----

#' @noRd
.plot_interactive <- function(tree, width = NULL, height = NULL, ...) {
  if (!requireNamespace("collapsibleTree", quietly = TRUE))
    stop("plot(style = 'interactive') requires the 'collapsibleTree' ",
         "package. Install it, or use style = 'dendrogram' / ",
         "'horizontal' / 'icicle'.", call. = FALSE)
  if (length(tree$nodes) == 0L)
    stop("Cannot plot an empty tree.", call. = FALSE)

  alphabet <- tree$alphabet
  pal      <- .pt_state_palette(alphabet)

  ## Build the parent/child rows that collapsibleTreeNetwork expects:
  ## - column 1 = parent name (NA for the root row)
  ## - column 2 = child name (unique across all rows)
  ## - optional tooltip + fill colour columns.
  ##
  ## We use the full pathway string as the unique node name; the root
  ## row is named "(root)".
  rows <- lapply(names(tree$nodes), function(ctx) {
    info     <- tree$nodes[[ctx]]
    is_root  <- identical(ctx, .ROOT)
    parent_ctx <- .pt_parent_ctx(ctx)
    parent_lbl <- if (is_root || is.na(parent_ctx))      NA_character_
                  else if (identical(parent_ctx, .ROOT)) "(root)"
                  else                                    parent_ctx
    node_id  <- if (is_root) "(root)" else ctx
    last     <- if (is_root) "(root)" else .pt_last_state(ctx)

    modal_idx   <- which.max(info$prob)
    modal_state <- alphabet[modal_idx]
    modal_prob  <- info$prob[modal_idx]

    dist_html <- paste(
      sprintf("&nbsp;&nbsp;%s: %.3f", alphabet, info$prob),
      collapse = "<br/>"
    )

    tooltip <- sprintf(
      paste0("<b>Pathway</b>: %s<br/>",
             "<b>Depth</b>: %d<br/>",
             "<b>Count</b>: %d<br/>",
             "<b>Modal next</b>: %s (%.3f)<br/>",
             "<b>Next-state distribution</b>:<br/>%s"),
      if (is_root) "(root)" else ctx,
      info$depth, as.integer(info$n),
      modal_state, modal_prob, dist_html
    )

    fill_colour <- if (is_root) "#222222" else unname(pal[[last]])

    ## Pass only columns collapsibleTreeNetwork needs.
    ## data.tree (used internally) treats every additional column as a
    ## node attribute and renames any whose name collides with its
    ## reserved-names list (count, depth, level, position, siblings,
    ## children, path, …). Renames generate noisy warnings and, on
    ## some data.tree versions, trip a downstream length() dispatch.
    ## Keep the row strictly minimal: parent, name, tooltip, fill.
    data.frame(
      parent       = parent_lbl,
      name         = node_id,
      tooltip_html = tooltip,
      fill_colour  = fill_colour,
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)

  ## Root must come first.
  ord <- order(is.na(df$parent), decreasing = TRUE)
  df  <- df[ord, , drop = FALSE]
  rownames(df) <- NULL

  collapsibleTree::collapsibleTreeNetwork(
    df,
    tooltipHtml = "tooltip_html",
    fill        = "fill_colour",
    collapsed   = FALSE,
    zoomable    = TRUE,
    width       = width,
    height      = height,
    ...
  )
}
