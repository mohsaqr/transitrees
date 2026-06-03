# ---- Tests for context_tree() ----

.simple_seqs <- function() {
  list(
    c("A", "B", "A", "B", "A"),
    c("B", "A", "B", "A", "B"),
    c("A", "A", "B", "B", "A"),
    c("B", "B", "A", "A", "B"),
    c("A", "B", "B", "A", "A"),
    c("B", "A", "A", "B", "B")
  )
}

# ---- structure ----

test_that("context_tree returns ctxtree with correct fields", {
  tree <- context_tree(.simple_seqs(), max_depth = 2L, min_count = 1L)
  expect_s3_class(tree, "transitrees")
  expect_named(tree, c("nodes", "edges", "alphabet",
                       "max_depth", "nmin",
                       "n_seq", "n_obs", "smoothing",
                       "pruned", "pruning", "data"))
  expect_setequal(tree$alphabet, c("A", "B"))
  expect_equal(tree$n_seq, 6L)
  expect_false(isTRUE(tree$pruned))
})

test_that("root node carries the marginal next-state distribution", {
  tree <- context_tree(.simple_seqs(), max_depth = 2L, min_count = 1L)
  root <- tree$nodes[["<root>"]]
  expect_false(is.null(root))
  expect_equal(sum(root$prob), 1, tolerance = 1e-12)
  expect_equal(length(root$prob), 2L)
  expect_equal(root$depth, 0L)
})

test_that("root node is retained even when nmin exceeds observations", {
  tree <- context_tree(.simple_seqs(), max_depth = 2L, min_count = 1000L)
  expect_true("<root>" %in% names(tree$nodes))
  expect_equal(length(tree$nodes), 1L)
  expect_equal(tree$max_depth, 0L)
})

test_that("node probabilities sum to 1 (after smoothing)", {
  tree <- context_tree(.simple_seqs(), max_depth = 3L, min_count = 1L)
  for (info in tree$nodes) {
    expect_equal(sum(info$prob), 1, tolerance = 1e-10)
  }
})

test_that("nmin filters out rare contexts", {
  tree_loose  <- context_tree(.simple_seqs(), max_depth = 2L, min_count = 1L)
  tree_strict <- context_tree(.simple_seqs(), max_depth = 2L, min_count = 5L)
  expect_lte(length(tree_strict$nodes), length(tree_loose$nodes))
})

test_that("depth respects max_depth", {
  tree <- context_tree(.simple_seqs(), max_depth = 2L, min_count = 1L)
  depths <- vapply(tree$nodes, function(x) x$depth, integer(1))
  expect_true(all(depths <= 2L))
})

# ---- input dispatch ----

test_that("context_tree accepts a character matrix", {
  m <- do.call(rbind, lapply(.simple_seqs(),
                              function(x) c(x, rep(NA, 6L - length(x)))))
  tree <- context_tree(m, max_depth = 2L, min_count = 1L)
  expect_s3_class(tree, "transitrees")
  expect_setequal(tree$alphabet, c("A", "B"))
})

test_that("context_tree accepts a wide data.frame", {
  df <- as.data.frame(do.call(rbind,
        lapply(.simple_seqs(),
               function(x) c(x, rep(NA, 6L - length(x))))),
        stringsAsFactors = FALSE)
  tree <- context_tree(df, max_depth = 2L, min_count = 1L)
  expect_s3_class(tree, "transitrees")
})

test_that("context_tree errors on unusable input", {
  expect_error(context_tree(matrix(1:6, 2, 3)),
               regexp = "data.*must be|wide data\\.frame")
  expect_error(context_tree(123L), regexp = "data.*must be|wide data\\.frame")
})

# ---- print/summary ----

test_that("print and summary dispatch correctly", {
  tree <- context_tree(.simple_seqs(), max_depth = 2L, min_count = 1L)
  expect_output(print(tree), "<transitrees>")
  s <- summary(tree)
  expect_s3_class(s, "summary.transitrees")
  expect_true(is.data.frame(s$table))
  expect_true(any(s$table$pathway == "(start)"))
  expect_output(print(s), "transitrees summary")
})

test_that("summary(tree)$table carries the canonical 7-column schema", {
  tree <- context_tree(.simple_seqs(), max_depth = 2L, min_count = 1L)
  s <- summary(tree)
  expect_named(s$table, c("pathway", "depth", "count",
                          "likely_next", "next_probability", "divergence",
                          "changes_prediction"))
  ## And it sorts by (depth, -count) — structural tree order.
  expect_equal(s$table$depth, sort(s$table$depth))
})

test_that("as.data.frame.transitrees is the canonical tidy view", {
  tree <- context_tree(.simple_seqs(), max_depth = 2L, min_count = 1L)
  df <- as.data.frame(tree)
  expect_s3_class(df, "data.frame")
  expect_named(df, c("pathway", "depth", "count",
                     "likely_next", "next_probability", "divergence",
                     "changes_prediction"))
  expect_equal(nrow(df), length(tree$nodes))
  expect_true(all(df$next_probability > 0 & df$next_probability <= 1))
  ## Identical to tree_pathways()
  expect_identical(df, tree_pathways(tree))
})

# ---- Nestimate netobject input ----

.fake_netobject <- function() {
  ## Mirrors the slots context_tree() reads from a real
  ## Nestimate::build_network() result: $data is a wide,
  ## trailing-NA-padded character frame; $nodes$label is the
  ## network's canonical node set (the default alphabet).
  d <- data.frame(
    T1 = c("A", "B", "A", "B", "A", "B"),
    T2 = c("B", "A", "A", "B", "B", "A"),
    T3 = c("A", "B", "B", "A", "B", "A"),
    T4 = c("B", "A", NA,  NA,  "A", "B"),
    stringsAsFactors = FALSE
  )
  structure(
    list(data  = d,
         nodes = data.frame(id = 1:2, label = c("A", "B"),
                            stringsAsFactors = FALSE)),
    class = c("netobject", "cograph_network")
  )
}

test_that("context_tree accepts a Nestimate netobject via $data", {
  no <- .fake_netobject()
  tn <- context_tree(no, max_depth = 2L, min_count = 1L)
  td <- context_tree(no$data, max_depth = 2L, min_count = 1L,
                      alphabet = no$nodes$label)
  drop_call <- function(x) x[setdiff(names(x), "call")]
  expect_s3_class(tn, "transitrees")
  expect_equal(drop_call(tn), drop_call(td))
})

test_that("netobject default alphabet comes from $nodes$label", {
  no <- .fake_netobject()
  ## Add a node never seen in $data: it must still enter the alphabet
  ## because the network node set is canonical.
  no$nodes <- data.frame(id = 1:3, label = c("A", "B", "C"),
                         stringsAsFactors = FALSE)
  tn <- context_tree(no, max_depth = 1L, min_count = 1L)
  expect_true("C" %in% tn$alphabet)
})

test_that("explicit alphabet overrides the netobject node set", {
  no <- .fake_netobject()
  tn <- context_tree(no, max_depth = 1L, min_count = 1L,
                      alphabet = c("A", "B", "Z"))
  expect_true("Z" %in% tn$alphabet)
})

test_that("netobject without a usable $data slot errors clearly", {
  bad <- structure(list(nodes = data.frame(label = "A")),
                    class = "netobject")
  expect_error(context_tree(bad), "carries no sequence data")
})

# ---- cograph network-family objects ----

test_that("a cograph_network carrying $data fits as is", {
  no <- .fake_netobject()
  ## same payload, cograph's class vector instead of netobject's
  cg <- structure(unclass(no), class = c("cograph_network", "list"))
  tc <- context_tree(cg, max_depth = 2L, min_count = 1L)
  td <- context_tree(no, max_depth = 2L, min_count = 1L)
  drop_call <- function(x) x[setdiff(names(x), "call")]
  expect_s3_class(tc, "transitrees")
  expect_equal(drop_call(tc), drop_call(td))
})

test_that("a pure-graph cograph object (no $data) errors with guidance", {
  ## nodes/edges/weights only — an aggregated transition network
  graph_only <- structure(
    list(nodes   = data.frame(id = 1:2, label = c("A", "B")),
         edges   = data.frame(from = 1L, to = 2L, weight = 3),
         weights = matrix(c(0, 3, 1, 0), 2),
         data    = NULL),
    class = c("cograph_network", "list"))
  expect_error(context_tree(graph_only),
               "carries no sequence data")
  expect_error(context_tree(graph_only),
               "cannot be recovered from edge weights")
})

# ---- sibling-family models accepted directly ----

test_that("a tna-style model (class 'tna', $data + $labels) fits direct", {
  ## tna shape: integer-coded $data, positional $labels (code k =
  ## state k), NO $nodes table. Must be accepted by class and decoded
  ## via $labels.
  tm <- structure(
    list(weights = matrix(0, 2, 2),
         inits   = c(plan = 0.5, discuss = 0.5),
         labels  = c("plan", "discuss"),
         data    = matrix(c(1L, 2L, 1L, 2L,
                            2L, 1L, 2L, 1L,
                            1L, 2L, 1L, NA), nrow = 4)),
    class = "tna")
  tr <- context_tree(tm, max_depth = 2L, min_count = 1L)
  expect_s3_class(tr, "transitrees")
  expect_setequal(tr$alphabet, c("plan", "discuss"))
  expect_false(any(c("1", "2") %in% tr$alphabet))   # decoded, not codes
})

test_that("an unclassed family-shaped list is detected structurally", {
  ## A future sibling that just follows the $data convention works
  ## with no code change here.
  obj <- list(data = data.frame(T1 = c("X","Y","X"),
                                T2 = c("Y","X","Y"),
                                stringsAsFactors = FALSE),
              labels = c("X", "Y"))
  tr <- context_tree(obj, max_depth = 1L, min_count = 1L)
  expect_s3_class(tr, "transitrees")
  expect_setequal(tr$alphabet, c("X", "Y"))
})

test_that("a plain list of character vectors is NOT mis-detected", {
  ## Regression: the structural detector must not swallow the normal
  ## ragged-list input form.
  lst <- list(c("A","B","A","B"), c("B","A","B","A"))
  tr  <- context_tree(lst, max_depth = 1L, min_count = 1L)
  expect_s3_class(tr, "transitrees")
  expect_setequal(tr$alphabet, c("A","B"))
})

test_that("sequences are extracted from a non-$data slot ($sequences)", {
  obj <- structure(
    list(sequences = data.frame(T1 = c("A","B","A","B"),
                                T2 = c("B","A","B","A"),
                                T3 = c("A","A","B","B"),
                                stringsAsFactors = FALSE),
         nodes = data.frame(label = c("A","B"))),
    class = "cograph_network")
  tr <- context_tree(obj, max_depth = 2L, min_count = 1L)
  expect_s3_class(tr, "transitrees")
  expect_setequal(tr$alphabet, c("A","B"))
})

test_that("sequences are extracted from an embedded netobject", {
  no  <- .fake_netobject()
  wrap <- structure(list(nodes = data.frame(label = c("A","B")),
                          wrapped = no),
                     class = c("cograph_network", "list"))
  tw <- context_tree(wrap, max_depth = 2L, min_count = 1L)
  td <- context_tree(no,   max_depth = 2L, min_count = 1L)
  drop_call <- function(x) x[setdiff(names(x), "call")]
  expect_equal(drop_call(tw), drop_call(td))
})

test_that("integer-coded sequence frame is decoded via $nodes labels", {
  ## How tna stores sequences (surfaced by cograph::as_cograph(<tna>)):
  ## an integer matrix in $data + an id/label table in $nodes. NA is
  ## trailing padding. Codes must map to labels, NOT be rejected as a
  ## "numeric matrix".
  obj <- structure(
    list(data  = matrix(c(1L, 2L, 1L, 2L,
                           2L, 1L, 2L, 1L,
                           1L, 1L, 2L, NA), nrow = 4),
         nodes = data.frame(id = 1:2, label = c("plan", "discuss"),
                            stringsAsFactors = FALSE)),
    class = c("cograph_network", "list"))
  tr <- context_tree(obj, max_depth = 2L, min_count = 1L)
  expect_s3_class(tr, "transitrees")
  expect_setequal(tr$alphabet, c("plan", "discuss"))
  ## decoded, not "1"/"2"
  expect_false(any(c("1", "2") %in% tr$alphabet))
})

test_that("numeric frame without a label table casts to character", {
  obj <- structure(
    list(data  = matrix(c(1L, 2L, 2L, 1L, 1L, 2L), nrow = 2),
         nodes = NULL),
    class = c("cograph_network", "list"))
  tr <- context_tree(obj, max_depth = 1L, min_count = 1L)
  expect_setequal(tr$alphabet, c("1", "2"))
})

test_that("a network object with an empty (0-row) $data errors", {
  empty_df <- structure(
    list(data  = data.frame(T1 = character(0)),
         nodes = data.frame(label = c("A", "B"))),
    class = c("netobject", "cograph_network"))
  expect_error(context_tree(empty_df), "carries no sequence data")
})

# ---- internal-NA handling: matrix and list paths agree ----

test_that("internal NAs are dropped, not turned into an 'NA' state", {
  ## Regression: the matrix/data.frame path used to keep interior NAs up
  ## to the last non-empty cell, so paste(traj) produced a literal "NA"
  ## context. The list path stripped them. Identical data must now yield
  ## identical trees regardless of container shape.
  m <- rbind(c("A", "B", "A", "B", "A"),
             c("A", NA,  "A", "B", "A"),
             c("B", "A", "B", "A", "B"))
  lst <- lapply(seq_len(nrow(m)), function(i) m[i, ])
  tr_m <- context_tree(m,   max_depth = 2L, min_count = 1L)
  tr_l <- context_tree(lst, max_depth = 2L, min_count = 1L)

  expect_false(any(grepl("NA", names(tr_m$nodes), fixed = TRUE)))
  expect_setequal(names(tr_m$nodes), names(tr_l$nodes))
  expect_equal(tr_m$n_obs, tr_l$n_obs)
  expect_equal(tr_m$n_seq, tr_l$n_seq)
})

# ---- integer-coded frame validation ----

test_that("0-based / out-of-range positional codes are rejected clearly", {
  ## Positional $labels are 1-based (code k = state k). A 0 would index
  ## labels[0] and silently drop, corrupting the frame.
  tm0 <- structure(
    list(labels = c("plan", "discuss"),
         data   = matrix(c(0L, 1L, 0L, 1L,
                           1L, 0L, 1L, 0L), nrow = 2, byrow = TRUE)),
    class = "tna")
  expect_error(context_tree(tm0, max_depth = 1L, min_count = 1L),
               "outside 1\\.\\.2")
  ## the equivalent 1-based frame still fits
  tm1 <- structure(
    list(labels = c("plan", "discuss"),
         data   = matrix(c(1L, 2L, 1L, 2L,
                           2L, 1L, 2L, 1L), nrow = 2, byrow = TRUE)),
    class = "tna")
  expect_setequal(context_tree(tm1, max_depth = 1L, min_count = 1L)$alphabet,
                  c("plan", "discuss"))
})

test_that("id-table codes absent from $nodes are rejected", {
  obj <- structure(
    list(data  = matrix(c(1L, 9L, 1L, 2L), nrow = 2),
         nodes = data.frame(id = 1:2, label = c("A", "B"),
                            stringsAsFactors = FALSE)),
    class = c("cograph_network", "list"))
  expect_error(context_tree(obj, max_depth = 1L, min_count = 1L),
               "absent from the \\$nodes id table")
})

# ---- argument validation ----

test_that("context_tree rejects malformed max_depth / nmin / weights", {
  m <- rbind(c("A", "B", "A", "B"), c("B", "A", "B", "A"))
  expect_error(context_tree(m, max_depth = -1L, min_count = 1L),
               "max_depth")
  expect_error(context_tree(m, max_depth = NA, min_count = 1L),
               "max_depth")
  expect_error(context_tree(m, max_depth = 2L, min_count = NA),
               "min_count")
  expect_error(context_tree(m, max_depth = 2L, min_count = 0L),
               "min_count")
  expect_error(context_tree(m, max_depth = 2L, min_count = 1L,
                            weights = c(1, NA)),
               "must not contain NA")
})
