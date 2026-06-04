# ---- Tests for prepare_input(): long -> wide sequence frame ----

.pd_long <- function() {
  base <- as.POSIXct("2020-01-01 09:00:00", tz = "UTC")
  data.frame(
    user  = c("a", "a", "a", "a", "b", "b", "b"),
    t     = base + c(0, 60, 3660, 3720, 0, 30, 90),  # a: 1h gap after 2nd
    state = c("X", "Y", "X", "Z", "Y", "X", "Y"),
    stringsAsFactors = FALSE)
}

test_that("time gap > time_threshold starts a new session", {
  w <- prepare_input(.pd_long(), actor = "user", time = "t", action = "state")
  ## user a splits into 2 sessions (gap of 1h > 900s); user b stays 1
  expect_equal(nrow(w), 3L)
  expect_equal(unname(w[["T1"]]), c("X", "X", "Y"))
  expect_equal(unname(w[["T2"]]), c("Y", "Z", "X"))
  expect_true(is.na(w[1L, "T3"]) && is.na(w[2L, "T3"]))
  expect_equal(w[3L, "T3"], "Y")
})

test_that("a large time_threshold keeps all of an actor's events in one session", {
  w <- prepare_input(.pd_long(), actor = "user", time = "t", action = "state",
                    time_threshold = 1e6)
  expect_equal(nrow(w), 2L)                 # one session per actor
  expect_equal(unname(w[["T1"]]), c("X", "Y"))
})

test_that("output is a wide character frame that context_tree() accepts", {
  w  <- prepare_input(.pd_long(), actor = "user", time = "t", action = "state")
  expect_s3_class(w, "data.frame")
  expect_true(all(grepl("^T[0-9]+$", names(w))))
  tr <- context_tree(w, max_depth = 2L, min_count = 1L)
  expect_s3_class(tr, "transitiontrees")
  expect_setequal(tr$alphabet, c("X", "Y", "Z"))
})

test_that("without time, one sequence per actor in `order`", {
  long <- data.frame(id = c("a","a","a","b","b"),
                     o  = c(2L,1L,3L,1L,2L),
                     s  = c("Y","X","Z","B","A"),
                     stringsAsFactors = FALSE)
  w <- prepare_input(long, actor = "id", order = "o", action = "s")
  expect_equal(nrow(w), 2L)
  expect_equal(unname(w[["T1"]]), c("X", "B"))   # ordered by `o`
  expect_equal(unname(w[["T2"]]), c("Y", "A"))
})

test_that("an explicit session column overrides time-gap splitting", {
  long <- data.frame(id = c("a","a","a","a"),
                     sess = c("m","m","n","n"),
                     s  = c("X","Y","X","Z"),
                     stringsAsFactors = FALSE)
  w <- prepare_input(long, actor = "id", session = "sess", action = "s")
  expect_equal(nrow(w), 2L)
  expect_equal(unname(w[["T1"]]), c("X", "X"))
  expect_equal(unname(w[["T2"]]), c("Y", "Z"))
})

test_that("numeric time is read as Unix seconds", {
  long <- data.frame(id = "a",
                     t  = c(0, 100, 2000),       # gap 1900s > 900 -> split
                     s  = c("X", "Y", "Z"),
                     stringsAsFactors = FALSE)
  w <- prepare_input(long, actor = "id", time = "t", action = "s")
  expect_equal(nrow(w), 2L)
})

test_that("invalid inputs error", {
  expect_error(prepare_input(.pd_long(), actor = "user", action = "nope"),
               "must name a column")
  expect_error(prepare_input(.pd_long(), actor = "nope", action = "state"),
               "not a column")
  expect_error(prepare_input(.pd_long(), action = "state", time = "t",
                            time_threshold = -1),
               "non-negative")
})

# ---- context_tree() long-format integration -----------------------------

test_that("context_tree() reshapes long data when action is named", {
  long <- data.frame(id = c("a","a","a","b","b"),
                     o  = c(1L,2L,3L,1L,2L),
                     s  = c("X","Y","X","Y","X"),
                     stringsAsFactors = FALSE)
  t1 <- context_tree(long, actor = "id", order = "o", action = "s",
                     max_depth = 2L, min_count = 1L)
  t2 <- context_tree(prepare_input(long, actor = "id", order = "o",
                                   action = "s"),
                     max_depth = 2L, min_count = 1L)
  expect_s3_class(t1, "transitiontrees")
  expect_identical(t1$nodes, t2$nodes)        # one-call == two-step
})

test_that("context_tree() errors if actor/time given without action", {
  long <- data.frame(id = c("a","a"), s = c("X","Y"),
                     stringsAsFactors = FALSE)
  expect_error(context_tree(long, actor = "id"), "action =")
  expect_error(context_tree(long, time = "id"), "action =")
})

test_that("context_tree() on already-wide data is unaffected", {
  m <- matrix(c("A","B","A","B","A","C"), nrow = 2, byrow = TRUE)
  expect_s3_class(context_tree(m, max_depth = 1L, min_count = 1L), "transitiontrees")
})
