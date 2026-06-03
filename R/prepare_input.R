# ---- prepare_input(): long-format -> wide sequence frame ----
#
# Base-R reimplementation of the long->sequence reshaping done by
# tna::prepare_data() / Nestimate::build_network(). pathtree stays pure
# base R, so the dplyr/tidyr machinery of those packages is replaced by
# order() / ave() / matrix-indexing, but the session-splitting logic is
# identical: within each actor, a new session begins when the gap to the
# previous timestamp exceeds `time_threshold` seconds.

#' Parse a time column to POSIXct.
#' @noRd
.pt_parse_time <- function(x, format = NULL, is_unix_time = FALSE,
                           unix_time_unit = "seconds") {
  if (inherits(x, c("POSIXct", "POSIXt"))) return(as.POSIXct(x))
  if (inherits(x, "Date"))                 return(as.POSIXct(x))
  if (is.numeric(x) || isTRUE(is_unix_time)) {
    div <- switch(unix_time_unit,
                  seconds = 1, milliseconds = 1e3, microseconds = 1e6,
                  stop("'unix_time_unit' must be 'seconds', ",
                       "'milliseconds', or 'microseconds'.", call. = FALSE))
    return(as.POSIXct(as.numeric(x) / div, origin = "1970-01-01",
                      tz = "UTC"))
  }
  x <- as.character(x)
  if (!is.null(format))
    return(as.POSIXct(x, format = format, tz = "UTC"))
  parsed <- as.POSIXct(x, tz = "UTC")          # ISO-8601 / standard
  if (all(is.na(parsed) == is.na(x)) == FALSE || anyNA(parsed[!is.na(x)]))
    stop("Could not parse 'time'; pass an explicit 'format' ",
         "(see strptime) or 'is_unix_time = TRUE'.", call. = FALSE)
  parsed
}

#' Reshape Long Event Data into a Wide Sequence Frame
#'
#' @description
#' Turns a long, one-row-per-event table into the wide, one-row-per-
#' sequence character frame that \code{\link{context_tree}()} consumes.
#' Events are grouped by \code{actor} and ordered by \code{time} (or
#' \code{order}); when \code{time} is given, each actor's events are
#' split into \strong{sessions} whenever the gap to the previous event
#' exceeds \code{time_threshold} seconds. This mirrors the timestamp /
#' session logic of \code{tna::prepare_data()} and
#' \code{Nestimate::build_network()}, in pure base R.
#'
#' @param data A long-format \code{data.frame}, one row per event.
#' @param actor Character. Column(s) naming the unit each sequence
#'   belongs to (e.g. a user id). Several columns are combined with
#'   \code{"-"}. \code{NULL} (default) treats the whole table as one
#'   actor.
#' @param time Character. Column holding the event timestamp, used both
#'   to order events and to split sessions by \code{time_threshold}.
#'   Numeric values are read as Unix time. \code{NULL} (default) uses
#'   \code{order} (or row order) instead and does no session splitting.
#' @param action Character. Column holding the event's state / code —
#'   the symbol that becomes a cell of the sequence. Required.
#' @param order Character. Optional column giving an explicit within-
#'   actor ordering (used when \code{time} is absent, or to break time
#'   ties). Defaults to row order.
#' @param session Character. Optional column giving an explicit session
#'   id within an actor. If supplied, sessions are taken from it directly
#'   and no time-gap splitting is done.
#' @param time_threshold Numeric. Seconds; a gap larger than this starts
#'   a new session. Default \code{900} (15 minutes), matching
#'   \code{tna::prepare_data()}.
#' @param format Character. Optional \code{\link{strptime}} format for a
#'   character \code{time} column.
#' @param is_unix_time Logical. Force \code{time} to be read as Unix
#'   time. Default \code{FALSE}.
#' @param unix_time_unit One of \code{"seconds"} (default),
#'   \code{"milliseconds"}, \code{"microseconds"}.
#'
#' @return A wide character \code{data.frame}, one row per sequence
#'   (session), columns \code{T1, T2, ...} holding the ordered states and
#'   trailing \code{NA}s past the end of each sequence. Row names are the
#'   session ids. Pass it straight to \code{\link{context_tree}()}.
#'
#' @examples
#' long <- data.frame(
#'   user  = c("a","a","a","a","b","b"),
#'   t     = as.POSIXct("2020-01-01 09:00:00", tz = "UTC") +
#'             c(0, 60, 3600, 3660, 0, 30),
#'   state = c("X","Y","X","Z","Y","X"),
#'   stringsAsFactors = FALSE)
#' ## one-hour gap splits user a into two sessions
#' wide <- prepare_input(long, actor = "user", time = "t", action = "state")
#' wide
#' context_tree(wide, max_depth = 2L, min_count = 1L)
#'
#' @seealso \code{\link{context_tree}}
#' @export
prepare_input <- function(data, actor = NULL, time = NULL, action = NULL,
                         order = NULL, session = NULL,
                         time_threshold = 900,
                         format = NULL, is_unix_time = FALSE,
                         unix_time_unit = c("seconds", "milliseconds",
                                            "microseconds")) {
  stopifnot(is.data.frame(data), nrow(data) > 0L)
  unix_time_unit <- match.arg(unix_time_unit)
  if (is.null(action) || !action %in% names(data))
    stop("'action' must name a column of 'data'.", call. = FALSE)
  for (nm in c(actor, time, order, session))
    if (!is.null(nm) && !all(nm %in% names(data)))
      stop("'", paste(nm, collapse = "', '"),
           "' is not a column of 'data'.", call. = FALSE)
  if (!is.numeric(time_threshold) || length(time_threshold) != 1L ||
      time_threshold < 0)
    stop("'time_threshold' must be a single non-negative number.",
         call. = FALSE)

  n   <- nrow(data)
  act <- as.character(data[[action]])

  ## actor identifier (default: a single actor)
  default_actor <- is.null(actor)
  actor_vec <- if (default_actor) rep("session", n)
    else if (length(actor) > 1L)
      do.call(paste, c(lapply(actor, function(a) data[[a]]), sep = "-"))
    else as.character(data[[actor]])

  ord <- if (is.null(order)) seq_len(n) else data[[order]]

  if (!is.null(session)) {
    ## Explicit session column (Nestimate-style): order within
    ## actor+session, take the session id directly.
    sess <- paste0(actor_vec, " | ", as.character(data[[session]]))
    o    <- order(actor_vec, sess, ord)
    session_id <- sess[o]
    act_s      <- act[o]
  } else if (!is.null(time)) {
    ## Timestamp logic: sort within actor by (time, order), then split
    ## into sessions on gaps > time_threshold (tna::prepare_data rule).
    tt <- .pt_parse_time(data[[time]], format, is_unix_time, unix_time_unit)
    o  <- order(actor_vec, tt, ord)
    actor_s <- actor_vec[o]; tt_s <- tt[o]; act_s <- act[o]
    gap <- c(NA_real_, as.numeric(difftime(tt_s[-1L], tt_s[-n],
                                           units = "secs")))
    new_actor   <- c(TRUE, actor_s[-1L] != actor_s[-n])
    new_session <- new_actor | is.na(gap) | gap > time_threshold
    session_nr  <- stats::ave(as.integer(new_session), actor_s,
                              FUN = cumsum)
    session_id  <- if (default_actor) paste0("session", session_nr)
      else paste0(actor_s, " session", session_nr)
  } else {
    ## No time and no session: one sequence per actor, ordered by `order`.
    o <- order(actor_vec, ord)
    session_id <- actor_vec[o]
    act_s      <- act[o]
  }

  ## Position within each session, then scatter into a wide matrix.
  seqpos   <- stats::ave(seq_along(session_id), session_id,
                         FUN = seq_along)
  ## rows sorted by session id as a string, matching tna::prepare_data's
  ## arrange(.session_id) (so "session10" precedes "session2").
  sessions <- sort(unique(session_id))
  max_len  <- max(seqpos)
  wide <- matrix(NA_character_, length(sessions), max_len,
                 dimnames = list(sessions, paste0("T", seq_len(max_len))))
  wide[cbind(match(session_id, sessions), seqpos)] <- act_s
  as.data.frame(wide, stringsAsFactors = FALSE)
}
