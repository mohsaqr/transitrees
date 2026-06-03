#' Student engagement trajectories
#'
#' An example set of categorical learning-engagement trajectories used in
#' the \pkg{pathtree} examples and the \dQuote{trajectories} vignette.
#' Each row is one learner, each column a time-step, and each cell the
#' engagement state at that step. Trailing \code{NA}s mark the end of a
#' trajectory. This wide character matrix is exactly the shape
#' \code{\link{context_tree}()} consumes.
#'
#' @format A character matrix with 138 rows (learners) and 15 columns
#'   (time-steps). Three states: \code{"Active"}, \code{"Average"},
#'   \code{"Disengaged"}.
#'
#' @source Example engagement-trajectory data from the \code{mohsaqr}
#'   package family.
#'
#' @examples
#' data(trajectories)
#' dim(trajectories)
#' tree <- context_tree(trajectories)
#' tree
"trajectories"
