#' transitiontrees: Context Trees and Variable-Order Markov Models
#'
#' Fits prediction suffix trees (Ron, Singer & Tishby 1996) and
#' variable-order Markov models from categorical sequence data.
#' Provides multiple pruning criteria, per-context Kullback-Leibler
#' diagnostics, smoothing, prediction, and tree visualisation.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import ggplot2
#' @importFrom stats predict simulate logLik nobs AIC BIC qchisq setNames
#' @importFrom utils tail head capture.output modifyList
#' @importFrom utils globalVariables
## usethis namespace: end
NULL

## ggraph evaluates `weight = <col>` via NSE inside a tbl_graph;
## declare the column names here to silence R CMD check.
utils::globalVariables(c("count", "leaf_weight"))
