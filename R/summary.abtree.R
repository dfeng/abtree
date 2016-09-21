#' summary.abtree
#'
#' @param obj of class 'abtree'
#'
#' @return a data frame summarizing the fitted tree
#' @export
#'
summary.abtree <- function(obj) {
  print(obj$tree)
  invisible(obj$tree)
}