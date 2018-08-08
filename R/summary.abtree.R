#' summary.abtree
#'
#' @param object of class 'abtree'
#' @param ... optional parameters
#'
#' @return a data frame summarizing the fitted tree
#' @export
#'
summary.abtree <- function(object, ...) {
  print(object$frame)
  invisible(object$frame)
}