#' print.abtree
#'
#' @param x object of class 'abtree'
#' @param ... optional arguments
#'
#' @return a data frame summarizing the fitted tree
#' @export
#'
print.abtree <- function(x, ...) {
  summary(x)
}