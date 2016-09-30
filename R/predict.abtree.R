#' predict.abtree
#'
#' @param obj an object of class 'abtree' returned by MakeTree
#' @param new.data a new data frame containing the variables used in MakeTree
#'
#' @return a vector of predicted optimal treatments for all of the observations
#' @export
#'
predict.abtree <- function(obj, new.data) {
  m <- ParseFormula(obj$formula, new.data)
  x.types <- vapply(m$x, class, character(1))
  y <- as.numeric(m$y)
  x <- data.matrix(m$x) # converts everything to numeric, which is what we want
  x[,x.types=="factor"] <- x[,x.types=="factor"]-1L # correcting for 0-index
  trt <- as.integer(m$trt)-1L # correcting for 0-index
  ord <- apply(x, 2, order)-1L # correcting for 0-index
  out <- rcpp_Predict(obj$cpp.ptr,
                      y, x, trt, obj$ncat)
  preds <- LETTERS[out$predict.trt+1L]
  # attr(preds, "table") <- matrix(unlist(out$test), ncol=18, byrow=T)[,15:18]
  preds
}