#' predict.abtree
#'
#' @param object an object of class 'abtree' returned by MakeTree
#' @param new.data a new data frame containing the variables used in MakeTree
#' @param type the type of data to return
#' @param ... optional arguments
#'
#' @return a vector of predicted optimal treatments for all of the observations
#' @export
#'
predict.abtree <- function(object, new.data, type="response", ...) {
  m <- ParseFormula(object$formula, new.data)
  # TODO: Handling missing data, via surrogate splits
  if (any(is.na(m$x))) {
    is_missing <- apply(m$x, 1, function(x) any(is.na(x)))
    m$x <- na.omit(m$x)
    warning("Missing values found in predictors. Rows containing missing values have been discarded.")
  }
  x.types <- vapply(m$x, class, character(1))
  y <- as.numeric(m$y)
  x <- data.matrix(m$x) # converts everything to numeric, which is what we want
  x[,x.types=="factor"] <- x[,x.types=="factor"]-1L # correcting for 0-index
  trt <- as.integer(m$trt)-1L # correcting for 0-index
  ord <- apply(x, 2, order)-1L # correcting for 0-index
  out <- rcpp_Predict(object$cpp.ptr,
                      y, x, trt, object$ncat, length(object$trt.levels))
  if (type=="response") {
    preds <- rep(NA, nrow(m$x))
    if (exists("is_missing")) {
      preds[!is_missing] <- object$trt.levels[out$predict.trt+1L]
    } else {
      preds <- object$trt.levels[out$predict.trt+1L]
    }
  } else if (type=="prob") {
    preds <- matrix(NA, nrow=nrow(m$x), ncol=length(object$trt.levels))
    if (exists("is_missing")) {
      preds[!is_missing,] <- out$predict.prob
    } else {
      preds <- out$predict.prob
    }
    colnames(preds) <- object$trt.levels
  }
  # attr(preds, "table") <- matrix(unlist(out$test), ncol=18, byrow=T)[,15:18]
  preds
}