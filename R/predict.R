#' predict.abtree
#'
#' @param object an object of class 'abtree' returned by MakeTree
#' @param newdata a new data frame containing the variables used in MakeTree
#'
#' @return a vector of predicted optimal treatments for all of the observations
#' @export
#'
predict.abtree <- function(object, newdata) {

  data <- ParseFormula(object$formula, newdata)
  data$trt <- TreatToAB(data$trt, object$treatment)

  x <- data.matrix(data$x)
  for (j in which(object$var.types[1,]>=1)) {
    x[,j] <- as.numeric(x[,j]-1)
  }

  # figure out ordering
  ord <- apply(x, 2, order)-1
  out <- rcpp_Predict(cbind(matrix(unlist(object$tree), ncol=14, byrow=T)),
                      data$y, x, data$trt, as.integer(object$var.types[2,]))

  z <- LETTERS[out$trt+1]
  # DF: I forget what this table is, but I'm going to comment it out for now
  # attr(z, "table") <- matrix(unlist(out$test), ncol=18, byrow=T)[,15:18]
  return(z)
}


