#' prune.abtree
#'
#' @param object an object of class 'abtree' returned by MakeTree
#' @param valid.data a new data frame containing the variables used in MakeTree
#'
#' @return a pruned tree (object of class 'abtree')
#' @importFrom rpart prune
#' @exportClass abtree
#' @export
prune.abtree <- function(object, valid.data) {
  # if (!exists("cp")) {
  #   cp <- object$cp.table$cp.param[which.max(object$cp.table$profit)]
  # }
  if (length(object$cp.table) == 0) {
    return(object)
  }
  data <- ParseFormula(object$formula, valid.data)
  data$trt <- TreatToAB(data$trt, object$treatment)
  # print(head(data$trt))

  x <- data.matrix(data$x)
  for (j in which(object$var.types[1,]>=1)) {
    x[,j] <- as.numeric(x[,j]-1)
  }

  out <- rcpp_Prune(cbind(matrix(unlist(object$tree), ncol=14, byrow=T)),
                    data$y, x, data$trt, as.integer(object$var.types[2,]),
                    cbind(matrix(unlist(object$cp.table), ncol=2, byrow=T)))

  object$tree <- out$tree
  object$cp.table <- out$cp.table
  return(object)
}
