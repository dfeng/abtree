#' prune.abtree
#'
#' @param object an object of class 'abtree' returned by MakeTree
#' @param valid.data a new data frame containing the variables used in MakeTree
#'
#' @return a pruned tree (object of class 'abtree')
#' @importFrom rpart prune
#' @exportClass abtree
#' @export
prune.abtree <- function(obj, valid.data) {
  if (length(obj$cp.table) == 0)
    return(obj)

  m <- ParseFormula(obj$formula, valid.data)
  x.types <- vapply(m$x, class, character(1))
  y <- as.numeric(m$y)
  x <- data.matrix(m$x) # converts everything to numeric, which is what we want
  x[,x.types=="factor"] <- x[,x.types=="factor"]-1L # correcting for 0-index
  trt <- as.integer(m$trt)-1L # correcting for 0-index

  if (any(obj$trt.levels != levels(m$trt)))
    stop("The treatment variable levels have changed!")
  if (!all.equal(obj$x.levels, sapply(m$x, function(t) levels(t))))
    stop("At least one predictor has different levels!")

  out <- rcpp_Prune(obj$cpp.ptr,
                    y, x, trt, obj$ncat, length(obj$trt.levels),
                    obj$cp.table)

  obj$cpp.tree <- out$cpp.prune.tree
  obj$frame <- FormatTree(obj)
  obj$cp.table <- out$cp.table
  obj
}
