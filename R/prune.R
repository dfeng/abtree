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

  if (obj$trt.levels != levels(m$trt))
    stop("The treatment variable levels have changed!")
  if (!all.equal(obj$x.levels, sapply(m$x, function(t) levels(t))))
    stop("At least one predictor has different levels!")

  out <- rcpp_Prune(cbind(matrix(unlist(obj$tree), ncol=14, byrow=T)),
                    m$y, m$x, m$trt, obj$ncat, length(obj$trt.levels),
                    cbind(matrix(unlist(obj$cp.table), ncol=2, byrow=T)))

  obj$tree <- out$tree
  obj$cp.table <- out$cp.table
  obj
}
