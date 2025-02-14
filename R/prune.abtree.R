#' prune.abtree
#'
#' @param tree an object of class 'abtree' returned by MakeTree
#' @param valid.data a new data frame containing the variables used in MakeTree
#' @param ... optional arguments
#'
#' @return a pruned tree (object of class 'abtree')
#' @importFrom rpart prune
#' @exportClass abtree
#' @export
prune.abtree <- function(tree, valid.data, ...) {
  if (length(tree$cp.table) == 0)
    return(tree)

  m <- ParseFormula(tree$formula, valid.data)
  if(any(is.na(m$x))) {
    m$x <- na.omit(m$x)
    warning("Missing values found in predictors. Rows containing missing values have been discarded.")
  }
  x.types <- vapply(m$x, class, character(1))
  y <- as.numeric(m$y)
  x <- data.matrix(m$x) # converts everything to numeric, which is what we want
  x[,x.types=="factor"] <- x[,x.types=="factor"]-1L # correcting for 0-index
  trt <- as.integer(m$trt)-1L # correcting for 0-index

  if (any(tree$trt.levels != levels(m$trt)))
    stop("The treatment variable levels have changed!")
  if (!all.equal(tree$x.levels, sapply(m$x, function(t) levels(t))))
    stop("At least one predictor has different levels!")

  out <- rcpp_Prune(tree$cpp.ptr,
                    y, x, trt, tree$ncat, tree$cp.table)

  tree$cpp.tree <- out$cpp.prune.tree
  tree$frame <- FormatTree(tree)
  tree$cp.table <- out$cp.table
  tree
}
