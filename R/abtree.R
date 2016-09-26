#' abtree
#'
#' @param formula an expression in the form of y ~ trt | cov1 + cov2 + ...
#' @param data the data frame where the variables reside
#' @param min.bucket the minimum number of observations in each leaf
#' @param min.split the minimum number of observations to return to the user
#'
#' @return tree: a data frame containing the nodes of the tree
#'         treatments: a mapping of A and B to the treatments passed in
#' @export
#'
#' @examples
#' \dontrun{
#' tree <- abtree(y ~ grp | hour + browser, data=x)
#' }
abtree <- function(formula, data, min.bucket=10, min.split=30,
                     max.depth=5) {
  m <- ParseFormula(formula, data)
  x.types <- vapply(m$x, class, character(1))

  # ========================  Conditions  ========================  # 
  # Question: why don't we just convert them all to factors?
  if (any(x.types == "character"))
    stop("At least one predictor is of type 'character'. Please convert to factor.")
  if (class(m$trt) != "factor")
    stop("The treatment variable should be of type 'factor'.")
  if (class(m$y) == "character")
    stop("The response variable cannot be of type 'character'. Please convert to factor.")

  x.levels <- sapply(m$x, function(t) levels(t))
  ncat <- sapply(x.levels, length)
  y <- as.numeric(m$y)
  x <- data.matrix(m$x) # converts everything to numeric, which is what we want
                        # since the Rcpp code will take a NumericMatrix
  x[,x.types=="factor"] <- x[,x.types=="factor"]-1L # correcting for 0-index
  trt <- as.integer(m$trt)-1L # correcting for 0-index
  ord <- apply(x, 2, order)-1L # correcting for 0-index
  trt.levels <- levels(m$trt)
  out <- rcpp_BuildTree(y, x, trt, ord,
                        ncat, length(trt.levels),
                        as.integer(min.bucket), # should these be here?
                        as.integer(min.split),  # or should there
                        as.integer(max.depth))  # be conditions?

  class(out) <- "abtree"
  out$trt.levels <- trt.levels
  # TODO: make this lapply(m$x, levels)
  out$cat.levels <- lapply(1:ncol(m$x), function(c) levels(m$x[,c]))
  out$formula <- formula
  out$x.types <- x.types
  out$x.levels <- x.levels
  out$ncat <- ncat
  out$y.name <- m$y.name
  out$trt.name <- m$trt.name
  out$frame <- FormatTree(out)
  out$cp.table <- matrix(unlist(out$cp.table), ncol=2, byrow=T)
  out
}