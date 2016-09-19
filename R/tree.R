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
  if (any(x.types == "character"))
    stop("At least one predictor is of type 'character'. Please convert to factor.")
  if (class(m$trt) == "character")
    stop("The treatment variable should be of type 'factor'.")
  if (class(m$y) == "character")
    stop("The response variable cannot be of type 'character'. Please convert to factor.")

  ncat <- sapply(m$x, function(t) length(levels(t)))
  y <- as.numeric(m$y)
  x <- data.matrix(m$x) # converts everything to numeric, which is what we want
                        # since the Rcpp code will take a NumericMatrix
  x[,x.types=="factor"] <- x[,x.types=="factor"]-1 # correcting for 0-index
  trt <- as.integer(m$trt)

  ord <- apply(x, 2, order)-1 # correcting for 0-index
  out <- rcpp_BuildTree(y, x, trt, ord,
                        ncat,
                        as.integer(min.bucket),
                        as.integer(min.split),
                        as.integer(max.depth))

  class(out) <- "abtree"
  out$treatments <- levels(m$trt)

  # split_vars <- colnames(data$x)[unique(sapply(out$tree, function(x) ifelse(x[12]!=-1, x[2], NA))+1)]
  # split_vars <- split_vars[!is.na(split_vars)]
  # out$cat_levels <- list()

  # for (s in split_vars) {
  #   if (var.types[1,s] >= 1)  {
  #     out$cat_levels[[s]] <- levels(data$x[,s])
  #   }
  # }
  # out$formula <- formula
  # out$var.types <- var.types
  # out$trt <- colnames(data$trt)
  return(out)
}

#' print.abtree
#'
#' @param object of class 'abtree'
#'
#' @return a data frame summarizing the fitted tree
#' @export
#'
print.abtree <- function(object) {
  summary(object)
}


#' summary.abtree
#'
#' @param object of class 'abtree'
#'
#' @return a data frame summarizing the fitted tree
#' @export
#'
summary.abtree <- function(object) {
  cat("\nabtree Summary:\n\n")
  print(FormatTree(object))
  invisible(FormatTree(object))
}


