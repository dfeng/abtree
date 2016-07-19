#' MakeTree
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
#' tree <- MakeTree(y ~ grp | hour + browser, data=x)
#' }
MakeTree <- function(formula, data, min.bucket=10, min.split=30,
                     max.depth=5) {
  data <- ParseFormula(formula, data)
  AB <- sort(unique(data$trt)) # temporary
  data$trt <- TreatToAB(data$trt, AB)
  # print(head(data$trt))

  # coerce character types to factors
  types <- sapply(data$x, class)
  if (any(types == "character"))
    stop("Error: At least one predictor is of type 'character'. Please convert to factor.")

  # todo: What if a level doesn't exist in the new data?
  # todo: Or worse... if the levels are ordered differently?
  # figure out class of each covariate
  # and temporarily convert ordinals to numerics
  x <- data.matrix(data$x)
  var.types <- sapply(data$x, GetVarType)
  for (j in which(var.types[1,]>=1)) {
    x[,j] <- as.numeric(x[,j]-1)
  }

  # figure out ordering
  ord <- apply(x, 2, order)-1
  # split the data into train, test
  # data <- lapply(data, SplitData, sample(length(data$trt)), floor(length(data$trt)*7/10))
  # init a tree
  #tree <- Node$new("root")
  #BuildTree(tree,  AB, min.bucket, min.split)
  out <- rcpp_BuildTree(data$y, x, data$trt, ord,
                        as.integer(var.types[2,]),
                        as.integer(min.bucket), as.integer(min.split),
                        as.integer(max.depth))
  class(out) <- "abtree"
  out$treatments <- AB
  split_vars <- colnames(data$x)[unique(sapply(out$tree, function(x) ifelse(x[12]!=-1, x[2], NA))+1)]
  split_vars <- split_vars[!is.na(split_vars)]
  out$cat_levels <- list()

  for (s in split_vars) {
    if (var.types[1,s] >= 1)  {
      out$cat_levels[[s]] <- levels(data$x[,s])
    }
  }
  out$formula <- formula
  out$var.types <- var.types
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


