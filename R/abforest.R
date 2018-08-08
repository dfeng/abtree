#' abforest
#'
#' @param formula an expression in the form of y ~ trt | cov1 + cov2 + ..., 
#'        where 'y' is a numeric/binary response variable, 
#'        'trt' must be a factor for the treatments, and cov1, cov2, ... are
#'        predictors
#' @param data the data frame where the variables reside
#' @param min.bucket the minimum number of observations in each leaf
#' @param min.split the minimum number of observations to return to the user
#' @param max.depth the maximum depth of each tree
#' @param n.tree the number of trees to grow
#' @param mtry the number of features to select at random to grow each tree
#'
#' @return out: a list of trees
#' 
#' @examples 
#' \dontrun{
#' single_B <- function(x, a=0, b=1) a+b*(x > 1)
#' single_A <- function(x) 0
#' f.single <- c(single_A, single_B)
#' n <- 2000
#' X <- seq(from=0, to=2, length.out=n)
#' funs <- f.single
#' A <- Vectorize(funs[[1]])(X)
#' B <- Vectorize(funs[[2]])(X)
#' A <- A + rnorm(n, sd=0.1)
#' B <- B + rnorm(n, sd=0.1)
#' 
#' data <- data.frame(
#'   y = c(A, B),
#'   x = rep(X, times=2),
#'   type = rep(c("A", "B"), each = n)
#' )
#' 
#' set.seed(5)
#' forest <- abforest(y ~ type | x, data=data, max.depth = 1000)
#' }
#' @export
abforest <- function(formula, data, min.bucket=10, min.split=30,
                   max.depth=4, n.tree=100,
                   mtry=floor(sqrt(length(all.vars(formula[[3]][[3]]))))) {
  # if (is.null(mtry)) mtry <- all.vars(formula[[3]][[3]])
  # getting the number of columns from formula+data

  n <- nrow(data)
  out <- list()
  out$trees <- list()

  # oob.estimates <- vector("list", n)
  for (i in 1:n.tree) {
    # print(i)
    bagged.index <- sample(n, replace=TRUE)
    # oobers <- (1:n)[-bagged.index]
    tree <- abtree(formula, data[bagged.index,], min.bucket=min.bucket, min.split=min.split,
                   max.depth=max.depth, mtry=mtry)
    # oob.pred <- predict(tree, data[oobers,])
    # for (i.oob in 1:length(oobers)) {
    #   oob <- oobers[i.oob]
    #   oob.estimates[[oob]] <- c(oob.estimates[[oob]], oob.pred[i.oob])
    # }
    out$trees[[i]] <- tree
  }
  # oob.pred <- lapply(oob.estimates, function(x) names(sort(table(x), decreasing=TRUE))[1])
  # oob.match <- table(match=(oob.pred == data[,tree$trt.name]), outcome=data[,tree$y.name])
  # full.response <- lapply(1:length(oob.estimates), function(i) rep(data[i,tree$y.name], length(oob.estimates[[i]])))
  # full.trt <- lapply(1:length(oob.estimates), function(i) rep(data[i,tree$trt.name], length(oob.estimates[[i]])))
  # oob.fullmatch <- table(match=(unlist(oob.estimates) == unlist(full.trt)), outcome=unlist(full.response))

  class(out) <- "abforest"
  # out$oob.estimates <- oob.estimates
  # out$oob.match <- oob.match
  # out$oob.fullmatch <- oob.fullmatch
  out
}
#' predict.abforest
#'
#' @param object an object of class 'abforest'
#' @param new.data a new data frame containing the variables used in MakeTree
#' @param type the type of data to return
#' @param ... optional arguments
#'
#' @return a vector of predicted optimal treatments for all of the observations
#' @export
predict.abforest <- function(object, new.data, type="response", ...) {
  n <- length(object$trees)
  m <- nrow(new.data)
  trt.levels <- object$trees[[1]]$trt.levels
  # uplift-type score / predicted probabilities
  if (type=="pred.response" | type=="raw.response") {
    preds <- array(, dim=c(m,length(trt.levels),n))
    for (i in 1:n) {
      preds[,,i] <- predict(object$trees[[i]], new.data, type="prob")
    }
    if (type=="raw.response") {
      return(preds)
    } else {
      resp <- apply(preds, 1:2, mean)
      colnames(resp) <- trt.levels
      return(resp)      
    }
  # returned values using only the response
  } else {
    preds <- matrix(, nrow=m, ncol=n)
    for (i in 1:n) {
      preds[,i] <- predict.abtree(object$trees[[i]], new.data, type="response")
    }
    if (type=="response") {
      return(apply(preds, 1, function(x) names(which.max(table(x)))))
    } else if (type=="prob") {
      return(t(apply(preds, 1, function(x) prop.table(table(factor(x, levels=trt.levels))))))
    } else if (type=="votes") {
      return(t(apply(preds, 1, function(x) table(factor(x, levels=trt.levels)))))
    }
  }
}