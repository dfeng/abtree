#' @export
abforest <- function(formula, data, min.bucket=10, min.split=30,
                   max.depth=4, n.tree=100,
                   mtry=floor(sqrt(length(all.vars(formula[[3]][[3]]))))) {
  # if (is.null(mtry)) mtry <- all.vars(formula[[3]][[3]])
  # getting the number of columns from formula+data

  n <- nrow(data)
  out <- list()
  for (i in 1:n.tree) {
    bagged.data <- data[sample(n, replace=TRUE),]
    tree <- abtree(formula, bagged.data, min.bucket=min.bucket, min.split=min.split,
                     max.depth=max.depth, mtry)
    out[[i]] <- tree
  }
  class(out) <- "abforest"
  out$trt.levels <- out[[1]]$trt.levels
  out
}
#' @export
predict.abforest <- function(abforest, new.data, type="response") {
  n <- length(abforest)
  m <- nrow(new.data)
  preds <- matrix(,nrow=n,ncol=m)
  for (i in 1:n) {
    preds[i,] <- predict.abtree(abforest[[i]], new.data, type="response")
  }
  if (type=="response") {
    return(apply(preds, 2, function(x) names(table(x))[which.max(table(x))]))
  } else if (type=="prob") {
    return(t(apply(preds, 2, function(x) prop.table(table(factor(x, levels=abforest$trt.levels))))))
  } else if (type=="votes") {
    return(t(apply(preds, 2, function(x) table(factor(x, levels=abforest$trt.levels)))))
  }
}