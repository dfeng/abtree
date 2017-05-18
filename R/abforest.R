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
  out
}
#' @export
predict.abforest <- function(abforest, new.data, type="response") {
  n <- length(abforest)
  m <- nrow(new.data)
  trt.levels <- abforest[[1]]$trt.levels
  # uplift-type score / predicted probabilities
  if (type=="pred.response") {
    preds <- array(,dim=c(m,length(trt.levels),n))
    for (i in 1:n) {
      preds[,,i] <- predict.abtree(abforest[[i]], new.data, type="prob")
    }
    return(apply(preds, 1:2, mean))
  # returned values using only the response
  } else {
    preds <- matrix(,nrow=m,ncol=n)
    for (i in 1:n) {
      preds[,i] <- predict.abtree(abforest[[i]], new.data, type="response")
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