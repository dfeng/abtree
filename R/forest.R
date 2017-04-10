#' @export
forest <- function(formula, data, min.bucket=10, min.split=30,
                   max.depth=4, n.tree=100) {
  n <- nrow(data)
  ltree <- list()
  for (i in 1:n.tree) {
    bagged.data <- data[sample(n, replace=TRUE),]
    tree <- abtree(formula, bagged.data, min.bucket=min.bucket, min.split=min.split,
                     max.depth=max.depth)
    ltree[[i]] <- tree
  }
  ltree
}
#' @export
pforest <- function(ltree, new.data) {
  n <- length(ltree)
  m <- nrow(new.data)
  preds <- matrix(,nrow=n,ncol=m)
  for (i in 1:n) {
    preds[i,] <- predict.abtree(ltree[[i]], new.data)
  }
  apply(preds, 2, function(x) names(table(x))[which.max(table(x))])
}