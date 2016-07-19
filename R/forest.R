#' @export
forest <- function(formula, data, min.bucket=10, min.split=30,
                   max.depth=5, n.tree=50) {
  n <- nrow(data)
  ltree <- list()
  for (i in 1:n.tree) {
    bagdata <- data[sample(n, replace=TRUE),]
    tree <- MakeTree(formula, bagdata, min.bucket=10, min.split=30,
                     max.depth=5)
    ltree[[i]] <- tree
  }
  ltree
}
#' @export
pforest <- function(ltree, newdata) {
  n <- length(ltree)
  m <- nrow(newdata)
  preds <- matrix(,nrow=n,ncol=m)
  for (i in 1:n) {
    preds[i,] <- predict.abtree(ltree[[i]], newdata)
  }
  apply(preds, 2, function(x) which.max(table(x))-1)
}