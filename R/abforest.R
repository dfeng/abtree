#' @export
abforest <- function(formula, data, min.bucket=10, min.split=30,
                   max.depth=4, n.tree=100,
                   mtry=floor(sqrt(length(all.vars(formula[[3]][[3]]))))) {
  # if (is.null(mtry)) mtry <- all.vars(formula[[3]][[3]])
  # getting the number of columns from formula+data

  n <- nrow(data)
  out <- list()
  out$trees <- list()

  oob.estimates <- vector("list", n)
  for (i in 1:n.tree) {
    bagged.index <- sample(n, replace=TRUE)
    oobers <- (1:n)[-bagged.index]
    tree <- abtree(formula, data[bagged.index,], min.bucket=min.bucket, min.split=min.split,
                   max.depth=max.depth, mtry=mtry)
    oob.pred <- predict(tree, data[oobers,])
    for (i.oob in 1:length(oobers)) {
      oob <- oobers[i.oob]
      oob.estimates[[oob]] <- c(oob.estimates[[oob]], oob.pred[i.oob])
    }
    out$trees[[i]] <- tree
  }
  oob.pred <- lapply(oob.estimates, function(x) names(sort(table(x), decreasing=TRUE))[1])
  oob.match <- table(match=(oob.pred == data[,tree$trt.name]), outcome=data[,tree$y.name])
  # full.response <- lapply(1:length(oob.estimates), function(i) rep(data[i,tree$y.name], length(oob.estimates[[i]])))
  # full.trt <- lapply(1:length(oob.estimates), function(i) rep(data[i,tree$trt.name], length(oob.estimates[[i]])))
  # oob.fullmatch <- table(match=(unlist(oob.estimates) == unlist(full.trt)), outcome=unlist(full.response))

  class(out) <- "abforest"
  out$oob.estimates <- oob.estimates
  out$oob.match <- oob.match
  # out$oob.fullmatch <- oob.fullmatch
  out
}
#' @export
predict.abforest <- function(abforest, new.data, type="response") {
  n <- length(abforest$trees)
  m <- nrow(new.data)
  trt.levels <- abforest$trees[[1]]$trt.levels
  # uplift-type score / predicted probabilities
  if (type=="pred.response") {
    preds <- array(, dim=c(m,length(trt.levels),n))
    for (i in 1:n) {
      preds[,,i] <- predict(abforest$trees[[i]], new.data, type="prob")
    }
    resp <- apply(preds, 1:2, mean)
    colnames(resp) <- trt.levels
    return(resp)
  # returned values using only the response
  } else {
    preds <- matrix(, nrow=m, ncol=n)
    for (i in 1:n) {
      preds[,i] <- predict.abtree(abforest$trees[[i]], new.data, type="response")
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