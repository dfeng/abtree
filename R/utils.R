# =================  #
# ===  Helpers  ===  #
# =================  #
ParseFormula <- function(formula, data) {
  rhs <- formula[[3]] # grp | hour + browser
  response <- formula[[2]] # y
  treat <- rhs[[2]] # grp
  covariates <- rhs[[3]] # hour + browser

  if (deparse(response) %in% colnames(data))
    y <- eval(response, data)
  else
    stop("Response variable not found in data.")
  if (deparse(treat) %in% colnames(data))
    trt <- eval(treat, data)
  else
    stop("Treatment variable not found in data.")

  cov.names <- all.vars(covariates) # c("hour", "browser")
  # TODO: proper try/catch error handling in case covariate note in data.
  x <- as.data.frame(data[,cov.names], stringsAsFactors=FALSE)
  list(y=y,x=x,trt=trt)
}

FormatTree <- function(object, child_id=FALSE, digits=3, scientific=FALSE) {
  tree <- rbind(matrix(unlist(object$tree), ncol=14, byrow=T))
  rownames(tree) <- tree[,1]
  is_leaf <- tree[,13] == -1
  colnames(tree)[1:14] <- c("id" ,"split_var", "split_value", "optimal_trt",
                      "n_A", "p_A", "n_B", "p_B",
                      "total_Q", "complexity", "branch", "pruned",
                      "childleft_id", "childright_id")
  if (child_id) {
    tree <- tree[,-c(1), drop=F]
  } else {
    tree <- tree[,-c(1, 13:14), drop=F]
  }
  tree <- as.data.frame(tree)
  tree$optimal_trt <- LETTERS[tree$optimal_trt+1]
  split_value <- rep(NA, nrow(tree))
  quantsplit <- rep(FALSE, nrow(tree))
  quantsplit[tree$split_var + 1 > 0] <- object$var.types[1,tree$split_var+1]==0
  split_value[quantsplit] <-
    format(tree$split_value[quantsplit], digits=digits, scientific=scientific)

  tree$var_type <- rep(NA, length(is_leaf))
  tree$var_type[!is_leaf] <- object$var.types[1,tree$split_var[!is_leaf]+1]
  tree$split_var[!is_leaf] <- colnames(object$var.types)[tree$split_var[!is_leaf]+1]
  tree$split_var[is_leaf] <- NA
  for (k in colnames(object$var.types)[object$var.types[1,]>=1] ) {
    var_splits <- which(tree$split_var == k)
    split_value[var_splits] <- object$cat_levels[[k]][tree$split_value[var_splits]+1]
  }
  split_value[is_leaf] <- NA
  tree$split_value <- split_value
  tree$p_A <- tree$p_A/tree$n_A
  tree$p_B <- tree$p_B/tree$n_B
  tree$split_var[is_leaf] <- "<leaf>"
  return(tree)

}

SplitData <- function(data, order, splitpos) {
  if (is.data.frame(data)) {
    data <- data[order,]
    l <- nrow(data)
    train <- data[1:splitpos,]
    test <- data[(splitpos+1):l,]
  } else {
    data <- data[order]
    l <- length(data)
    train <- data[1:splitpos]
    test <- data[(splitpos+1):l]
  }
  return(list(train=train, test=test))
}