# =================  #
# ===  Helpers  ===  #
# =================  #

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
ParseFormula <- function(formula, data) {
  rhs <- formula[[3]]
  lhs <- formula[[2]]
  condition <- rhs[[2]]
  rhs <- rhs[[3]]
  out <- list()
  if (deparse(lhs) %in% colnames(data)) {
    y <- eval(lhs, data)
    # names(y) <- deparse(lhs)
    out[["y"]] <- y
  }
  if (deparse(condition) %in% colnames(data)) {
    trt <- eval(condition, data)
    # names(trt) <- deparse(condition)
    out[["trt"]] <- trt
  }
  rhs <- all.names(rhs)
  rhs <- rhs[rhs != "+"]
  x <- as.data.frame(data[,rhs], stringsAsFactors=FALSE)
  out[["x"]] <- x
  # return(list(y=y,x=x,trt=trt))
  return(out)
}

# convert treatment to 0,1 (A,B)
TreatToAB <- function(x, AB) { ifelse(x==AB[2], 1, 0) }

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

# input: a variable
# output: integer vector of length 2
#         type: integer (0: quantitative, 1: nominal, 2: ordinal)
#         ncat: number of nominal categories; 0 for quantitative/ordinal vars
GetVarType <- function(y) {
  type <- 0
  ncat <- 0
  if ("ordered" %in% class(y)) {
    type <- 2
  } else if("factor" %in% class(y)) {
    type <- 1
    ncat <- length(levels(y))
  }
  return(c(type=type, ncat=ncat))
}
