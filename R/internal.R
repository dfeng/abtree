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
  # TODO: proper try/catch error handling in case covariate not in data.
  x <- as.data.frame(data[,cov.names], stringsAsFactors=FALSE)
  list(y=y, x=x, trt=trt, y.name=deparse(response), trt.name=deparse(treat))
}

# TODO: make the columns actually data.frames inside a data.frame,
#       like how rpart does
FormatTree <- function(obj) {
  ntrt <- length(obj$trt.levels)
  tree <- as.data.frame(matrix(unlist(obj$cpp.tree), ncol=(10+ntrt*2), byrow=T))
  colnames(tree) <- c("id",
                      "split_var", "split_value","optimal_trt",
                      "opt_Q", "complexity", "branch", "level",
                      "childleft_id", "childright_id",
                      paste(
                        rep(c("mu", "n"), each=ntrt),
                        rep(LETTERS[1:ntrt], times=2),
                        sep="_")
                      )
  # 
  # following rpart frame output
  # 
  rownames(tree) <- tree$id
  tree <- tree[,-1]

  tree[tree == -1] <- NA # replace -1 with NA
  tree$split_value[is.na(tree$split_var)] <- NA
  tree$optimal_trt <- LETTERS[tree$optimal_trt+1L]
  tree$var_type <- obj$x.types[tree$split_var+1L]
  tree$split_var <- names(obj$x.types)[tree$split_var+1L]
  tree$split_factor <- sapply(1:nrow(tree), function(i) {
    ifelse(tree$var_type[i] == "factor",
      obj$x.levels[[tree$split_var[i]]][tree$split_value[i]+1L],
      NA)
  })
  tree$split_var[is.na(tree$split_var)] <- "<leaf>"
  
  ## maybe add back later
  tree$complexity <- NULL
  tree$branch <- NULL
  tree
}