#' abtree
#'
#' @param formula an expression in the form of y ~ trt | cov1 + cov2 + ..., 
#'        where 'y' is a numeric/binary response variable, 
#'        'trt' must be a factor for the treatments, and cov1, cov2, ... are
#'        predictors
#' @param data the data frame where the variables reside
#' @param min.bucket the minimum number of observations in each leaf
#' @param min.split the minimum number of observations to return to the user
#'
#' @return tree: a data frame containing the nodes of the tree
#'         treatments: a mapping of A and B to the treatments passed in
#' @export
#'
#' @examples
#' # fitting the model on the entire 'nsw' dataset
#' tree <- abtree(outcome ~ treat | age + educ + race + marr + nodegr + log.re75 + u75,
#'                data = nsw)
#' plot(tree, xpd=TRUE)
#' rec_treats <- predict(tree, new.data = nsw) # recommended treatments for each persons
#' 
#' # it is better to divide the dataset into a training set, a validation set 
#' # (for pruning), and a test set for assessing model on unseen data
#' 
#' 
#' set.seed(5)
#' s <- sample(1:nrow(nsw))
#' train <- nsw[s[1:500], ] # use 500 obs for training
#' valid <- nsw[s[501:650],] # use 150 obs for pruning the tree
#' test <- nsw[s[651:nrow(nsw)],] # use the remainder for model assessment
#' 
#' # fitting the model on the entire dataset
#' tree <- abtree(outcome ~ treat | age + educ + race + marr + nodegr + log.re75 + u75,
#'                data = train)
#' ptree <- prune(tree, valid) # prune step
#' rec_treats <- predict(tree, new.data = test) # recommended treatments 
#' 
#' 
#' # model assessment
#' 
#' # 1. divide test set into two groups: individuals who received our 
#' #    recommended treatments and individuals who did not receive
#' #    our recommended treatments
#' 
#' match_grp <- test$treat == rec_treats
#' 
#' # 2. we want the average outcome to be better for the "match" group; check:
#' mean(test$outcome[match_grp])
#' mean(test$outcome[!match_grp])
#' 
#' # 3. indeed, the individuals who received our recommended treatment had a 9% higher
#' #    chance of improving their income after the program than those who did not 

abtree <- function(formula, data, min.bucket=10, min.split=30,
                     max.depth=5, mtry=NULL) {
  m <- ParseFormula(formula, data)
  if(any(is.na(m$x))) {
    m$x <- na.omit(m$x)
    warning("Missing values found in predictors. Rows containing missing values have been discarded.")
  }
  x.types <- vapply(m$x, class, character(1))

  # ========================  Conditions  ========================  # 
  # Question: why don't we just convert them all to factors?
  if (any(x.types == "character"))
    stop("At least one predictor is of type 'character'. Please convert to factor.")
  if (class(m$trt) != "factor")
    stop("The treatment variable should be of type 'factor'.")
  if (class(m$y) == "character")
    stop("The response variable cannot be of type 'character'. Please convert to factor.")

  # Random Forest parameters
  if (is.null(mtry)) mtry <- ncol(m$x)

  x.levels <- sapply(m$x, levels)
  ncat     <- sapply(x.levels, length)
  y        <- as.numeric(m$y)
  x        <- data.matrix(m$x) # converts everything to numeric, which is what we want
                               # since the Rcpp code will take a NumericMatrix
  x[,x.types == "factor"] <- x[,x.types == "factor"] - 1L # correcting for 0-index
  trt        <- as.integer(m$trt) - 1L # correcting for 0-index
  ord        <- apply(x, 2, order) - 1L # correcting for 0-index
  trt.levels <- levels(m$trt)
  out <- rcpp_BuildTree(y, x, trt, ord,
                        ncat, length(trt.levels),
                        as.integer(min.bucket), # should these (as.integer)
                        as.integer(min.split),  # be here? or should they
                        as.integer(max.depth),  # be conditions?
                        as.integer(mtry))

  class(out)     <- "abtree"
  out$formula    <- formula
  out$x.types    <- x.types
  out$x.levels   <- x.levels
  out$trt.levels <- trt.levels
  out$ncat       <- ncat
  out$y.name     <- m$y.name
  out$trt.name   <- m$trt.name
  out$frame      <- FormatTree(out)
  # if (length(out$cp.table) == 0) {
  #   out$cp.table <- NULL
  # } else {
  #   out$cp.table <- matrix(unlist(out$cp.table), ncol=2, byrow=T)
  # }
  out
}