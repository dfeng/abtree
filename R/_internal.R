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
  list(y=y, x=x, trt=trt, y.name=response, trt.name=treat)
}