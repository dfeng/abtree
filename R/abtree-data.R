#' @title National Supported Work Study Experimental Data
#' @name nsw
#' @description This dataset, analyzed in LaLonde (1986), comes from the 
#' National Supported Work Study Demonstration (NSW) that ran from 1975-1978. 
#' This was a temporary employment program that experimented with providing 
#' disadvantaged workers some work experience and other forms of assistance, with
#' the hope of improving their employability at the end of the program. Applicants
#' were randomly assigned to a treatment group that receives benefits or a
#' control group that received nothing. This dataset is adapted from 
#' \code{\link[FindIt]{LaLonde}}.
#' 
#' @docType data
#' @usage nsw
#' @format A data frame with 722 observations and 10 variables, described below:
#' \describe{
#'  \item{outcome}{whether earnings in 1978 are larger than in 1975}
#'  \item{treat}{whether the individual received the treatment}
#'  \item{age}{age in years}
#'  \item{educ}{education in years}
#'  \item{race}{subject's race}
#' \item{marr}{1 = married, 0 = not married}
#' \item{nodegr}{1 = no high school degree, 0 = has high school degree}
#' \item{log.re75}{log of earnings in 1975 }
#' \item{u75}{1 = unemployed in 1975, 0 = employed in 1975}
#' \item{wts.extrap}{extrapolation weights to the 1978 Panel Study for Income Dynamics dataset}
#' }
#' @references 
#' \enumerate{ 
#' \item LaLonde, R.J. (1986). Evaluating the econometric evaulations of training programs with experimental data. \emph{American Economic Review}, \bold{76}(4), 604-620.
#' \item Egami, N, Ratkovic, M., and Imai, K. (2015). FindIt: Finding Heterogeneous Treatment Effects. \emph{R package version 0.5.} \url{https://CRAN.R-project.org/package=FindIt} 
#' }
#' @source The \code{LaLonde} dataset in the 'FindIt' R package.
#' @keywords datasets
#' @examples 
#' # fitting the model on the entire 'nsw' dataset
#' tree <- abtree(outcome ~ treat | age + educ + race + marr + nodegr + log.re75 + u75,
#'                data = nsw)
#' plot(tree, xpd=TRUE)
#' rec_treats <- predict(tree, new.data = nsw) # recommended treatments for each person
#' 
#' 
#' # it is better to divide the dataset into a training set, a validation set 
#' # (for pruning), and a test set for assessing model on unseen data
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
#'
#' 
NULL