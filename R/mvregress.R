
#' @title  mvregress
#' @description Performs multivariate normal regression for one to several response variables and one to several covariate variables.
#' @param X matrix of covariates with each column representing a variable, and each row representing an observation
#' @param Y matrix of response variables with each column representing a variable, and each row representing an observation
#' @param ... extra parameters passed to lm(Y~X, ...)
#' @examples 
#' y1 = seq(0.1,1,0.1)
#' y2 = y1*y1+1
#' 
#' Y = cbind(y1,y2)
#' 
#' x1 = log(1:10)
#' x2 = x1^2-x1
#' X = cbind(x1,x2)
#'  
#' m = mvregress(X,Y)
#' @export
mvregress <- function(X, Y, ...) {
  # TODO: replace lm with likelihood optimization algorithm
  mod = lm(formula = Y~X, ...)
  sm = summary(mod)
  k=ncol(X)
  N=nrow(mod$residuals)
  p=ncol(mod$residuals)
  LL = logLikMVR(mod)
  nmod = lm(Y~1)
  LRT = likelihoodRatioTest(logLikMVR(nmod), LL, k-1)
  pval=LRT$pvalue
  LRTstatistic = LRT$statistic
  
  results = list(model=mod, nullmodel=nmod, summary=sm, logLik=LL, LRTstatistic=LRTstatistic, LRpvalue=pval)
  return(results)
}
