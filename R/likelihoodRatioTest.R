#' @title likelihoodRatioTest
#' @description Calculates the likelihood ratio test statistic and associated p-value under the chi squared approximation with degrees of freedom df.
#' @param logLik_H0 log likelihood of the null model
#' @param logLik_H1 log likelihood of the alternative model
#' @param df degrees of freedom for the chisquared random variable (usually the difference in degrees of freedom between the alternative and the null models)
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
#' LL_H1 = logLikMVR(lm(Y~X)) # here degrees of freedom is 3
#' LL_H0 = logLikMVR(lm(Y~1)) # here degrees of freedom is 1
#' likelihoodRatioTest(LL_H0, LL_H1, df=3-1)
#' @export
likelihoodRatioTest <- function(logLik_H0, logLik_H1, df) {
  # test statistic is -2*log(Likelihood H0 / Likelihood H1) ~ chisq(df)
  statistic = -2 * (logLik_H0 - logLik_H1)
  return(list(statistic=statistic, pvalue = pchisq(q=statistic, df=df, lower.tail = F)))
}
