#' @export
calcWilks <- function(E,H) {
  lambda = det(E)/det(E+H)
  return(lambda)
}

#' @title wilksLambdaTest
#' @description Performs Wilks likelihood ratio test on the multivariate regression model lm_mod versus the null hypothesis model lm_nullmod. Note that if n_adj=0 (the default) then the return values pval and pvalAdj will be the same. Otherwise, pval will be the unadjusted pvalue (full degrees of freedom), and pvalAdj will be the adjusted pvalue (reduced degrees of freedom). Data matrices should be named X for the covariates, and Y for the responses.
#' @param lm_mod fitted model (eg: lm(Y~X)) corresponding to the alternative hypothesis.
#' @param lm_nullmod fitted model (eg: lm(Y~1)) corresponding to the null hypothesis, default is lm(Y~1) when lm_nullmod==NULL.
#' @param n_adj default=0, subtracted from the error degrees of freedom (for example, may be used to correct for previous testing done on the data).
#' @examples 
#' x1 = rnorm(10)
#' x2 = rnorm(10)
#' y1 = rnorm(10)
#' y2 = rnorm(10)
#' Y = cbind(y1,y2)
#' X = cbind(x1,x2)
#' wilksLambdaTest(lm(Y~X), n_adj=3)
#' @export
wilksLambdaTest <- function(lm_mod, lm_nullmod=NULL, n_adj=0) {
  if(is.null(lm_nullmod)) {
    lm_nullmod = lm(lm_mod$model$Y~1)
  }
  
  dof  = ncol(lm_mod$model$X) - ncol(lm_nullmod$model$X)
  if(is.null(lm_nullmod$model$X)) {
    dof = ncol(lm_mod$model$X)
  }
  
  testAdj = likelihoodRatioTest(logLikMVR(lm_nullmod, n_adj), 
                             logLikMVR(lm_mod, n_adj), 
                             df=dof)
  
  test = likelihoodRatioTest(logLikMVR(lm_nullmod, 0), 
                             logLikMVR(lm_mod, 0), 
                             df=dof)
  return(list(pval=test$pvalue, pvalAdj=testAdj$pvalue, lambda=test$statistic, lambdaAdj=testAdj$statistic, chisq_df=dof, n_adj=n_adj))
}
