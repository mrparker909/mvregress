#' @export
calcWilks <- function(E,H) {
  lambda = det(E)/det(E+H)
  return(lambda)
}

#' @title wilksLambdaTest
#' @description Performs Wilks likelihood ratio test on the multivariate regression model lm_mod versus the null hypothesis model lm_nullmod. Note that if n_adj=0 (the default) then the return values pval and pvalAdj will be the same. Otherwise, pval will be the unadjusted pvalue (full degrees of freedom), and pvalAdj will be the adjusted pvalue (reduced degrees of freedom).
#' @param lm_mod fitted model (eg: lm(Y~X)) corresponding to the alternative hypothesis.
#' @param lm_nullmod fitted model (eg: lm(Y~1)) corresponding to the null hypothesis, default is lm(Y~1) when lm_nullmod==NULL.
#' @param n_adj default=0, subtracted from the error degrees of freedom (for example, may be used to correct for previous testing done on the data).
#' @examples 
#' 
#' lam = NULL
#' pval = NULL
#' for(i in 1:1000) {
#' x1 = rnorm(100)
#' x2 = rnorm(100)
#' y1 = rnorm(100)
#' y2 = rnorm(100)
#' A = cbind(y1,y2)
#' B = cbind(x1,x2)
#' wi = wilksLambdaTest(lm(A~B), n_adj=0)
#' lam = c(lam, wi$lambda)
#' pval = c(pval, wi$pval)
#' }
#' 
#' hist(lam, freq=F, breaks=100)
#' curve(dchisq(x, wi$chisq_df), 0, max(lam), lwd=2, xlab = "", ylab = "", add = T)
#' hist(pval, freq=F, breaks=100)
#' 
#' lam = NULL
#' pval = NULL
#' for(i in 1:1000) {
#' x1 = rnorm(100)
#' x2 = rnorm(100)
#' y1 = rnorm(100)
#' y2 = rnorm(100)
#' y3 = rnorm(100)
#' y4 = rnorm(100)
#' y5 = rnorm(100)
#' y6 = rnorm(100)
#' A = cbind(y1,y2,y3,y4,y5,y6)
#' B = cbind(x1,x2)
#' wi = wilksLambdaTest(lm(A~B), n_adj=0)
#' lam = c(lam, wi$lambda)
#' pval = c(pval, wi$pval)
#' }
#' 
#' hist(lam, freq=F, breaks=100)
#' curve(dchisq(x, wi$chisq_df), 0, max(lam), lwd=2, xlab = "", ylab = "", add = T)
#' hist(pval, freq=F, breaks=100)
#' @export
wilksLambdaTest <- function(lm_mod, lm_nullmod=NULL, n_adj=0) {
  
  Ymod=lm_mod$model[[names(lm_mod$model)[1]]]
  Xmod=lm_mod$model[[names(lm_mod$model)[2]]]
  if(is.null(dim(Xmod))) { Xmod = cbind(Xmod) }
  
  p=ncol(Ymod)
  
  if(is.null(lm_nullmod)) {
    lm_nullmod = lm(Ymod~1)
  }
  Ynullmod=lm_nullmod$model[[names(lm_nullmod$model)[1]]]
  Xnullmod=lm_nullmod$model[[names(lm_nullmod$model)[2]]]
  
  # dof  = ncol(Xmod) - ncol(Xnullmod)
  df1 = ncol(Xmod)+1
  df2 = ncol(Xnullmod) + 1
  if(is.null(Xnullmod)) {
    #dof = ncol(Xmod)
    df2 = 1
  }
  
  testAdj = likelihoodRatioTest(logLikMVR(lm_nullmod, n_adj), 
                             logLikMVR(lm_mod, n_adj), 
                             df=(df1-df2)*p) #df=dof*p)
  
  test = likelihoodRatioTest(logLikMVR(lm_nullmod, 0), 
                             logLikMVR(lm_mod, 0), 
                             df=(df1-df2)*p) #df=dof*p)
  return(list(pval=test$pvalue, pvalAdj=testAdj$pvalue, lambda=test$statistic, lambdaAdj=testAdj$statistic, chisq_df=(df1-df2)*p, n_adj=n_adj))
}
