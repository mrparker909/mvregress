#' @title logLikMVR
#' @description Calculates the log likelihood of the classic multivariate regression model.
#' @param lm_mod a linear model fit using lm.
#' @return log likelihood
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
#' logLikMVR(lm(Y~X))
#' logLikMVR(lm(Y~1))
logLikMVR <- function(lm_mod) {
  E  <- lm_mod$residuals
  n  <- nrow(E)
  p  <- ncol(E)
  S  <- cov(E)
  Si <- solve(S)
  LL <- -(0.5)*( n*p * log(2*pi) - n * log(det(S)) - sum(diag(E %*% Si %*% t(E))))
  return(LL)
}
