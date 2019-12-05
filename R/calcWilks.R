#' @export
calcWilks <- function(E,H) {
  lambda = det(E)/det(E+H)
  return(lambda)
}
