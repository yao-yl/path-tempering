#' The link function
#'
#' transform temperature lambda to a.
#'
#' @export
#' @param  a transform temperature in [0,2].
#' @param  a_lower The lower bound of the quadrature. When a < a_lower, the inverse
#' temperature  is 0 and the sample is from the base.
#' @param a_upper The upper bound of the quadrature. When a > a_upper, the inverse
#' temperature  is 1 and the sample is from the target.
#' @return temperature lambda in [0,1].


g_lambda=function(a, a_lower=0.1, a_upper=0.8){
	a_reflected = 1 - abs(1 - a)
	a_scaled <- (a_reflected - a_lower)/(a_upper - a_lower)
	lambda <- ifelse(a_scaled <= 0, 0, ifelse(a_scaled < 1, 0.5 + 0.25*(3*(2*a_scaled - 1) - (2*a_scaled - 1)^3), 1))
	return(lambda)
}
