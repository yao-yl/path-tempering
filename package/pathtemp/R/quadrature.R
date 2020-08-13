#' Numerical quadrature using trapezoid rule
#'
#' Use numerical quadrature using trapezoid rule to compute the log normalization
#' constant from its potinwise gradient.
#'
#' @export
#' @param a a vector, all sampled transformed parameters.
#' @param u a vector, pointwise gradient at those points.
#' @param  a_lower The lower bound of the quadrature. When a < a_lower, the invse
#' temperture  is 0 and the sample is from the base.
#' @param a_upper The upper bound of the quadrature. When a > a_upper, the invse
#' temperture  is 1 and the sample is from the target.

#'
#' @details The function use numerical quadrature using trapezoid rule to compute
#' the normalization constant from its potinwise gradient \code{u}.
#'
#' @return The pointwise   log normalization constant at all a.
#'
#
path_quadrature <- function(a, u, a_lower=0.1, a_upper=0.8){
	if(length(a)!=length(u))
		stop("wrong dimension input.")
	keep <- (a > a_lower) & (a < a_upper)
	if (sum(keep) > 0){
		a_uniq <- sort(unique(a[keep]))
		N_uniq <- length(a_uniq)
		u_bar <- rep(NA, N_uniq)
		for (i in 1:N_uniq){
			ok <- a == a_uniq[i]
			u_bar[i] <- mean(u[ok])
		}
		width <- c(a_uniq, a_upper) - c(a_lower, a_uniq)
		u_bar_avg <- (c(u_bar, u_bar[N_uniq]) + c(u_bar[1], u_bar))/ 2
		log_z <- c(0, cumsum(width*u_bar_avg))
		return(list(a=c(0, a_lower, a_uniq, a_upper, 1), log_z=c(log_z[1], log_z, log_z[length(log_z)])))
	}
	else
		return(list(a=c(0, 1), log_z=c(0, 0)))
}
