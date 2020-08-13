#' Extract pointwise gradient for path sampling
#'
#' Compute the pointwise gradient for path sampling, including the gradient of
#' log normalization constant and the log marginal density.
#'
#'
#' @export
#' ' @param  a_lower The lower bound of the quadrature. When a < a_lower, the invse
#' temperature  is 0 and the sample is from the base.
#' @param a_upper The upper bound of the quadrature. When a > a_upper, the invse
#' temperature  is 1 and the sample is from the target.
#' @param b The current parameterization of the pseudo prior.
#' @param K_logit a vector, pointwise gradient at those points.


#'
#' @details The function use numerical quadrature using trapezoid rule to compute
#' the normalization constant from its potinwise gradient \code{u}.
#'
#' @return The pointwise  log normalization constant at all a.
#'
#


path_gradients <- function(a_lower, a_upper, b, K_logit, mu_logit, sigma_logit, K_gaussian, mu_gaussian, sigma_gaussian, a, log_lik){
	N <- length(a)
	a_reflected = 1 - abs(1 - a)
	a_scaled <- (a_reflected - a_lower)/(a_upper - a_lower)
	lambda <- ifelse(a_scaled <= 0, 0, ifelse(a_scaled < 1, 0.5 + 0.25*(3*(2*a_scaled - 1) - (2*a_scaled - 1)^3), 1))
	grad_lambda = ifelse(a_scaled <= 0, 0, ifelse(a_scaled < 1, (2*0.25*(3 - 3*(2*a_scaled - 1)^2)/(a_upper - a_lower)), 0))
	u_prior <- rep(NA, N)
	kernel_logit <- array(NA, c(N, K_logit))
	kernel_gaussian <- array(NA, c(N, K_gaussian))
	grad_logit <- array(NA, c(N, K_logit))
	grad_gaussian <- array(NA, c(N, K_gaussian))
	b_linear <- b[1]
	b_logit <- b[1 + (1:K_logit)]
	b_gaussian <- b[1 + K_logit + (1:K_gaussian)]
	for (i in 1:N){
		kernel_logit[i,] <- 1 / (1 + exp(-(lambda[i] - mu_logit) / sigma_logit));
		kernel_gaussian[i,] = exp(-.5*(lambda[i] - mu_gaussian) * (lambda[i] - mu_gaussian) / (sigma_gaussian * sigma_gaussian))
		grad_logit[i,] <-  kernel_logit[i,] * kernel_logit[i,] * exp(-(lambda[i] - mu_logit) / sigma_logit) / sigma_logit
		grad_gaussian[i,] = kernel_gaussian[i,] * (mu_gaussian - lambda[i]) / (sigma_gaussian * sigma_gaussian)
		u_prior[i] = grad_lambda[i] * (b_linear + sum(grad_logit[i,]*b_logit) + sum(grad_gaussian[i,]*b_gaussian))
	}
	u_lik <- rep(NA, N)
	for (i in 1:N){
		u_lik[i] <- grad_lambda[i] * (log_lik[i] )
	}
	u_post <- u_prior + u_lik
	return(list(u_prior=u_prior, u_lik=u_lik, u_post=u_post))
}
