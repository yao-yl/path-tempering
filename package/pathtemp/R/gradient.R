#' Extract pointwise gradient for path sampling
#'
#' Compute the pointwise gradient for path sampling, including the gradient of
#' log normalization constant and the log marginal density.
#'
#'
#' @export
#' @param  a_lower The lower bound of the quadrature. When a < a_lower, the inverse
#' temperature  is 0 and the sample is from the base.
#' @param a_upper The upper bound of the quadrature. When a > a_upper, the inverse
#' temperature  is 1 and the sample is from the target.
#' @param b The current parameterization of the pseudo prior.
#' @param K_logit An integer, the length of the logit kernels.
#' @param K_gaussian An integer, the length of the Gaussian kernels.
#' @param sigma_logit The width of logit kernels.
#' @param mu_logit The location of logit kernels.
#' @param mu_gaussian The location of Gaussian kernels.
#' @param sigma_gaussian The width of Gaussian kernels.
#' @param a All sampled transformed parameters
#' @param log_lik the pointwise log density difference between two models
#'
#' @details The function use numerical quadrature using trapezoid rule to compute
#' the normalization constant from its potinwise gradient \code{u}. For the gradient of
#' the log normalization constant, the gradient is the sum of the log deritive of lambda,
#' and the pointwise log density difference between two models. For the gradient of
#' the log marginal density, the gradient is the sum of the previous two terms and the
#' gradient of log pseudo prior.
#'
#' For current implementations, all these gradients terms have closed form.
#'
#' @return Names list of u_prior, u_lik, and u_post: the gradient of log pseudo prior,
#' log normalization constant, and the log marginal density. These list can be used in
#' nest adaptation in \code{\link{path_quadrature}}.
#'
#' @seealso \code{\link{path_quadrature}}
#'
#'
#


path_gradients <- function(a_lower=0.1, a_upper=0.8, b, K_logit, mu_logit, sigma_logit,
													 K_gaussian, mu_gaussian, sigma_gaussian,
													 a, log_lik){
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
