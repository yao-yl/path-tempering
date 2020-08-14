#' Fit one adaptation of path sampling
#'
#' Run path sampling with given parameters and pseudo priors.
#'
#' @export
#' @param  sampling_model The stan model generated from \code{\link{code_temperature_augmented}}.
#' @param  data_list The list of data used in the original stan model.
#' @param   iter The number of iterations for the augmented stan program.
#' @param   max_treedepth The max tree depth  for the augmented stan program. A smaller
#' max_treedepth may increase exploration. efficiency in early adaptations.
#' @param thin The number of thining for the augmented stan program.
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
#' @param all_a (optional) All sampled transformed tempeture from previous adaptations.
#' If not NULL, they will be remixed and used for normalization constant estimation.
#' @param all_log_lik (optional) All sampled pointwise log density difference between two
#' models from previous adaptations.  If not NULL, they will be remixed and used for
#'  normalization constant estimation.
#' @param N_grid The number of internal interpolations in parametric estimation
#' The default is 100.
#' @param visualize_progress whether to visualize the progress.
#' @param output_file The name of the output file that saves local variables, which
#' otherwise would be saved in stan. The default is \code{stan-output_path.txt}.
#'
#' @details The function fits one adaptations of path sampling with given
#' parameters and  pseudo priors that are specified by \code{b}. It runs the stan program,
#' run path sampling to compute the log normalization constants and log marginal,
#' and these two estimations with extra parametric regularization.
#'
#'
#'
#' @examples
#' \dontrun{
#' #' Make sure to pre-compile the stan optimization file
#' cat("Compiling Stan optimizer...")
#' update_model <- stan_model("solve_tempering.stan")
#' cat("done.\n")
#' }
#'
#'
#' @return Names list of fit_main: the stan fit object; b: the updated pseudo prior
#' that approximate the marginal density after this adaptations;  all_a: all sampled
#' inverse temperature including the current adaptations; all_log_lik: all saved
#' log likelihood, including the current adaptations; path_post_z: the estimation of the
#' log normalization constants, path_post_a: the location of a in path_post_z.
#'
#' @seealso \code{\link{path_gradients}}, \code{\link{path_quadrature}},\code{\link{path_sample}}
#'
#'
#

path_fit <- function(sampling_model, data_list=list(), a_lower, a_upper, b, K_logit,
										 mu_logit, sigma_logit, K_gaussian,
										 mu_gaussian, sigma_gaussian, all_a=NULL, all_log_lik=NULL,
										 N_grid=100, iter=2000, max_treedepth=8,thin=2,
										 visualize_progress=FALSE, output_file='stan-output_path.txt'){
	if(!exists("update_model"))
		stop("Make sure to pre-compile the stan optimization file. \n update_model <- stan_model(\"solve_tempering.stan\").")
	sink(output_file, append=FALSE)
	fit_main=sampling(sampling_model, data=c(list(a_lower=a_lower, a_upper=a_upper,
																								b=b, K_logit=K_logit, mu_logit=mu_logit,
																								sigma_logit=sigma_logit, K_gaussian=K_gaussian,
																								mu_gaussian=mu_gaussian,
																								sigma_gaussian=sigma_gaussian), data_list),
										chains=1, iter=iter, thin=thin, control=list(max_treedepth=max_treedepth)) # only support 1 chain
	sink()
	sim=extract(fit_main, permuted=FALSE)
	extract_lp_maxrix=extract_lp(output_file)
	a=sim[,,"a"]
	a=a[-length(a)]
	if( max (abs( extract_lp_maxrix$a - a))>0.01)
		stop("Mismatch, report a bug.")
	log_lik=extract_lp_maxrix$log_lik
	N_sims <- length(a)
	all_a <- c(all_a, a)
	all_log_lik <- c(all_log_lik, log_lik)
	gradients <- path_gradients(a_lower, a_upper, b, K_logit, mu_logit, sigma_logit, K_gaussian, mu_gaussian, sigma_gaussian, all_a, all_log_lik)
	all_a_reflected = 1 - abs(1 - all_a)
	path_post <- path_quadrature(all_a_reflected, gradients$u_post, a_lower, a_upper)
	fit_measure <- max(abs(path_post$log_z))
	path_lik <- path_quadrature(all_a_reflected, gradients$u_lik, a_lower, a_upper)
	approx_grid <- approx(path_lik$a, path_lik$log_z, seq(0, 1, length=N_grid))
	a_grid <- approx_grid$x
	log_p_grid <- approx_grid$y
	scael_max=max(max(abs(log_p_grid))/30, 1)
	path_fit <- optimizing(update_model, data=list(N=N_grid, a=a_grid, log_p=log_p_grid/scael_max, K_logit=K_logit, mu_logit=mu_logit, sigma_logit=sigma_logit, K_gaussian=K_gaussian, mu_gaussian=mu_gaussian, sigma_gaussian=sigma_gaussian, a_lower=a_lower, a_upper=a_upper), as_vector=FALSE, iter=6000)
	curve_fit <- path_fit$par$fit*scael_max
	b <- path_fit$par$b*scael_max
	a_reflected = 1 - abs(1 - a)
	in_target <- a_reflected > a_upper
	in_prior=a_reflected<a_lower
	in_prior_rate=mean(in_prior)
	if(in_prior_rate>0.8){
		in_prior_rate=in_prior_rate*0.95
		b[1]=in_prior_rate*	(log_mean_exp(log_lik[in_prior])) + (1-in_prior_rate)*b[1]
	}
	if(visualize_progress==TRUE){
		hist(a_reflected, xlim=c(0, 1), xaxs="i", breaks=seq(0, 1, .01), probability = TRUE,
				 main="Post. draws of a_reflected", xlab="a_reflected", cex.main=.9, cex.axis=.9, cex.lab=.9)
		abline(v=a_lower, col="gray")
		abline(v=a_upper, col="gray")
		plot(path_lik$a, path_lik$log_z, xlim=c(0, 1), xaxs="i", bty="l", pch=20, cex=.2, col="blue", main="Path est of log norml. const.", xlab="a", ylab="log z", cex.main=.9, cex.axis=.9, cex.lab=.9)
		lines(a_grid, curve_fit, lwd=.5, col="red")
		abline(v=a_lower, col="gray")
		abline(v=a_upper, col="gray")
		plot(path_post$a, path_post$log_z, xlim=c(0, 1), xaxs="i", bty="l", pch=20, cex=.2, col="blue", main="Path est of log marginal", xlab="a", ylab="log p", cex.main=.9, cex.axis=.9, cex.lab=.9)
		abline(v=a_lower, col="gray")
		abline(v=a_upper, col="gray")
	}
	return(list(fit_main=fit_main, b=b, all_a=all_a, all_log_lik=all_log_lik, fit_measure=fit_measure, path_post_a=path_lik$a, path_post_z=path_lik$log_z))
}

# internal ----------------------------------------------------------------
log_mean_exp=function(v)
{
	max(v)+ log(mean(exp(v- max(v) )))
}


