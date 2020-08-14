#' Adaptive path sampling
#'
#' Run path sampling with adaptations.
#'
#' @export
#' @param  sampling_model The stan model generated from \code{\link{code_temperature_augmented}}.
#' @param  data_list The list of data used in the original stan model.
#' @param  N_loop The max adaptations. The default is 10.
#' @param  max_sample_store The max sample to store from previous adaptaions. The default is
#' 4000.
#' @param  iter The number of iterations for the augmented stan program.
#' @param  iter_final The number of iterations during the final adaptation.
#' @param   max_treedepth The max tree depth  for the augmented stan program. A smaller
#' max_treedepth may increase exploration. efficiency in early adaptations.
#' @param thin The number of thining for the augmented stan program.
#' @param  a_lower The lower bound of the quadrature. When a < a_lower, the inverse
#' temperature  is 0 and the sample is from the base.
#' @param a_upper The upper bound of the quadrature. When a > a_upper, the inverse
#' temperature  is 1 and the sample is from the target.
#' @param K_logit An integer, the length of the logit kernels. The default is 20.
#' @param K_gaussian An integer, the length of the Gaussian kernels. The default is 20.
#' @param N_grid The number of internal interpolations in parametric estimation
#' The default is 100.
#' @param visualize_progress whether to visualize the progress.
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
#' library(rstan)
#' rstan_options(auto_write = TRUE)
#' sampling_model=stan_model(file_new)
#' path_sample_fit=path_sample(sampling_model=sampling_model,
#'                             iter_final=6000, iter=2000,
#'                             data=list(gap=10), N_loop = 6,
#'                              visualize_progress = TRUE)
#' sim_cauchy=extract(path_sample_fit$fit_main)
#' in_target= sim_cauchy$lambda==1
#' in_prior = sim_cauchy$lambda==0
#' # sample from the target
#' hist(sim_cauchy$theta[in_target])
#' # sample from the base
#' hist(sim_cauchy$theta[in_prior])
# the joint "path"
#' plot(sim_cauchy$a, sim_cauchy$theta)
# the normalization constant
#' plot(g_lambda(path_sample_fit$path_post_a), path_sample_fit$path_post_z)
#'}
#'
#'
#' @return a \code{\link{path_fit}} object.
#'
#' @seealso \code{\link{path_quadrature}}, \code{\link{path_gradients}}, \code{\link{path_fit}}.
#'
#'
#


path_sample=function(sampling_model,data_list, visualize_progress=FALSE,
										 N_loop =10, max_sample_store=4000,
										 a_lower =.1, a_upper = .8,  K_logit = 20, K_gaussian=20,  N_grid=100,
										 thin=2, iter=2000,  max_treedepth=10,  iter_final=4000){

	# initialization:
	b <- rep(0, 1 + K_logit + K_gaussian)
	mu_logit <-  a_lower + (a_upper - a_lower)*seq(1,2*K_logit-1,2)/(K_logit)
	sigma_logit <- 2*rep((a_upper - a_lower)/K_logit, K_logit)
	mu_gaussian <- mu_logit
	sigma_gaussian <- sigma_logit
	all_a=NULL
	all_log_lik=NULL

	if(visualize_progress==TRUE)
		par(mfrow=c(N_loop,3), mar=c(2,2,2,1), mgp=c(1.5,.5,0), tck=-.01)
	for (i in 1:N_loop){
		fit <- path_fit(sampling_model=sampling_model, data_list=data_list,
										a_lower, a_upper, b, K_logit, mu_logit, sigma_logit, K_gaussian, mu_gaussian, sigma_gaussian,
										all_a, all_log_lik, N_grid=N_grid,iter=ifelse(i==N_loop, iter_final, iter),
										max_treedepth=max_treedepth, thin=thin, visualize_progress=visualize_progress)
		b <- -fit$b
		all_a <- fit$all_a
		all_log_lik <- fit$all_log_lik
		if(length(all_a)> max_sample_store )
		{
			discard=length(all_a)- max_sample_store
			all_a=all_a[-c(1:discard)]
			all_log_lik=all_log_lik[-c(1:discard)]
		}
		print( paste("------ adaptation=", i,"---------fit measure = ", fit$fit_measure,"---------"))
	}
	return(fit)
}




