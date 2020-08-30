
## not related to other code
# an example why we need to rescael the gradient in log z est.

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



#run LDA model and obtain some results:
#save(a, log_lik_sum, file="LDA_data.RData")
load( file="LDA_data.RData")

N_sims <- length(a)
all_a <- c(all_a, a)
all_log_lik <- c(all_log_lik, log_lik_sum)
gradients <- path_gradients(a_lower, a_upper, b, K_logit, mu_logit, sigma_logit, K_gaussian, mu_gaussian, sigma_gaussian, all_a, all_log_lik)
all_a_reflected = 1 - abs(1 - all_a)
path_post <- path_sampling(all_a_reflected, gradients$u_post, a_lower, a_upper)
fit_measure <- max(abs(path_post$log_z))
path_lik <- path_sampling(all_a_reflected, gradients$u_lik, a_lower, a_upper)
approx_grid <- approx(path_lik$a, path_lik$log_z, seq(0, 1, length=N_grid))
a_grid <- approx_grid$x
log_p_grid <- approx_grid$y
# So far it is fine


path_fit <- optimizing(update_model, data=list(N=N_grid, a=a_grid, log_p=log_p_grid, K_logit=K_logit, mu_logit=mu_logit, sigma_logit=sigma_logit, K_gaussian=K_gaussian, mu_gaussian=mu_gaussian, sigma_gaussian=sigma_gaussian, a_lower=a_lower, a_upper=a_upper), as_vector=FALSE, iter=1e5, verbose=T)
curve_fit <- path_fit$par$fit

## the optimization is supposed to fit log_p_grid, but no, it does not.

qplot(a_grid, log_p_grid, main = "the samples are fine")
qplot(a_grid, curve_fit,  main = "but the optimzer does not fit")
print(path_fit$value)  ## likely the value is too large. may need to 


## convert the input to unit scale 
path_fit <- optimizing(update_model, data=list(N=N_grid, a=a_grid, log_p=log_p_grid/1000, K_logit=K_logit, mu_logit=mu_logit, sigma_logit=sigma_logit, K_gaussian=K_gaussian, mu_gaussian=mu_gaussian, sigma_gaussian=sigma_gaussian, a_lower=a_lower, a_upper=a_upper), as_vector=FALSE, iter=1000, verbose=T)
curve_fit <- path_fit$par$fit
qplot(a_grid, log_p_grid/1000, main = "the samples are fine")
qplot(a_grid, curve_fit,  main = "it is fine")
path_fit$par$b
