# The second experiment is running on a cluster with SLURM_ARRAY_TASK_ID in 1:30
args <-  Sys.getenv("SLURM_ARRAY_TASK_ID")
print(Sys.getenv("SLURM_ARRAY_JOB_ID"))
arrayid <- as.integer(args[1])
set.seed(as.integer(arrayid))

library(rstan)
rstan_options(auto_write = TRUE)
update_model <- stan_model("solve_tempering.stan")
comp=stan_model(file="HS.stan")
probit=stan_model(file="probit.stan")
logit=stan_model(file="logit.stan")

inv_logit=function (x)
{
	1/(1 + exp(-x))
}

path_sampling <- function(a, u, a_lower=0, a_upper=1){
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
log_mean_exp=function(v)
{
	max(v)+ log(mean(exp(v- max(v) )))
}


path_fit <- function(graph_filename, a_lower, a_upper, b, K_logit,
										 mu_logit, sigma_logit, K_gaussian,
										 mu_gaussian, sigma_gaussian, all_a, all_log_lik,
										 N_grid=100, iter=2000, chains=3,max_treedepth=10,thin=2,hs_data){
	lda_fit <-sampling(comp, data=c(list(a_lower=0, a_upper=1, b=b, K_logit=K_logit, mu_logit=mu_logit, sigma_logit=sigma_logit, K_gaussian=K_gaussian, mu_gaussian=mu_gaussian, sigma_gaussian=sigma_gaussian), hs_data), chains=4,	iter=iter,control = list(max_treedepth=max_treedepth))
	sim=extract(lda_fit)
	log_lik_sum=sim$log_lik_sum
	a <- sim$a
	beta=sim$beta
	
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
	scael_max=max(max(abs(log_p_grid))/30, 1)
	path_fit <- optimizing(update_model, data=list(N=N_grid, a=a_grid, log_p=log_p_grid/scael_max, K_logit=K_logit, mu_logit=mu_logit, sigma_logit=sigma_logit, K_gaussian=K_gaussian, mu_gaussian=mu_gaussian, sigma_gaussian=sigma_gaussian, a_lower=a_lower, a_upper=a_upper), as_vector=FALSE, iter=6000)
	curve_fit <- path_fit$par$fit*scael_max
	b <- path_fit$par$b*scael_max
	a_reflected = 1 - abs(1 - a)
	in_target <- a_reflected > a_upper
	log_test=sim$log_lik_test
	
	in_prior=a_reflected<a_lower
	in_prior_rate=mean(in_prior)
	if(in_prior_rate>0.8){
		in_prior_rate=in_prior_rate*0.95
		b[1]=in_prior_rate*	(log_mean_exp(log_lik_sum[in_prior])) + (1-in_prior_rate)*b[1]
	}
	beta_t= beta[in_target,]
	beta_p=  beta[in_prior,]
	time=sum(get_elapsed_time(lda_fit)[,2])  # count only sampling time, same for invididual fit
	return(list(lda_fit=lda_fit, b=b, all_a=all_a, all_log_lik=all_log_lik, fit_measure=fit_measure, path_post_a=path_lik$a, path_post_z=path_lik$log_z,  curve_fit=curve_fit, N_sims=N_sims, beta_t=beta_t, beta_p=beta_p, time=time))
}


compare_eff=function(hs_data, a_lower, a_upper){
	K_logit=K_gaussian=2
	b <- rep(0, 1 + K_logit + K_gaussian)
	mu_logit <-  a_lower + (a_upper - a_lower)*seq(1,2*K_logit-1,2)/(2*K_logit-1)
	sigma_logit <- 2*rep((a_upper - a_lower)/K_logit, K_logit)
	mu_gaussian <- mu_logit
	sigma_gaussian <- sigma_logit
	N_loop <- 1   
	N_grid=10 
	N_sims_v=c()
	time_v=c()
	max_treedepth=10
	curve_fit_list=b_list=c()
	set.seed(1)
	for (i in 1:N_loop){
		print(i)
		iter=2000
		thin=1
		all_a=NULL
		all_log_lik=NULL
		fit <- path_fit(paste("path_", i, ".pdf", sep=""), a_lower, a_upper, b, K_logit, mu_logit, sigma_logit, K_gaussian, mu_gaussian, sigma_gaussian, all_a, all_log_lik, N_grid=N_grid, chains=chains,iter=iter, max_treedepth=max_treedepth, thin=thin, hs_data=hs_data)
		b <- -fit$b
		all_a <- fit$all_a
		all_in_target =  fit$beta_t 
		all_in_prior =  fit$beta_p
		all_log_lik <- fit$all_log_lik
		curve_fit_list[[i]]=fit$curve_fit
		b_list[[i]]=fit$b
		N_sims_v[i]=fit$N_sims
		time_v[i]=fit$time
		#print( paste("====================iter=", i,"=========================", fit$fit_measure,"==========="))
	}
	
	probit_fit=sampling(probit, data=hs_data,	iter=2000, chains=4, control = list(max_treedepth=max_treedepth) )
	beta_probit=extract(probit_fit)$beta
	logit_fit=sampling(logit, data=hs_data,	iter=2000,chains=4, control = list(max_treedepth=max_treedepth) )
	beta_logit=extract(logit_fit)$beta
	time_sep=sum(get_elapsed_time(probit_fit)[,2])+sum(get_elapsed_time(logit_fit)[,2])
	return(c(min(apply(beta_probit, 2, ess_tail)/time_sep),
					 min(apply(all_in_target, 2, ess_tail)/sum(time_v)),
					 min(apply(beta_logit, 2, ess_tail)/time_sep),
					 min(apply(all_in_prior, 2, ess_tail)/sum(time_v))))
}








flag=0
d_v=40*c(2:8)
# prop=c(0.4,0.5, 0.6)
prop=c(0.1)
eff_mat2=array(NA, c(1,length(d_v), 4))
for(i.sim in 1:length(d_v)){
	for(j in 1:length(prop)){
		d <- d_v[i.sim]
		n <- 40
		n_test=0
		set.seed(as.integer(arrayid))
		x_full <- matrix(rnorm((n+n_test)*d, 0, 0.5), n+n_test,d)
		beta=rep(0,d)
		beta[1:3]=3:1/3   #HS--2
		#beta[1:10]=10:1/10
		temp_common=rnorm(n+n_test,0,0.4)
		temp_common=temp_common-mean(temp_common)
		for( i in 1: (d*prop[j]))
			x_full[, i]= x_full[, i]+temp_common
		f_full=x_full%*%beta
		f_full=  inv_logit (x_full%*%beta)
		y_full=rbinom(n= n+n_test, size=1, prob=f_full)
		x=x_full[1:n,]
		y=y_full[1:n]
		tau0 <- 1/(d-1) * 2/sqrt(n) # should be a reasonable scale for tau (see Piironen&Vehtari 2017, EJS paper)
		hs_data <- list(n=n, d=d, x=x, y=as.vector(y),
										scale_icept=1, scale_global=tau0,
										slab_scale=5, slab_df=4)
		eff_mat2[j,i.sim,]=compare_eff(hs_data=hs_data, a_lower=0.5, a_upper =0.7)
	}
	save(eff_mat2,file=paste(arrayid, ".RData", sep=""))
}



