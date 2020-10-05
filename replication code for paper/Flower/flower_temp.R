library(rstan)
library(MASS)

######################################
#Helper Functions & Configuration Data
######################################

fold_a <- function(a){
  return(1-abs(1-a))
}

smooth_est <- function(x,b){
  b_lin <- b[1]
  b_log <- b[2:(1+K_logit)]
  b_gaus <- b[(K_logit+2):length(b)]
  logit_ks <- matrix(nrow=K_logit,ncol=length(x))
  gauss_ks <- matrix(nrow=K_gaussian,ncol=length(x))
  for (i in 1:K_logit){
    logit_ks[i,] <- 1 / (1 + exp(-(x - mu_logit[i]) / sigma_logit[i]))
  }
  for (i in 1:K_gaussian){
    gauss_ks[i,] <- exp(-0.5*(x - mu_gaussian[i])^2 / (sigma_gaussian[i]^2));
  }
  return(as.numeric(exp(b_lin*x + b_log%*%logit_ks + b_gaus%*%gauss_ks)))
} 

lambda_list <- seq(0,1,0.1) # discrete (inverse) tempeture ladder
n_lambda <- length(lambda_list)

log_temp_dist <- function(lambda,theta,omega,A){
  if (is.null(dim(theta))){
    log_target <- -1*(sqrt(sum(theta^2))-10-A*cos(omega*atan2(theta[2],theta[1])))/2
  } else {
    log_target <- -1*(sqrt(rowSums(theta^2))-10-A*cos(omega*atan2(theta[,2],theta[,1])))/2 
  }
  log_psi <- sum(dnorm(theta,0,10,log=TRUE))
  return(lambda*log_target+(1-lambda)*log_psi)
}

update_probability=function(theta,r,z,omega,A)
{
  q_k <- rep(0, n_lambda)
  for( i in 1:n_lambda)
  {
    q_k[i] <- log_temp_dist(lambda_list[i],theta,omega,A)
    if (is.na(q_k[i])){
      q_k[i] <- 0.001
    }
  }
  update_prob=q_k+log(r)-log(z)
  update_prob=exp(update_prob)
  update_prob=update_prob/sum(update_prob)
  return(update_prob)
}

tempeture_sample=function(theta,r,z,omega,A)
{
  update_prob=update_probability(theta,r,z,omega,A)
  lambda_index=which(rmultinom(1, size=1, prob=update_prob)==1)
  return(lambda_index)
}

#Configuration parameters for simulated tempering algorithm
a_lower <- 0.1
a_upper <- 0.8
K_logit <- 20
K_gaussian <- K_logit
mu_logit <- a_lower + (a_upper - a_lower)*seq(1,2*K_logit-1,2)/(K_logit)
mu_gaussian <- mu_logit
sigma_logit <- 2*rep((a_upper - a_lower)/K_logit, K_logit)
sigma_gaussian <- sigma_logit
conf_temp <- list(a_lower = a_lower,
                  a_upper = a_upper,
                  K_logit = K_logit,
                  K_gaussian = K_gaussian,
                  mu_logit = mu_logit,
                  sigma_logit = sigma_logit,
                  mu_gaussian = mu_logit,
                  sigma_gaussian = sigma_logit)

#Tempered model and model for parametric prior approximation
temp_mod <- stan_model("flower.stan")
temp_mod_disc <- stan_model("flower_disc.stan")
smooth_mod <- stan_model("solve_tempering.stan")

#####################################
#Continuous Tempering & Path Sampling
#####################################

#Function for simulated tempering algorithm
cont_temper_path <- function(conf_temp, omega=6, A=6, scale_ratio = 4, temp_mod, smooth_mod, max_iter = 10, epsilon = 0.01, samp_iter = 1000, N_grid = 100, ad=0.8){

  b <- rep(0, 1 + K_logit + K_gaussian)
  temper_data <- list(a_lower=conf_temp$a_lower,
                      a_upper=conf_temp$a_upper,
                      b=b,
                      K_logit=conf_temp$K_logit,
                      mu_logit=conf_temp$mu_logit,
                      sigma_logit=conf_temp$sigma_logit,
                      K_gaussian=conf_temp$K_gaussian,
                      mu_gaussian=conf_temp$mu_gaussian,
                      sigma_gaussian=conf_temp$sigma_gaussian,
                      omega = omega,
                      A = A)
  
  fit_measure <- Inf
  iter <- 1
  draws_temp_cum <- list(a=c(),log_posterior=c(),log_psi=c())
  
  history <- vector(mode="list",length=max_iter)
  
  while(iter <= max_iter){
    print(paste("Running tempering iteration",iter))
    fit_tempered <- sampling(temp_mod,data=temper_data,iter=samp_iter,chains=3,control=list(adapt_delta=ad),refresh=0)
    fit_draws_temp <- extract(fit_tempered,pars=c("a","log_posterior","log_psi"))
    draws_temp_cum <- lapply(1:3,function(var_num) c(fit_draws_temp[[var_num]],draws_temp_cum[[var_num]]))
    
    gradients <- path_gradients(conf_temp, b, draws_temp_cum[[1]], draws_temp_cum[[2]], draws_temp_cum[[3]])
    path_posterior <- path_sampling(fold_a(draws_temp_cum[[1]]), gradients$u_post, conf_temp$a_lower, conf_temp$a_upper)
    #fit_measure <- max(abs(path_posterior$log_z))
    path_likelihood <- path_sampling(fold_a(draws_temp_cum[[1]]), gradients$u_lik, conf_temp$a_lower, conf_temp$a_upper)
    
    interp_grid <- approx(path_likelihood$a, path_likelihood$log_z, seq(0, 1, length=N_grid))
    smoothing_data <- list(N=N_grid,
                           a=interp_grid$x,
                           log_p=interp_grid$y,
                           K_logit=conf_temp$K_logit,
                           mu_logit=conf_temp$mu_logit,
                           sigma_logit=conf_temp$sigma_logit, 
                           K_gaussian=conf_temp$K_gaussian,
                           mu_gaussian=conf_temp$mu_gaussian,
                           sigma_gaussian=conf_temp$sigma_gaussian, 
                           a_lower=conf_temp$a_lower,
                           a_upper=conf_temp$a_upper)
    smooth_path_fit <- optimizing(smooth_mod, data=smoothing_data, as_vector=FALSE)
    b <- -smooth_path_fit$par$b
    temper_data$b <- b
    
    interp_grid_prior <- approx(path_posterior$a, path_posterior$log_z, seq(0, 1, length=N_grid))
    prior_data <- smoothing_data
    prior_data$a <- interp_grid_prior$x
    prior_data$log_p <- interp_grid_prior$y
    prior_coefs <- optimizing(smooth_mod, data=prior_data, as_vector=FALSE)$par$b
    history[[iter]] <- list(fit = fit_tempered, a_cum = fold_a(draws_temp_cum[[1]]), a = fold_a(fit_draws_temp$a), theta = extract(fit_tempered,pars="theta")$theta, prior_a = path_posterior, z_est = path_likelihood, coefs = prior_coefs)
    
    iter <- iter+1
  }
  
  # print("Running final fit.")
  # final_fit <- sampling(temp_mod,data=temper_data,iter=10*samp_iter,chains=4,refresh=0)
  # final_draws <- extract(final_fit,pars=c("a","theta"))
  # theta_final <- final_draws$theta[fold_a(final_draws$a) > conf_temp$a_upper,]
  # return(list(theta_final=theta_final,history=history, model=final_fit))
  return(list(history=history))
}

path_gradients <- function(temp_conf, b, a, log_posterior, log_psi){
  N <- length(a)
  a_reflected = 1 - abs(1 - a)
  a_scaled <- (a_reflected - temp_conf$a_lower)/(temp_conf$a_upper - temp_conf$a_lower)
  lambda <- ifelse(a_scaled <= 0, 0, ifelse(a_scaled < 1, 0.5 + 0.25*(3*(2*a_scaled - 1) - (2*a_scaled - 1)^3), 1))
  grad_lambda = ifelse(a_scaled <= 0, 0, ifelse(a_scaled < 1, (2*0.25*(3 - 3*(2*a_scaled - 1)^2)/(temp_conf$a_upper - temp_conf$a_lower)), 0))
  u_prior <- rep(NA, N)
  kernel_logit <- array(NA, c(N, temp_conf$K_logit))
  kernel_gaussian <- array(NA, c(N, temp_conf$K_gaussian))
  grad_logit <- array(NA, c(N, temp_conf$K_logit))
  grad_gaussian <- array(NA, c(N, temp_conf$K_gaussian))
  b_linear <- b[1]
  b_logit <- b[1 + (1:temp_conf$K_logit)]
  b_gaussian <- b[1 + temp_conf$K_logit + (1:temp_conf$K_gaussian)]
  for (i in 1:N){
    kernel_logit[i,] <- 1 / (1 + exp(-(lambda[i] - temp_conf$mu_logit) / temp_conf$sigma_logit));
    kernel_gaussian[i,] = exp(-.5*(lambda[i] - temp_conf$mu_gaussian) * (lambda[i] - temp_conf$mu_gaussian) / (temp_conf$sigma_gaussian * temp_conf$sigma_gaussian))
    grad_logit[i,] <-  kernel_logit[i,] * kernel_logit[i,] * exp(-(lambda[i] - temp_conf$mu_logit) / temp_conf$sigma_logit) / temp_conf$sigma_logit
    grad_gaussian[i,] = kernel_gaussian[i,] * (temp_conf$mu_gaussian - lambda[i]) / (temp_conf$sigma_gaussian * temp_conf$sigma_gaussian)
    u_prior[i] = grad_lambda[i] * (b_linear + sum(grad_logit[i,]*b_logit) + sum(grad_gaussian[i,]*b_gaussian))
  }
  u_lik <- rep(NA, N)
  for (i in 1:N){
    u_lik[i] <- grad_lambda[i] * (log_posterior[i] - log_psi[i])
  }
  u_post <- u_prior + u_lik
  return(list(u_prior=u_prior, u_lik=u_lik, u_post=u_post))
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

############################################
#Continuous Tempering & Empirical Estimation
############################################

cont_temper_imp <- function(conf_temp, omega = 6, A = 6, scale_ratio = 4, temp_mod, smooth_mod, max_iter = 10, epsilon = 0.01, samp_iter = 1000, N_grid = 100, ad=0.8){

  b <- rep(0, 1 + K_logit + K_gaussian)
  temper_data <- list(a_lower=conf_temp$a_lower,
                      a_upper=conf_temp$a_upper,
                      b=b,
                      K_logit=conf_temp$K_logit,
                      mu_logit=conf_temp$mu_logit,
                      sigma_logit=conf_temp$sigma_logit,
                      K_gaussian=conf_temp$K_gaussian,
                      mu_gaussian=conf_temp$mu_gaussian,
                      sigma_gaussian=conf_temp$sigma_gaussian,
                      omega = omega,
                      A = A)
  
  fit_measure <- Inf
  iter <- 1
  draws_temp_cum <- list(a=c(),log_posterior=c(),log_psi=c())
  
  history <- vector(mode="list",length=max_iter)
  
  while(iter <= max_iter){
    print(paste("Running tempering iteration",iter))
    fit_tempered <- sampling(temp_mod,data=temper_data,iter=samp_iter,chains=3,control=list(adapt_delta=ad),refresh=0)
    fit_draws_temp <- extract(fit_tempered,pars=c("a","log_posterior","log_psi"))
    draws_temp_cum <- lapply(1:3,function(var_num) c(fit_draws_temp[[var_num]],draws_temp_cum[[var_num]]))
    
    p_a_est <- density(fold_a(draws_temp_cum[[1]]))
    interp_grid <- approx(p_a_est$x, p_a_est$y, seq(0, 1, length=N_grid))
    log_p_a_grid <- log(interp_grid$y)
    incl <- which(!is.na(log_p_a_grid))
    
    smoothing_data <- list(N=length(incl),
                           a=interp_grid$x[incl],
                           log_p=log_p_a_grid[incl],
                           K_logit=conf_temp$K_logit,
                           mu_logit=conf_temp$mu_logit,
                           sigma_logit=conf_temp$sigma_logit, 
                           K_gaussian=conf_temp$K_gaussian,
                           mu_gaussian=conf_temp$mu_gaussian,
                           sigma_gaussian=conf_temp$sigma_gaussian, 
                           a_lower=conf_temp$a_lower,
                           a_upper=conf_temp$a_upper)
    smooth_path_fit <- optimizing(smooth_mod, data=smoothing_data, as_vector=FALSE)
    prior_coefs <- smooth_path_fit$par$b
    b_old <- b
    b <- (prior_coefs-b)
    temper_data$b <- b
    
    sampled_a <- fold_a(draws_temp_cum[[1]])
    z_est <- list(a=sampled_a,log_z=numeric(length=length(sampled_a)))
    z_est$log_z <- log(smooth_est(sampled_a,-b_old)) + log(smooth_est(sampled_a,prior_coefs))
    
    history[[iter]] <- list(a_cum = fold_a(draws_temp_cum[[1]]), 
                            a = fold_a(fit_draws_temp$a), 
                            theta = extract(fit_tempered,pars="theta")$theta, 
                            z_est = z_est, 
                            coefs = prior_coefs)
    
    iter <- iter+1
  }
  
  print("Running final fit.")
  final_fit <- sampling(temp_mod,data=temper_data,iter=10*samp_iter,chains=4,refresh=0)
  final_draws <- extract(final_fit,pars=c("a","theta"))
  theta_final <- final_draws$theta[fold_a(final_draws$a) > conf_temp$a_upper,]
  return(list(theta_final=theta_final,history=history, model=final_fit))
}

##########################################
#Discrete Tempering & Empirical Estimation
##########################################

disc_temper_imp <- function(omega = 6, A = 6, scale_ratio = 4, temp_mod_disc, max_iter = 10, epsilon = 0.01, samp_iter = 150, ad=0.8){

  temper_data <- list(omega = omega,
                      A = A)
  
  Z=rep(1/n_lambda,n_lambda)
  r=rep(1/n_lambda,n_lambda)
  fit_measure <- Inf
  iter <- 1
  #draws_temp_cum <- list(a=c(),log_posterior=c(),log_psi=c())
  
  history <- vector(mode="list",length=max_iter)
  
  while(iter <= max_iter){
    print(paste("Running tempering iteration",iter))
    
    theta_sample <- matrix(nrow=samp_iter,ncol=2)
    lambda_sample <- rep(NA,samp_iter)
    lambda_index <- which(rmultinom(1, size=1, prob=rep(1/n_lambda,n_lambda))==1)
    temper_data$lambda <- lambda_list[lambda_index]
    c <- rep(0, n_lambda)
    
    for( i in 1:samp_iter){
      theta_sample_stan <- sampling(temp_mod_disc,data=temper_data,iter=100,chains=1,refresh = 0)
      theta_sample[i,] <- extract(theta_sample_stan)[["theta"]][50,]
      lambda_index <- tempeture_sample(theta_sample[i,], r, Z,omega,A)
      lambda_sample[i] <- lambda_index
      c[lambda_index] <- c[lambda_index]+  1/samp_iter
      temper_data$lambda <- lambda_list[lambda_index]
    }
    
    if (min(c)==0)
    {
      log_c=log(c+exp(-10))
      c=exp(log_c)
      c=c/sum(c)
    }
    
    Z <- Z*r[1]/r * c/c[1]
    z_est <- list(a = lambda_list, log_z = log(Z)-log(Z[1]))
    history[[iter]] <- list(z_est = z_est,
                            a = lambda_list[lambda_sample],
                            theta = theta_sample)
    iter <- iter+1
  }
  
  return(list(history=history))
}

##########################################
#Discrete Tempering & Rao-Blackwellization
##########################################

disc_temper_rb <- function(omega = 6, A = 6, scale_ratio = 4, temp_mod_disc, max_iter = 10, epsilon = 0.01, samp_iter = 150, ad=0.8){

    temper_data <- list(omega = omega,
                        A = A)
  
  Z=rep(1/n_lambda,n_lambda)
  r=rep(1/n_lambda,n_lambda)
  fit_measure <- Inf
  iter <- 1
  #draws_temp_cum <- list(a=c(),log_posterior=c(),log_psi=c())
  
  history <- vector(mode="list",length=max_iter)
  
  while(iter <= max_iter){
    print(paste("Running tempering iteration",iter))
    
    theta_sample <- matrix(nrow=samp_iter,ncol=2)
    lambda_sample <- rep(NA,samp_iter)
    lambda_index <- which(rmultinom(1, size=1, prob=rep(1/n_lambda,n_lambda))==1)
    temper_data$lambda <- lambda_list[lambda_index]
    c <- rep(0, n_lambda)
    
    for( i in 1:samp_iter){
      theta_sample_stan <- sampling(temp_mod_disc,data=temper_data,iter=100,chains=1,refresh = 0)
      theta_sample[i,] <- extract(theta_sample_stan)[["theta"]][50,]
      lambda_index <- tempeture_sample(theta_sample[i,], r, Z,omega,A)
      lambda_sample[i] <- lambda_index
      c <- c + update_probability(theta_sample[i,],r,Z,omega,A)/samp_iter
      temper_data$lambda <- lambda_list[lambda_index]
    }
    
    if (min(c)==0)
    {
      log_c <- log(c+exp(-10))
      c <- exp(log_c)
      c <- c/sum(c)
    }
    
    Z <- Z*r[1]/r * c/c[1]
    z_est <- list(a = lambda_list, log_z = log(Z)-log(Z[1]))
    #r <- (Z)/sum(Z)
    history[[iter]] <- list(z_est = z_est,
                            a = lambda_list[lambda_sample],
                            theta = theta_sample)
    iter <- iter+1
  }
  
  return(list(history=history))
}
