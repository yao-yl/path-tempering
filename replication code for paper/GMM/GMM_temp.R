library(rstan)
library(cowplot)
library(MASS)

######################################
#Helper Functions & Configuration Data
######################################

fold_a <- function(a){
  return(1-abs(1-a))
}

gen_simplex <- function(n){
  S <- diag(n)
  v <- ((1+sqrt(n+1))/n)*rep(1,n)
  S <- rbind(S,v)
  S <- t(apply(S,1,function(x){
    w <-x- colMeans(S)
    return(w/sqrt(sum(w^2)))
  }))
  return(S)
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

assign_count <- function(assigns,n){
  counts <- sapply(1:n,function(k){
    return(sum(assigns == k))
  })
  return(counts)
}

lambda_list <- seq(0,1,0.1) # discrete (inverse) tempeture ladder
n_lambda <- length(lambda_list)

log_temp_dist <- function(lambda,theta,centers,scale){
  M <- nrow(centers)
  target <- 0
  for (i in 1:M){
    target <- target + (1/M)*prod(dnorm(theta,centers[i,],scale))
  }
  log_target <- log(target)
  log_psi <- sum(dnorm(theta,0,1.5,log=TRUE))
  return(lambda*log_target+(1-lambda)*log_psi)
}

update_probability=function(theta,r,z,centers,scale)
{
  q_k=rep(0, n_lambda)
  for( i in 1:n_lambda)
  {
    q_k[i]=log_temp_dist(lambda_list[i],theta,centers,scale)
    if (is.na(q_k[i])){
      q_k[i] <- 0.001
    }
  }
  update_prob=q_k+log(r)-log(z)
  update_prob=exp(update_prob)
  update_prob=update_prob/sum(update_prob)
  return(update_prob)
}

tempeture_sample=function(theta,r,z,centers,scale)
{
  update_prob=update_probability(theta,r,z,centers,scale)
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
temp_mod <- stan_model("GMM.stan")
temp_mod_disc <- stan_model("GMM_disc.stan")
smooth_mod <- stan_model("solve_tempering.stan")

#####################################
#Continuous Tempering & Path Sampling
#####################################

#Function for simulated tempering algorithm
cont_temper_path <- function(conf_temp, centers, scale_ratio = 4, temp_mod, smooth_mod, max_iter = 10, epsilon = 0.01, samp_iter = 1000, N_grid = 100, ad=0.8){
  sim_dist <- min(dist(centers))

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
                      M = nrow(centers),
                      D = ncol(centers),
                      mode_centers = centers,
                      mode_scale = sim_dist/scale_ratio)
  
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
    history[[iter]] <- list(fit = fit_tempered,a_cum = fold_a(draws_temp_cum[[1]]), a = fold_a(fit_draws_temp$a), theta = extract(fit_tempered,pars="theta")$theta, prior_a = path_posterior, z_est = path_likelihood, coefs = prior_coefs)
    
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

cont_temper_imp <- function(conf_temp, centers, scale_ratio = 4, temp_mod, smooth_mod, max_iter = 10, epsilon = 0.01, samp_iter = 1000, N_grid = 100, ad=0.8){
  sim_dist <- min(dist(centers))
  
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
                      M = nrow(centers),
                      D = ncol(centers),
                      mode_centers = centers,
                      mode_scale = sim_dist/scale_ratio)
  
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

disc_temper_imp <- function(centers, scale_ratio = 4, temp_mod_disc, max_iter = 30, epsilon = 0.01, samp_iter = 150, ad=0.8){
  sim_dist <- min(dist(centers))
  scale <- sim_dist/scale_ratio
  
  temper_data <- list(M = nrow(centers),
                      D = ncol(centers),
                      mode_centers = centers,
                      mode_scale = scale)
  
  Z=rep(1/n_lambda,n_lambda)
  r=rep(1/n_lambda,n_lambda)
  fit_measure <- Inf
  iter <- 1
  #draws_temp_cum <- list(a=c(),log_posterior=c(),log_psi=c())
  
  history <- vector(mode="list",length=max_iter)
  
  while(iter <= max_iter){
    print(paste("Running tempering iteration",iter))
    
    theta_sample <- matrix(nrow=samp_iter,ncol=ncol(centers))
    lambda_sample <- rep(NA,samp_iter)
    lambda_index <- which(rmultinom(1, size=1, prob=rep(1/n_lambda,n_lambda))==1)
    temper_data$lambda <- lambda_list[lambda_index]
    c <- rep(0, n_lambda)
    
    for( i in 1:samp_iter){
      theta_sample_stan <- sampling(temp_mod_disc,data=temper_data,iter=100,chains=1,refresh = 0)
      theta_sample[i,] <- extract(theta_sample_stan)[["theta"]][50,]
      lambda_index <- tempeture_sample(theta_sample[i,], r, Z,centers,scale)
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

disc_temper_rb <- function(centers, scale_ratio = 4, temp_mod_disc, max_iter = 30, epsilon = 0.01, samp_iter = 150, ad=0.8){
  sim_dist <- min(dist(centers))
  scale <- sim_dist/scale_ratio
  
  temper_data <- list(M = nrow(centers),
                      D = ncol(centers),
                      mode_centers = centers,
                      mode_scale = scale)
  
  Z=rep(1/n_lambda,n_lambda)
  r=rep(1/n_lambda,n_lambda)
  fit_measure <- Inf
  iter <- 1
  #draws_temp_cum <- list(a=c(),log_posterior=c(),log_psi=c())
  
  history <- vector(mode="list",length=max_iter)
  
  while(iter <= max_iter){
    print(paste("Running tempering iteration",iter))
    
    theta_sample <- matrix(nrow=samp_iter,ncol=ncol(centers))
    lambda_sample <- rep(NA,samp_iter)
    lambda_index <- which(rmultinom(1, size=1, prob=rep(1/n_lambda,n_lambda))==1)
    temper_data$lambda <- lambda_list[lambda_index]
    c <- rep(0, n_lambda)
    
    for( i in 1:samp_iter){
      theta_sample_stan <- sampling(temp_mod_disc,data=temper_data,iter=100,chains=1,refresh = 0)
      theta_sample[i,] <- extract(theta_sample_stan)[["theta"]][50,]
      lambda_index <- tempeture_sample(theta_sample[i,], r, Z,centers,scale)
      lambda_sample[i] <- lambda_index
      c <- c + update_probability(theta_sample[i,],r,Z,centers,scale)/samp_iter
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

###################
#Plotting functions
###################

gen_centers <- function(M,dims){
  return(lapply(dims,function(d) return (matrix(rnorm(M*d),ncol=d))))
}

gen_plots_cont_path <- function(centers_list){
  dims <- lapply(centers_list,function(centers) return(ncol(centers)))
  N <- length(dims)
  temp_runs <- as.list(vector(length=N))
  for (i in 1:N){
    centers <- centers_list[[i]]
    temp_runs[[i]] <- cont_temper_path(conf_temp,centers,6,temp_mod,smooth_mod,samp_iter = 4000,max_iter=30) 
    plot_prior_a(temp_runs[[i]]$history,label=paste("dim",dims[i],sep=''),ts=seq(1,30,3))
    plot_norm_const(temp_runs[[i]]$history,label=paste("dim",dims[i],sep=''),ts=seq(1,30,3))
    #plot_mode_props(temp_runs[[i]]$history,label=paste("dim",dims[i],sep=''),ts=seq(1,30,1),centers=centers)
    if (dims[i] == 2){
      plot_cross(temp_runs[[i]]$history,label="dim2")
      plot_thetas(temp_runs[[i]]$history,label="dim2")
    }
  }
  return(temp_runs)
}

gen_plots_cont_imp <- function(centers_list){
  dims <- lapply(centers_list,function(centers) return(ncol(centers)))
  N <- length(dims)
  temp_runs <- as.list(vector(length=N))
  for (i in 1:N){
    centers <- centers_list[[i]]
    temp_runs[[i]] <- cont_temper_imp(conf_temp,centers,6,temp_mod,smooth_mod,samp_iter = 4000, max_iter=30) 
    plot_prior_a(temp_runs[[i]]$history,label=paste("dim",dims[i],sep=''),ts=seq(1,30,3))
    plot_norm_const(temp_runs[[i]]$history,label=paste("dim",dims[i],sep=''),ts=seq(1,30,3))
    #plot_mode_props(temp_runs[[i]]$history,label=paste("dim",dims[i],sep=''),centers=centers)
    if (dims[i] == 2){
      plot_cross(temp_runs[[i]]$history,label="dim2")
      plot_thetas(temp_runs[[i]]$history,label="dim2")
    }
  }
  return(temp_runs)
}

gen_plots_disc_imp <- function(centers_list){
  dims <- lapply(centers_list,function(centers) return(ncol(centers)))
  N <- length(dims)
  temp_runs <- as.list(vector(length=N))
  for (i in 1:N){
    centers <- centers_list[[i]]
    temp_runs[[i]] <- disc_temper_imp(centers,6,temp_mod_disc,samp_iter = 1000,max_iter=30) 
    plot_prior_a_disc(temp_runs[[i]]$history,label=paste("dim",dims[i],sep=''),ts=seq(1,30,3))
    plot_norm_const(temp_runs[[i]]$history,label=paste("dim",dims[i],sep=''),ts=seq(1,30,3))
    #plot_mode_props_disc(temp_runs[[i]]$history,label=paste("dim",dims[i],sep=''),centers=centers)
    if (dims[i] == 2){
      plot_cross(temp_runs[[i]]$history,label="dim2")
      plot_thetas(temp_runs[[i]]$history,label="dim2")
    }
  }
  return(temp_runs)
}

gen_plots_disc_rb <- function(centers_list){
  dims <- lapply(centers_list,function(centers) return(ncol(centers)))
  N <- length(dims)
  temp_runs <- as.list(vector(length=N))
  for (i in 1:N){
    centers <- centers_list[[i]]
    temp_runs[[i]] <- disc_temper_rb(centers,6,temp_mod_disc,samp_iter = 1000,max_iter=30) 
    plot_prior_a_disc(temp_runs[[i]]$history,label=paste("dim",dims[i],sep=''),ts=seq(1,30,3))
    plot_norm_const(temp_runs[[i]]$history,label=paste("dim",dims[i],sep=''),ts=seq(1,30,3))
    #plot_mode_props_disc(temp_runs[[i]]$history,label=paste("dim",dims[i],sep=''),centers=centers)
    if (dims[i] == 2){
      plot_cross(temp_runs[[i]]$history,label="dim2")
      plot_thetas(temp_runs[[i]]$history,label="dim2")
    }
  }
  return(temp_runs)
}

#Plotting functions
plot_prior_a <- function(history,ts = c(1,3,5,7,10),label){
  funcs <- lapply(ts,function(t){
    b <- history[[t]]$coefs
    func <- function(x) smooth_est(x,b)
    vfunc <- Vectorize(func)
    z_pr <- integrate(vfunc,lower=0,upper=1)$value
    return(function(x) vfunc(x)/z_pr)
  })
  
  plots <- as.list(vector(length=length(ts)))
  iter <- 1
  for (t in ts){
    plots[[iter]] <- ggplot(data=data.frame(x=history[[t]]$a_cum),aes(x)) +
      geom_histogram(aes(y=stat(20*count/sum(count))),breaks=seq(0,1,0.05)) + 
      stat_function(fun = funcs[[iter]],color='red') +
      xlab('a') + 
      ylab('')
    #ggsave(paste("plot_",label,"_",iter,".png"),plots[[iter]])
    iter <- iter + 1
  }
  grid_plot <- plot_grid(plotlist=plots,labels=ts)
  ggsave(paste("priors_",label,".png",sep=''),grid_plot,width=9,height=(1/3)*9)
}

plot_prior_a_disc <- function(history,ts = c(1,3,5,7,10),label){
  plots <- as.list(vector(length=length(ts)))
  iter <- 1
  for (t in ts){
    plots[[iter]] <- ggplot(data=data.frame(x=history[[t]]$a),aes(x)) +
      geom_histogram() + 
      xlab('a') + 
      ylab('')
    #ggsave(paste("plot_",label,"_",iter,".png"),plots[[iter]])
    iter <- iter + 1
  }
  grid_plot <- plot_grid(plotlist=plots,labels=ts)
  ggsave(paste("priors_",label,".png",sep=''),grid_plot,width=9,height=(1/3)*9)
}
  
plot_norm_const <- function(history,ts = c(1,3,5,7,10),label){
  plots <- as.list(vector(length=length(ts)))
  iter <- 1
  for (t in ts){
    plots[[iter]] <- ggplot(data=data.frame(x=history[[t]]$z_est$a,y=history[[t]]$z_est$log_z),aes(x,y)) + 
      geom_line() +
      geom_point() + 
      xlab('a') +
      ylab('log(z)')
    iter <- iter + 1
  }
  grid_plot <- plot_grid(plotlist=plots,labels=ts)
  ggsave(paste("zs_",label,".png",sep=''),grid_plot,width=9,height=(1/3)*9)
}

plot_thetas <- function(history,ts = c(1,3,5,7,10),label){
  plots <- as.list(vector(length=length(ts)))
  iter <- 1
  for (t in ts){
    as <- history[[t]]$a
    incl <- (as >= 0.8)
    plots[[iter]] <- ggplot(data=data.frame(x=history[[t]]$theta[incl,1],y=history[[t]]$theta[incl,2]),aes(x,y)) + 
      geom_point() + 
      xlab('theta_1') +
      ylab('theta_2') + 
      theme(aspect.ratio = 1)
    iter <- iter + 1
  }
  grid_plot <- plot_grid(plotlist=plots,labels=ts)
  ggsave(paste("thetas_",label,".png",sep=''),grid_plot,width=9,height=6)
}

plot_cross <- function(history,t=10,k=5,label){
  a_seq <- seq(0,1,length.out=k+1)
  a_means <- vector(length=k)
  iter <- 1
  plots <- as.list(vector(length=k))
  for (i in 1:k){
    a_means[i] <- (a_seq[i+1]+a_seq[i])/2
    ths <- history[[t]]$theta
    as <- history[[t]]$a
    incl <- ((as > a_seq[i]) & (as <= a_seq[i+1]))
    plots[[iter]] <- ggplot(data=data.frame(x=history[[t]]$theta[incl,1],y=history[[t]]$theta[incl,2]),aes(x,y)) + 
      geom_point() + 
      xlab('theta_1') +
      ylab('theta_2') + 
      theme(aspect.ratio = 1)
    iter <- iter + 1
  }
  grid_plot <- plot_grid(plotlist=plots,labels=paste("a~",a_means))
  ggsave(paste("thetas_cross_",label,".png",sep=''),grid_plot,width=9,height=6)
}

plot_mode_props <- function(history,ts = seq(1,10,1),label,centers){
  plots <- as.list(vector(length=length(ts)))
  iter <- 1
  props <- matrix(nrow=length(ts),ncol=nrow(centers))
  for (t in ts){
    as <- history[[t]]$a
    incl <- (as >= 0.8)
    thetas <- history[[t]]$theta[incl,]
    assigns <- apply(thetas,1,function(th){
      return(which.min(colSums((t(centers)-th)^2)))
    })
    new <- as.numeric(assign_count(assigns,nrow(centers)))/length(assigns)
    if (length(new) == 0){
      print('oops')
    }
    props[iter,] <- new
    iter <- iter + 1
  }
  plot <- ggplot() + 
    geom_hline(yintercept=(1/nrow(centers)),color='red') + 
    xlab('Iteration') + 
    ylab('Mode proportions') + 
    ylim(0,1)
  for (i in 1:(ncol(props))){
    plot <- plot + geom_line(data=data.frame(x=ts,y=props[,i]),aes(x,y))
  }
  ggsave(paste("props_",label,".png",sep=''),plot,width=9,height=6)
}

plot_mode_props_disc <- function(history,ts = seq(1,10,1),label,centers){
  plots <- as.list(vector(length=length(ts)))
  iter <- 1
  props <- matrix(nrow=length(ts),ncol=(ncol(centers)+1))
  for (t in ts){
    as <- history[[t]]$a
    incl <- (as == 1)
    if (sum(incl)>1){
      thetas <- history[[t]]$theta[incl,]
      assigns <- apply(thetas,1,function(th){
        return(which.min(colSums((t(centers)-th)^2)))
      })
      new <- as.numeric(assign_count(assigns,nrow(centers)))/length(assigns)
      if (length(new) == 0){
        print('oops')
      }
      props[iter,] <- new  
    } else {
      props[iter,] <- rep(0,nrow(centers))
    }
    iter <- iter + 1
  }
  plot <- ggplot() + 
    geom_hline(yintercept=(1/(ncol(centers)+1)),color='red') + 
    xlab('Iteration') + 
    ylab('Mode proportions') + 
    ylim(0,1)
  for (i in 1:(ncol(props))){
    plot <- plot + geom_line(data=data.frame(x=ts,y=props[,i]),aes(x,y))
  }
  ggsave(paste("props_",label,".png",sep=''),plot,width=9,height=6)
}
