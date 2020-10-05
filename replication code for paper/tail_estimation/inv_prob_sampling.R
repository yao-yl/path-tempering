library(rstan)
library(ggplot2)
options(mc.cores = parallel::detectCores())

#Compile models for 8 schools and smoothing the marginal estimate
smooth_mod <- stan_model("./fit_prior.stan")
main_mod <- stan_model("./8schools.stan")

#Function for updaing tempering configuration list with new lower and upper bounds for inverse temperature
update_conf <- function(lbound,ubound){
  divs <- 5
  lower <- lbound
  upper <- ubound
  K_logit <- round(divs*(upper-lower))
  K_gaussian <- K_logit
  if (K_logit > 0){
    mu_logit <- lower + (upper - lower)*seq(1,K_logit,1)/(K_logit)
    mu_gaussian <- mu_logit
    sigma_logit <- 2*rep((upper - lower)/K_logit, K_logit)
    sigma_gaussian <- sigma_logit
  } else {
    mu_logit <- numeric(length=0)
    mu_gaussian <- mu_logit
    sigma_logit <- numeric(length=0)
    sigma_gaussian <- sigma_logit
  }
  conf_temp <- list(K_logit = K_logit,
                    K_gaussian = K_gaussian,
                    mu_logit = mu_logit,
                    sigma_logit = sigma_logit,
                    mu_gaussian = mu_logit,
                    sigma_gaussian = sigma_logit)
  return(conf_temp)
}

#Compute gradients for path sampling for the eight schools model
calc_U <- function(params){
  tau <- params[1]
  mu <- params[2]
  theta <- params[3:10]
  return((1/tau^3)*sum((theta-mu)^2)-(8/tau))
}

#Simple window-based estimation of optimal prior from gradients
op_windowed <- function(us,num_samps = 10){
  N <- length(us)
  eus <- numeric(length = N)
  eus[1] <- sum(us[1:num_samps]^2)
  for (i in 2:N){
    upper <- min(i + num_samps - 1,N)
    eus[i] <- eus[i-1] - us[i-1]^2 + us[upper]^2
  }
  pr_vals <- sqrt((1/num_samps)*eus)
  return(pr_vals)
}

#Data for 8 schools model
y <- c(28,8,-3,7,-1,1,18,12)
sds <- c(15,10,16,11,9,11,10,18)
N <- length(y)

lbound <- 0.1
ubound <- 100

#Function for computing implicit divide-and-conquer marginal estimate until convergence
idc_sample <- function(lbound,ubound,iter = 15, samps_per_iter = 4000){
  conf_temp <- update_conf(lbound,ubound)
  b <- rep(0,1+conf_temp$K_gaussian+conf_temp$K_logit)
  b_marg_sum <- rep(0,1+conf_temp$K_gaussian+conf_temp$K_logit)
  b_marg_avg <- b_marg_sum
  smooth_errs <- numeric(length=iter)
  
  params <- matrix(nrow=0,ncol=10)
  i <- 1
  sampling <- TRUE
  max_path_diff <- Inf
  all_b <- b
  
  while(sampling & (i < iter)){
    print(paste("Running iteration",i))
    print("-- Sampling parameters")
    main_dat <- append(list(y=y,sds=sds,N=N,b=b,lbound=lbound,ubound=ubound),conf_temp) 
    main_fit <- sampling(main_mod,data=main_dat,control=list(adapt_delta=0.9),iter = samps_per_iter,refresh=0) 
    new_params <- Reduce(cbind,extract(main_fit,pars=c('tau','mu','theta')))
    params <- rbind(params,new_params)
    params <- params[order(params[,1]),]
    
    hist(params[,1],main=paste("Iteration",i),breaks=round((ubound-lbound)*5))
    
    U_vals <- apply(params,1,calc_U)
    S <- nrow(params)
    tau_diffs <- mapply(`-`,params[2:S,1],params[1:(S-1),1])
    U_means <- mapply(mean,U_vals[2:S],U_vals[1:(S-1)])
    
    path_marg <- cumsum(tau_diffs*U_means)
    #path_marg <- 100*(path_marg -  min(path_marg)+1)
    x_marg <- params[2:S,1]
    
    if(i > 1){
      min_old <- min(path_marg_old)
      log_rats <- path_marg - approx(x=x_marg_old,y=path_marg_old,xout = x_marg,yleft=min_old,yright=min_old,ties="ordered")$y
      k_hat <- psis(log_ratios=log_rats,r_eff=NA)$diagnostics$pareto_k
      print(k_hat)
      if(k_hat < 0.5){
        sampling <- FALSE
      }
    }
    
    print(paste("--- Worst case path sampling differece:",round(max_path_diff,3)))
    
    int_space <- 0.2
    x_int <- seq(lbound,ubound,by=int_space)
    path_int <- approx(x=x_marg,y=path_marg,xout=x_int,ties="ordered")
    min_path <- min(path_marg)
    path_int$y[is.na(path_int$y)] <- min_path
    
    print("-- Smoothing true marginal estimate")
    smooth_path_marg <- optimizing(smooth_mod, data=append(list(x=x_int,log_p=path_int$y,N=length(x_int)),conf_temp), as_vector=FALSE)
    
    prior_est_exp <- op_windowed(U_vals,num_samps = 5)
    prior_int <- approx(x=params[,1],y=log(prior_est_exp),xout=x_int,ties="ordered")
    prior_int$y[is.na(prior_int$y)] <- max(prior_int$y,na.rm=TRUE)
    
    print("-- Smooothing target marginal estimate")
    smooth_prior <- optimizing(smooth_mod, data=append(list(x=x_int,log_p=prior_int$y,N=length(x_int)),conf_temp), as_vector=FALSE)
    
    blin <- smooth_path_marg$par$b_linear
    blog <- smooth_path_marg$par$b_logit
    bgaus <- smooth_path_marg$par$b_gaussian
    b_marg <- c(blin,blog,bgaus)
    b_marg_avg_old <- b_marg_avg
    b_marg_sum <- b_marg_sum + c(blin,blog,bgaus)
    b_marg_avg <- (1/i)*b_marg_sum
    
    blin_pr <- smooth_prior$par$b_linear
    blog_pr <- smooth_prior$par$b_logit
    bgaus_pr <- smooth_prior$par$b_gaussian
    b_pr <- c(blin_pr,blog_pr,bgaus_pr)
    b <- b_marg - b_pr
    all_b <- rbind(all_b,b)
    
    smooth_errs[i] <- smooth_path_marg$par$sigma_err
    
    path_marg_old <- path_marg
    x_marg_old <- x_marg
    i <- i+1
  }

  return(list(path_marg = path_marg, x_marg = x_marg, x_diffs= tau_diffs, Us = U_means, all_b = all_b, b_marg = b_marg, conf_temp = conf_temp))
}

#Return function that evaluates to parametrically smoothed marginal estimate
smooth_marg_func <- function(b_marg,conf_temp){
  smooth_func <- function(x){
    smooth_est(x,b_marg,conf_temp)
  }
  return(smooth_func)
}

#Return path sampling estimate of cdf at each of the sampled points
cdf_est <- function(path_marg,x_diffs){
  S <- length(x_diffs)
  path_cum <- cumsum(x_diffs*exp(path_marg))
  path_cum <- path_cum/path_cum[S]
  return(path_cum)
}

#Estimate quantile for a given probability from cdf estimates
find_q <- function(x_marg,path_cum,p){
  disc <- Inf
  cum_func <- function(x){approxfun(x=x_marg,y=path_cum)(x)}
  gr <- seq(min(x_marg),max(x_marg),length.out=10)
  while(disc > min(p,1-p)/1000){
    ps <- sapply(gr,cum_func)
    p_ind <- max(which(ps < 1-p))
    disc <- abs(ps[p_ind] - (1-p))
    gr <- seq(gr[p_ind],gr[p_ind+1],length.out=10)
  }
  return(gr[p_ind])
}

#Compute parametric smoothed estimate of marginal at a given point
smooth_est <- function(x,b,conf_temp){
  b_lin <- b[1]
  b_log <- b[2:(1+conf_temp$K_logit)]
  b_gaus <- b[(conf_temp$K_logit+2):length(b)]
  logit_ks <- matrix(nrow=conf_temp$K_logit,ncol=length(x))
  gauss_ks <- matrix(nrow=conf_temp$K_gaussian,ncol=length(x))
  for (i in 1:conf_temp$K_logit){
    logit_ks[i,] <- 1 / (1 + exp(-(x - conf_temp$mu_logit[i]) / conf_temp$sigma_logit[i]))
  }
  for (i in 1:conf_temp$K_gaussian){
    gauss_ks[i,] <- exp(-0.5*(x - conf_temp$mu_gaussian[i])^2 / (conf_temp$sigma_gaussian[i]^2));
  }
  return(as.numeric(b_lin*x + b_log%*%logit_ks + b_gaus%*%gauss_ks))
} 

#Compute derivative of parametric estimate at a given point
smooth_est_deriv <- function(x,b,conf_temp){
  N <- length(x)
  deriv <- numeric(length=N)
  kernel_logit <- array(NA, c(N, temp_conf$K_logit))
  kernel_gaussian <- array(NA, c(N, temp_conf$K_gaussian))
  grad_logit <- array(NA, c(N, temp_conf$K_logit))
  grad_gaussian <- array(NA, c(N, temp_conf$K_gaussian))
  b_linear <- b[1]
  b_logit <- b[1 + (1:conf_temp$K_logit)]
  b_gaussian <- b[1 + conf_temp$K_logit + (1:conf_temp$K_gaussian)]
  for (i in 1:N){
    kernel_logit[i,] <- 1 / (1 + exp(-(x[i] - temp_conf$mu_logit) / temp_conf$sigma_logit));
    kernel_gaussian[i,] = exp(-.5*(x[i] - temp_conf$mu_gaussian) * (x[i] - temp_conf$mu_gaussian) / (temp_conf$sigma_gaussian * temp_conf$sigma_gaussian))
    grad_logit[i,] <-  kernel_logit[i,] * kernel_logit[i,] * exp(-(x[i] - temp_conf$mu_logit) / temp_conf$sigma_logit) / temp_conf$sigma_logit
    grad_gaussian[i,] = kernel_gaussian[i,] * (temp_conf$mu_gaussian - x[i]) / (temp_conf$sigma_gaussian * temp_conf$sigma_gaussian)
    deriv[i] = b_linear + sum(grad_logit[i,]*b_logit) + sum(grad_gaussian[i,]*b_gaussian)
  }
  return(deriv)
}

#Spline smooth the path sampling marginal estimate from output of IDC algorithm and extrapolate to the left so that the marginal estimate is defined at zero
spline_smooth <- function(idc,lbound){
  ssf <- smooth.spline(x=idc$x_marg,y=idc$path_marg)
  x_left <- seq(0,lbound,0.002)
  s_left <- rep(ssf$y[1],length(x_left))
  ssp <- vector(mode="list")
  ssp$y <- c(s_left,ssft$y)
  ssp$x <- c(x_left,ssf$x)
  return(ssp)
}
