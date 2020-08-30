library(rstan)
options(mc.cores = parallel::detectCores())
main_model <- stan_model('./8schools_dc.stan')

# Recursive divide-and-conquer sampling algorithm
## Generates samples over a given range such that for each lambda, theta ~ p(theta|lambda) and
## such that there is a sufficient coverage of lambda samples across the entire range to
## facilitate path samplinga
dc_sample <- function(low,upp){
  
  # Fit 8 schools over specified range of tau
  fit <- sampling(main_model,data=list(y=y,sds=sds,N=length(y),lbound=low,ubound=upp),control=list(adapt_delta=0.9),iter=max(200*(upp-low),500))
  
  # Extract and order parameters
  params <- Reduce(cbind,extract(fit,pars=c('tau','mu','theta')))
  params <- params[order(params[,1]),]
  xs <- params[,1] #xs = tau in 8 schools
  
  # Divide range into sub-intervals over which samples appear sufficiently uniform
  # (This step helps to cut down on the number of recursive resamples needed for insufficiently covered sub-intervals)
  regions <- calc_regions(xs,low,upp)
  redos <- (regions$r_sizes < 200*regions$r_lens)
  
  # Resample sub-intervals with insufficient coverage
  x_samps_2 <- Reduce(rbind,lapply(which(redos),function(i){ dc_sample(regions$end_ps[i],regions$end_ps[i+1]) })) 
  # Collect samples from sub-intervals with sufficient coverage
  x_samps_1 <- Reduce(rbind,lapply(which(!redos),function(i) { return(params[regions$r_is[i,1]:regions$r_is[i,2],]) })) 
  
  #upper_lim <- regions$end_ps[length(regions$end_ps)]
  
  # Combine all samples together
  if (is.null(x_samps_2)){
    return(x_samps_1)
  } else {
    return(rbind(x_samps_2,x_samps_1))
  }
}

# Calculate sub-intervals of sampled range over which samples are sufficiently uniform
calc_regions <- function(xs,low,upp){
  
  run <- TRUE
  
  #Index of first sample in current sub-interval/region
  start <- 1
  
  #Upper and lower limits for frequencies determining cutoff for uniformity
  lim_l <- 50
  lim_u <- 200
  
  #Defining quantities for sub-intervals: end-points, lengths, number of elements, starting and ending indices
  end_ps <- c(low)
  r_lens <- c()
  r_sizes <- c()
  r_is <- matrix(ncol=2,nrow=0)
  
  while(run){
    
    #Sub-interval divided into bins of length len
    len <- xs[start+100]-xs[start]
    testing <- TRUE
    bins <- 0
    
    #Step through bins, counting the number of elements in each
    while(testing) {
      bins <- bins + 1
      
      #Number of elements in next bin 
      bin_num <- sum( (xs > xs[start] + bins*len) & (xs <= xs[start] + (bins+1)*len) )
      
      #If too different from other bins, declare this the cutoff for this sub-interval
      if ((bin_num  < lim_l) || (bin_num > lim_u)){
        
        testing <- FALSE
        
        end_p <- min(xs[start] + bins*len,upp)
        r_lens <- c(r_lens,end_p-end_ps[length(end_ps)])
        end_ps <- c(end_ps,end_p)
        r_sizes <- c(r_sizes,sum(xs[start:length(xs)] < end_p))
        end <- max(which(xs <= end_p))
        r_is <- rbind(r_is,c(start,end))
        start <- min(which(xs > end_p))
      }
    }
    
    # Once we have fewer than 100 samples left, declare all remaining samples the last sub-interval
    if (start+100 > length(xs)){
      
      if (start < Inf){
        
        r_lens <- c(r_lens,xs[length(xs)]-end_ps[length(end_ps)])
        r_is <- rbind(r_is,c(start,length(xs)))
        end_ps <- c(end_ps,xs[length(xs)])
        r_sizes <- c(r_sizes,length(xs)-start+1)
      }
      
      run <- FALSE
    }
    
  }
  
  return(list(end_ps=end_ps,r_sizes=r_sizes,r_is=r_is,r_lens=r_lens))
}

# Calculate U functions for path sampling
calc_U <- function(params){
  
  tau <- params[1]
  mu <- params[2]
  theta <- params[3:10]
  return((1/tau^3)*sum((theta-mu)^2)-(8/tau)-(2.2/tau)+(20/(tau^2)))
}

# Using a divide-and-conquer search, find the estimate of Q(p) where Q is the quantile function
find_p <- function(path_cum,x_marg,q){
  
  # Discrepancy between cdf of current estimate and p
  disc <- Inf
  
  # Compute estimate of cdf from path sampling estimates
  cum_func <- function(x){approxfun(x=x_marg,y=path_cum)(x)}
  
  # Grid of points on which to evaluate cdf estimate
  gr <- seq(min(x_marg),max(x_marg),length.out=10)
  
  # Find point for which F(q) is nearest to p
  # And repeat search on smaller interval containing q if needed
  while(disc > q/10){
    
    qs <- sapply(gr,cum_func)
    q_ind <- max(which(qs < 1-q))
    disc <- abs(qs[q_ind] - (1-q))
    gr <- seq(gr[q_ind],gr[q_ind+1],length.out=10)
  }
  
  return(gr[q_ind])
}

# Calculate path sampling estimate of cdf and find quantile estimate
est_p_from_params <- function(params,q){
  
  # Compute U functions
  U_vals <- apply(params,1,calc_U)
  S <- nrow(params)
  tau_diffs <- mapply(`-`,params[2:S,1],params[1:(S-1),1])
  U_means <- mapply(mean,U_vals[2:S],U_vals[1:(S-1)])
  
  # Compute PS marginal density and cdf estimates
  path_marg <- cumsum(tau_diffs*U_means)
  path_cum <- cumsum(tau_diffs*exp(path_marg))
  path_cum <- path_cum/path_cum[S-1]
  x_marg <- params[2:S,1]
  
  # Compute quantile estimate
  p_est <- find_p(path_cum,x_marg,q)
  
  return(list(p_est=p_est,path_cum=path_cum))
}

# Evaluate the estimated cdf at a point
eval_cdf_est <- function(path_cum,x_marg,q){
  
  cum_func <- function(x){approxfun(x=x_marg,y=path_cum)(x)}
  return(cum_func(q))
}

# Estimate a given quantile of tau
est_p <- function(q){
  
  # Initialize
  p_old <- 1.5
  disc <- Inf
  
  # Region over which to begin estimation
  low <- 1
  upp <- 20
  
  samps <- matrix(ncol=10,nrow=0)
  
  # Continue estimation over successive regions until estimated cdf value stabilizes
  while(disc > q/10){
    
    # Combine previous samples with new samples and re-sort
    samps <- rbind(dc_sample(low,upp),samps)
    samps <- samps[order(samps[,1]),]
    x_marg <- samps[2:nrow(samps),1]
    
    # Estimate cdf and quantile
    est_new <- est_p_from_params(samps,q)
    p_new <- est_new$p_est
    path_cum <- est_new$path_cum
    
    # Combine cdf at previous quantile estimate and compute difference with target probability
    new_q <- eval_cdf_est(path_cum,x_marg,p_old)
    disc <- abs((1-q) - new_q)
    
    # Prepare for next iteration
    p_old <- est_new$p_est
    low <- upp
    upp <- upp+20
  }
  
  return(p_new)
}
