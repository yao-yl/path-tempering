###The true normalizing constant###
###################################
setwd("~/Desktop/consinuous_tempering/code")

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

a=seq(0,2,0.002)
a_reflected = 1 - abs(1 - a)
a_lower=0.1
a_upper=0.8
a_scaled <- (a_reflected - a_lower)/(a_upper - a_lower)
lambda <- ifelse(a_scaled <= 0, 0, ifelse(a_scaled < 1, 0.5 + 0.25*(3*(2*a_scaled - 1) - (2*a_scaled - 1)^3), 1))


g=function(a){
  a_reflected = 1 - abs(1 - a)
  a_scaled <- (a_reflected - a_lower)/(a_upper - a_lower)
  lambda <- ifelse(a_scaled <= 0, 0, ifelse(a_scaled < 1, 0.5 + 0.25*(3*(2*a_scaled - 1) - (2*a_scaled - 1)^3), 1))
  return(lambda)
}

z=function(l, n, k, a, b){
  return((beta(k+1, n-k+1)*(n+1) )^(-l) /beta(alpha, beta) * beta(l*k+a, l*(n-k)+b)) 
}

alpha=9
beta=0.75
k=115
n=550


alpha=2
beta=1
k=60
n=80


 
a.grid =seq(0,1,0.001)
lambda.grid.true=seq(0,1,0.001)
za.grid.true=a.grid
for( i in 1:length(lambda.grid.true))
  za.grid.true[i]=z(g(a.grid[i]),n,k,alpha, beta)
plot(a.grid, log(za.grid.true),type='l')


z.grid.true=lambda.grid.true
for( i in 1:length(lambda.grid.true))
  z.grid.true[i]=z(lambda.grid.true[i],n,k,alpha, beta)



cols <- c( "blue","red" )
cols2 <- sapply(cols, function(i) {
  c2r <- col2rgb(i) / 255
  c2r <- rgb(c2r[1], c2r[2], c2r[3], alpha=0.5)
})


theta.grid=seq(0,1,0.005)
base_density.grid=dbeta(theta.grid,alpha, beta)
base_density.grid[201]=1e5





density_joint=function(lambda,theta, alpha, beta, n,k){
return(dbinom(x=k,prob=theta,size=n)^lambda*dbeta(theta,alpha, beta)/z(lambda,n,k,alpha,beta))
}
lambda.grid=seq(0,1,0.005)


library(RColorBrewer)

col.grid <- colorRampPalette( brewer.pal(9,'Blues'))(64)
col.grid[1]='white' 
prob.grid=matrix(NA,length(lambda.grid),length(theta.grid))
for(j in 1:length(theta.grid))
  for(i in 1:length( lambda.grid ) )
{
    prob.grid[i,j]=density_joint(theta=theta.grid[j], lambda=lambda.grid[i],alpha, beta, n,k)
}

density_joint_a_theta=function(a,theta, alpha, beta, n,k){
  return(dbinom(x=k,prob=theta,size=n)^(g(a))*dbeta(theta,alpha, beta)/z(g(a),n,k,alpha,beta))
}
prob.grid_az=matrix(NA,length(a.grid),length(theta.grid))
for(j in 1:length(theta.grid))
  for(i in 1:length( a.grid ) )
  {
    prob.grid_az[i,j]=density_joint_a_theta(theta=theta.grid[j], a=a.grid[i],alpha, beta, n,k)
  }



####################################################
#Continuous Simulated Tempering with Path Sampling##   
###################################################
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



path_gradients <- function(a_lower, a_upper, b, K_logit, mu_logit, sigma_logit, K_gaussian, mu_gaussian, sigma_gaussian, a, log_bin){
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
    u_lik[i] <- grad_lambda[i] * (log_bin[i] )
  }
  u_post <- u_prior + u_lik
  return(list(u_prior=u_prior, u_lik=u_lik, u_post=u_post))
}


beta_bin_code='
data {
  real a_lower;
  real a_upper;
  int K_logit;
  vector[K_logit] mu_logit;
  vector[K_logit] sigma_logit;
  int K_gaussian;
  vector[K_gaussian] mu_gaussian;
  vector[K_gaussian] sigma_gaussian;
  vector[1 + K_logit + K_gaussian] b;
  int k;
  int n;
  real alpha;
  real beta;
}
transformed data {
  real b_linear;
  vector[K_logit] b_logit;
  vector[K_gaussian] b_gaussian;
  b_linear = b[1];
  b_logit = segment(b, 2, K_logit);
  b_gaussian = segment(b, 2 + K_logit, K_gaussian);
}
parameters {
  real<lower=0,upper=1> theta;
  real<lower=0,upper=2> a;
}
transformed parameters {
  real log_posterior;
  real log_psi;
  real a_reflected;    
  real a_scaled;
  real lambda;
  vector[K_logit] kernel_logit;
  vector[K_gaussian] kernel_gaussian;
  a_reflected = 1 - fabs(1 - a);
  a_scaled = (a_reflected - a_lower)/(a_upper - a_lower);
  if (a_scaled <= 0)
    lambda = 0;
  else if (a_scaled < 1)
    lambda = 0.5 + 0.25*(3*(2*a_scaled - 1) - (2*a_scaled - 1)^3);
  else
    lambda = 1;
  kernel_logit = 1 ./ (1 + exp(-(lambda - mu_logit) ./ sigma_logit));
  kernel_gaussian = exp(-.5*(lambda - mu_gaussian) .* (lambda - mu_gaussian) ./ (sigma_gaussian .* sigma_gaussian));
  log_posterior  =  binomial_lpmf(k|n,theta)  ;
  log_psi= beta_lpdf( theta|alpha,beta)   ;
}
model {
  target += lambda*b_linear;
  target += dot_product(kernel_logit, b_logit) + dot_product(kernel_gaussian, b_gaussian);
  target += lambda*log_posterior + log_psi;
}
'
compile_tempering <-  stan(model_code = beta_bin_code ,data=list(a_lower=0.1,a_upper=0.8,K_logit=2,mu_logit=c(0.2,0.6),sigma_logit=c(0.1,0.1),K_gaussian=2, mu_gaussian=c(0.2,0.6),sigma_gaussian=c(0.1,0.1),b=c(0,0,0,0,0), k=k,n=n,alpha=alpha, beta=beta),  iter=4,chains=1)

path_fit <- function( graph_filename, a_lower, a_upper, b, K_logit, mu_logit, sigma_logit, K_gaussian, mu_gaussian, sigma_gaussian, all_a, all_log_bin,N_grid, iter){
  fit <- stan(model_code = beta_bin_code, data=list(a_lower=a_lower, a_upper=a_upper, b=b, K_logit=K_logit, mu_logit=mu_logit, sigma_logit=sigma_logit, K_gaussian=K_gaussian, mu_gaussian=mu_gaussian, sigma_gaussian=sigma_gaussian, k=k,n=n,alpha=alpha, beta=beta), chains=3, iter=iter,fit=compile_tempering)
  print(fit, pars=c("a", "theta"))
  sim=extract(fit)
  a_chains <- extract(fit, permuted=FALSE)[,,"a"]
  theta_chains <- extract(fit, permuted=FALSE)[,,"theta"]
  log_posterior=sim$log_posterior
  theta=sim$theta
  a <- sim$a
  N_sims <- length(a)
  all_a <- c(all_a, a)
  all_log_bin <- c(all_log_bin, log_posterior)
  gradients <- path_gradients(a_lower, a_upper, b, K_logit, mu_logit, sigma_logit, K_gaussian, mu_gaussian, sigma_gaussian, all_a, all_log_bin)
  all_a_reflected = 1 - abs(1 - all_a)
  path_post <- path_sampling(all_a_reflected, gradients$u_post, a_lower, a_upper)
  fit_measure <- max(abs(path_post$log_z))
  path_lik <- path_sampling(all_a_reflected, gradients$u_lik, a_lower, a_upper)
  approx_grid <- approx(path_lik$a, path_lik$log_z, seq(0, 1, length=N_grid))
  a_grid <- approx_grid$x
  log_p_grid <- approx_grid$y
  path_fit <- optimizing(update_model, data=list(N=N_grid, a=a_grid, log_p=log_p_grid, K_logit=K_logit, mu_logit=mu_logit, sigma_logit=sigma_logit, K_gaussian=K_gaussian, mu_gaussian=mu_gaussian, sigma_gaussian=sigma_gaussian, a_lower=a_lower, a_upper=a_upper), as_vector=FALSE)
  curve_fit <- path_fit$par$fit
  b <- path_fit$par$b
  a_reflected = 1 - abs(1 - a)
  in_target <- a_reflected > a_upper
  if (sum(in_target) > 0  ) {
    beta_target <- theta[in_target]
  }else
    beta_target <- 0
  MeanlogPosDen=mean(dbeta(beta_target, k+alpha, n-k+beta, log = T ))
  
  # pdf(graph_filename, height=8, width=7)
  # par(mfrow=c(3,3), mar=c(3,3,2,1), mgp=c(1.5,.5,0), tck=-.01)
  # for (i in 1:3){
  #   plot(a_chains[,i], theta_chains[,i], xlim=c(0,2), ylim=c(-15,30), xaxs="i", bty="l", pch=20, cex=.2, main=paste("Chain", i), xlab="a", ylab="beta_1", cex.main=.9, cex.axis=.9, cex.lab=.9)
  # }
  # plot(a, theta, xlim=c(0,2), ylim=c(0,1), xaxs="i", bty="l", pch=20, cex=.2, main="Joint post draws from Stan", xlab="a", ylab="beta", cex.main=.9, cex.axis=.9, cex.lab=.9)
  # abline(v=a_lower, col="gray")
  # abline(v=a_upper, col="gray")
  # abline(v=2-a_lower, col="gray")
  # abline(v=2-a_upper, col="gray")
  # if (sum(in_target) > 0){
  #   plot(density(beta_target))
  # }
  # hist(a_reflected, xlim=c(0, 1), xaxs="i", breaks=seq(0, 1, .01), main="Post draws of a_reflected", xlab="a_reflected", cex.main=.9, cex.axis=.9, cex.lab=.9)
  # abline(v=a_lower, col="gray")
  # abline(v=a_upper, col="gray")
  # 
  # plot(path_lik$a, path_lik$log_z, xlim=c(0, 1), xaxs="i", bty="l", pch=20, cex=.2, col="blue", main="Path est of p(y|a_reflected)", xlab="a_reflected", ylab="log p(y|a_reflected)", cex.main=.9, cex.axis=.9, cex.lab=.9)
  # lines(a_grid, curve_fit, lwd=.5, col="red")
  # abline(v=a_lower, col="gray")
  # abline(v=a_upper, col="gray")
  # plot(path_post$a, path_post$log_z, xlim=c(0, 1), xaxs="i", bty="l", pch=20, cex=.2, col="blue", main="Path est of p(a_reflected)", xlab="a_reflected", ylab="log p(a_reflected", cex.main=.9, cex.axis=.9, cex.lab=.9)
  # abline(v=a_lower, col="gray")
  # abline(v=a_upper, col="gray")
  # dev.off()
  return(list(fit=fit, b=b, all_a_reflected=all_a_reflected, all_log_bin=all_log_bin, fit_measure=fit_measure, path_post_a=path_lik$a, path_post_z=path_lik$log_z , MeanlogPosDen=MeanlogPosDen,a_chains_1=a_chains[,1], theta_chains_1=theta_chains[,1] )  )
}





a_lower <- .1
a_upper <- .8
K_logit <- 20
mu_logit <-  a_lower + (a_upper - a_lower)*seq(1,2*K_logit-1,2)/(K_logit)
sigma_logit <- 2*rep((a_upper - a_lower)/K_logit, K_logit)
K_gaussian <- K_logit
mu_gaussian <- mu_logit
sigma_gaussian <- sigma_logit
b <- rep(0, 1 + K_logit + K_gaussian)
all_a <- NULL
all_log_bin <- NULL
update_model <- stan_model("~/Desktop/tempering/solve_tempering.stan")
iter=2000
N_loop <- 20
a_dynamic=list()
z_dynamic=list()
a_chains_1=matrix(NA,N_loop,iter/2)
theta_chains_1=matrix(NA,N_loop,iter/2)
MLPD=c()
set.seed(1)
for (i in 1:N_loop){
  print(i)
  fit <- path_fit(paste("path_", i, ".pdf", sep=""), a_lower, a_upper, b, K_logit, mu_logit, sigma_logit, K_gaussian, mu_gaussian, sigma_gaussian, all_a, all_log_bin, N_grid=100,  iter=iter)
  b <- -fit$b
  all_a <- fit$all_a
  all_log_bin <- fit$all_log_bin
  a_dynamic[[i]]=fit$path_post_a
  z_dynamic[[i]]=fit$path_post_z
  MLPD[i]=fit$MeanlogPosDen
  a_chains_1[i,]=fit$a_chains_1
  theta_chains_1[i,]=fit$theta_chains_1
  # if(fit$fit_measure < .4) {
  #   i <- i + 1
  #   fit <- path_fit(paste("path_bin_", i, ".pdf", sep=""), a_lower, a_upper, b, K_logit, mu_logit, sigma_logit, K_gaussian, mu_gaussian, sigma_gaussian, all_a, all_log_bin, N_grid=100, iter=1e4)
  #   break
  # }
}

L2_error=matrix(NA, 4,20) 
for( i in 1:10)
  L2_error[1,i]=  mean((     approx(g(a_dynamic[[i]]), z_dynamic[[i]], lambda.grid.true)$y-log(z.grid.true)     )^2    )
for( i in 11:20)
L2_error[1,i]=L2_error[1,10]

###############################################################
####################################################
#2. Continuous Simulated Tempering with emperical estination ##   
###################################################
beta_bin_code='
data {
real a_lower;
real a_upper;
int K_logit;
vector[K_logit] mu_logit;
vector[K_logit] sigma_logit;
int K_gaussian;
vector[K_gaussian] mu_gaussian;
vector[K_gaussian] sigma_gaussian;
vector[1 + K_logit + K_gaussian] b;
int k;
int n;
real alpha;
real beta;
}
transformed data {
real b_linear;
vector[K_logit] b_logit;
vector[K_gaussian] b_gaussian;
b_linear = b[1];
b_logit = segment(b, 2, K_logit);
b_gaussian = segment(b, 2 + K_logit, K_gaussian);
}
parameters {
real<lower=0,upper=1> theta;
real<lower=0,upper=2> a;
}
transformed parameters {
real log_posterior;
real log_psi;
real a_reflected;    
real a_scaled;
real lambda;
vector[K_logit] kernel_logit;
vector[K_gaussian] kernel_gaussian;
a_reflected = 1 - fabs(1 - a);
a_scaled = (a_reflected - a_lower)/(a_upper - a_lower);
if (a_scaled <= 0)
lambda = 0;
else if (a_scaled < 1)
lambda = 0.5 + 0.25*(3*(2*a_scaled - 1) - (2*a_scaled - 1)^3);
else
lambda = 1;
kernel_logit = 1 ./ (1 + exp(-(lambda - mu_logit) ./ sigma_logit));
kernel_gaussian = exp(-.5*(lambda - mu_gaussian) .* (lambda - mu_gaussian) ./ (sigma_gaussian .* sigma_gaussian));
log_posterior  =  binomial_lpmf(k|n,theta)  ;
log_psi= beta_lpdf( theta|alpha,beta)   ;
}
model {
target += lambda*b_linear;
target += dot_product(kernel_logit, b_logit) + dot_product(kernel_gaussian, b_gaussian);
target += lambda*log_posterior + log_psi;
}
'
compile_tempering <-  stan(model_code = beta_bin_code ,data=list(a_lower=0.1,a_upper=0.8,K_logit=2,mu_logit=c(0.2,0.6),sigma_logit=c(0.1,0.1),K_gaussian=2, mu_gaussian=c(0.2,0.6),sigma_gaussian=c(0.1,0.1),b=c(0,0,0,0,0), k=k,n=n,alpha=alpha, beta=beta),  iter=4,chains=1)

emperical_fit <- function( a_lower, a_upper, c, K_logit, mu_logit, sigma_logit, K_gaussian, mu_gaussian, sigma_gaussian ,N_grid, iter,n.iter){
  c_old=c
  fit <- stan("mixture_tempering.stan", data=list(a_lower=a_lower, a_upper=a_upper, b=-c_old, K_logit=K_logit, mu_logit=mu_logit, sigma_logit=sigma_logit, K_gaussian=K_gaussian, mu_gaussian=mu_gaussian, sigma_gaussian=sigma_gaussian, k=k,n=n,alpha=alpha, beta=beta), chains=2, iter=iter,fit=compile_tempering)
  print(fit, pars=c("a", "theta"))
  sim=extract(fit)
  a_chains <- extract(fit, permuted=FALSE)[,,"a"]
  theta_chains <- extract(fit, permuted=FALSE)[,,"theta"]
  #log_posterior=sim$log_posterior
  theta=sim$theta
  a <- sim$a
  a_reflected = 1 - abs(1 - a);
  N_sims <- length(a_reflected)
  p_a_density=density(a_reflected)   ## emperical estimation of the density
  approx_grid <- approx(p_a_density$x, p_a_density$y, seq(0, 1, length=N_grid))
  a_grid <- approx_grid$x
  log_p_grid <- log(approx_grid$y)
  log_p_grid[is.na(log_p_grid)]=-10
  fit_measure=max(abs(log_p_grid))
  index_ok=which(!is.na(log_p_grid))
  kernel_fit <- optimizing(update_model, data=list(N=length(index_ok), a=a_grid[index_ok], log_p=log_p_grid[index_ok], K_logit=K_logit, mu_logit=mu_logit, sigma_logit=sigma_logit, K_gaussian=K_gaussian, mu_gaussian=mu_gaussian, sigma_gaussian=sigma_gaussian, a_lower=a_lower, a_upper=a_upper), as_vector=FALSE)
  curve_fit <- kernel_fit$par$fit
  b <- kernel_fit$par$b
  lambda.grid=g(a_grid)
  logC.grid=rep(0,length(lambda.grid))
  predict_log_density=function(b_old){
  b_linear = b_old[1]
  b_logit = b_old[1: K_logit+1]
  b_gaussian = b_old[1: K_gaussian+K_logit+1] 
  for(i.lambda in 1:length(lambda.grid)){
  kernel_logit = 1 / (1 + exp(-(lambda.grid[i.lambda] - mu_logit) / sigma_logit))
  kernel_gaussian = exp(-.5*(lambda.grid[i.lambda] - mu_gaussian)^2  / (sigma_gaussian^2))
  logC.grid[i.lambda]=lambda.grid[i.lambda]*b_linear+t(kernel_logit)%*%(b_logit)+t(kernel_gaussian)%*%b_gaussian;                      
  }
  return(logC.grid)
  }
  logC.grid=predict_log_density(c_old)
  logP.grid=predict_log_density(b)
  logP_density_grid=log_p_grid
  log_z_est=logC.grid+logP.grid
  a_reflected = 1 - abs(1 - a)
  in_target <- a_reflected > a_upper
  if (sum(in_target) > 0  ) {
    beta_target <- theta[in_target]
  }else
    beta_target <- 0
  MeanlogPosDen=mean(dbeta(beta_target, k+alpha, n-k+beta, log = T ))
  return(list(fit=fit, b=b, fit_measure=fit_measure, a.grid=a_grid, lambda.grid=lambda.grid, logP.grid=logP.grid, logZ.grid=log_z_est,MeanlogPosDen=MeanlogPosDen, sim_theta=theta,sim_a = a,logP_density_grid=logP_density_grid,a_chains_1=a_chains[,1], theta_chains_1=theta_chains[,1]))
}

a_lower <- .1
a_upper <- .8

K_logit <- 20
mu_logit <-  a_lower + (a_upper - a_lower)*seq(1,2*K_logit-1,2)/(K_logit)
sigma_logit <- 2*rep((a_upper - a_lower)/K_logit, K_logit)
K_gaussian <- K_logit
mu_gaussian <- mu_logit
sigma_gaussian <- sigma_logit


update_model <- stan_model("solve_tempering.stan")

iter=4000
N_loop <- 20
 
a_chains_2=matrix(NA,N_loop,iter/2)
theta_chains_2=matrix(NA,N_loop,iter/2)

N_grid=100
MLPD=c()
logP.grid.iter=matrix(NA,N_grid,N_loop)
logZ.grid.iter=matrix(NA,N_grid,N_loop)
a.grid.iter=seq(0, 1, length=N_grid)
lambda.grid.iter=g(a.grid.iter)
c <- rep(0, 1 + K_logit + K_gaussian)
b <- rep(NA, 1 + K_logit + K_gaussian)

set.seed(1)
for (i in 1:N_loop){
  print(i)
  c_old=c
  fit <- emperical_fit( a_lower, a_upper, c, K_logit, mu_logit, sigma_logit, K_gaussian, mu_gaussian, sigma_gaussian, N_grid=N_grid, iter=iter,n.iter=i)
  b <- fit$b
  c=c_old+b
  MLPD[i]=fit$MeanlogPosDen
  logP.grid.iter[,i]=fit$logP.grid
  logZ.grid.iter[,i]=fit$logZ.grid
  a_chains_2[i,]=fit$a_chains_1
  theta_chains_2[i,]=fit$theta_chains_1
  
  #theta.sim[[i]]=fit$sim_theta
  #a.sim[[i]]=fit$sim_a
  #fit_measure=fit$fit_measure
  # if( !is.na(fit_measure)   &  fit_measure < 0.5) {
  #   i <- i + 1
  #   fit <- path_fit(paste("path_bin_", i, ".pdf", sep=""), a_lower, a_upper, b, K_logit, mu_logit, sigma_logit, K_gaussian, mu_gaussian, sigma_gaussian, all_a, all_log_bin, N_grid=100, iter=1e4)
  #   MLPD[i]=fit$MeanlogPosDen
  #   logP.grid.iter[,i]=fit$logP.grid
  #   logZ.grid.iter[,i]=fit$logZ.grid
  #   theta.sim[[i]]=fit$sim_theta
  #   a.sim[[i]]=fit$sim_a
  #   break
  # }
}

tempsave=logZ.grid.iter

for(i in 1:ncol(logZ.grid.iter))  # normalize such that the log z(0)=0
  logZ.grid.iter[,i]=logZ.grid.iter[,i]-logZ.grid.iter[1,i]

  
for( i in 1:20)
  L2_error[2,i]=  mean((     approx( lambda.grid.iter , logZ.grid.iter[,i], lambda.grid.true)$y-log(z.grid.true)     )^2    )




#############################################################################
###################################################
###  R-B #########################################
##################################################

lambda_list=seq(0,1,0.1) # discrete (inverse) tempeture ladder
n_lambda=length(lambda_list)
log_posterior=function(lambda,theta){
  log_link= dbinom(x=k,size=n,prob=theta,log=T)
  log_psi= dbeta(x=theta, alpha, beta,log=T)      
  return( lambda*log_link+log_psi)
}

stan_tempering_discrete='
data {
real<lower=0,upper=1> lambda;
int k;
int n;
real alpha;
real beta;
}
parameters {
real<lower=0,upper=1> theta;
}
transformed parameters {
real log_posterior;
real log_psi;
log_posterior  =  binomial_lpmf(k|n,theta)  ;
log_psi= beta_lpdf( theta|alpha,beta)   ;
}
model {
target += lambda*log_posterior + log_psi;
}
'
compile_RB=stan(model_code = stan_tempering_discrete,data=list(lambda=0,k=k,n=n,alpha=alpha, beta=beta ),iter=1,chains = 1)


update_probability=function(theta,r,z)
{
  q_k=rep(0, n_lambda)
  for( i in 1:n_lambda)
  {
    q_k[i]=log_posterior(lambda_list[i],theta)
  }
  update_prob=q_k+log(r)-log(z)
  update_prob=exp(update_prob)
  update_prob=update_prob/sum(update_prob)
  return(update_prob)
}

tempeture_sample=function(theta,r,z)
{
  update_prob=update_probability(theta,r,z)
  lambda_index=which(rmultinom(1, size=1, prob=update_prob)==1)
  return(lambda_index)
}
N=150


Z=rep(1/n_lambda,n_lambda)
r=rep(1/n_lambda,n_lambda)
N_RB_loop=20
Z.RB.iter=matrix(NA, length(lambda_list),N_RB_loop)


lambda_chain_3=matrix(NA, N_RB_loop, N)
theta_chain_3=matrix(NA, N_RB_loop, N)

set.seed(2)
for(iter_index in 1:N_RB_loop){
  theta_sample=rep(NA,N)
  lambda_sample=rep(NA,N)
  lambda_index=which(rmultinom(1, size=1, prob=rep(1/n_lambda,n_lambda))==1)
  c=rep(0, n_lambda)
  for( i in 1:N){
    theta_sample_stan=stan(model_code = stan_tempering_discrete,data=list(lambda=lambda_list[lambda_index],k=k,n=n,alpha=alpha, beta=beta ),iter=100,chains=1, fit=compile_RB,refresh=-1)
    theta_sample[i]=extract(theta_sample_stan)[["theta"]][50]
    lambda_index=tempeture_sample(theta_sample[i], r, Z)
    lambda_sample[i]=lambda_index
    c=c+update_probability(theta_sample[i],r,Z)/N
  }
  Z=Z*r[1]/r * c/c[1]
  Z.RB.iter[,iter_index]=log(Z)-log(Z[1])
  lambda_chain_3[iter_index,]=lambda_list[lambda_sample]
  theta_chain_3[iter_index,]=theta_sample
}

for( i in 1:20)
  L2_error[3,i]=  mean((     approx( lambda_list , Z.RB.iter[,i], lambda.grid.true)$y-log(z.grid.true)     )^2    )

####################################################
####### discrete + emperical #######################
####################################################
#######################################################
lambda_list=seq(0,1,0.1) # discrete (inverse) tempeture ladder
n_lambda=length(lambda_list)
log_posterior=function(lambda,theta){
  log_link= dbinom(x=k,size=n,prob=theta,log=T)
  log_psi= dbeta(x=theta, alpha, beta,log=T)      
  return( lambda*log_link+log_psi)
}

stan_tempering_discrete='
data {
real<lower=0,upper=1> lambda;
int k;
int n;
real alpha;
real beta;
}
parameters {
real<lower=0,upper=1> theta;
}
transformed parameters {
real log_posterior;
real log_psi;
log_posterior  =  binomial_lpmf(k|n,theta)  ;
log_psi= beta_lpdf( theta|alpha,beta)   ;
}
model {
target += lambda*log_posterior + log_psi;
}
'
compile_RB=stan(model_code = stan_tempering_discrete,data=list(lambda=0,k=k,n=n,alpha=alpha, beta=beta ),iter=1,chains = 1)


update_probability=function(theta,r,z)
{
  q_k=rep(0, n_lambda)
  for( i in 1:n_lambda)
  {
    q_k[i]=log_posterior(lambda_list[i],theta)
  }
  update_prob=q_k+log(r)-log(z)
  update_prob=exp(update_prob)
  update_prob=update_prob/sum(update_prob)
  return(update_prob)
}

tempeture_sample=function(theta,r,z)
{
  update_prob=update_probability(theta,r,z)
  lambda_index=which(rmultinom(1, size=1, prob=update_prob)==1)
  return(lambda_index)
}
N=150
Z=rep(1/n_lambda,n_lambda)
r=rep(1/n_lambda,n_lambda)

N_ST_loop=20
Z.ST.iter=matrix(NA, length(lambda_list), N_ST_loop)
lambda_chain_4=matrix(NA, N_ST_loop, N)
theta_chain_4=matrix(NA, N_ST_loop, N)

set.seed(1)
for(iter_index in 1:N_ST_loop){
  theta_sample=rep(NA,N)
  lambda_sample=rep(NA,N)
  lambda_index=which(rmultinom(1, size=1, prob=rep(1/n_lambda,n_lambda))==1)
  c=rep(0, n_lambda)
  for( i in 1:N){
    theta_sample_stan=stan(model_code = stan_tempering_discrete,data=list(lambda=lambda_list[lambda_index],k=k,n=n,alpha=alpha, beta=beta ),iter=100,chains=1, fit=compile_RB,refresh = -1)
    theta_sample[i]=extract(theta_sample_stan)[["theta"]][50]
    lambda_index=tempeture_sample(theta_sample[i], r, Z)
    lambda_sample[i]=lambda_index
    c[lambda_index]=c[lambda_index]+  1/N
  }
  if (min(c)==0)
  {
    log_c=log(c+exp(-10))
    c=exp(log_c)
    c=c/sum(c)
  }
  Z=Z*r[1]/r * c/c[1]
  Z.ST.iter[,iter_index]=log(Z)-log(Z[1])
  lambda_chain_4[iter_index,]=lambda_list[lambda_sample]
  theta_chain_4[iter_index,]=theta_sample
}
for( i in 1:20)
  L2_error[4,i]=  mean((     approx( lambda_list , Z.ST.iter[,i], lambda.grid.true)$y-log(z.grid.true)     )^2    )
 
 
#################################################
#Figures:  Est. of Normalising constant.
pdf("easy_z.pdf",width=6.5, height=2.5)
cex.label=0.5
par(mfrow=c(1,4),oma=c(2,2,2,0.5),mar=c(1,1,1,0) ,mgp=c(1.5,0.25,0), lwd=0.5,tck=-0.01, cex.axis=0.5, cex.lab=0.9, cex.main=0.9)
col.grid <- colorRampPalette( brewer.pal(9,'Blues'))(N_loop)
plot( g(a_dynamic[[1]]), z_dynamic[[1]],type='l',ylim=c(-4,0.5),col=col.grid[15],xlab ="",ylab="",axes=F)
lines( c(0,1),c(0,0),col='grey')
lines( lambda.grid.true, log(z.grid.true), type='l',col=2,lwd=1.5)
text(0.4,-2.4,"Ture", cex=cex.label,col=2)
text(0.6,-2.3,"iter 1", cex=cex.label,col='blue')
text(0.7,-.1,"initial guess", cex=cex.label,col='blue')
axis(1,lwd=0.5)
axis(2,lwd=0.5)
box(lwd=0.5)
mtext(1, text=expression(lambda),line=1,cex=0.8)
mtext(2, text="Log normalizing constant",line=1,cex=0.6)
mtext(3, text="with path sampling",line=0,cex=0.6)
mtext(3, text="continuous tempering",line=0.8,cex=0.6)

#######2.
col.grid <- colorRampPalette( brewer.pal(9,'Blues'))(N_loop)
plot( lambda.grid.iter,logZ.grid.iter[,1],type='l',ylim=c(-4,0.4),col=col.grid[2], xlab ="",ylab="",axes=F)
lines( c(0,1),c(0,0),col='grey')
for(i in 2:N_loop)
  lines( lambda.grid.iter, logZ.grid.iter[,i],type='l',col=col.grid[i])
lines( lambda.grid.true, log(z.grid.true), type='l',col=2,lwd=1.5)
text(0.4,-2.4,"Ture", cex=cex.label,col=2)
text(0.6,-1.3,"iter 1", cex=cex.label,col='blue')
text(0.45,-2.05,"iter 20", cex=cex.label,col='blue')
text(0.7,-.1,"initial guess", cex=cex.label,col='blue')
axis(1,lwd=0.5)
axis(2,lwd=0.5)
box(lwd=0.5)
mtext(1, text=expression(lambda),line=1,cex=0.6)
mtext(3, text="with emperical estimation",line=0,cex=0.6)
mtext(3, text="continuous tempering",line=0.8,cex=0.6)

 
######3.
col.grid <- colorRampPalette( brewer.pal(9,'Blues'))(N_RB_loop+1)
plot( lambda_list,Z.RB.iter[,1],type='l',ylim=c(-4,0.5),col=col.grid[2], xlab ="",ylab="",axes=F)
points(lambda_list,Z.RB.iter[,1],pch=19,col=col.grid[2],cex=0.5)
lines( c(0,1),c(0,0),col='grey')
for(i in 2:N_RB_loop){
  lines( lambda_list, Z.RB.iter[,i],type='l',col=col.grid[i+1])
  points(lambda_list,Z.RB.iter[,i],pch=19,col=col.grid[i+1],cex=0.5)
  
}
lines( lambda.grid.true, log(z.grid.true), type='l',col=2,lwd=1.5)
text(0.35,-2.2,"Ture", cex=cex.label,col=2)
text(0.6,-2.2,"iter 1", cex=cex.label,col='blue')
text(0.55,-2.5,"iter 5", cex=cex.label,col='blue')
text(0.7,-.1,"initial guess", cex=cex.label,col='blue')
axis(1,lwd=0.5)
axis(2,lwd=0.5)
box(lwd=0.5)
mtext(1, text=expression(lambda),line=1,cex=0.6)
mtext(3, text="Rao-Blackwellized",line=0,cex=0.6)
mtext(3, text="discrete tempering",line=0.8,cex=0.6)

#######4.
col.grid <- colorRampPalette( brewer.pal(9,'Blues'))(30+1)  #N_ST_loop+1
plot( lambda_list,Z.ST.iter[,1],type='l',ylim=c(-4,0.5),col=col.grid[2], xlab ="",ylab="",axes=F)
lines( c(0,1),c(0,0),col='grey')
for(i in 1:N_ST_loop){
  lines( lambda_list, Z.ST.iter[,i],type='l',col=col.grid[i+1])
  #points(lambda_list,Z.ST.iter[,i],pch=19,col=col.grid[i+1],cex=0.5)
}
for(i in c(20)){
  lines( lambda_list, Z.ST.iter[,i],type='l',col='blue' )
  points(lambda_list,Z.ST.iter[,i],pch=19,col='blue',cex=0.5)
}
lines( lambda.grid.true, log(z.grid.true), type='l',col=2,lwd=1.5)
text(0.44,-2.6,"Ture", cex=cex.label,col=2)
text(0.7,-2.8,"iter 20", cex=cex.label,col='blue')
text(0.8,-.1,"initial guess", cex=cex.label,col='blue')
axis(1,lwd=0.5)
axis(2,lwd=0.5)
box(lwd=0.5)
mtext(1, text=expression(lambda),line=1,cex=0.6)
mtext(3, text="discrete tempering",line=0.8,cex=0.6)
mtext(3, text="with emperical estimation",line=0,cex=0.6)
dev.off()




