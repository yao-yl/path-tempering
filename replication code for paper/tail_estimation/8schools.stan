data {
  int<lower=1> N;
  real y[N];
  real<lower=0> sds[N];
  
  int<lower=0> K_logit;
  vector[K_logit] mu_logit;
  vector[K_logit] sigma_logit;
  int<lower=0> K_gaussian;
  vector[K_gaussian] mu_gaussian;
  vector[K_gaussian] sigma_gaussian;
  vector[1+K_logit+K_gaussian] b;
  
  //real<lower=0> a;
  real<lower=0> lbound;
  real<lower=0> ubound;
}

transformed data{
  real b_linear;
  vector[K_logit] b_logit;
  vector[K_gaussian] b_gaussian;
  
  b_linear = b[1];
  if (K_logit > 0){
    b_logit = segment(b, 2, K_logit);  
  } 
  if (K_gaussian > 0){
    b_gaussian = segment(b, 2 + K_logit, K_gaussian);  
  }
}

parameters {
  //vector[N] theta_un;
  vector[N] theta;
  real<lower=lbound,upper=ubound> tau;
  real mu;
}

transformed parameters {
  vector[K_logit] kernel_logit;
  vector[K_gaussian] kernel_gaussian;
  //vector[N] theta;
  
  //theta = tau*theta_un + mu;
  kernel_logit = 1 ./ (1 + exp(-(tau - mu_logit) ./ sigma_logit));
  kernel_gaussian = exp(-.5*(tau - mu_gaussian) .* (tau - mu_gaussian) ./ (sigma_gaussian .* sigma_gaussian));
}

model {
  target += normal_lpdf(y|theta,sds);
  //target += normal_lpdf(theta_un|0,1);
  target += normal_lpdf(theta|mu,tau);
  target += -1*tau*b_linear;
  target += -1*(dot_product(kernel_logit, b_logit) + dot_product(kernel_gaussian, b_gaussian));
  //target += normal_lpdf(tau|0,30);
}
