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
  int<lower=1> M;
  int<lower=1> D;
  real<lower=0> mode_scale;
  vector[D] mode_centers[M];
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
  vector[D] theta;
  real<lower=0,upper=2> a;
}

transformed parameters {
  real log_posterior;
  real log_psi;
  real a_reflected;   
  real a_scaled;
  real lambda;
  real mode_densities[M];
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
  
  for (i in 1:M){
    mode_densities[i] = normal_lpdf(theta | mode_centers[i], mode_scale);
  }
  log_posterior = log_sum_exp(mode_densities);
  log_psi = normal_lpdf(theta | rep_vector(0,D),1.5);
}

model {
  target += lambda*b_linear;
  target += dot_product(kernel_logit, b_logit) + dot_product(kernel_gaussian, b_gaussian);
  target += (lambda)*log_posterior + (1-lambda)*log_psi;
}
