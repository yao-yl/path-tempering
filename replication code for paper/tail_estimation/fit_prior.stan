data {
  int N;
  vector[N] log_p;
  vector[N] x;
  int K_logit;
  vector[K_logit] mu_logit;
  vector[K_logit] sigma_logit;
  int K_gaussian;
  vector[K_gaussian] mu_gaussian;
  vector[K_gaussian] sigma_gaussian;
}
  
parameters {
  real intercept;
  vector[1 + K_logit + K_gaussian] b;
  real<lower=0> sigma_err;
}

transformed parameters {
  real b_linear;
  vector[K_logit] b_logit;
  vector[K_gaussian] b_gaussian;
  vector[N] fit;
  vector[K_logit] kernel_logit;
  vector[K_gaussian] kernel_gaussian;
  b_linear = b[1];
  b_logit = segment(b, 2, K_logit);
  b_gaussian = segment(b, 2 + K_logit, K_gaussian);
  for (j in 1:N){
    kernel_logit = 1 ./ (1 + exp(-(x[j] - mu_logit) ./ sigma_logit));  
    kernel_gaussian = exp(-.5*(x[j] - mu_gaussian) .* (x[j] - mu_gaussian) ./ (sigma_gaussian .* sigma_gaussian));
    fit[j] = x[j]*b_linear + dot_product(kernel_logit, b_logit) + dot_product(kernel_gaussian, b_gaussian);
  }
}

model {
  log_p ~ normal(fit, sigma_err);
}
