data {
  real a_lower;
  real a_upper;
  int N;
  vector[N] a;
  vector[N] log_p;
  int K_logit;
  vector[K_logit] mu_logit;
  vector[K_logit] sigma_logit;
  int K_gaussian;
  vector[K_gaussian] mu_gaussian;
  vector[K_gaussian] sigma_gaussian;
}
transformed data {
  real a_scaled;
  vector[N] lambda;
  for (i in 1:N){
    a_scaled = (a[i] - a_lower)/(a_upper - a_lower);
  if (a_scaled <= 0)
    lambda[i] = 0;
  else if (a_scaled < 1)
    lambda[i] = (0.5 + 0.25*(3*(2*a_scaled - 1) - (2*a_scaled - 1)^3));
  else
    lambda[i] = 1;
  }
}
  
parameters {
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
    kernel_logit = 1 ./ (1 + exp(-(lambda[j] - mu_logit) ./ sigma_logit));  
    kernel_gaussian = exp(-.5*(lambda[j] - mu_gaussian) .* (lambda[j] - mu_gaussian) ./ (sigma_gaussian .* sigma_gaussian));
    fit[j] = lambda[j]*b_linear + dot_product(kernel_logit, b_logit) + dot_product(kernel_gaussian, b_gaussian);
  }
}
model {
  log_p ~ normal(fit, sigma_err);
  b_gaussian~normal(0,6);
}
