data {
  real<lower=0,upper=1> lambda;
  int<lower=1> M;
  int<lower=1> D;
  real<lower=0> mode_scale;
  vector[D] mode_centers[M];
}

parameters {
  vector[D] theta;
}

transformed parameters {
  real log_posterior;
  real log_psi;
  real mode_densities[M];
  for (i in 1:M){
    mode_densities[i] = normal_lpdf(theta | mode_centers[i], mode_scale);
  }
  log_posterior = log_sum_exp(mode_densities);
  log_psi = normal_lpdf(theta | rep_vector(0,D),1.5);
}

model {
  target += (lambda)*log_posterior + (1-lambda)*log_psi;
}
