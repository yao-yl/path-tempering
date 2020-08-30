data {
  int<lower=1> N;
  real y[N];
  real<lower=0> sds[N];
  real<lower=0> lbound;
  real<lower=0> ubound;
}

parameters {
  vector[N] theta_un;
  real<lower=lbound,upper=ubound> tau;
  real mu;
}

transformed parameters {
  vector[N] theta;
  theta = tau*theta_un + mu;
}

model {
  target += normal_lpdf(y|theta,sds);
  target += normal_lpdf(theta_un|0,1);
  target += inv_gamma_lpdf(tau|1.2,20);
}
