functions {
  real flower_lpdf(vector x, real omega, real A){
    real rate = (-1*(sqrt(sum(square(x))) - 10 -A*cos(omega*atan2(x[2],x[1]))))/2;
    return(rate);
  }
}

data {
  real<lower=0,upper=1> lambda;
  real<lower=0> omega;
  real<lower=0> A;
}

parameters {
  vector[2] theta;
}

transformed parameters {
  real log_posterior;
  real log_psi;
  log_posterior = flower_lpdf(theta | omega, A);
  log_psi = normal_lpdf(theta | rep_vector(0,2),10);
}

model {
  target += (lambda)*log_posterior + (1-lambda)*log_psi;
}
