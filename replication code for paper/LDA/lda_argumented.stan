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
	//original:
  int<lower=2> K;               // num topics
  int<lower=2> V;               // num words
  int<lower=1> M;               // num docs
  int<lower=1> N;               // total word instances
  int<lower=1,upper=V> w[N];    // word n
  int<lower=1,upper=M> doc[N];  // doc ID for word n
  vector<lower=0>[K] alpha;     // topic prior
  vector<lower=0>[V] beta;      // word prior
}
transformed data {
	//argument:
	real b_linear;
	vector[K_logit] b_logit;
	vector[K_gaussian] b_gaussian;
	b_linear = b[1];
	b_logit = segment(b, 2, K_logit);
	b_gaussian = segment(b, 2 + K_logit, K_gaussian);
}
parameters {
	real<lower=0,upper=2> a;
  positive_ordered[K] theta_first;
  simplex[K] theta_ex_first[M-1];   // topic dist for doc m
  simplex[V] phi[K];     // word dist for topic k
}
transformed parameters {
	//argument transformed:
	real a_reflected;
	real a_scaled;
	real lambda;
	vector[K_logit] kernel_logit;
	vector[K_gaussian] kernel_gaussian;
	//original:
  simplex[K] theta_first_transform = theta_first / sum(theta_first);
  simplex[K] theta[M];
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

  theta[1]=theta_first_transform;
  theta[2:M]=theta_ex_first;
}
model {
	real log_q;//lp of the alternative model
   for(k in 1:K)
     theta_first[k]~gamma(alpha[k], 1);
  for (m in 1:(M-1))
    theta_ex_first[m] ~ dirichlet(alpha);  // prior
  for (k in 1:K)
    phi[k] ~ dirichlet(beta);     // prior
	log_q=target();
	target+=-target();
   for(k in 1:K)
     theta_first[k]~gamma(alpha[k], 1);
  for (m in 1:(M-1))
    theta_ex_first[m] ~ dirichlet(alpha);  // prior
  for (k in 1:K)
    phi[k] ~ dirichlet(beta);     // prior
  for(k in 1:K)
     theta_first[k]~gamma(alpha[k], 1);
  for (m in 1:(M-1))
    theta_ex_first[m] ~ dirichlet(alpha);  // prior
  for (k in 1:K)
    phi[k] ~ dirichlet(beta);     // prior
  for (n in 1:N) {
    real gamma[K];
    for (k in 1:K)
      gamma[k] = log(theta[doc[n], k]) + log(phi[k, w[n]]);
    target += log_sum_exp(gamma);  // likelihood;
  }
	target+=target()*(lambda)-target();
	target += (1-lambda)*log_q;
	target += lambda*b_linear;
	target += dot_product(kernel_logit, b_logit) + dot_product(kernel_gaussian, b_gaussian);
}



generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
     real gamma[K];
     for (k in 1:K)
       gamma[k] = log(theta[doc[n], k]) + log(phi[k, w[n]]);
     log_lik[n] = log_sum_exp(gamma); 
   }
}
