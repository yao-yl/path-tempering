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
	int<lower=2> V;               // num words
	int<lower=1> M;               // num docs
	int<lower=1> N;               // total word instances
	int<lower=1> N_test;
	int<lower=1> K;               // num topics
	int<lower=1,upper=V> w[N];    // word n
	int<lower=1,upper=V> w_test[N_test];    // word in test set
	int<lower=1,upper=M> doc[N];  // doc ID for word n
	int<lower=1,upper=M> doc_test[N_test];  // doc ID for test words
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
  simplex[K] theta[M];   // topic dist for doc m
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
   real log_lik_sum=0;

  
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
  for (n in 1:N) {
    real gamma[K];
    for (k in 1:K)
      gamma[k] = log(theta[doc[n], k]) + log(phi[k, w[n]]);
    log_lik_sum += log_sum_exp(gamma);  // likelihood;
  }
}

model {
	for (m in 1:M)
		theta[m] ~ dirichlet(alpha);  // prior
	for (k in 1:K)
		phi[k] ~ dirichlet(beta);     // prior
	target+=log_lik_sum*lambda;
	target += lambda*b_linear;
	target += dot_product(kernel_logit, b_logit) + dot_product(kernel_gaussian, b_gaussian);
	}

generated quantities{
	real log_lik_test[N_test];
	for (n in 1:N_test) {
		real temp[K];
		for (k in 1:K)
			temp[k] = log(theta[doc_test[n], k]) + log(phi[k, w_test[n]]);
		log_lik_test[n] = log_sum_exp(temp);
	}
}