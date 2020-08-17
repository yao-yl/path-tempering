
data {
		//augmented:
	real a_lower;
	real a_upper;
	int K_logit;
	vector[K_logit] mu_logit;
	vector[K_logit] sigma_logit;
	int K_gaussian;
	vector[K_gaussian] mu_gaussian;
	vector[K_gaussian] sigma_gaussian;
	vector[1 + K_logit + K_gaussian] b;
	
	int<lower=0> n;				      // number of observations
	int<lower=0> d;             // number of predictors
	//int<lower=0> n_test;				      // number of observations
	int<lower=0,upper=1> y[n];	// outputs
	matrix[n,d] x;				      // inputs
	//int<lower=0,upper=1> y_test[n_test];	// outputs
//	matrix[n_test,d] x_test;				      // inputs
	real<lower=0> scale_icept;	// prior std for the intercept
	real<lower=0> scale_global;	// scale for the half-t prior for tau
	real<lower=0> slab_scale;
	real<lower=0> slab_df;
}
transformed data {
	//augmented:
	real b_linear;
	vector[K_logit] b_logit;
	vector[K_gaussian] b_gaussian;
	//original:
	b_linear = b[1];
	b_logit = segment(b, 2, K_logit);
	b_gaussian = segment(b, 2 + K_logit, K_gaussian);
}
parameters {
	real<lower=0,upper=2> a;
	real beta0; // intercept
	vector[d] z; // auxiliary parameter
	real<lower=0> tau;			// global shrinkage parameter
	vector<lower=0>[d] lam;	// local shrinkage parameter
	real<lower=0> caux; // auxiliary
}
transformed parameters {
	real log_p;
	real log_q;	
	real log_lik_sum;
	real<lower=0> c;
	vector[d] beta;				// regression coefficients
	vector[n] f;				// latent values
	vector[n] f_phi;				// latent values
	vector<lower=0>[d] lambda_tilde;
 	//augmented transformed:
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
	c = slab_scale * sqrt(caux);
	lambda_tilde = sqrt( c^2 * square(lam) ./ (c^2 + tau^2* square(lam)) );
	beta = z .* lambda_tilde*tau;
	f = beta0 + x*beta;
	f_phi=Phi(f);
	log_p =  bernoulli_lpmf(y|f_phi);
	log_q = bernoulli_logit_lpmf(y|f);
	log_lik_sum=log_p-log_q;
}

model {
	z ~ normal(0,1);
	lam ~ cauchy(0,1);
	tau ~ cauchy(0, scale_global);
	caux ~ inv_gamma(0.5*slab_df, 0.5*slab_df);
	beta0 ~ normal(0,scale_icept);
	target+=log_p*lambda + log_q * (1-lambda);
	target += lambda*b_linear;
	target += dot_product(kernel_logit, b_logit) + dot_product(kernel_gaussian, b_gaussian);
}


 


