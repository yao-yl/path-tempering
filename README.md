#  Adaptive Path Sampling in Metastable Posterior Distributions

Code for paper  *Adaptive Path Sampling in Metastable Posterior Distributions* by  Yuling Yao, Collin  Cademartori, Aki Vehtari, Andrew Gelman.

## Usage:
`source("tempering_function.R")` will load all functions.  

`code_tempeture_argument()` constructs a geometric bridge betwwen two densitis. It can often be driven by two motivations:
1. In metastable sampling, we want to start from a easy to sampled base measrument, which can guide the target sampling. 
2. In model expansion, we want to work with multiple models. We may also want to fit two models at the same time. 

In either case, just add another module `alternative model` to your stan code. For example, you want to sample a bimdoal distribution with a hot start:

```
...
model{
theta ~ cauchy(-10, 0.5);
theta ~ cauchy(10, 0.5);
}
alternative model{
theta ~ normal (0,5);
}
```

Or maybe you are fitting a logit link and a probit link at the same time, which only requries to add an `alternative model` as if you have new models.

```
data{
int y[n];
real x[n];
}
parameter{
real beta;
}
model{
y ~  binomial_logit(1, x* beta);
}
alternative model{
y ~  bernoulli(Phi(x* beta);
}
```
For exmple, run the code 
```file_new=code_tempeture_argument(stan_file="example/cauchy_mixture.stan")```
will generate a new working stan model.

To run this example, use
```
library(rstan)
rstan_options(auto_write = TRUE)
# currently only supports 1 chain
source("tempering_function.R") 
# generate the new stan file, please check if it is OK.
file_new=code_tempeture_argument(stan_file="example/cauchy_mixture.stan")
# compile the new working model
sampling_model=stan_model(file_new)
path_sample_fit=path_sample(sampling_model=sampling_model, data=list(gap=10), N_loop = 6, visualize_progress = TRUE )

```

It is easy to monitor any any in target samples:
```
sim_cauchy=extract(path_sample_fit$fit_main)
in_target= sim_cauchy$lambda==1
hist(sim_cauchy$theta[in_target], breaks = 30)
```

