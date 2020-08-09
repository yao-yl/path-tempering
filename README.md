#  Adaptive Path Sampling in Metastable Posterior Distributions

Code for paper  *Adaptive Path Sampling in Metastable Posterior Distributions* by  Yuling Yao, Collin  Cademartori, Aki Vehtari, Andrew Gelman.

## Usage:
`source("tempering_function.R")` will load all functions.  

`code_tempeture_argument()` constructs a geometric bridge between two densitis. It can often be driven by three motivations:
1. In metastable sampling, we want to start from a easy to sampled base measrument, which can guide the target sampling. 
2. In model expansion, we want to work with multiple models. We may also want to fit two models at the same time. 
3. We want to compute the marginal likelihood, or Bayes factor of two models.

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
y ~  bernoulli(Phi(x* beta));
}
```

For the third use, simply put the prior part into the alternative model:
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
beta ~ normal(0,1);
}
alternative model{
beta ~ normal(0,1);
}
```


## Exmple:
Let's first constuct a simple bimodal density:
```
data {
  real gap;
}
parameters {
  real theta;
}
model{
gap ~ cauchy(theta, 0.5);
-gap ~ cauchy(theta, 0.5);
}
```
This will not work in stan with large `gap`. It exhibits posteiror bimodality  and a large R hat if running multiple chains. 

A User only need to specify a base model,  say N$(0,5)$,  and list it in an  \texttt{alternative model} block as if it is a regular model.   
 
```
...
model{ // keep the original model  
 y~ cauchy(theta,0.2);   
 -y~ cauchy(theta,0.2);   
}
alternative model{ // add a new block 
 theta~ normal(0,5);   
}
```

After saving this stan code to a  `cauchy.stan`, we run the function `code_tempeture_argument()` that automatically constructs a tempered path between the orginal model and the alternative model:

```
source("tempering_function.R")
file_new=code_tempeture_argument(stan_file="example/cauchy_mixture.stan")
> output:
> A new stan file has been created: cauchy_augmented.stan.
```

We have automated path sampling and its adaptation into a function `path_sample()`.
To run this example, call
```
library(rstan)
rstan_options(auto_write = TRUE)
# currently only supports 1 chain
# generate the new stan file, please check if it is OK.
# compile the new working model
sampling_model=stan_model(file_new)
path_sample_fit=path_sample(sampling_model=sampling_model, data=list(gap=10), N_loop = 6, visualize_progress = TRUE, iter_final=6000 )
```

The returned value `path_sample_fit` provides access to the draw  the posterior draws $\theta$ from the target density and base density, the 
joint path and the normalization constant:
```
sim_cauchy=extract(path_sample_fit$fit_main)
in_target= sim_cauchy$lambda==1
in_prior = sim_cauchy$lambda==0
# sample from the target 
hist(sim_cauchy$theta[in_target])
# sample from the base 
hist(sim_cauchy$theta[in_prior])
# the joint "path"
plot(sim_cauchy$a, sim_cauchy$theta)
# the normalization constant
plot(g_lambda(path_sample_fit$path_post_a), path_sample_fit$path_post_z)
```
![cauchy exampl output](/example/img/Cauchy.jpg)

 

