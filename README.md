#   Adaptive Path Sampling and Continuous Tempering

This package implements adaptive path sampling. It iteratively reduces the gap between the proposal and the target density, and provide a reliable  normalizing constant  estimation with practical diagnostic using importance sampling theory. 
By equipping simulated tempering with a continuous temperature, path tempering enables efficient sampling from multimodal densities.    

Reference:
*Adaptive Path Sampling in Metastable Posterior Distributions* by  Yuling Yao, Collin  Cademartori, Aki Vehtari, Andrew Gelman.
([Replication code]<https://github.com/yao-yl/path-tempering/tree/master/replication%20code%20for%20paper>)

## Installation in R
```R
library(devtools)
devtools::install_github("yao-yl/path-tempering/package/pathtemp",  upgrade="never")
``` 

## Semaless integration in Stan
The new `alternative model` block in Stan enables more than one model to be fit at the same time. To fully express a continous expansoon of these two models, a geometric bridge is contrtucted between two posterior densitis. Such design can often be driven by three motivations:
1. In metastable sampling, we want to start from an easy-to-sample base density, which guides the target sampling. 
2. In model expansion, we want to work with multiple models. We may also want to fit two models at the same time for both computation efficiency and model comparison. 
3. We want to compute the marginal likelihood, or Bayes factor of two models.

In any of these cases, just add another module `alternative model` to your stan code. For example, you want to sample a bimdoal distribution with a hot start:

```stan
parameter{
real theta;
}
model{
theta ~ cauchy(-10, 0.5);
theta ~ cauchy(10, 0.5);
}
alternative model{
theta ~ normal (0,5);
}
```

Or maybe you are fitting a logit link and a probit link regression at the same time, which only requires to add an `alternative model` as if you have new models.

```stan
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
```stan
data{
int y[n];
real x[n];
}
parameter{
real beta;
}
model{
target +=  binomial_logit_ipmf(y| 1, x* beta);
beta ~ normal(0,1);
}
alternative model{
beta ~ normal(0,1);
}
```
Notably,  the prior is allowed to be unnormalized (when we use the sampling statement `~`). the normalization constant of the prior appear twice and has no effect for the marginal likelihood computation.


## Example:
### Multimodal sampling
Let's first construct a simple bimodal density:
```stan
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
This will not work in Stan with large `gap`. It exhibits posterior bimodality  and a large R-hat if running multiple chains. 

A user only needs to specify a base model,  which can be the prior, or some better-shaped posterior approximation.   Just write it in an  \texttt{alternative model} block as if it is a regular model.    
 
```stan
...
model{ // keep the original model  
 gap~ cauchy(theta,0.2);   
 -gap~ cauchy(theta,0.2);   
}
alternative model{ // add a new block 
 theta~ normal(0,5);   
}
```

After saving this stan code to a  `cauchy.stan`, we run the function `code_temperature_augment()` that automatically constructs a tempered path between the original model and the alternative model:

```
library(pathtemp)
file_new=code_temperature_augment(stan_file="example/cauchy_mixture.stan")
> output:
> A new stan file has been created: cauchy_augmented.stan.
```

We have automated path sampling and its adaptation into a function `path_sample()`.
To run this example, call
```
library(rstan)
rstan_options(auto_write = TRUE)
# compile stan optimizer
update_model <- stan_model("solve_tempering.stan")# need to compile it first
# currently only supports 1 chain
# generate the new stan file, please check if it is OK.
# compile the new working model
sampling_model=stan_model(file_new)
path_sample_fit=path_sample(sampling_model=sampling_model, data=list(gap=10), N_loop = 6, visualize_progress = TRUE, iter_final=6000 )
```
Here `solve_tempering.stan` needs to be downloaded first. 

The returned value `path_sample_fit` provides access to the posterior draws from the target and base density, the joint path,  and the log normalization constant:
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
Here is the output:
![cauchy exampl output](/example/img/Cauchy.jpg)

 
### Fitting logit and probit models together
In this example, we fit two models (probit and logiistic regressions)in one joint sampling.  A path between them effectively expands the model continuously such both individual model are special cases of the augmented model. The computation efficiency is enhanced as we are fitting one slightly larger model rather than fitting two models. In addition,  the log normalization constant tells which models fits the data more, which is related to but not the same as log Bayes factor in model comparisons. 


Let's consider a logit and probit link regression.  It is  avaialbe in `example/logit.stan`.
```stan
data {
   int n;
   int y[n];
   real x[n];
}

parameters {
   real beta;
}
model {
   for (i in 1:n)
     y[i]~ bernoulli_logit(beta * x[i]);
}
alternative model {
   for (i in 1:n)
     y[i]~  bernoulli(Phi(beta * x[i]));
}

```

It is recommended to write `bernoulli_lpdf` if the goal is to compute the marginal likelihood.  

We generate some x and y to fit this model
```R
file_new=code_temperature_augmented(stan_file="logit/logit.stan")
sampling_model=stan_model(file_new)
update_model <- stan_model("solve_tempering.stan")
set.seed(100)
n=30
x=runif(n, 0,1)
y=rbinom(n=n, size=1, prob=exp(3*x)/max( exp(3*x) +1) )
path_sample_fit=path_sample(sampling_model=sampling_model, data=list(n=n, x=x, y=y), N_loop = 6, visualize_progress = TRUE,iter_final=6000)
```

Here is the output:
![logit exampl output](/example/img/logit.jpg)
The log z here is NOT the Bayes factor! But it can be interpreted as the posterior density of lambda: hence telling which model is more supported by the the data.  
Not surprisingly, the logit model fits the data better as that is how we generate y.


