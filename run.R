library(rstan)
rstan_options(auto_write = TRUE)
# currently only supports 1 chain
source("tempering_function.R") 

# generate the new stan file, please check if it is OK.
file_new=code_tempeture_argument(stan_file="example/cauchy_mixture.stan")

# compile the new working model
sampling_model=stan_model(file_new)
 

path_sample_fit=path_sample(sampling_model=sampling_model, data=list(gap=10), N_loop = 6, visualize_progress = TRUE,iter_final=6000)


# collect final samples
sim_cauchy=extract(path_sample_fit$fit_main)
in_target= sim_cauchy$lambda==1
sim_cauchy=extract(path_sample_fit$fit_main)
in_prior = sim_cauchy$lambda==0


# sample from the target 
hist(sim_cauchy$theta[in_target], breaks = 30, prob=T)
# sample from the base 
hist(sim_cauchy$theta[in_prior], breaks = 30)
# the "path"
plot(sim_cauchy$a, sim_cauchy$theta)
# the normalization constant
plot(g_lambda(path_sample_fit$path_post_a), path_sample_fit$path_post_z)


