library(rstan)
rstan_options(auto_write = TRUE)
# currently only supports 1 chain
source("tempering_function.R") 

# generate the new stan file, please check if it is OK.
file_new=code_tempeture_argument(stan_file="example/cauchy_mixture.stan")

# compile the new working model
sampling_model=stan_model(file_new)
 


path_sample_fit=path_sample(sampling_model=sampling_model, data=list(gap=10), N_loop = 6, visualize_progress = TRUE )


sim_cauchy=extract(path_sample_fit$fit_main)
in_target= sim_cauchy$lambda==1
hist(sim_cauchy$theta[in_target], breaks = 30)
 