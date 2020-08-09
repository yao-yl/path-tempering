add_line=function(string_vec, add_line_id, add_on){
	if(length(string_vec)<add_line_id)
		stop("wrong line id")
	temp1=string_vec[1:add_line_id]
	temp2=string_vec[-c(1:add_line_id)]
	return(c(temp1,add_on, temp2))
}


block_start=function(stan_code){
	grep("{",stan_code, fixed = TRUE)  
	# "Fixed = TRUE" disables regex
}
block_end=function(stan_code) {
	temp=grep("\\}",stan_code)  
	temp[!temp %in% grep(" \\}",stan_code)]   
}

arg_string_data=c(
	"\t//argument:",
	"\treal a_lower;",
	"\treal a_upper;",
	"\tint K_logit;",
	"\tvector[K_logit] mu_logit;",
	"\tvector[K_logit] sigma_logit;",
	"\tint K_gaussian;",
	"\tvector[K_gaussian] mu_gaussian;",
	"\tvector[K_gaussian] sigma_gaussian;",
	"\tvector[1 + K_logit + K_gaussian] b;",
	"\t//original:"
	)

arg_string_trans_data_pre=c("\t//argument:",
"\treal b_linear;",
"\tvector[K_logit] b_logit;",
"\tvector[K_gaussian] b_gaussian;",
"\t//original:"
)


arg_string_trans_data_suffix=c(
	"\tb_linear = b[1];",
	"\tb_logit = segment(b, 2, K_logit);",
	"\tb_gaussian = segment(b, 2 + K_logit, K_gaussian);"
)


arg_string_trans_parameter_pre=c("\t//argument transformed:",
	"\treal a_reflected;",    
	"\treal a_scaled;",
	"\treal lambda;",
	"\tvector[K_logit] kernel_logit;",
	"\tvector[K_gaussian] kernel_gaussian;",
	"\t//original:")
	
	
arg_string_trans_parameter_suffix=c(	
	"\ta_reflected = 1 - fabs(1 - a);",
	"\ta_scaled = (a_reflected - a_lower)/(a_upper - a_lower);",
	"\tif (a_scaled <= 0)",
	"\t	lambda = 0;",
	"\telse if (a_scaled < 1)",
		"\t	lambda = 0.5 + 0.25*(3*(2*a_scaled - 1) - (2*a_scaled - 1)^3);",
	"\telse",
	"\t	lambda = 1;",
	"\tkernel_logit = 1 ./ (1 + exp(-(lambda - mu_logit) ./ sigma_logit));",
	"\tkernel_gaussian = exp(-.5*(lambda - mu_gaussian) .* (lambda - mu_gaussian) ./ (sigma_gaussian .* sigma_gaussian));"
)

arg_string_model_define=c(
	"\treal log_q;//lp of the alternative model",
	"\treal base_jacobian;//lp of the jacobian",
  "\treal log_p;//lp of the main model",
	"\tbase_jacobian = target();"
)



arg_string_model=c(
	"\tlog_p = target()-base_jacobian;",
	"\tprint(\" log_p = \", log_p,\" log_q = \", log_q, \" a = \", a);",
	"\ttarget += (lambda-1)*log_p;",
	"\ttarget += (1-lambda)*log_q;",
	"\ttarget += lambda*b_linear;",
	"\ttarget += dot_product(kernel_logit, b_logit) + dot_product(kernel_gaussian, b_gaussian);"
)

code_tempeture_argument=function(stan_file=NULL){
	stan_code=readLines(stan_file)
	if(sum(grep(" a;| a\\[| a ", stan_code))+ sum(grep(" lambda;| lambda\\[| lambda ", stan_code))>0 ) 
		stop("Please reserve parameter name a and lambda")
	# locate data block
	line_data_start=which(grepl("data",stan_code)&!grepl("transformed",stan_code))
	if(length(line_data_start)>1)
		stop("must contain one data block!")
	if(length(line_data_start)==1){
	line_data_start=block_start(stan_code)[min (which(block_start(stan_code) >= line_data_start))]
	stan_code_new= add_line(string_vec=stan_code, add_line_id=line_data_start , arg_string_data)
 	}else{
	stan_code_new= c( "data{", arg_string_data, "}", stan_code)
 	}
	
	# locate transformed data block
	line_transformed_start=which(grepl("transformed",stan_code_new) &grepl("data",stan_code_new) )
	if(length(line_transformed_start)>1){stop("wrong format of transformed data")
	} else if(length(line_transformed_start)==1) {
		line_transformed_start=block_start(stan_code_new) [min (which(block_start(stan_code_new) >= line_transformed_start))]
		stan_code_new= add_line(string_vec=stan_code_new, add_line_id=line_transformed_start , arg_string_trans_data_pre)
	  line_transformed_end=block_end(stan_code_new)[min (which(block_end(stan_code_new) >= line_transformed_start))]
		stan_code_new= add_line(string_vec=stan_code_new, add_line_id=line_transformed_end-1,  add_on=arg_string_trans_data_suffix)
	}else{
		line_data_start=which(grepl("data",stan_code_new)&!grepl("transformed",stan_code_new))
		line_transformed_start=block_end(stan_code_new)[min (which(block_end(stan_code_new) >= line_data_start))]
		stan_code_new= add_line(string_vec=stan_code_new, add_line_id=line_transformed_start , c("transformed data {", arg_string_trans_data_pre,arg_string_trans_data_suffix, "}"))
	}
	
	# locate parameters
	line_parameter_start=which(!grepl("transformed",stan_code_new) & grepl("parameters",stan_code_new) )
	if(length(line_parameter_start)!=1)
		stop("must contain one parameter block!")
	line_parameter_start=block_start(stan_code_new)[ (min(which(block_start(stan_code_new) >= line_parameter_start))  ) ]
	stan_code_new= add_line(string_vec=stan_code_new, add_line_id=line_parameter_start , "\treal<lower=0,upper=2> a;")
	
	# transformed parameters
	line_transformed_start=which(grepl("transformed",stan_code_new) &grepl("parameters",stan_code_new) )
	if(length(line_transformed_start)>1){stop("wrong format of transformed parameters")
	} else if(length(line_transformed_start)==1) {
		line_transformed_start=block_start(stan_code_new) [min (which(block_start(stan_code_new) >= line_transformed_start))]
		stan_code_new= add_line(string_vec=stan_code_new, add_line_id=line_transformed_start , arg_string_trans_parameter_pre)
		line_transformed_end=block_end(stan_code_new)[min (which(block_end(stan_code_new) >= line_transformed_start))]
		stan_code_new= add_line(string_vec=stan_code_new, add_line_id=line_transformed_end-1,  add_on=arg_string_trans_parameter_suffix)
	}else{
		line_transformed_start=block_end(stan_code_new)[min (which(block_end(stan_code_new) >= line_parameter_start))]
		stan_code_new= add_line(string_vec=stan_code_new, add_line_id=line_transformed_start , c("transformed parameters {", arg_string_trans_parameter_pre,arg_string_trans_parameter_suffix,	 "}"))			
	}
	line_model_alternative_start=which(grepl("alternative",stan_code_new) & grepl("model",stan_code_new) )
	if(length(line_model_alternative_start)!=1)
		stop("must contain an alternative model")
	flag_1=block_start(stan_code_new)[min (which(block_start(stan_code_new) >= line_model_alternative_start))]
	flag_2=block_end(stan_code_new)[min (which(block_end(stan_code_new) >= line_model_alternative_start))]
	if(flag_2-flag_1<2) 
		stop("must contain an alternative model")
	alt_model=stan_code_new[c((flag_1+1):(flag_2-1))]
	stan_code_new=stan_code_new[-c(flag_1:flag_2)]
	
	line_model_start=which(!grepl("alternative",stan_code_new) & grepl("model",stan_code_new) )
	if(length(line_model_start)!=1)
		stop("must contain model block!")
	
	
	line_model_start=block_start(stan_code_new)[min (which(block_start(stan_code_new) >= line_model_start))]
	alt_model=c(alt_model, "\tlog_q = target()-base_jacobian;",
							"\ttarget += -log_q;")
	stan_code_new= add_line(string_vec=stan_code_new, add_line_id=line_model_start,  
													add_on=	c(arg_string_model_define, alt_model))
	line_model_end=block_end(stan_code_new)[min (which(block_end(stan_code_new) >= line_model_start))]
	stan_code_new= add_line(string_vec=stan_code_new, add_line_id=line_model_end-1,  add_on=arg_string_model)
	
	
	#generated quantities
	line_generated_start=which(grepl("generated",stan_code_new) &grepl("quantities",stan_code_new) )
	if(length(line_transformed_start)>1){stop("wrong format of generated quantities")
	} else if(length(line_generated_start)==1) {
		line_generated_start=block_start(stan_code_new) [min (which(block_start(stan_code_new) >= line_generated_start))]
		stan_code_new= add_line(string_vec=stan_code_new, add_line_id=line_generated_start , "\tprint(\"----\");")
	}else{
		line_model_end=block_end(stan_code_new)[min (which(block_end(stan_code_new) >= line_model_end))]
		stan_code_new= add_line(string_vec=stan_code_new, add_line_id=line_model_end , 
														add_on = c("generated quantities {", "\tprint(\"----\");", "}"))
		
	}
	
	
	tempering_file=sub(".stan","_augmented.stan" ,stan_file)	
	writeLines(stan_code_new, tempering_file)
	cat(paste("A new stan file has been created:",tempering_file, "\n"  ))
	return(tempering_file)
}

path_sampling <- function(a, u, a_lower=0, a_upper=1){
	keep <- (a > a_lower) & (a < a_upper)
	if (sum(keep) > 0){
		a_uniq <- sort(unique(a[keep]))
		N_uniq <- length(a_uniq)
		u_bar <- rep(NA, N_uniq)
		for (i in 1:N_uniq){
			ok <- a == a_uniq[i]
			u_bar[i] <- mean(u[ok])
		}
		width <- c(a_uniq, a_upper) - c(a_lower, a_uniq)
		u_bar_avg <- (c(u_bar, u_bar[N_uniq]) + c(u_bar[1], u_bar))/ 2
		log_z <- c(0, cumsum(width*u_bar_avg))
		return(list(a=c(0, a_lower, a_uniq, a_upper, 1), log_z=c(log_z[1], log_z, log_z[length(log_z)])))
	}
	else
		return(list(a=c(0, 1), log_z=c(0, 0)))
}






path_gradients <- function(a_lower, a_upper, b, K_logit, mu_logit, sigma_logit, K_gaussian, mu_gaussian, sigma_gaussian, a, log_lik){
	N <- length(a)
	a_reflected = 1 - abs(1 - a)
	a_scaled <- (a_reflected - a_lower)/(a_upper - a_lower)
	lambda <- ifelse(a_scaled <= 0, 0, ifelse(a_scaled < 1, 0.5 + 0.25*(3*(2*a_scaled - 1) - (2*a_scaled - 1)^3), 1))
	grad_lambda = ifelse(a_scaled <= 0, 0, ifelse(a_scaled < 1, (2*0.25*(3 - 3*(2*a_scaled - 1)^2)/(a_upper - a_lower)), 0))
	u_prior <- rep(NA, N)
	kernel_logit <- array(NA, c(N, K_logit))
	kernel_gaussian <- array(NA, c(N, K_gaussian))
	grad_logit <- array(NA, c(N, K_logit))
	grad_gaussian <- array(NA, c(N, K_gaussian))
	b_linear <- b[1]
	b_logit <- b[1 + (1:K_logit)]
	b_gaussian <- b[1 + K_logit + (1:K_gaussian)]
	for (i in 1:N){
		kernel_logit[i,] <- 1 / (1 + exp(-(lambda[i] - mu_logit) / sigma_logit));
		kernel_gaussian[i,] = exp(-.5*(lambda[i] - mu_gaussian) * (lambda[i] - mu_gaussian) / (sigma_gaussian * sigma_gaussian))
		grad_logit[i,] <-  kernel_logit[i,] * kernel_logit[i,] * exp(-(lambda[i] - mu_logit) / sigma_logit) / sigma_logit
		grad_gaussian[i,] = kernel_gaussian[i,] * (mu_gaussian - lambda[i]) / (sigma_gaussian * sigma_gaussian)
		u_prior[i] = grad_lambda[i] * (b_linear + sum(grad_logit[i,]*b_logit) + sum(grad_gaussian[i,]*b_gaussian))
	}
	u_lik <- rep(NA, N)
	for (i in 1:N){
		u_lik[i] <- grad_lambda[i] * (log_lik[i] )
	}
	u_post <- u_prior + u_lik
	return(list(u_prior=u_prior, u_lik=u_lik, u_post=u_post))
}




extract_sim_string=function(line_string){
	match_a=gregexpr("a = ", line_string)[[1]][1]
	a_string=substring(line_string, first = match_a+4)
	match_p=gregexpr("log_p = ", line_string)[[1]][1]
	p_string=substring(line_string, first = match_p+8)
	match_space=gregexpr(" ", p_string)[[1]][1]
	p_string=substring(p_string, first=1, last =match_space-1)
	
	match_q=gregexpr("log_q = ", line_string)[[1]][1]
	q_string=substring(line_string, first = match_q+8)
	match_space=gregexpr(" ", q_string)[[1]][1]
	q_string=substring(q_string, first=1, last =match_space-1)
	return( c(as.numeric(p_string),  as.numeric(q_string),   as.numeric(a_string) ))
}				



path_initialize=function(){
	a_lower <- .1
	a_upper <- .8
	K_logit <- 20
	mu_logit <-  a_lower + (a_upper - a_lower)*seq(1,2*K_logit-1,2)/(K_logit)
	sigma_logit <- 2*rep((a_upper - a_lower)/K_logit, K_logit)
	K_gaussian <- K_logit
	mu_gaussian <- mu_logit
	sigma_gaussian <- sigma_logit
	all_a=NULL
	all_log_lik=NULL
	return(list(a_lower=a_lower, a_upper=a_upper, b=b, K_logit=K_logit, 
							mu_logit=mu_logit, sigma_logit=sigma_logit, 
							K_gaussian=K_gaussian, mu_gaussian=mu_gaussian, 
							sigma_gaussian=sigma_gaussian, all_a=all_a,all_log_lik=all_log_lik)) 
}



extract_lp=function(){
	stan_out=readLines('stan-output_path.txt')
	stan_out=stan_out[-(1:grep("(Sampling)",stan_out)[1])]  
	stan_out=stan_out[- which(stan_out=="") ]  
	stan_out=stan_out[-  grep("Iteration",stan_out)  ]
	flag=grep("----", stan_out)
	flag=flag[-length(flag)]
	extract_sim_string_vec=t(sapply(stan_out[flag+1] ,extract_sim_string))
	rownames(extract_sim_string_vec)=NULL
	log_lik=extract_sim_string_vec[,1]-extract_sim_string_vec[,2] #log_p-log_q
	return(list(log_lik=log_lik,a= extract_sim_string_vec[,3]))
}



path_fit <- function(sampling_model, data_list=list(), a_lower, a_upper, b, K_logit, 
										 mu_logit, sigma_logit, K_gaussian, 
										 mu_gaussian, sigma_gaussian, all_a, all_log_lik,
										 N_grid=100, iter=2000, max_treedepth=8,thin=2,
										 visualize_progress=FALSE){
	sink('stan-output_path.txt', append=FALSE)
	fit_main=sampling(sampling_model, data=c(list(a_lower=a_lower, a_upper=a_upper, 
																			b=b, K_logit=K_logit, mu_logit=mu_logit, 
																			sigma_logit=sigma_logit, K_gaussian=K_gaussian, 
																			mu_gaussian=mu_gaussian, 
																			sigma_gaussian=sigma_gaussian), data_list),
								chains=1, iter=iter, thin=thin, control=list(max_treedepth=max_treedepth)) # only support 1 chain
	sink()
	sim=extract(fit_main, permuted=FALSE)
	extract_lp_maxrix=extract_lp()
	a=sim[,,"a"]
	a=a[-length(a)]
	if( max (abs( extract_lp_maxrix$a - a))>0.01)
		stop("Mismatch, report a bug.")
	log_lik=extract_lp_maxrix$log_lik
	N_sims <- length(a)
	all_a <- c(all_a, a)
	all_log_lik <- c(all_log_lik, log_lik)
	gradients <- path_gradients(a_lower, a_upper, b, K_logit, mu_logit, sigma_logit, K_gaussian, mu_gaussian, sigma_gaussian, all_a, all_log_lik)
	all_a_reflected = 1 - abs(1 - all_a)
	path_post <- path_sampling(all_a_reflected, gradients$u_post, a_lower, a_upper)
	fit_measure <- max(abs(path_post$log_z))
	path_lik <- path_sampling(all_a_reflected, gradients$u_lik, a_lower, a_upper)
	approx_grid <- approx(path_lik$a, path_lik$log_z, seq(0, 1, length=N_grid))
	a_grid <- approx_grid$x
	log_p_grid <- approx_grid$y
	scael_max=max(max(abs(log_p_grid))/30, 1)
	path_fit <- optimizing(update_model, data=list(N=N_grid, a=a_grid, log_p=log_p_grid/scael_max, K_logit=K_logit, mu_logit=mu_logit, sigma_logit=sigma_logit, K_gaussian=K_gaussian, mu_gaussian=mu_gaussian, sigma_gaussian=sigma_gaussian, a_lower=a_lower, a_upper=a_upper), as_vector=FALSE, iter=6000)
	curve_fit <- path_fit$par$fit*scael_max
	b <- path_fit$par$b*scael_max
	a_reflected = 1 - abs(1 - a)
	in_target <- a_reflected > a_upper
	in_prior=a_reflected<a_lower
	in_prior_rate=mean(in_prior)
	if(in_prior_rate>0.8){
		in_prior_rate=in_prior_rate*0.95
		b[1]=in_prior_rate*	(log_mean_exp(log_lik_sum[in_prior])) + (1-in_prior_rate)*b[1]
	}
	if(visualize_progress==TRUE){
	hist(a_reflected, xlim=c(0, 1), xaxs="i", breaks=seq(0, 1, .01), probability = TRUE,
			 main="Post. draws of a_reflected", xlab="a_reflected", cex.main=.9, cex.axis=.9, cex.lab=.9)
	abline(v=a_lower, col="gray")
	abline(v=a_upper, col="gray")
	plot(path_lik$a, path_lik$log_z, xlim=c(0, 1), xaxs="i", bty="l", pch=20, cex=.2, col="blue", main="Path est of log norml. const.", xlab="a", ylab="log z", cex.main=.9, cex.axis=.9, cex.lab=.9)
	lines(a_grid, curve_fit, lwd=.5, col="red")
	abline(v=a_lower, col="gray")
	abline(v=a_upper, col="gray")
	plot(path_post$a, path_post$log_z, xlim=c(0, 1), xaxs="i", bty="l", pch=20, cex=.2, col="blue", main="Path est of log marginal", xlab="a", ylab="log p", cex.main=.9, cex.axis=.9, cex.lab=.9)
	abline(v=a_lower, col="gray")
	abline(v=a_upper, col="gray")
	}
	return(list(fit_main=fit_main, b=b, all_a=all_a, all_log_lik=all_log_lik, fit_measure=fit_measure, path_post_a=path_lik$a, path_post_z=path_lik$log_z))
}




path_sample=function(sampling_model,data_list, visualize_progress=FALSE,
										 N_loop =10, max_sample_store=4000,
										 a_lower =.1, a_upper = .8,  K_logit = 20, K_gaussian=20,  N_grid=100,   
										 thin=2, iter=2000,  max_treedepth=10,  iter_final=4000){
	
	# initialization:
	b <- rep(0, 1 + K_logit + K_gaussian)
	mu_logit <-  a_lower + (a_upper - a_lower)*seq(1,2*K_logit-1,2)/(K_logit)
	sigma_logit <- 2*rep((a_upper - a_lower)/K_logit, K_logit)
	mu_gaussian <- mu_logit
	sigma_gaussian <- sigma_logit
	all_a=NULL
	all_log_lik=NULL
	
	if(visualize_progress==TRUE)
		par(mfrow=c(N_loop,3), mar=c(2,2,2,1), mgp=c(1.5,.5,0), tck=-.01)
	for (i in 1:N_loop){
		fit <- path_fit(sampling_model=sampling_model, data_list=data_list,  
										a_lower, a_upper, b, K_logit, mu_logit, sigma_logit, K_gaussian, mu_gaussian, sigma_gaussian, 
										all_a, all_log_lik, N_grid=N_grid,iter=ifelse(i==N_loop, iter_final, iter), 
										max_treedepth=max_treedepth, thin=thin, visualize_progress=visualize_progress)
		b <- -fit$b
		all_a <- fit$all_a
		all_log_lik <- fit$all_log_lik
		if(length(all_a)> max_sample_store )
		{
			discard=length(all_a)- max_sample_store
			all_a=all_a[-c(1:discard)]
			all_log_lik=all_log_lik[-c(1:discard)]
		}
		print( paste("------ adaptation=", i,"---------fit measure = ", fit$fit_measure,"---------"))
	}
	return(fit)
}



g_lambda=function(a){
	a_reflected = 1 - abs(1 - a)
	a_scaled <- (a_reflected - a_lower)/(a_upper - a_lower)
	lambda <- ifelse(a_scaled <= 0, 0, ifelse(a_scaled < 1, 0.5 + 0.25*(3*(2*a_scaled - 1) - (2*a_scaled - 1)^3), 1))
	return(lambda)
}


cat("Compiling Stan optimizer...")
update_model <- stan_model("solve_tempering.stan")
cat("done.\n")
