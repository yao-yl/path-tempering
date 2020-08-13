#' Augmented a Stan file with alternative model blocks
#'
#' Convenience function for decroating a Stan file with alternative model 
#' blocks and contruct a geometric bridge in between.
#'
#' @export
#' @param stan_file The path to the \code{Stan} program to use (\pkg{rstan} 
#' package). It should be a character string file name.
#' 
#' @details The function constructs a geometric bridge in between two posterior 
#' densities that is propotional to the joint posterior in two model blocks. 
#' 
#' @return The new \code{Stan} file path.
#' 
#' @examples
#' \dontrun{
#' file_new=code_tempeture_augment(stan_file="example/cauchy_mixture.stan")
#' }



code_tempeture_augmented=function(stan_file=NULL){
	stan_code=readLines(stan_file)
	if(sum(grep(" a;| a\\[| a ", stan_code))+ sum(grep(" lambda;| lambda\\[| lambda ", stan_code))>0 ) 
		stop("Please reserve parameter name a and lambda")
	# locate data block
	line_data_start=which(grepl("data",stan_code)&!grepl("transformed",stan_code))
	if(length(line_data_start)>1)
		stop("must contain one data block!")
	if(length(line_data_start)==1){
		line_data_start=block_start(stan_code)[min (which(block_start(stan_code) >= line_data_start))]
		stan_code_new= add_line(string_vec=stan_code, add_line_id=line_data_start , aug_string_data)
	}else{
		stan_code_new= c( "data{", aug_string_data, "}", stan_code)
	}
	# locate transformed data block
	line_transformed_start=which(grepl("transformed",stan_code_new) &grepl("data",stan_code_new) )
	if(length(line_transformed_start)>1){stop("wrong format of transformed data")
	} else if(length(line_transformed_start)==1) {
		line_transformed_start=block_start(stan_code_new) [min (which(block_start(stan_code_new) >= line_transformed_start))]
		stan_code_new= add_line(string_vec=stan_code_new, add_line_id=line_transformed_start , aug_string_trans_data_pre)
		line_transformed_end=block_end(stan_code_new)[min (which(block_end(stan_code_new) >= line_transformed_start))]
		stan_code_new= add_line(string_vec=stan_code_new, add_line_id=line_transformed_end-1,  add_on=aug_string_trans_data_suffix)
	}else{
		line_data_start=which(grepl("data",stan_code_new)&!grepl("transformed",stan_code_new))
		line_transformed_start=block_end(stan_code_new)[min (which(block_end(stan_code_new) >= line_data_start))]
		stan_code_new= add_line(string_vec=stan_code_new, add_line_id=line_transformed_start , c("transformed data {", aug_string_trans_data_pre,aug_string_trans_data_suffix, "}"))
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
		stan_code_new= add_line(string_vec=stan_code_new, add_line_id=line_transformed_start , aug_string_trans_parameter_pre)
		line_transformed_end=block_end(stan_code_new)[min (which(block_end(stan_code_new) >= line_transformed_start))]
		stan_code_new= add_line(string_vec=stan_code_new, add_line_id=line_transformed_end-1,  add_on=aug_string_trans_parameter_suffix)
	}else{
		line_transformed_start=block_end(stan_code_new)[min (which(block_end(stan_code_new) >= line_parameter_start))]
		stan_code_new= add_line(string_vec=stan_code_new, add_line_id=line_transformed_start , c("transformed parameters {", aug_string_trans_parameter_pre,aug_string_trans_parameter_suffix,	 "}"))			
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
													add_on=	c(aug_string_model_define, alt_model))
	line_model_end=block_end(stan_code_new)[min (which(block_end(stan_code_new) >= line_model_start))]
	stan_code_new= add_line(string_vec=stan_code_new, add_line_id=line_model_end-1,  add_on=aug_string_model)
	
	
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

# internal ----------------------------------------------------------------
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

aug_string_data=c(
	"\t//augmented:",
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

aug_string_trans_data_pre=c("\t//augmented:",
														"\treal b_linear;",
														"\tvector[K_logit] b_logit;",
														"\tvector[K_gaussian] b_gaussian;",
														"\t//original:"
)


aug_string_trans_data_suffix=c(
	"\tb_linear = b[1];",
	"\tb_logit = segment(b, 2, K_logit);",
	"\tb_gaussian = segment(b, 2 + K_logit, K_gaussian);"
)


aug_string_trans_parameter_pre=c("\t//augmented transformed:",
																 "\treal a_reflected;",    
																 "\treal a_scaled;",
																 "\treal lambda;",
																 "\tvector[K_logit] kernel_logit;",
																 "\tvector[K_gaussian] kernel_gaussian;",
																 "\t//original:")


aug_string_trans_parameter_suffix=c(	
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

aug_string_model_define=c(
	"\treal log_q;//lp of the alternative model",
	"\treal base_jacobian;//lp of the jacobian",
	"\treal log_p;//lp of the main model",
	"\tbase_jacobian = target();"
)



aug_string_model=c(
	"\tlog_p = target()-base_jacobian;",
	"\tprint(\" log_p = \", log_p,\" log_q = \", log_q, \" a = \", a);",
	"\ttarget += (lambda-1)*log_p;",
	"\ttarget += (1-lambda)*log_q;",
	"\ttarget += lambda*b_linear;",
	"\ttarget += dot_product(kernel_logit, b_logit) + dot_product(kernel_gaussian, b_gaussian);"
)

