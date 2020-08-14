#' Extract path gradient from the augmented stan sampling results
#'
#'
#' @export
#' @param output_file The name of the output file that saves local variables, which
#' otherwise would be saved in stan. The default is \code{stan-output_path.txt}.
#'
#' @details Currently Stan cannot save local variables, while the path sampling gradient
#' often depends on partial sum of the target lp. This is an convenient function
#' that prints the need value during each leapfrog step and collect them afterwards.
#'
#' @return The pointwise log density difference between two model blocks.
#'


extract_lp=function(output_file='stan-output_path.txt'){
	stan_out=readLines(output_file)
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


# internal ----------------------------------------------------------------
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
