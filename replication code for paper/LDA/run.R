remotes::install_github("MansMeg/posteriordb", subdir = "rpackage")
library(posteriordb) 
ghpdb <- pdb_github()
po <- posterior("prideprejustice_chapter-ldaK5", ghpdb) 
sd <- stan_data(po)
sc <- stan_code(po)
set.seed(2020) 
subsample=sample(1:sd$N, round(0.6*sd$N))
train_sample=sort(subsample[1:round(0.3*sd$N)])
test_sample=sort(subsample[-(1:round(0.3*sd$N))])
w_train=sd$w[train_sample]
w_test=sd$w[test_sample]
doc_train=sd$doc[train_sample]
doc_test=sd$doc[test_sample]
N=length(w_train)
N_test=length(w_test)
sd_train=sd
sd_train$w=w_train
sd_train$w_test=w_test
sd_train$N=N
sd_train$doc=doc_train
sd_train$doc_test=doc_test
sd_train$N_test=N_test
sd_train$K=as.integer(5) ## varies in the experiment from 3 to 15.
save(sd_train, file="sd_train_chapter.RData")



### too slow, need to subsample.
lda_fit=stan(file="lda.stan", data=c(list(a_lower=a_lower, a_upper=a_upper, b=b, K_logit=K_logit, mu_logit=mu_logit, sigma_logit=sigma_logit, K_gaussian=K_gaussian, mu_gaussian=mu_gaussian, sigma_gaussian=sigma_gaussian), sd_train), iter=10, thin =1, chains = 1)

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



log_mean_exp=function(v)
{
	max(v)+ log(mean(exp(v- max(v) )))
}



path_fit <- function(graph_filename, a_lower, a_upper, b, K_logit, 
										 mu_logit, sigma_logit, K_gaussian, 
										 mu_gaussian, sigma_gaussian, all_a, all_log_lik,
										 N_grid=100, iter=2000, chains=3,max_treedepth=5,thin=2){
	lda_fit <- stan(file="lda.stan", chains = chains,thin = thin,
									data=c(list(a_lower=a_lower, a_upper=a_upper, b=b, K_logit=K_logit, 
															mu_logit=mu_logit, sigma_logit=sigma_logit, 
															K_gaussian=K_gaussian, mu_gaussian=mu_gaussian, 
															sigma_gaussian=sigma_gaussian), sd_train),
									iter=iter,control = list(max_treedepth=max_treedepth))
	sim=extract(lda_fit)
	log_lik_sum=sim$log_lik_sum
	a <- sim$a
	N_sims <- length(a)
	all_a <- c(all_a, a)
	all_log_lik <- c(all_log_lik, log_lik_sum)
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
	log_test=sim$log_lik_test
	
	in_prior=a_reflected<a_lower
	in_prior_rate=mean(in_prior)
	if(in_prior_rate>0.8){
		in_prior_rate=in_prior_rate*0.95
		b[1]=in_prior_rate*	(log_mean_exp(log_lik_sum[in_prior])) + (1-in_prior_rate)*b[1]
	}
	test_lpd=mean(apply(log_test[in_target,], 2, log_mean_exp))
	pdf(graph_filename, height=3, width=7)
	par(mfrow=c(1,2), mar=c(3,3,2,1), mgp=c(1.5,.5,0), tck=-.01)
	hist(a_reflected, xlim=c(0, 1), xaxs="i", breaks=seq(0, 1, .01), main="Post draws of a_reflected", xlab="a_reflected", cex.main=.9, cex.axis=.9, cex.lab=.9)
	abline(v=a_lower, col="gray")
	abline(v=a_upper, col="gray")
	plot(path_lik$a, path_lik$log_z, xlim=c(0, 1), xaxs="i", bty="l", pch=20, cex=.2, col="blue", main="Path est of p(a_reflected)", xlab="a_reflected", ylab="log p(a_reflected)", cex.main=.9, cex.axis=.9, cex.lab=.9)
	lines(a_grid, curve_fit, lwd=.5, col="red")
	abline(v=a_lower, col="gray")
	abline(v=a_upper, col="gray")
	dev.off()
	
	return(list(b=b, all_a=all_a, all_log_lik=all_log_lik, fit_measure=fit_measure, path_post_a=path_lik$a, path_post_z=path_lik$log_z, test_lpd=test_lpd, curve_fit=curve_fit))
}

a_lower <- .1
a_upper <- .8
K_logit <- 20
mu_logit <-  a_lower + (a_upper - a_lower)*seq(1,2*K_logit-1,2)/(K_logit)
sigma_logit <- 2*rep((a_upper - a_lower)/K_logit, K_logit)
K_gaussian <- K_logit
mu_gaussian <- mu_logit
sigma_gaussian <- sigma_logit
#initial
b <- rep(0, 1 + K_logit + K_gaussian)

## better init.
b[1] =70000
b[23:25]=210000
b[30:31]=205000
b[35:36]=200000
b_list=list()

all_a <- NULL
all_log_lik <- NULL
N_grid=400
update_model <- stan_model("/Users/YulingYao/Documents/GitHub/tempering/code/beta_bin/solve_tempering.stan")
iter=4000
N_loop <- 10
chains=3
curve_fit_list=c()
test_lpd=c()
set.seed(1)
for (i in 1:N_loop){
	print(i)
	max_treedepth=ifelse(i<=8,5,7)
	iter=ifelse(i<=8,2000,4000)
	thin=ifelse(i<=5,4,2)
	fit <- path_fit(paste("path_", i, ".pdf", sep=""), a_lower, a_upper, b, K_logit, mu_logit, sigma_logit, K_gaussian, mu_gaussian, sigma_gaussian, all_a, all_log_lik, N_grid=N_grid, chains=chains,iter=iter, max_treedepth=max_treedepth, thin=thin)
	b <- -fit$b
	all_a <- fit$all_a
	all_log_lik <- fit$all_log_lik
	curve_fit_list[[i]]=fit$curve_fit
	b_list[[i]]=fit$b
	test_lpd[i]=fit$test_lpd
	print( paste("====================iter=", i,"=========================", fit$fit_measure,"==========="))
}

print(test_lpd)



save(fit, all_a, all_log_lik, curve_fit_list, b_list, test_lpd, a_grid, file="LDA_out.RData")


all_a_reflected = 1 - abs(1 - all_a)
final_a_re=1 - abs(1 -all_a[14250:(14250-3000)] )

pdf(graph_filename, height=3, width=7)
par(mfrow=c(1,3), mar=c(3,3,2,1), mgp=c(1.5,.5,0), tck=-.01)
hist(all_a_reflected, xlim=c(0, 1), xaxs="i",  main="Post draws of a_reflected", xlab="a_reflected", cex.main=.9, cex.axis=.9, cex.lab=.9)
hist(final_a_re, xlim=c(0, 1), xaxs="i",  main="Post draws of a_reflected", xlab="a_reflected", cex.main=.9, cex.axis=.9, cex.lab=.9)
abline(v=c(a_lower,a_upper), col="gray")




#### results 
pdf("~/Desktop/lda_tempering.pdf", width = 7.5, height = 2)
par(mfrow=c(1,5),oma=c(1.3,1,1,0),mar=c(1,1,1,1) ,mgp=c(1.5,0.25,0), lwd=0.5,tck=-0.01, cex.axis=0.7, cex.lab=0.6, cex.main=0.9)
hist(all_a_reflected, xlim=c(0, 1), xaxs="i", breaks=seq(0,1,length.out = 30),  axes = FALSE, main="", xlab="", ylab="", col="#B97C7C", border = "#A25050")
axis(1, at=c(0,0.5,1), padj=-0.5, lwd=0.5)
axis(2,  las=2, at=c(0,4000,8000), lwd=0.5)
abline(v=c(0.1,0.8), col='gray', lty=2)
mtext(3, line=0, text="dist. of temperature\n all adaptations" ,cex=0.7)
mtext(1, line=1, text="temperature a",cex=0.7)


hist(final_a_re, xlim=c(0, 1), xaxs="i", breaks=seq(0,1,length.out = 30),    axes = FALSE,xlab="", ylab="", col="#B97C7C", border = "#A25050", main="")
axis(1,  at=c(0,0.5,1), padj=-0.5, lwd=0.5)
axis(2,  las=2, at=c(0,1000,2000), lwd=0.5)
abline(v=c(0.1,0.8), col='gray', lty=2)
mtext(3, line=0, text="dist. of temperature\n last adaptations" ,cex=0.7)
mtext(1, line=1, text="temperature a",cex=0.7)


a_grid=seq(0, 1, length=N_grid)
plot(g_lambda(fit$path_post_a), fit$path_post_z, xlim=c(0, 1), xaxs="i",  ylim=c(-91110, 0),
		 pch=16, col=alpha("darkred", alpha=0.5),cex=0.9 , xlab="", ylab="",
		 axes = FALSE, main="", xpd=T)
lines(g_lambda(fit$path_post_a[order(fit$path_post_a)]), fit$path_post_z[order(fit$path_post_a)], col="#A25050",lwd=1)
mtext(1, line=1, text="temperature",cex=0.7)
mtext(3, line=0.5, text="est. log normaliz. const." ,cex=0.7)
axis(2, las=2, lwd=0.5)
box(bty='l')
axis(1,  at=c(0, 0.5, 1), padj=-0.5, lwd=0.5)
	for(i in 1:9)
		lines(g_lambda (a_grid), curve_fit_list[[i]], lwd=.5, col="darkorange")
lines(g_lambda (a_grid), curve_fit_list[[10]], lwd=.5, col="#A25050")
legend("bottomleft", col=c("darkred", "darkorange"), lwd=0.8, legend = c("last adaptaion", "first 9 adaptatons"), bty='n')
dev.off()


pdf("~/Desktop/lda2.pdf", width = 6.5, height = 2)
 par(mfrow=c(1,4),oma=c(1.3,2,1,0),mar=c(1,2,1,1) ,mgp=c(1.5,0.25,0), lwd=0.5,tck=-0.01, cex.axis=0.7, cex.lab=0.6, cex.main=0.9)
sub_idx=sample(1:length(all_a_reflected),5000)
xx=c(all_a_reflected[sub_idx])
yy=c(all_log_lik[sub_idx])
plot(xx, yy,pch=16, col=alpha("darkred", alpha=0.4),cex=0.9, xaxs="i",  xlim=c(0,1), xlab="", ylab="",axes = FALSE, main="")
mtext(1, line=1, text="a",cex=0.7)
mtext(2,  text="log\n lik",cex=0.7, line=2, las=2)
mtext(3, line=0.5, text="log likelihood" ,cex=0.7)
axis(2, las=2, at=seq(-66000,-72000, by=-2000), lwd=0.5)
abline(h=seq(-66000,-72000, by=-2000), col='grey', lty=2)
box(bty='l')
axis(1,  at=c(0, 0.5, 0.8, 1), padj=-0.5, ylim=c(0,2.4), lwd=0.5)
abline(v=c(0.1,0.8), col='grey', lty=2)
 


sub_idx=sample(1:length(all_a_reflected),1000)
xx=c(0,0.1, all_a_reflected[sub_idx])
yy=c(0,0,gradient2$u_lik[sub_idx])
plot(xx,  yy,  pch=16, col=alpha("darkred", alpha=0.6),cex=0.9, xaxs="i", xlim=c(0,1), xlab="", ylab="",axes = FALSE, main="")
lines(xx[order(xx)],  yy[order(xx)],  pch=16, col="#A25050",lwd=0.9)
mtext(1, line=1, text="a",cex=0.7)
mtext(2, line=1, text="u",cex=0.7, las=2)
mtext(3, line=0.5, text="computed gradient" ,cex=0.7)
axis(2, las=2, at=c(0,-40000,-80000,-120000), lwd=0.5)
abline(h=c(0,-40000,-80000,-120000), col='grey', lty=2)
axis(1,  at=c(0, 0.5, 0.8, 1), padj=-0.5, lwd=0.5)
box(bty='l')
abline(v=c(0.1,0.8), col='grey', lty=2)

yy=c(0,0, sapply(all_a_scaled, grad_lambda_fun)[sub_idx])
plot(xx,  yy,  pch=16, col=alpha("darkred", alpha=0.6),cex=0.9, xaxs="i", xlim=c(0,1), xlab="", ylab="",axes = FALSE, main="", ylim=c(0,2.4))
lines(seq(0,1,length.out = 40), sapply((seq(0,1,length.out = 40) - a_lower)/(a_upper - a_lower), grad_lambda_fun) ,  pch=16, col="#A25050",lwd=0.9)
mtext(1, line=1, text="a",cex=0.7)
mtext(2, line=1, text="d lambda\n da",cex=0.7, las=2)
mtext(3, line=0.5, text="derivative of lambda" ,cex=0.7)
box(bty='l')
axis(2, las=2, at=c(0,1,2), lwd=0.5)
abline(h=c(0,1,2), col='grey', lty=2)
axis(1,  at=c(0, 0.5, 0.8, 1), padj=-0.5, lwd=0.5)
abline(v=c(0.1,0.8), col='grey', lty=2)





sub_idx=sample(1:length(all_a_reflected))
xx=c(0,0.1, all_a_reflected[sub_idx])
yy=c(0,0,gradient2$u_lik[sub_idx])
plot(xx,  yy,  pch=16, col=alpha("darkred", alpha=0.6),cex=0.9, xaxs="i", xlim=c(0,1), ylim=c(-140000,-120000),
		 xlab="", ylab="",axes = FALSE, main="")
lines(xx[order(xx)],  yy[order(xx)],  pch=16, col="#A25050",lwd=0.9)
mtext(1, line=1, text="a",cex=0.7)
mtext(3, line=0, text="computed gradient \n zoom in" ,cex=0.7)
axis(2, las=2, at=c(0,-140000,-130000,-120000), lwd=0.5)
abline(h=c(0,-130000,-120000,-120000), col='grey', lty=2)
axis(1,  at=c(0, 0.5, 0.8, 1), padj=-0.5, lwd=0.5)
box(bty='l')
abline(v=c(0.1,0.8), col='grey', lty=2)

dev.off()

 

 
