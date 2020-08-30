library(rstan)
rstan_options(auto_write = TRUE)
source("utlitity.R")

d <- 100
n <- 40
n_test=100
set.seed(2020)
x_full <- matrix(rnorm((n+n_test)*d, 0, 1), n+n_test,d)
beta=rep(0,d)
beta[1:3]=3:1/3   #HS--2
#beta[1:10]=10:1/10

temp_common=rnorm(n+n_test,0,0.5)   # cor=  0.25/ 1.25
for( i in 1:10)
	x_full[, i]= x_full[, i]+temp_common
f_full=x_full%*%beta
print(mean(f_full))
inv_logit=function (x)
{
	1/(1 + exp(-x))
}
f_full=  inv_logit (x_full%*%beta)
y_full=rbinom(n= n+n_test, size=1, prob=f_full)

x=x_full[1:n, ]
y=y_full[1:n]

# plot(x[,1], y)
# cor(x[,1],x[,4])

tau0 <- 1/(d-1) * 2/sqrt(n) # should be a reasonable scale for tau (see Piironen&Vehtari 2017, EJS paper)
hs_data <- list(n=n, d=d, x=x, y=as.vector(y),
								n_test=n_test, x_test=x_test, y_test=as.vector(y_test),
								scale_icept=1, scale_global=tau0,
								slab_scale=5, slab_df=4)




update_model <- stan_model("solve_tempering.stan")
comp=stan_model(file="HS.stan")
probit=stan_model(file="probit.stan")
logit=stan_model(file="logit.stan")


a_length=c(0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45)  # a_min, a_max=1-a_min
eff_mat=matrix(NA, length(a_length), 4)    # tail eff per second
for(i.sim in 1:length(a_length)){
	a_temp=a_length[i.sim]
	eff_mat[i.sim,]=compare_eff(hs_data=hs_data, a_lower =a_temp, a_upper = 1- a_temp)
}





pdf("~/desktop/HS_vary_a.pdf", width=4, height=1.9)
par(mfrow=c(1,2),oma=c(1,1,1,0),mar=c(1,1,1,1) ,mgp=c(1.5,0.25,0), lwd=0.5,tck=-0.01, cex.axis=0.7, cex.lab=0.6, cex.main=0.9)
L=length(a_length)
plot(a_length, rep(mean(eff_mat[,1]), L), axes = FALSE, main="", xlab="", ylab="", lty=1, type='l', lwd=2, col="darkgreen",  xpd=TRUE,  yaxs='i', ylim=c(0,50))
abline(h=c(10, 20, 30,40), col='gray', lty=2)
points(a_length, eff_mat[,2], pch=19, col="#A25050", cex=0.6)
lines(a_length, eff_mat[,2],   col="#A25050",lwd=1.5)
axis(1, lwd=0.5, at=c(0.1,0.2,0.3,0.4), padj=-0.8)
axis(2, lwd=0.5, las=2, at=c(0,20, 40,60))
mtext(1, line=0.5, text=expression(a[min]),cex=0.7)
mtext(3, line=0.5, text="min ESS /second \n probit model" ,cex=0.7)
text(0.35,15,labels = "path sampling", cex=0.7,col="#A25050")
text(0.3,3,labels = "individual fit", cex=0.6,col='darkgreen', xpd=T)
box(lwd=0.5, bty='l')

plot(a_length, rep(mean(eff_mat[,3]), L), axes = FALSE, main="", xlab="", ylab="", lty=1, type='l', lwd=2, col="darkgreen",  xpd=TRUE,  yaxs='i', ylim=c(0,50))
abline(h=c(10, 20, 30,40), col='gray', lty=2)
points(a_length, eff_mat[,4], pch=19, col="#A25050", cex=0.6)
lines(a_length, eff_mat[,4],   col="#A25050",lwd=1.5, xpd=T)
axis(1, lwd=0.5, at=c(0.1,0.2,0.3,0.4), padj=-0.8)
axis(2, lwd=0.5, las=2, at=c(0,20, 40,60))
mtext(1, line=0.5, text=expression(a[min]),cex=0.7)
mtext(3, line=0.5, text="min ESS /second  \n logit model" ,cex=0.7)
box(lwd=0.5, bty='l')
text(0.35,20,labels = "path sampling", cex=0.7,col="#A25050")
text(0.25,5,labels = "individual fit", cex=0.7,col='darkgreen', xpd=T)
dev.off()
