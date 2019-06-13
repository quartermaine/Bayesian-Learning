
# Exercise 1-Time series models in Stan -----------------------------------

# a) ------------------------------------------------------------------

set.seed(123456)

require(rstan)
require(RColorBrewer)
pal<-brewer.pal(8,"Set2")

# set prior values
mu=10
Sigma2=2
phi_grid=seq(-1,1,by=0.05);phi_grid=round(phi_grid,2)
Taf=200

# AR(1)-process simulation 
ARSim<-function(mu,Sigma2,phi,nIter){
  xt=double(Taf)
  xt[1]=mu
  for(i in 2:nIter){
  xt[i]=mu+phi*(xt[i-1]-mu)+rnorm(1,0,sqrt(Sigma2))
  }
  L=list(xt)
  names(L)[1] <- paste("AR(1)-phi=",phi,sep="")
  return(L)
}



results<-sapply(phi_grid, function(phi) ARSim(mu,Sigma2=Sigma2,phi=phi,nIter=Taf))
nPhi<-length(results)


indx<-which(phi_grid%in%c(-1.00,-0.95,-0.90,-0.05,0.00,0.05,0.90,0.95,1.00) )
indx

# plot AR(1)-process
par(mfrow=c(round(length(indx)/3),3))
for(i in indx){
  plot(results[[i]],main=names(results)[i],type="l",
       col=sample(pal,1),panel.first=grid(25,25),lwd=2)
}

# plot AR(1)-process as Time series
par(mfrow=c(round(length(indx)/3),3))
for(i in indx){
  plot.ts(results[[i]],main=names(results)[i],type="l",
          col=sample(pal,1),panel.first=grid(25,25),lwd=2)
}


# b) -------------------------------------------------------------

MCMC_model='
data{
int<lower=0> N;
vector[N] y;
}
parameters{
real mu;
real<lower=0> sigma;
real phi;
}
model{
for(n in 2:N)
y[n] ~ normal(mu + phi*(y[n-1]-mu), sqrt(sigma));

}'

ar1=ARSim(mu,Sigma2,0.3,Taf)
ar2=ARSim(mu,Sigma2,0.95,Taf)

#N1=length(ar1)

data1<-list(N = Taf, 
           y=ar1[[1]] )

data2<-list(N=Taf,
               y=ar2[[1]])

fit.mod1<-stan(model_code = MCMC_model, # Stan model
              data = data1,              # named list of data
              chains = 1,              # number of Markov chains
              warmup = 1000,           # number of warmup iterations per chain
              iter = 2000,             # total number of iterations per chain
              cores = 2,               # number of cores (could use one per chain)
              refresh = 0 )


fit.mod2<-stan(model_code = MCMC_model, # Stan model
               data = data2,              # named list of data
               chains = 1,              # number of Markov chains
               warmup = 1000,           # number of warmup iterations per chain
               iter = 2000,             # total number of iterations per chain
               cores = 2,               # number of cores (could use one per chain)
               refresh = 0 )




model_stats1<-extract(fit.mod1)

model_stats2<-extract(fit.mod2)

# i)
## 95% credible intervals

q1_mu=quantile(model_stats1$mu,probs=c(0.025,0.975)) # for mu
q1_sigma=quantile(model_stats1$sigma,probs=c(0.025,0.975)) # for sigma
q1_phi=quantile(model_stats1$phi,probs=c(0.025,0.975)) # for phi


q2_mu=quantile(model_stats2$mu,probs=c(0.025,0.975)) # for mu
q2_sigma=quantile(model_stats2$sigma,probs=c(0.025,0.975)) # for sigma
q2_phi=quantile(model_stats2$phi,probs=c(0.025,0.975)) # for phi

library(knitr)
# 
dataFrame1=data.frame( rbind(q1_mu,q1_sigma,q1_phi) )
colnames(dataFrame1)<-c("2.5%","97.5%")
rownames(dataFrame1)<-NULL
rownames(dataFrame1)<-c("mu","sigma","phi")

#
kable(dataFrame1,caption = "Table for phi=0.3")
#
dataFrame2=data.frame( rbind(q2_mu,q2_sigma,q2_phi) )
colnames(dataFrame2)<-c("2.5%","97.5%")
rownames(dataFrame2)<-NULL
rownames(dataFrame2)<-c("mu","sigma","phi")
#
kable(dataFrame2,caption = "Table for phi=0.95")


# ii)

# convergence plots 

par(mfrow=c(1,2))
plot(model_stats1$mu,col=sample(pal,1),
     main=expression(paste("Convergence plot for ",mu," with ",phi,"=0.3")),
     xlab=expression(mu),type="l")
plot(model_stats1$sigma,col=sample(pal,1),
     main=expression(paste("Convergence plot for ",sigma," with ",phi,"=0.3")),
     xlab=expression(mu),type="l")


par(mfrow=c(1,2))
plot(model_stats2$mu,col=sample(pal,1),
     main=expression(paste("Convergence plot for ",mu," with ",phi,"=0.95")),
     xlab=expression(mu),type="l")
plot(model_stats2$sigma,col=sample(pal,1),
     main=expression(paste("Convergence plot for ",sigma," with ",phi,"=0.95")),
     xlab=expression(mu),type="l")



# posterior plots
par(mfrow=c(1,2))
plot(model_stats1$mu,model_stats1$phi,col=sample(pal,1),
     main=expression(paste("Plot joint posterior ",phi,"=0.3")),
     xlab=expression(mu),ylab=expression(phi))
plot(model_stats2$mu,model_stats2$phi,col=sample(pal,1),
     main=expression(paste("Plot joint posterior ",phi,"=0.95")),
     xlab=expression(mu),ylab=expression(phi))




# c) ----------------------------------------------------------------------

campy<-read.table("campy.dat",head=T)


poisson_model='
data{
int<lower=0> N;
int y[N];
}
parameters{
real mu ;
real<lower=0> sigma;
real<lower=-1, upper=1> phi;
vector[N] xt;
}
model{
mu ~ normal(0,100); // Normal with mean 0, st.dev. 100
sigma ~ scaled_inv_chi_square(1,2); // Scaled-inv-chi2 with nu 1, sigma 2
for(n in 2:N){
xt[n]~normal(mu + phi*(xt[n-1]-mu), sqrt(sigma));
y[n]~poisson(exp(xt[n]));
}
}'


dataP=list(N=length(campy$c),y=campy$c)

fit.poisson=stan(model_code = poisson_model ,
                 data = dataP,              # named list of data
                 chains = 1,              # number of Markov chains
                 warmup = 1000,           # number of warmup iterations per chain
                 iter = 2000,             # total number of iterations per chain
                 cores = 2,               # number of cores (could use one per chain)
                 refresh = 0 )

p<-extract(fit.poisson)
theta_post=exp(p$xt)
theta_post_mean=apply(theta_post,2,mean)
theta_post_mean

CI<-apply(theta_post,2,function(x) quantile(x,prob=c(0.025,0.975)))


plot(campy$c,main="Plot of post mean and 95% CI",
     pch=20,ylab="campy")
lines(theta_post_mean,col="red",lwd=2)
lines(CI[2,],col="green",lwd=2)
lines(CI[1,],col="blue",lwd=2)
legend("topleft",legend=c("95% CI","Post Mean","5% CI"),col=c("green","red","blue"),
       lty=1,lwd=2)

#library(shinystan)
#launch_shinystan(fit.poisson)


# d) ----------------------------------------------------------------------



poisson_model2='
data{
int<lower=0> N;
int y[N];
}
parameters{
real mu ;
real<lower=0> sigma;
real<lower=-1, upper=1> phi;
vector[N] xt;
}
model{
mu ~ normal(0,100); // Normal with mean 0, st.dev. 100
sigma ~ scaled_inv_chi_square(100,0.5); // Scaled-inv-chi2 with nu 100, sigma 0.5
for(n in 2:N){
xt[n]~normal(mu + phi*(xt[n-1]-mu), sqrt(sigma));
y[n]~poisson(exp(xt[n]));
}
}'



fit.poisson2=stan(model_code = poisson_model2 ,
                 data = dataP,              # named list of data
                 chains = 1,              # number of Markov chains
                 warmup = 1000,           # number of warmup iterations per chain
                 iter = 2000,             # total number of iterations per chain
                 cores = 2,               # number of cores (could use one per chain)
                 refresh = 0 )

p2<-extract(fit.poisson2)
theta_post2=exp(p2$xt)
theta_post_mean2=apply(theta_post2,2,mean)
theta_post_mean2


CI2<-apply(theta_post2,2,function(x) quantile(x,prob=c(0.025,0.975)))


plot(campy$c,main="Plot of post mean and 95% CI",
     pch=20,ylab="campy")
lines(theta_post_mean2,col="red",lwd=2)
lines(CI2[2,],col="green",lwd=2)
lines(CI2[1,],col="blue",lwd=2)
legend("topleft",legend=c("95% CI","Post Mean","5% CI"),col=c("green","red","blue"),
       lty=1,lwd=2)









