##################################
### Question 1                   #
##################################

# a 

# import dataset
templink<-read.delim("TempLinkoping.txt")
# import library
library(mvtnorm)

prior_m0=c(-10,100,-100)
prior_v=4
prior_sigma_sq0=1
prior_omega=0.5*diag(3)

# inverse function
NormalNonInfoPrior<-function(v0,sigma_sq0){
  
  
  PostDraws<-(v0*sigma_sq0)/rchisq(1,v0)
  
  return(PostDraws)
}


# simulation function
conjugate_prior_sim<-function(N,v0,sigma_sq0,mu0,omega0){
  temp_pred<-matrix(0,N,dim(templink)[1])
  betas<-matrix(0,N,length(mu0))
  for(i in 1:N){
    sigmaSq_prior<-NormalNonInfoPrior(v0,sigma_sq0)
    betas[i,]<-rmvnorm(1,mean=mu0,sigma=sigmaSq_prior*solve(omega0))
    temp_pred[i,]<-betas[i,1]+betas[i,2]*templink[,1]+betas[i,3]*templink[,1]^2+rnorm(1,0,sigmaSq_prior)
  }
  return(list("temp_pred"=temp_pred,"beta_preds"=betas))
}

# 


# sample of betas
prior_sample<-conjugate_prior_sim(20,prior_v,prior_sigma_sq0,prior_m0,prior_omega)
# plot the 
plot(templink[,1],templink[,2],panel.first=grid(25,25),ylab = "temp",xlab="time",
     main="Plot of the data and curves",pch="o")
for (j in 1:dim(prior_sample$temp_pred)[1]) lines(templink[,1], prior_sample$temp_pred[j,],col=sample(colors(), 1),lwd=2)


# b


posterior_sim<-function(nsample,X,y,v0,mu0,sigma_sq0,omega0){
  betas_post=matrix(0,nsample,length(mu0))
  sigma_post=rep(0,nsample)
  beta_hat<-solve(t(X)%*%X)%*%t(X)%*%y
  
  for (i in 1:nsample){
  
    
    #
    mu_n<-solve( t(X)%*%X+omega0 )%*%(t(X)%*%X%*%beta_hat+omega0%*%mu0 )
    #
    omega_n<-t(X)%*%X+omega0
    #
    v_n<-v0+dim(templink)[1]
    #
    sigma_m<-( v0*sigma_sq0+(t(y)%*%y+t(mu0)%*%omega0%*%mu0-t(mu_n)%*%omega_n%*%mu_n) )/v_n
    #
    sim_sigma<-NormalNonInfoPrior(v_n,sigma_m)
    
    sigma_post[i]<-sim_sigma
    
    sim_betas<-rmvnorm(1,mean=mu_n,sigma=sim_sigma[1]*solve(omega_n))
    
    betas_post[i,]<-sim_betas
  }
  
  return(list("sigma_post"=sigma_post,"betas_post"=betas_post))
  
}

Xt<-cbind(1,templink[,1],templink[,1]^2)
yt=templink[,2]

post_sample<-posterior_sim(100,Xt,yt,prior_v,prior_m0,prior_sigma_sq0,prior_omega)

betas_pr<-post_sample$betas_post



Dt<-matrix(0,dim(betas_pr)[1],length(templink[,1]) )

for (i in 1:dim(betas_pr)[1]){

  Dt[i,]<-betas_pr[i,1]+betas_pr[i,2]*templink[,1]+betas_pr[i,3]*(templink[,1]^2)

}

plot(templink[,1],templink[,2],panel.first=grid(25,25),ylab = "temp",xlab="time",
     main="Plot of the data and curves",pch="o")
for (j in 1:dim(Dt)[1]) lines(templink[,1], Dt[j,],col=sample(colors(), 1),lwd=2)




posterior_CI<-apply(Dt,2,function(y) quantile(y, probs = c(0.025, 0.975)))

plot(templink[,1],templink[,2],panel.first=grid(25,25),ylab = "temp",xlab="time",
     main="Plot of the data and curves",pch="o")
lines(templink[,1],posterior_CI[1,],col="blue",lwd=2)
lines(templink[,1],posterior_CI[2,],col="orange",lwd=2)


# c
#calcutate derivative of time function and set to 0
timehat=betas_pr[,2]/(-2*betas_pr[,3])
timehat


# d



##################################
### Question 2                   #
##################################
# Loading data from file
df<-read.table("WomenWork.dat",header=T)

# a

library(glmnet)

glmModel <-glm(Work ~ 0 + ., data = df, family = binomial)


cat("The model coefficients are :\n")
glmModel$coefficients
cat("The models missclassification error is : ",mean(ifelse(glmModel$fitted.values>0.5,1,0)==df$Work))

summary(glmModel)

# b
chooseCov <- c(2:9)
tau <- 10; # Prior scaling factor such that Prior Covariance = (tau^2)*I

# Data from the read.table function is a data frame. Let's convert y and X to vector and matrix.
y <- as.vector(df[,which(names(df)=="Work")]); 
X <- as.matrix(df[,-which(names(df)=="Work")]);
covNames <- names(df)[-1];
#X <- X[,chooseCov]; # Here we pick out the chosen covariates.
nPara <- dim(X)[2];

# Setting up the prior
mu <- as.vector(rep(0,nPara)) # Prior mean vector
Sigma <- tau^2*diag(nPara); # Prior covariance matrix

# Defining the functions that returns the log posterior (Logistic and Probit models). Note that the first input argument of

# this function must be the one that we optimize on, i.e. the regression coefficients.

LogPostLogistic <- function(betaVect,y,X,mu,Sigma){
  
  nPara <- length(betaVect);
  linPred <- X%*%betaVect;
  
  # evaluating the log-likelihood                                    
  logLik <- sum( linPred*y -log(1 + exp(linPred)));
  if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
  
  # evaluating the prior
  logPrior <- dmvnorm(betaVect, matrix(0,nPara,1), Sigma, log=TRUE);
  
  # add the log prior and log-likelihood together to get log posterior
  return(logLik + logPrior)
}

# Different starting values. Ideally, any random starting value gives you the same optimum (i.e. optimum is unique)
initVal <- as.vector(rep(0,dim(X)[2]));
# Or a random starting vector: as.vector(rnorm(dim(X)[2]))
# Set as OLS estimate: as.vector(solve(crossprod(X,X))%*%t(X)%*%y); # Initial values by OLS

logPost = LogPostLogistic;

OptimResults<-optim(initVal,logPost,gr=NULL,y,X,mu,Sigma,method=c("BFGS"),
                    control=list(fnscale=-1),hessian=TRUE)

# Printing the results to the screen
postMode <- OptimResults$par
postCov <- -solve(OptimResults$hessian) # Posterior covariance matrix is -inv(Hessian)
names(postMode) <- covNames # Naming the coefficient by covariates
approxPostStd <- sqrt(diag(postCov)) # Computing approximate standard deviations.
names(approxPostStd) <- covNames # Naming the coefficient by covariates
print('The posterior mode is:')
print(postMode)
print('The approximate posterior standard deviation is:')
print(approxPostStd)

# Plotting some of the marginal posteriors
par(mfrow = c(2,2))
for (k in 1:4){
  betaGrid <- seq(0, postMode[k] + 4*approxPostStd[k], length = 1000)
  plot(betaGrid, dnorm(x = betaGrid, mean = postMode[k], sd = approxPostStd[k]), 
       type = "l", lwd = 2, main = names(postMode)[k], 
       ylab = '', xlab = expression(beta))
}

credInter=quantile(df$NSmallChild, c(0.025, 0.975))

credInter


# c




pred_dist_logreg<-function(nSample,mode,Cov,query){
  
  sample_postMode<-rmvnorm(nSample,mean=mode,sigma=Cov)

  y_preds<-query%*%t(sample_postMode)
  
  
  y_logpreds<-ifelse(exp(y_preds)/(1+exp(y_preds))>runif(nSample,0,1),1,0)
  
  return(y_logpreds)
  
}
  
query<-c(1,10,8,10,1,40,1,1)

pr<-pred_dist_logreg(10000,postMode,postCov,query)

plot(density(pr),main="Density of the predictive \ndistirbution for the query",
     panel.first=grid(25,25),col="mediumvioletred",lwd=2)
polygon(density(pr), density = 18, angle = 45,col="lightslateblue")

hist(pr,main="Histogram of the predictive \ndistirbution for the query",
     panel.first=grid(25,25),border="dodgerblue",col="gray",lwd=2)

