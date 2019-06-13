
# Exercise 1 -Normal model, mixture of normal model with semi-conj --------

# a)

#  i) 


rainfall<-read.table('rainfall.dat',head=F)

hist(rainfall[,1])
plot(density(rainfall[,1]))


# inverse chi-sqare function
rInvChi2<-function(v0,sigma_sq0){
  # returns one sample from Inv-chi square for given (df,sigma)
  inverse.chi<-(v0*sigma_sq0)/rchisq(1,v0)
  return(inverse.chi)
}


GibbsSampler<-function(nIter,mu0,tau0_squared,Sigma0_squared,v0,y){
  n=length(y)
  muPost<-rep(0,nIter)
  SigmaPost<-rep(0,nIter)
  #ySample<-rep(0,nIter)
  
  numerator= ( n/Sigma0_squared )+( 1/tau0_squared)
  w=(n/Sigma0_squared)/numerator
  mu_n=w*mean(y)+(1-w)*mu0
  tau_n_squared=1/numerator
  
  for (i in 1:nIter){
  
    muPost[i]<-rnorm(n = 1,mean = mu_n,sd = tau_n_squared)
    SigmaPost[i]<-rInvChi2(v0+n, (v0*Sigma0_squared+sum( (y-muPost[i])^2)) / (v0+n) )
    #ySample[i]<-rnorm(n=1,muPost[i],SigmaPost[i])
  }
  
  return(list('muPost'=muPost,'SigmaPost'=SigmaPost)) #,'samplePost'=ySample))
  
}

mu0=mean(rainfall[,1])
tau0=1
v0=1
Sigma2=var(rainfall[,1])

gibbs_sample=GibbsSampler(1000,mu0,tau0,Sigma2,v0,rainfall[,1] )


#   ii)
plot(gibbs_sample$muPost,type = "l",main = expression(paste("Trace plot for ", mu)),
     col="lightblue",pch="o")
plot(gibbs_sample$SigmaPost,type = 'l',main = expression(paste("Trace plot for ",sigma)),
     col="lightblue")


plot(dplyr::cummean(gibbs_sample$muPost),type="l",col="blue")
abline(h=mean(rainfall[,1]),col="red")

plot(dplyr::cummean(gibbs_sample$SigmaPost),type="l",col="blue")
abline(h=var(rainfall[,1]),col="red")



# b)
x<-as.matrix(rainfall[,1])

# Model options
nComp <- 2    # Number of mixture components

# Prior options
alpha <- 10*rep(1,nComp) # Dirichlet(alpha)
muPrior <- rep(0,nComp) # Prior mean of mu
tau2Prior <- rep(10,nComp) # Prior std of mu
sigma2_0 <- rep(var(x),nComp) # s20 (best guess of sigma2)
nu0 <- rep(4,nComp) # degrees of freedom for prior on sigma2

# MCMC options
nIter <- 1000 # Number of Gibbs sampling draws



# Plotting options
plotFit <- TRUE
lineColors <- c("blue", "green", "magenta", 'yellow')
sleepTime <- 0.1 # Adding sleep time between iterations for plotting
################   END USER INPUT ###############

###### Defining a function that simulates from the 
rScaledInvChi2 <- function(n, df, scale){
  return((df*scale)/rchisq(n,df=df))
}

####### Defining a function that simulates from a Dirichlet distribution
rDirichlet <- function(param){
  nCat <- length(param)
  piDraws <- matrix(NA,nCat,1)
  for (j in 1:nCat){
    piDraws[j] <- rgamma(1,param[j],1)
  }
  piDraws = piDraws/sum(piDraws) # Diving every column of piDraws by the sum of the elements in that column.
  return(piDraws)
}

# Simple function that converts between two different representations of the mixture allocation
S2alloc <- function(S){
  n <- dim(S)[1]
  alloc <- rep(0,n)
  for (i in 1:n){
    alloc[i] <- which(S[i,] == 1)
  }
  return(alloc)
}

# Initial value for the MCMC
nObs <- length(x)
S <- t(rmultinom(nObs, size = 1 , prob = rep(1/nComp,nComp))) # nObs-by-nComp matrix with component allocations.
mu <- quantile(x, probs = seq(0,1,length = nComp))
sigma2 <- rep(var(x),nComp)
probObsInComp <- rep(NA, nComp)

# Setting up the plot
xGrid <- seq(min(x)-1*apply(x,2,sd),max(x)+1*apply(x,2,sd),length = 100)
xGridMin <- min(xGrid)
xGridMax <- max(xGrid)
mixDensMean <- rep(0,length(xGrid))
effIterCount <- 0
ylim <- c(0,2*max(hist(x)$density))


for (k in 1:nIter){
  message(paste('Iteration number:',k))
  alloc <- S2alloc(S) # Just a function that converts between different representations of the group allocations
  nAlloc <- colSums(S)
  print(nAlloc)
  # Update components probabilities
  pi <- rDirichlet(alpha + nAlloc)
  
  # Update mu's
  for (j in 1:nComp){
    precPrior <- 1/tau2Prior[j]
    precData <- nAlloc[j]/sigma2[j]
    precPost <- precPrior + precData
    wPrior <- precPrior/precPost
    muPost <- wPrior*muPrior + (1-wPrior)*mean(x[alloc == j])
    tau2Post <- 1/precPost
    mu[j] <- rnorm(1, mean = muPost, sd = sqrt(tau2Post))
  }
  
  # Update sigma2's
  for (j in 1:nComp){
    sigma2[j] <- rScaledInvChi2(1, df = nu0[j] + nAlloc[j], scale = (nu0[j]*sigma2_0[j] + sum((x[alloc == j] - mu[j])^2))/(nu0[j] + nAlloc[j]))
  }
  
  # Update allocation
  for (i in 1:nObs){
    for (j in 1:nComp){
      probObsInComp[j] <- pi[j]*dnorm(x[i], mean = mu[j], sd = sqrt(sigma2[j]))
    }
    S[i,] <- t(rmultinom(1, size = 1 , prob = probObsInComp/sum(probObsInComp)))
  }
  
  # Printing the fitted density against data histogram
  if (plotFit && (k%%1 ==0)){
    effIterCount <- effIterCount + 1
    hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), main = paste("Iteration number",k), ylim = ylim)
    mixDens <- rep(0,length(xGrid))
    components <- c()
    for (j in 1:nComp){
      compDens <- dnorm(xGrid,mu[j],sd = sqrt(sigma2[j]))
      mixDens <- mixDens + pi[j]*compDens
      lines(xGrid, compDens, type = "l", lwd = 2, col = lineColors[j])
      components[j] <- paste("Component ",j)
    }
    mixDensMean <- ((effIterCount-1)*mixDensMean + mixDens)/effIterCount

    lines(xGrid, mixDens, type = "l", lty = 2, lwd = 3, col = 'red')
    legend("topright", box.lty = 1, legend = c("Data histogram",components, 'Mixture'),
           col = c("black",lineColors[1:nComp], 'red'), lwd = 2)
    Sys.sleep(sleepTime)
  }
  
}

hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), main = "Final fitted density")
lines(xGrid, mixDensMean, type = "l", lwd = 2, lty = 4, col = "red")
lines(xGrid, dnorm(xGrid, mean = mean(x), sd = apply(x,2,sd)), type = "l", lwd = 2, col = "blue")
legend("topright", box.lty = 1, legend = c("Data histogram","Mixture density","Normal density"), col=c("black","red","blue"), lwd = 2)

# c)



hist(rainfall[,1], breaks = 20, freq = FALSE, main = " ",xlim = c(xGridMin,xGridMax))
lines(xGrid, mixDensMean, type = "l", lwd = 2, lty = 4, col = "red")
lines(seq(-100,400,0.01), dnorm(seq(-100,400,0.01), mean = mean(gibbs_sample$muPost), 
                   sd = sqrt(mean(gibbs_sample$SigmaPost))), type = "l", lwd = 2, col = "blue")
legend("topright", box.lty = 1, legend = c("Data histogram","Mixture density","Normal density"), col=c("black","red","blue"), lwd = 2)



# Exercise 2-Metropolis Random walk for Poisson regression ----------------


# a)

ebay_data=read.table('eBayNumberOfBidderData.dat',head=T)
dim(ebay_data)

glmModel <-glm(nBids ~ 0 + ., data = ebay_data, family = poisson(link = "log"))

# create a boolen for the significant coefficients a=0.05
coeffs_toselect<-summary(glmModel)$coefficients[-1,4]<0.05
# select sig. variables
sign_coeffs<- names(coeffs_toselect)[coeffs_toselect == TRUE] 

cat('The significant coefficients are :',sign_coeffs)

# b)

library("mvtnorm")
# Data from the read.table function is a data frame. Let's convert y and X to vector and matrix.
y <- as.vector(ebay_data[,which(names(ebay_data)=="nBids")]);length(y)
X <- as.matrix(ebay_data[,-which(names(ebay_data)=='nBids')]);dim(x)
covNames <- names(X)
nPara <- dim(X)[2];

tau <- 10;  
Sigma <- tau^2*solve(t(X)%*%X);



LogPostPoison <- function(betaVect,y,X,Sigma){
  
  nPara <- length(betaVect);
  linPred <- X%*%betaVect;
  
  # evaluating the log-likelihood                                    
  logLik <- sum( linPred*y -exp(linPred) );
  if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
  
  # evaluating the prior
  logPrior <- dmvnorm(betaVect, matrix(0,nPara,1), Sigma, log=TRUE);

  # add the log prior and log-likelihood together to get log posterior
  return(logLik + logPrior)
}

initVal <- as.vector(rep(0,dim(X)[2])); 


OptimResults<-optim(initVal,LogPostPoison,gr=NULL,y,X,Sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

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


library(MASS)
require(RColorBrewer)
pal<-brewer.pal(9, "Paired")

samp<-mvrnorm(1000,postMode,postCov)

par(mfrow=c(2,4))
for(i in 2:9){
  hist(samp[,i],freq=F,
       main=paste(colnames(X)[i]),
       xlab='param',ylab='value',lwd=2)
  lines(density(samp[,i]),col =sample(pal,1))
}



#  c)

library(mvtnorm)

RWMSampler<-function(logPostFunc,n.sim,cx,betaVect,SigmaP,...){
  
  initSample=mvrnorm(n = 1,betaVect,SigmaP)
 
  mat_betas<-matrix(nrow=n.sim,ncol=length(initSample))

  mat_betas[1,]<-initSample
  
  for(i in 2:n.sim){
  
    propSample<-as.vector( mvrnorm(n=1,mat_betas[i-1,],cx*SigmaP) )
  
    r=exp( logPostFunc(propSample,...)-logPostFunc(as.vector(mat_betas[i-1,]),...) )
    
    if( runif(1)<min(1,r) ){
      mat_betas[i,]<-propSample
    }else{
      mat_betas[i,]<-mat_betas[i-1,]
    }
    
  }
  return(mat_betas)
  
}


betaVect=as.vector(rep(0,dim(X)[2]));y=y;X=X;SigmaP=postCov;n.sim=5000;cx=0.5
res<-RWMSampler(LogPostPoison,n.sim,cx,betaVect,SigmaP,y,X,Sigma)

# require(RColorBrewer)
# pal<-brewer.pal(9, "Paired")

par(mfrow=c(3,3))
for(i in 1:9){
  plot(res[,i],type="l",
       main=paste('Covergence of ',colnames(X)[i]),
       col =sample(pal,1),
       xlab='index',ylab='value',lwd=2)
}


par(mfrow=c(2,4))
for(i in 2:9){
  hist(res[,i],freq=F,breaks=30,
       main=paste('Covergence of ',colnames(X)[i]),
       xlab='index',ylab='value',lwd=2)
  lines(density(res[,i]),col =sample(pal,1),lwd=2)
}






