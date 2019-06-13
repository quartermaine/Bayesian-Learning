NormalNonInfoPrior<-function(NDraws,Data){
  
  #####################################################################################
  # PURPOSE: 	Generates samples from the joint posterior distribution of the parameters in the
  #
  #		MODEL: x1,....xn iid Normal(mu,sigma^2) model. 
  #		PRIOR: p(mu,sigma^2) propto 1/sigma^2				
  #
  #
  # INPUT:	NDraw: 		(Scalar)	The number of posterior draws
  #		Data:		(n-by-1)	Data vector of n observations
  #				
  # OUTPUT:	PostDraws:	(NDraws-by-2)	Matrix of posterior draws. 
  #						First column holds draws of mu		
  #						Second column holds draws of sigma^2
  #
  # AUTHOR:	Mattias Villani, Sveriges Riksbank and Stockholm University.
  #		E-mail: mattias.villani@riksbank.se.
  #
  # DATE:		2005-04-13
  #
  ######################################################################################
  
  n<-length(Data)
  Datamean<-mean(Data)
  tau2<-sum(Data-3.5)^2/n
  PostDraws=matrix(0,NDraws,2)
  PostDraws[,2]<-(n*tau2)/rchisq(NDraws,n)
  PostDraws[,1]<-Datamean+rnorm(NDraws,0,1)*sqrt(PostDraws[,2]/n)
  
  return(PostDraws)
}

PostDraws<-NormalNonInfoPrior(10000,log(data)) # Generating 1000 draws from the joint posterior density of mu and sigma^2
hist(PostDraws[,2],prob=T,border="blue")# Plotting the histogram of sigma-draws
lines(density(PostDraws[,2], adjust=2), col="red", lwd=2)








