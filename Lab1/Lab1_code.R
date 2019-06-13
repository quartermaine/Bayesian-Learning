BetaPrior <- function(a,b,s){
  xGrid <- seq(0.001, 0.999, by=0.001)
  n=20
  posterior = dbeta(xGrid, a+s, b+(n-s))
  maxDensity <- max(posterior) # Use to make the y-axis high enough
  plot(xGrid, posterior, type = 'l', lwd = 3, col = "blue", xlim <- c(0,1), ylim <- c(0, maxDensity), xlab = "theta", 
       ylab = 'Density', main = paste0("Beta,(" ,a+s,",",b+(n-s), ")density",sep=" "))
}

BetaPrior(2,2,14)

##########################################
## Exercise 1
##########################################

## a
set.seed(123456)

mean_beta=double(length(1:10000))
sd_beta=double(length(1:10000))
# Vector of beta parameters : a,b
ab <- c(16, 8)

# Simulate 1000 draws from the beta posterior: p_sim
for(i in 1:10000){
  beta_sim <- rbeta(i, ab[1], ab[2])
  mean_beta[i]<-mean(beta_sim)
  sd_beta[i]<-sd(beta_sim)
  
}

beta_true_mean=ab[1]/(ab[1]+ab[2])

beta_true_sd=sqrt((ab[1]*ab[2])/((ab[1]+ab[2])^2*(ab[1]+ab[2]+1)))

plot(mean_beta,type="l",col="cornflowerblue")
abline(h=beta_true_mean, col="red",lwd=2)

plot(sd_beta,type="l",col="darkslategray3")
abline(h=beta_true_sd,col="red",lwd=2)


## b

set.seed(123456)

n=10000
# Vector of beta parameters : a,b
ab <- c(16, 8)

# Simulate 1000 draws from the beta posterior: p_sim
p_sim <- rbeta(n, ab[1], ab[2])

# Construct a histogram of the simulated values
hist(p_sim,prob=T)
lines(density(p_sim),lwd=2,col="red")

# Compute the probability that P is smaller than 0.4
sum(p_sim < 0.4) / n
#the exact value is
pbeta(0.40, ab[1], ab[2])

quantile(p_sim, c(0.05, 0.95))

### c

phi=log(p_sim/(1-p_sim))

hist(phi,prob=T,border="blue")
lines(density(phi),lwd=2,col="red")


##########################################
## Exercise 2
##########################################


###  a


data<-c(14,25, 45, 25, 30, 33, 19, 50, 34, 67)
mean(log(data))
var(log(data))


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
  tau2<-(sum((Data-3.5)^2))/n
  PostDraws=matrix(0,NDraws,2)
  PostDraws[,2]<-(n*tau2)/rchisq(NDraws,n)
  PostDraws[,1]<-Datamean+rnorm(NDraws,0,1)*sqrt(PostDraws[,2]/n)
  
  return(PostDraws)
}

# library(LaplacesDemon)
# n<-10000
# d<-tau2<-(sum((log(data)-3.5)^2))/n
# chi_inv<-rinvchisq(n,d,scale=1/d)


PostDraws<-NormalNonInfoPrior(10000,log(data)) # Generating 1000 draws from the joint posterior density of mu and sigma^2
hist(PostDraws[,2],prob=T,border="blue",xlim=c(0,1.5))# Plotting the histogram of sigma-draws
lines(density(PostDraws[,2], adjust=2), col="red", lwd=2)

mean(PostDraws[,2])
var(PostDraws[,2])


# calculate the theoretical inverse chi-square

# make a grid
seq=seq(0,2.5,0.01)
#tao2=sum((log(salary)-3.5)^2)/10
# inverse chi-square function
invchi<-function(n,x){
  y=2^(-n/2)/gamma(n/2)*x^(-n/2+1)*exp(-1/(2*x))
  return(y)
}
# inverse scaled chi-squared function
invchiscaled<-function(n,x,s){
  y=(n/2)^(-n/2)/gamma(n/2)*s^n*x^(-n/2+1)*exp(-n*s^2/(2*x))
  return(y)
}
#values=invchis(10,seq,tao2)
# make a sample 
invchival=invchi(10,seq)
#ress=as.data.frame(cbind(seq,values))
plot(seq,invchival,type="l", panel.first=grid(25,25),
     col="aquamarine4",lwd=2,
     main=expression(paste("Plot of the theoretical Inv- ", chi^2)))


###  b

G=2*pnorm(sqrt(PostDraws[,2])/sqrt(2),mean=0,sd=1)-1
hist(G,col="lightblue",freq = F)
lines(density(G),col="red",lwd=2)



###  c

#tail credible intervalfor G
interval=quantile(G, c(0.025, 0.975))
interval

# hist(G,col="gray",freq = F)
# lines(density(G),col=col1,lwd=2)
# abline(v = interval[1], col="purple", lwd=2, lty=2)
# abline(v=interval[2],col="purple",lwd=2,lty=2)



library(WVPlots)

V2=G; V1=1:length(V2)


Z <- data.frame(V1,V2)


p <- WVPlots::ShadedDensity(frame = Z,
                            xvar="V2",
                            threshold =interval[1],
                            title = "Your title",tail="right")
p


library(dplyr)
library(ggplot2)

dens <- density(G)

data <- tibble(x = dens$x, y = dens$y) %>% 
  mutate(variable = case_when(
    (x >= interval[1] & x <= interval[2]) ~ "On",
    ( x<=interval[1]) ~ "Off",(x >=interval[2])~"onOff",
    TRUE ~ NA_character_))
#> Warning: package 'bindrcpp' was built under R version 3.4.4

ggplot(data, aes(x, y)) + geom_line(size=2) +
  geom_area(data = filter(data, variable == 'On'), fill = 'lightslateblue') + 
  geom_area(data = filter(data, variable == 'Off'), fill = 'orange')+
  geom_area(data = filter(data, variable == 'onOff'), fill = 'orange')


#Highest Posterior Density interval
dx <- density(G)
mat<-cbind(dx$x,dx$y)
mat<-mat[order(mat[,2], decreasing = TRUE),]
dn <- cumsum(mat[,2])/sum(mat[,2])
bi <- which(dn>=0.95)[1]
yy<-mat[bi,2]
yy


data1 <- tibble(x = dx$x, y = dx$y) %>% 
  mutate(variable = case_when(
    (y >= yy) ~ "On",
    ( y<yy) ~ "Off",
    TRUE ~ NA_character_))
#> Warning: package 'bindrcpp' was built under R version 3.4.4

ggplot(data1, aes(x, y)) + geom_line(size=2) +
  geom_area(data = filter(data1, variable == 'On'), fill = 'skyblue1') + 
  geom_area(data = filter(data1, variable == 'Off'), fill = 'slateblue')



#compute Highest Posterior Density interval with coda
library(coda)
HPDinterval(as.mcmc(G), prob=0.95)


##########################################
## Exercise 3
##########################################

###  a

dat<-c(-2.44,2.14,2.54,1.83,2.02,2.33,-2.79,2.23,2.07,2.02)

posterior<-function(k){
  return (exp(k*(sum(cos(dat-2.39))-1))/besselI(x = k, nu=0)^length(dat))
  #( exp(k*sum((cos(dat-2.39)-1)))/besselI(k,0))^(length(dat) )
}

steps=seq(0.1,8,0.01)
res<-double(length(seq(0.1,8,0.01)))


for(k in 1:length(steps)){
  res[k]<-posterior(steps[k])
}

#normalize the density so integral is 1
norm_constant=sum(res)*0.01

#normalized sample 
norm_res=res/norm_constant

par(mfrow=c(1,2))
plot(res,type="l")
plot(norm_res,type="l")



###  b

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(uniqv)]
}

getmode(norm_res)

#without function

indx_mode=which.max(norm_res)

mode=norm_res[indx_mode]

plot(norm_res,type="l")
abline(v=indx_mode,lwd=2,lty="dashed",col="red")



