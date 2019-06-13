# Rstan -Example 1: Eight Schools -----------------------------------------

library(rstan)

Stan_model='data {
  int<lower=0> J;         // number of schools 
  real y[J];              // estimated treatment effects
  real<lower=0> sigma[J]; // standard error of effect estimates 
}
parameters {
  real mu;                // population treatment effect
  real<lower=0> tau;      // standard deviation in treatment effects
  vector[J] eta;          // unscaled deviation from mu by school
}
transformed parameters {
  vector[J] theta = mu + tau * eta;        // school treatment effects
}
model {
  target += normal_lpdf(eta | 0, 1);       // prior log-density
  target += normal_lpdf(y | theta, sigma); // log-likelihood
}'


schools_data <- list(J = 8, 
                    y = c(28,  8, -3,  7, -1,  1, 18, 12),
                    sigma = c(15, 10, 16, 11,  9, 11, 10, 18))

fit = stan( model_code = Stan_model,  # Stan model
            data = schools_data,    # named list of data
            chains = 4,             # number of Markov chains
            warmup = 1000,          # number of warmup iterations per chain
            iter = 2000,            # total number of iterations per chain
            cores = 2,              # number of cores (could use one per chain)
            refresh = 0 )

# Print the fitted model
print(fit, pars=c("theta", "mu", "tau", "lp__"), probs=c(.1,.5,.9))


par(mfrow = c(1,1))
plot(fit)
# Do automatic traceplots of all chains
traceplot(fit, pars = c("mu", "tau"), inc_warmup = TRUE, nrow = 2)

# Pair plot
pairs(fit, pars = c("mu", "tau", "lp__"), las = 1)


# Rstan -Example 2: iid Normal  -------------------------------------------

y1 = c(4,5,6,4,0,2,5,3,8,6,10,8)
N1 = length(y1)

StanModel = '
data {
int<lower=0> N; // Number of observations
int<lower=0> y[N]; // Number of flowers
}
parameters {
real mu;
real<lower=0> sigma2;
}
model {
mu ~ normal(0,100); // Normal with mean 0, st.dev. 100
sigma2 ~ scaled_inv_chi_square(1,2); // Scaled-inv-chi2 with nu 1, sigma 2
for(i in 1:N)
y[i] ~ normal(mu,sqrt(sigma2));
}'

# make data list
dat <- list(N= N1, 
            y = y1)
# fit stan model
fit1<-stan(model_code = StanModel,
           data=dat,
           chains=2,
           warmup = 1000,
           iter=2000,
           cores=2,
           refresh=0)

# print fit
print(fit1)
# plot fit
plot(fit1)
# plot traceplot
traceplot(fit1, pars = c("mu", "sigma2"), inc_warmup = TRUE, nrow = 2)


# Rstan -Example 3: Fitting a nonlinear model -----------------------------

y_i=(a1*exp(-b1*x_i)+a_2*exp(-b_2*x_i))*e_i , where e_i~N(0,sigma2)

exponential_model='data { int N; 
       vector[N] x; 
       vector[N] y; } 

parameters { vector[2] log_a; 
             ordered[2] log_b; 
             real<lower=0> sigma; }

transformed parameters { vector<lower=0>[2] a; 
                         vector<lower=0>[2] b; 
                         a <- exp(log_a); 
                         b <- exp(log_b); } 

model { vector[N] ypred; 
        ypred <- a[1]*exp(-b[1]*x) + a[2]*exp(-b[2]*x); 
        y ~ lognormal(log(ypred), sigma); }'


# Set up the true parameter values 
a <- c(.8, 1); b <- c(2, .1); sigma <- .2

# Simulate data 
x <- (1:1000)/100 
N <- length(x) 
ypred <- a[1]*exp(-b[1]*x) + a[2]*exp(-b[2]*x) 
y <- ypred*exp(rnorm(N, 0, sigma))


plot(x,y,xlab = 'x',ylab='y')
lines(x,ypred,col="red",lwd=2)
legend("topright",legend=c("data","prediction"),
       col=c("black","red"),lty=c(NA,1),pch=c(1,NA))


# Fit the model 
fittedModel <- stan(model_code=exponential_model, 
                    data=list(N=N, x=x, y=y), 
                    iter=1000, 
                    chains=4)

print(fittedModel, pars=c("a", "b", "sigma"))

plot(fittedModel)

traceplot(fittedModel)

f.mod<-extract(fittedModel)

f.mod$a



# Example -Multilevel Normal ----------------------------------------------

StanMod = '
data {
int<lower=0> N; // Number of observations
int<lower=0> y[N]; // Number of flowers
int<lower=0> P; // Number of plants
}
transformed data {
int<lower=0> M; // Number of months
M = N / P;
}
parameters {
real mu;
real<lower=0> sigma2;
real mup[P];
real sigmap2[P];
}
model {
mu ~ normal(0,100); // Normal with mean 0, st.dev. 100
sigma2 ~ scaled_inv_chi_square(1,2); // Scaled-inv-chi2 with nu 1, sigma 2
for(p in 1:P){
mup[p] ~ normal(mu,sqrt(sigma2));
for(m in 1:M)
y[M*(p-1)+m] ~ normal(mup[p],sqrt(sigmap2[p]));
}
}'


y1 = c(4,5,6,4,0,2,5,3,8,6,10,8)
N1 = length(y1)
P=3

my_data = list(N=N1, y=y1, P=P)
burnin = 1000
niter = 2000
fit = stan(model_code=StanMod,data=my_data,
           warmup=burnin,iter=niter,chains=4)
# Print the fitted model
print(fit,digits_summary=3)



