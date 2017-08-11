# setwd("SV")
rm(list = ls())
library(rjags)
library(coda)
set.seed(1345221)

save_on = TRUE


ada=1000
iter=10000
th=1
cha= 2 #2



phi <- 0.98
sigma <- 0.2
sigma2 <- sigma^2
# beta = 0.05;
beta <- 0.5;
mu <- 2*log(beta);

param <- c(mu, phi, sigma2)

a1 = mu
P1 <- (sigma^2)/(1-phi^2)


################# simultated data ####
T <- 1000
h_true <- rep(NaN,T)

h_true[1] <- a1 + sqrt(P1)*rnorm(1)
for (t in c(2:T)){
  h_true[t] = mu + phi*(h_true[t-1]-mu) + sigma*rnorm(1)
}
y <- exp(h_true/2)*rnorm(T)
# plot(y,type='l')
# lines(h_true, type='l', col='red')

################# full DA ####

# inits <- function()(list(mu = 0, phi = 0.97, sigma = 0.15, h = var(y)*rep(1,T)))
inits <- function()(list(mu = 0, phi_star = (0.97+1)/2, sigma2_star = 1/(0.15^2), h = log(var(y))*rep(1,T)))
params <- c('mu','sigma2','phi','h')
data <- list(y=y,T=T)

sv_model_string <- "
model{
# define the observation process
for (t in 1:T){
y[t] ~ dnorm(0, 1.0/exp(h[t]))
}

# define the state process
for (t in 2:T){
h[t] ~ dnorm(mean_h[t],sigma2_star)
mean_h[t] <- mu + phi*(h[t-1] - mu) 
}

h[1] ~ dnorm(mean_h[1],sigma2_star)
mean_h[1] <- mu + phi*(h0 - mu) 

h0 ~ dnorm(mean_h0,1/P1)
P1 <- sigma2/(1-phi^2)
mean_h0 <- mu

# define the priors for parameters
mu ~ dnorm(0,0.1)

# (phi+1)/2 ~ betapdf((phi+1)/2, 20, 1.5)
phi <- 2*phi_star - 1
phi_star ~ dbeta(20, 1.5)
# (phi+1)/2 ~ dbeta(20, 1.5)

# 1/sigma2 ~ gampdf(1./s2, 5/2, 0.05/2);
sigma2 <- 1/sigma2_star
#   sigma2_star ~ dgamma(5/2, 2/0.05) # shape and rate
sigma2_star ~ dgamma(5/2, 0.05/2) # shape and rate
# 1/sigma2  ~  dgamma(5/2, 2/0.05)
}
"


### JAGS #####
tstart = proc.time()
sv_model_DA <- jags.model(textConnection(sv_model_string), 
                       data=data, inits=inits, n.chains = cha, n.adapt = ada)
time_init_DA = proc.time()-tstart

tstart = proc.time()
output_sv_DA <- coda.samples(sv_model_DA, params, n.iter = iter, thin=th)
time_sample_DA = proc.time()-tstart


if (save_on) {
  save(file="SV_DA.RData",y, h_true, param, output_sv_DA, sv_model_DA, time_init_DA, time_sample_DA)
}

quit()