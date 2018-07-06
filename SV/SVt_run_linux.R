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
nu <- 5;
# beta = 0.05;
beta <- 0.5;
mu <- 2*log(beta);
rho <- (nu-2)/nu

param <- c(mu, phi, sigma2, nu)

a1 = mu
P1 <- (sigma^2)/(1-phi^2)


################# simultated data ####
T <- 4000
h_true <- rep(NaN,T)

h_true[1] <- a1 + sqrt(P1)*rnorm(1)
for (t in c(2:T)){
  h_true[t] = mu + phi*(h_true[t-1]-mu) + sigma*rnorm(1)
}

y <- sqrt(rho)*exp(h_true/2)*rt(T,nu)
#y <- exp(h_true/2)*rt(T,nu)
# plot(y,type='l')
# lines(h_true, type='l', col='red')

################# full DA ####

# inits <- function()(list(mu = 0, phi = 0.97, sigma = 0.15, h = var(y)*rep(1,T)))
inits <- function()(list(mu = 0, phi_star = (0.97+1)/2, sigma2_star = 1/(0.15^2), nu = 8, h = log(var(y))*rep(1,T)))
params <- c('mu','sigma2','phi','nu','h')
data <- list(y=y,T=T)

svt_model_string <- "
model{
# define the observation process
for (t in 1:T){
# y[t] ~ dnorm(0, 1.0/exp(h[t]))
y[t] ~ dt(0, 1.0/(rho*exp(h[t])),nu)
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

rho <- (nu-2)/nu
nu ~ dexp(0.01)
}
"


### JAGS #####
tstart = proc.time()
svt_model_DA <- jags.model(textConnection(svt_model_string), 
                       data=data, inits=inits, n.chains = cha, n.adapt = ada)
time_init_DA = proc.time()-tstart

tstart = proc.time()
output_svt_DA <- coda.samples(svt_model_DA, params, n.iter = iter, thin=th)
time_sample_DA = proc.time()-tstart

# save selected output

mat1_DA = as.matrix(output_svt_DA[1]) 
mat2_DA = as.matrix(output_svt_DA[2]) 
mat_names_DA <- colnames(mat1_DA)

ESS_DA = lapply(output_svt_DA,effectiveSize)
ESS1_DA = as.matrix(ESS_DA[[1]])
ESS2_DA = as.matrix(ESS_DA[[2]])

theta1_DA = mat1_DA[,T+(1:4)] 
theta2_DA = mat2_DA[,T+(1:4)] 

H_short1_DA = mat1_DA[,seq(200,T,by=200)]
H_short2_DA = mat2_DA[,seq(200,T,by=200)]

mean_H1_DA = colMeans(mat1_DA[,1:T])
mean_H2_DA = colMeans(mat2_DA[,1:T])


if (save_on) {
  save(file=paste("SVt_DA_T",toString(T),"_selected.RData",sep=""),
       y, h_true, param, 
       svt_model_DA, time_init_DA, time_sample_DA, mat_names_DA,
       ESS1_DA, ESS2_DA, theta1_DA, theta2_DA, 
       H_short1_DA, H_short2_DA, mean_H1_DA, mean_H2_DA)
}
 
# if (save_on) {
#   save(file=paste("SV_DA_T",toString(T),".RData",sep=""),
#        y, h_true, param, output_sv_DA, sv_model_DA, time_init_DA, time_sample_DA)
# }

quit()