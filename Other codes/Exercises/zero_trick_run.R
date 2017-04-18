# Use the "zeros trick" to model a normal distribution with unknown mean mu and 
# unknown standard deviation sigma, including predicting a new observation.

library(rjags) 
library(coda)

ada=100
iter=1000
th=1
cha=1

# Data
z <- rep(0,8)
y <- c(-1, -0.3, 0.1, 0.2, 0.7, 1.2, 1.7, NA)


data <-list(y=y, z=z)

# Initial values
params <- c('mu','sigma','y[8]')
inits <- function(){list(mu = 0, sigma = 1, y = c(NA, NA, NA, NA, NA, NA, NA, 0))} 

# Run the MCMC:
mod <- jags.model('zero_trick_bugs.R',data,inits,n.chains=cha,n.adapt=ada)
output <- coda.samples(mod,params,n.iter=iter,thin=th)
summary(output)