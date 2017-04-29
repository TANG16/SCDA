# Set working directory: ####
# setwd("Abadi")

save_on = TRUE
# Load required packages and fix the random seed
library(rjags)
library(coda)
library(lattice)
set.seed(134522)

# MCMC details: ####

ada=100
iter=10000
th=1
cha=2

# Load data: ###
source("Abadi_Data.R")
# Set initial parameter values: ####
source("Abadi_StartingVals.R")

# Initialise the model: ####
tstart = proc.time()
mod <- jags.model('Abadi_Bugs.R',data,inits,n.chains=cha,n.adapt=ada)
temp = proc.time()-tstart
time_init <- temp

# Run the MCMC simulations: ####
tstart = proc.time()
output1 <- coda.samples(mod,params,n.iter=iter,thin=th)
temp = proc.time()-tstart
time_sample <- temp

if (save_on) {
  save(mod, output1, time_init, time_sample, file = paste("Results/Abadi_iter",toString(iter),"linux.RData",sep=""))
}

quit()