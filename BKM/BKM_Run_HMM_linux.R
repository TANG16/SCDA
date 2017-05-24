# Load required packages and fix the random seed
# setwd("BKM")
# install.packages(rjags)
# install.packages(coda)
# install.packages(lattice)
# rm(list=ls())

library(rjags)
library(coda)
library(lattice)
set.seed(134522)

save_on = TRUE
# MCMC details: ####

# ada=500
# iter=10000
# th=1
# cha=4

ada=100
iter=2000
th=1
cha=3

cat(sprintf("adapt = %i",ada),"\n")
cat(sprintf("iter = %i",iter),"\n")
cat(sprintf("thinning = %i",th),"\n")
cat(sprintf("chains = %i",cha),"\n")
# cat(sprintf("%.4f",pi),"\n")


# Read data ###
source("BKM_Data_HMM.R")
# Set parameters and inital values
source("BKM_StartingVals_HMM.R")


cat("V1: explicit pdf formulae\n")
# cat("V2: built-in pdf functions \n")


cat("Initialise the model:\n")
# Run the MCMC: ###
tstart=proc.time()
mod <- jags.model('BKM_Bugs_HMM.R',data,inits,n.chains=cha,n.adapt=ada)
# mod <- jags.model('BKM_Bugs_HMM_v2.R',data,inits,n.chains=cha,n.adapt=ada)
temp=proc.time()-tstart
time_HMM_init <- temp 
if (save_on) {
  save(mod, time_HMM_init, file = paste("Results/BKM_HMM_model_ada",toString(ada),"_linux.RData",sep=""))
  # save(mod, time_HMM_init, file = paste("Results/BKM_HMM_v2_model_ada",toString(ada),"_linux.RData",sep=""))
}


cat("Run the MCMC simulations:\n")
tstart=proc.time()
output1 <- coda.samples(mod,params,n.iter=iter,thin=th)
temp=proc.time()-tstart
time_HMM_sample <- temp 
if (save_on) {
  save(output1, time_HMM_sample, mod, time_HMM_init, file = paste("Results/BKM_HMM_iter",toString(iter),"_ada",toString(ada),"_linux.RData",sep=""))
  # save(output1, time_HMM_sample, mod, time_HMM_init, file = paste("Results/BKM_HMM_v2_iter",toString(iter),"_ada",toString(ada),"_linux.RData",sep=""))
}

quit()