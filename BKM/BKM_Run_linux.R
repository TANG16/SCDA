# setwd("BKM")
# rm(list=ls())

# Load required packages and fix the random seed
library(rjags)
library(coda)
library(lattice)
set.seed(134522)

save_on = TRUE
scaled_on = FALSE

# MCMC details: ####

# ada=1000
# iter=3000
# th=1
# cha=2
ada=100
iter=2000
th=1
cha=3


# Read data
if (scaled_on){
  source("BKM_Data_scaled.R")
  # Set parameters and inital values
  source("BKM_StartingVals_scaled.R") 
}else {
  source("BKM_Data.R")
  # Set parameters and inital values
  source("BKM_StartingVals.R")
}

cat("Initialise the model:\n")
# Compile the model: ####
if (scaled_on){
  tstart=proc.time()
  mod <- jags.model('BKM_Bugs_scaled.R',data,inits,n.chains=cha,n.adapt=ada)
  temp=proc.time()-tstart
  time_init <- temp 
  if (save_on) {
    save(mod, time_init, file = "Results/BKM_model_scaled_linux.RData")
  } 
}else {
  tstart=proc.time()
  mod <- jags.model('BKM_Bugs.R',data,inits,n.chains=cha,n.adapt=ada)
  temp=proc.time()-tstart
  time_init <- temp # ada = 100 --> PC: 1.23; ada = 1000 --> PC: 6.29
  if (save_on) {
    save(mod, time_init, file = "Results/BKM_model_linux.RData")
  }
}


cat("Run the MCMC simulations:\n")
# Run the MCMC: #### 
tstart=proc.time()
output1 <- coda.samples(mod,params,n.iter=iter,thin=th)
temp=proc.time()-tstart
time_sample <- temp # PC org:   5.23

if (save_on) {
  if (scaled_on){
    save(output1, time_sample, mod, time_init, file = paste("Results/BKM_iter",toString(iter),"_ada",toString(ada),"_scaled_linux.RData",sep=""))
}else {
    save(output1, time_sample, mod, time_init, file = paste("Results/BKM_iter",toString(iter),"_ada",toString(ada),"_linux.RData",sep=""))
  } 
}
quit()