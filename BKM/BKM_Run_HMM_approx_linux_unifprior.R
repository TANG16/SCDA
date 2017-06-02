# Load required packages and fix the random seed
# setwd("BKM")
# install.packages(rjags)
# install.packages(coda)
# install.packages(lattice)

library(rjags)
library(coda)
library(lattice)
set.seed(134522)

save_on = TRUE
# MCMC details:  ####

ada=1000
iter=2000
th=1
cha=3

cat("BKM HMM approximation","\n")
cat(sprintf("adapt = %i",ada),"\n")
cat(sprintf("iter = %i",iter),"\n")
cat(sprintf("thinning = %i",th),"\n")
cat(sprintf("chains = %i",cha),"\n")


# Read data ###
source("BKM_Data_HMM_approx_unifprior.R")
# Set parameters and inital values
source("BKM_StartingVals_HMM_approx_unifprior.R")
cat(sprintf("scale = %i",sc),"\n")
cat(sprintf("bin size = %i",bin_size),"\n")

cat("Initialise the model:\n")
tstart=proc.time()
mod <- jags.model('BKM_Bugs_HMM_approx_unifprior.R',data,inits,n.chains=cha,n.adapt=ada)
temp=proc.time()-tstart
time_HMM_init <- temp  
if (save_on) {
  save(mod, time_HMM_init, file = paste("Results/BKM_HMM_approx_model_ada",toString(ada),"_linux_sc",toString(sc),"_bin",toString(bin_size),"_unifprior.RData",sep=""))
}
 
cat("Run the MCMC simulations:\n")
  
tstart=proc.time()
output1 <- coda.samples(mod,params,n.iter=iter,thin=th)
temp=proc.time()-tstart
time_HMM_sample <- temp  
if (save_on) {
  save(output1, time_HMM_sample, mod, time_HMM_init, file = paste("Results/BKM_HMM_approx_iter",toString(iter),"_ada",toString(ada),"_linux_sc",toString(sc),"_bin",toString(bin_size),"_unifprior.RData",sep=""))
}

quit()