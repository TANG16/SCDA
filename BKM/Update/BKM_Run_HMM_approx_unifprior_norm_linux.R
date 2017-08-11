# Load required packages and fix the random seed
# setwd("BKM/Update")
# rm(list=ls())
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
iter=10000
# ada=10000
# iter=100000
th=1
cha=1

cat("BKM HMM adaptive approximation","\n")
cat(sprintf("adapt = %i",ada),"\n")
cat(sprintf("iter = %i",iter),"\n")
cat(sprintf("thinning = %i",th),"\n")
cat(sprintf("chains = %i",cha),"\n")


# source("BKM_Data_HMM_approx_unifprior_norm.R")
source("BKM_Data_HMM_approx_unifprior_norm_fun.R")

##########################################
cat("\n *** N_bin = 10 *** \n")

# Read data ###
data_10 <- data_fun(10)

# Set parameters and inital values
source("BKM_StartingVals_HMM_approx_unifprior_norm.R")
cat(sprintf("N_bin = %i",data_10$N_bin),"\n")


cat("Initialise the model:\n")
tstart=proc.time()
mod_10 <- jags.model('BKM_Bugs_HMM_approx_unifprior_sum_norm.R',data_10,inits,n.chains=cha,n.adapt=ada)
time_HMMnorm_init_10 <- proc.time()-tstart

cat("Run the MCMC simulations:\n")
tstart=proc.time()
outputHMMnorm_10 <- coda.samples(mod_10,params,n.iter=iter,thin=th)
time_HMMnorm_sample_10 <- proc.time()-tstart


##########################################

cat("\n *** N_bin = 20 *** \n")

# Read data ###
data_20 <- data_fun(20)

cat("Initialise the model:\n")
tstart=proc.time()
mod_20 <- jags.model('BKM_Bugs_HMM_approx_unifprior_sum_norm.R',data_20,inits,n.chains=cha,n.adapt=ada)
time_HMMnorm_init_20 <- proc.time()-tstart  

cat("Run the MCMC simulations:\n")
tstart=proc.time()
outputHMMnorm_20 <- coda.samples(mod_20,params,n.iter=iter,thin=th)
time_HMMnorm_sample_20 <- proc.time()-tstart

##########################################

cat("\n *** N_bin = 50 *** \n")

# Read data ###
data_50 <- data_fun(50)

cat("Initialise the model:\n")
tstart=proc.time()
mod_50 <- jags.model('BKM_Bugs_HMM_approx_unifprior_sum_norm.R',data_50,inits,n.chains=cha,n.adapt=ada)
time_HMMnorm_init_50 <- proc.time()-tstart  

cat("Run the MCMC simulations:\n")
tstart=proc.time()
outputHMMnorm_50 <- coda.samples(mod_50,params,n.iter=iter,thin=th)
time_HMMnorm_sample_50 <- proc.time()-tstart


##########################################

cat("\n *** N_bin = 5 *** \n")

# Read data ###
data_5 <- data_fun(5)

cat("Initialise the model:\n")
tstart=proc.time()
mod_5 <- jags.model('BKM_Bugs_HMM_approx_unifprior_sum_norm.R',data_5,inits,n.chains=cha,n.adapt=ada)
time_HMMnorm_init_5 <- proc.time()-tstart  

cat("Run the MCMC simulations:\n")
tstart=proc.time()
outputHMMnorm_5 <- coda.samples(mod_5,params,n.iter=iter,thin=th)
time_HMMnorm_sample_5 <- proc.time()-tstart


##########################################

if (save_on) {
  save(outputHMMnorm_10, outputHMMnorm_20, outputHMMnorm_50, outputHMMnorm_5,
       time_HMMnorm_sample_10, time_HMMnorm_sample_20, time_HMMnorm_sample_50, time_HMMnorm_sample_5,
       mod_10, mod_20, mod_50, mod_5,
       time_HMMnorm_init_10, time_HMMnorm_init_20, time_HMMnorm_init_50, time_HMMnorm_init_5, 
       file = paste("BKM_HMM_approx_iter",toString(iter),"_ada",toString(ada),"_COMB_unifprior_norm.RData",sep="")
       )
}
  
quit()

