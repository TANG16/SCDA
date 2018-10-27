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
iter=10000
# ada=10000
# iter=10000
th=1
cha=1

cat("BKM HMM approximation","\n")
cat(sprintf("adapt = %i",ada),"\n")
cat(sprintf("iter = %i",iter),"\n")
cat(sprintf("thinning = %i",th),"\n")
cat(sprintf("chains = %i",cha),"\n")


# Read data ###
# source("BKM_Data_HMM_approx_unifprior.R")
source("BKM_Data_HMM_approx_unifprior_fun.R")


# Set parameters and inital values
source("BKM_StartingVals_HMM_approx_unifprior.R")
# cat(sprintf("scale = %i",sc),"\n")
# cat(sprintf("bin size = %i",bin_size),"\n")


##########################################

cat("\n *** N_bin = 29, bin_size = 29 *** \n")

data_29_29 <- data_fun_fixed()

cat("Initialise the model:\n")
tstart=proc.time()
mod_29_29<- jags.model('BKM_Bugs_HMM_approx_unifprior.R',data_29_29,inits,n.chains=cha,n.adapt=ada)
temp=proc.time()-tstart
time_HMM_init_29_29 <- temp   

cat("Run the MCMC simulations:\n")
tstart=proc.time()
output_HMM_29_29 <- coda.samples(mod_29_29,params,n.iter=iter,thin=th)
temp=proc.time()-tstart
time_HMM_sample_29_29 <- temp  

##########################################

cat("\n *** N_bin = 15, bin_size = 55 *** \n")

data_15_55 <- data_fun_fixed(N_bin = 15,bin_size = 55)

cat("Initialise the model:\n")
tstart=proc.time()
mod_15_55 <- jags.model('BKM_Bugs_HMM_approx_unifprior.R',data_15_55,inits,n.chains=cha,n.adapt=ada)
temp=proc.time()-tstart
time_HMM_init_15_55 <- temp   

cat("Run the MCMC simulations:\n")
tstart=proc.time()
output_HMM_15_55 <- coda.samples(mod_15_55,params,n.iter=iter,thin=th)
temp=proc.time()-tstart
time_HMM_sample_15_55 <- temp  

##########################################

cat("\n *** N_bin = 5, bin_size = 169 *** \n")

data_5_169 <- data_fun_fixed(N_bin = 5, bin_size = 169)

cat("Initialise the model:\n")
tstart=proc.time()
mod_5_169 <- jags.model('BKM_Bugs_HMM_approx_unifprior.R',data_5_169,inits,n.chains=cha,n.adapt=ada)
temp=proc.time()-tstart
time_HMM_init_5_169 <- temp   

cat("Run the MCMC simulations:\n")
tstart=proc.time()
output_HMM_5_169 <- coda.samples(mod_5_169,params,n.iter=iter,thin=th)
temp=proc.time()-tstart
time_HMM_sample_5_169 <- temp  


##########################################

if (save_on) {
  file_save = paste("C:/Users/ab507t/Desktop/BKM_HMM_approx_iter",toString(iter),"_ada",toString(ada),"_COMB_unifprior.RData",sep="")
  
  save(output_HMM_29_29, output_HMM_15_55, output_HMM_5_169,
       time_HMM_sample_29_29, time_HMM_sample_15_55, time_HMM_sample_5_169,
       mod_29_29, mod_15_55, mod_5_169, 
       time_HMM_init_29_29, time_HMM_init_15_55, time_HMM_init_5_169,
       file = file_save)
}


quit()
