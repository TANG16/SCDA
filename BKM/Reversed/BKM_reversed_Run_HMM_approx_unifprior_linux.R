# Load required packages and fix the random seed
# setwd("BKM/Reversed")
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
iter=100000
# ada=10000
# iter=10000
th=1
cha=2

cat("BKM reversed HMM approximation","\n")
cat(sprintf("adapt = %i",ada),"\n")
cat(sprintf("iter = %i",iter),"\n")
cat(sprintf("thinning = %i",th),"\n")
cat(sprintf("chains = %i",cha),"\n")


# Read data ###
# source("BKM_Data_HMM_approx_unifprior.R")
source("BKM_reversed_Data_HMM_approx_unifprior_fun.R")


# Set parameters and inital values
source("BKM_reversed_StartingVals_HMM_approx_unifprior.R")
# cat(sprintf("scale = %i",sc),"\n")
# cat(sprintf("bin size = %i",bin_size),"\n")


##########################################

# cat("\n *** N_bin = 50, bin_size = 45 *** \n")
# 
# data_50_45 <- data_fun_fixed( N_bin = 50, bin_size = 45 )
# 
# cat("Initialise the model:\n")
# tstart=proc.time()
# mod_50_45<- jags.model('BKM_reversed_Bugs_HMM_approx_unifprior.R',data_50_45,inits,n.chains=cha,n.adapt=ada)
# temp=proc.time()-tstart
# time_HMM_init_50_45 <- temp   
# 
# cat("Run the MCMC simulations:\n")
# tstart=proc.time()
# output_HMM_50_45 <- coda.samples(mod_50_45,params,n.iter=iter,thin=th)
# temp=proc.time()-tstart
# time_HMM_sample_50_45 <- temp  

##########################################


# if (save_on) {
#   save(output_HMM_50_45, 
#        time_HMM_sample_50_45,  
#        mod_50_45, 
#        time_HMM_init_50_45, 
#        file = paste("BKM_HMM_approx_iter",toString(iter),"_ada",toString(ada),"_reversed.RData",sep=""))
# }

##########################################

cat("\n *** N_bin = 30, bin_size = 75 *** \n")

data_30_75 <- data_fun_fixed()

cat("Initialise the model:\n")
tstart=proc.time()
mod_30_75<- jags.model('BKM_reversed_Bugs_HMM_approx_unifprior.R',data_30_75,inits,n.chains=cha,n.adapt=ada)
temp=proc.time()-tstart
time_HMM_init_30_75 <- temp   

cat("Run the MCMC simulations:\n")
tstart=proc.time()
output_HMM_30_75 <- coda.samples(mod_30_75,params,n.iter=iter,thin=th)
temp=proc.time()-tstart
time_HMM_sample_30_75 <- temp  

##########################################


if (save_on) {
  save(output_HMM_30_75, 
       time_HMM_sample_30_75,  
       mod_30_75, 
       time_HMM_init_30_75, 
       file = paste("BKM_HMM_approx_iter",toString(iter),"_ada",toString(ada),"_reversed.RData",sep=""))
}

quit()
