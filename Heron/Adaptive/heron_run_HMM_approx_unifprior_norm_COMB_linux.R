# setwd("Heron/Adaptive")
# rm(list=ls())
library(rjags)
library(coda)
set.seed(1345221)


# ada=10000
# iter=500000
# th=1000
# cha=1
save_on = TRUE

ada= 1000
iter=30000
th=1
cha=2 #2

cat("Heron HMM\n")

cat(sprintf("adapt = %i",ada),"\n")
cat(sprintf("iter = %i",iter),"\n")
cat(sprintf("thinning = %i",th),"\n")
cat(sprintf("chains = %i",cha),"\n")

# Read data ###
# source("heron_data_HMM_approx_unifprior_norm.R")
source("heron_data_HMM_approx_unifprior_norm_fun.R")
data_10_5 <- data_fun(10, 50)



# source("heron_startingvals_HMM_approx_unifprior_norm.R")
source("heron_startingvals_HMM_approx_unifprior_norm_fun.R")

startingvals <- startingvals_fun(data_10_5)
inits <- startingvals$inits
params <- startingvals$params


##########################################
cat("\n *** N_bin1 = 10, N_bin3 = 5 *** \n")

cat("Initialise the model:\n")
tstart = proc.time()
mod_10_5 <- jags.model('heron_jags_HMM_approx_unifprior_norm.R',data_10_5,
                       inits,n.chains=cha,n.adapt=ada)
time_init_HMM_10_5 = proc.time()-tstart

cat("Run the MCMC simulations:\n")
tstart = proc.time()
output_10_5 <- coda.samples(mod_10_5,params,n.iter=iter,thin=th)
time_sample_HMM_10_5 = proc.time()-tstart

##########################################
cat("\n *** N_bin1 = 50, N_bin3 = 40 *** \n")
data_50_40 <- data_fun(50, 40)

cat("Initialise the model:\n")
tstart = proc.time()
mod_50_40 <- jags.model('heron_jags_HMM_approx_unifprior_norm.R',data_50_40,
                        inits,n.chains=cha,n.adapt=ada)
time_init_HMM_50_40 = proc.time()-tstart

cat("Run the MCMC simulations:\n")
tstart = proc.time()
output_50_40 <- coda.samples(mod_50_40,params,n.iter=iter,thin=th)
time_sample_HMM_50_40 = proc.time()-tstart

##########################################
cat("\n *** N_bin1 = 100, N_bin3 = 70 *** \n")
# Read data ###

data_100_70 <- data_fun(100,70)

cat("Initialise the model:\n")
tstart = proc.time()
mod_100_70 <- jags.model('heron_jags_HMM_approx_unifprior_norm.R',data_100_70,
                         inits,n.chains=cha,n.adapt=ada)
time_init_HMM_100_70 = proc.time()-tstart

cat("Run the MCMC simulations:\n")
tstart = proc.time()
output_100_70 <- coda.samples(mod_100_70,params,n.iter=iter,thin=th)
time_sample_HMM_100_70 = proc.time()-tstart


##########################################

if (save_on) {
  save(output_10_5, output_50_40, output_100_70,
       time_sample_HMM_10_5, time_sample_HMM_50_40, time_sample_HMM_100_70,
       mod_10_5, mod_50_40, mod_100_70, 
       time_init_HMM_10_5, time_init_HMM_50_40, time_init_HMM_100_70,
       file = paste("Heron_HMM_approx_iter",toString(iter),"_ada",toString(ada),"_linux_COMB_unifprior_norm.RData",sep=""))
}


quit()
