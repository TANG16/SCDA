# setwd("Heron")
rm(list=ls())
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

source("heron_data_HMM_approx_unifprior.R")
source("heron_startingvals_HMM_approx_unifprior.R")

cat(sprintf("bin1 size = %i",bin_size1),"\n")
cat(sprintf("bin3 size = %i",bin_size3),"\n")

cat(sprintf("N bin min1 = %i",N_bin_min1),"\n")
cat(sprintf("N bin min3 = %i",N_bin_min3),"\n")

cat("Initialise the model:\n")
tstart = proc.time()
mod <- jags.model('heron_jags_HMM_approx_unifprior.R',data,inits,n.chains=cha,n.adapt=ada)
time_init_HMM = proc.time()-tstart

if (save_on) {
  save(mod, time_init_HMM, file = paste("Results/Heron_HMM_approx_model_ada",toString(ada),"_linux_bin1",toString(bin_size1),"_bin3",toString(bin_size3),"_unifprior.RData",sep=""))
}

cat("Run the MCMC simulations:\n")
tstart = proc.time()
output1 <- coda.samples(mod,params,n.iter=iter,thin=th)
time_sample_HMM = proc.time()-tstart
 
if (save_on) {
  save(output1, time_sample_HMM, mod, time_init_HMM, file = paste("Results/Heron_HMM_approx_iter",toString(iter),"_ada",toString(ada),"_linux_bin1",toString(bin_size1),"_bin3",toString(bin_size3),"_unifprior.RData",sep=""))
}
 

quit()
