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

ada=100
iter=5000
th=1
cha=1

cat("Heron HMM\n")

cat(sprintf("adapt = %i",ada),"\n")
cat(sprintf("iter = %i",iter),"\n")
cat(sprintf("thinning = %i",th),"\n")
cat(sprintf("chains = %i",cha),"\n")

source("heron_data_HMM.R")
source("heron_startingvals_HMM.R")

cat("Initialise the model:\n")
tstart = proc.time()
mod <- jags.model('heron_jags_HMM.R',data,inits,n.chains=cha,n.adapt=ada)
time_init_HMM = proc.time()-tstart

if (save_on) {
  save(mod, time_init_HMM, file = paste("Results/Heron_HMM_model_ada",toString(ada),"_linux.RData",sep=""))
}

cat("Run the MCMC simulations:\n")
tstart = proc.time()
output1 <- coda.samples(mod,params,n.iter=iter,thin=th)
time_sample_HMM = proc.time()-tstart
 
if (save_on) {
  save(output1, time_sample_HMM, file = paste("Results/Heron_HMM_iter",toString(iter),"_ada",toString(ada),"_linux.RData",sep=""))
}
 

quit()