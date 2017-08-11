# setwd("Heron")
rm(list = ls())
library(rjags)
library(coda)
set.seed(1345221)

save_on = TRUE

# ada=10000
# iter=500000
# th=1000
# cha=1

ada=1000
iter=30000
th=1
cha= 1 #2


cat("Heron DA\n")

cat(sprintf("adapt = %i",ada),"\n")
cat(sprintf("iter = %i",iter),"\n")
cat(sprintf("thinning = %i",th),"\n")
cat(sprintf("chains = %i",cha),"\n")


source("heron_data.R")
# source("heron_startingvals.R")
load(file = "X2_X4_inits_from_HMM.RData")

cat("Initialise the model:\n")
tstart = proc.time()
mod_DA <- jags.model('heron_jags.R',data,inits,n.chains=cha,n.adapt=ada)
time_init_DA = proc.time()-tstart


cat("Run the MCMC simulations:\n")
tstart = proc.time()
output_DA <- coda.samples(mod_DA,params,n.iter=iter,thin=th)
time_sample_DA = proc.time()-tstart
 
if (save_on) {
  save(output_DA, time_sample_DA, mod_DA, time_init_DA, 
       file = paste("Results/Heron_DA_iter",toString(iter),"_ada",toString(ada),"_HMMinit_linux.RData",sep=""))
}

quit()