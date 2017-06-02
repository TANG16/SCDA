rm(list = ls())
library(rjags)
library(coda)
set.seed(1345221)

save_on = TRUE

# ada=10000
# iter=500000
# th=1000
# cha=1

ada=500
iter=10000
th=1
cha=2


cat("Heron DA\n")

cat(sprintf("adapt = %i",ada),"\n")
cat(sprintf("iter = %i",iter),"\n")
cat(sprintf("thinning = %i",th),"\n")
cat(sprintf("chains = %i",cha),"\n")


source("heron_data.R")
source("heron_startingvals.R")


cat("Initialise the model:\n")
tstart = proc.time()
mod <- jags.model('heron_jags.R',data,inits,n.chains=cha,n.adapt=ada)
time_init = proc.time()-tstart

if (save_on) {
  save(mod, time_init, file = paste("Results/Heron_DA_model_ada",toString(ada),"_linux.RData",sep=""))
}


cat("Run the MCMC simulations:\n")
tstart = proc.time()
output1 <- coda.samples(mod,params,n.iter=iter,thin=th)
time_sample = proc.time()-tstart
 
if (save_on) {
  save(output1, time_sample, mod, time_init, file = paste("Results/Heron_DA_iter",toString(iter),"_ada",toString(ada),"_linux.RData",sep=""))
}

quit()