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
# MCMC details: ####

# ada=500
# iter=10000
# th=1
# cha=4

ada=100
iter=2000
th=1
cha=3

cat(sprintf("adapt = %i",ada),"\n")
cat(sprintf("iter = %i",iter),"\n")
cat(sprintf("thinning = %i",th),"\n")
cat(sprintf("chains = %i",cha),"\n")
# cat(sprintf("%.4f",pi),"\n")


# Read data ###
source("BKM_Data_HMM.R")
# Set parameters and inital values
source("BKM_StartingVals_HMM.R")


cat("Initialise the model:\n")
# Run the MCMC: ###
# if (file.exists("BKM_HMM_model_ada500_linux.RData")) {
#   load("BKM_HMM_model_ada500_linux.RData")
#   tstart=proc.time()
#   mod$recompile()
#   temp=proc.time()-tstart
#   time_HMM_init_recomp <- temp
#   if (save_on) {
#     save(mod, time_HMM_init_recomp, file = "BKM_HMM_model_ada500_linux_recompile.RData")
#   }
# } else {
  tstart=proc.time()
  mod <- jags.model('BKM_Bugs_HMM.R',data,inits,n.chains=cha,n.adapt=ada)
  temp=proc.time()-tstart
  time_HMM_init <- temp # ada=100 PC: 357.39 ~ 6min --> 188.28 ~3 min with 0!
  # ada = 1000 --> PC: 1798.68 
  # ada = 100 --> laptop: 659; linux: 490.122 ~ 8 min
  if (save_on) {
    save(mod, time_HMM_init, file = paste("Results/BKM_HMM_model_ada",toString(ada),"_linux.RData",sep=""))
  }
# }

cat("Run the MCMC simulations:\n")
  
tstart=proc.time()
output1 <- coda.samples(mod,params,n.iter=iter,thin=th)
temp=proc.time()-tstart
time_HMM_sample <- temp # PC:   1843.02 ~31 min
# ada = 1000 --> 1741.71
# ada = 100 --> laptop: 6449 ~ 108 min, linux 4868 ~ 81 min
if (save_on) {
  save(output1, time_HMM_sample, mod, time_HMM_init, file = paste("Results/BKM_HMM_iter",toString(iter),"_ada",toString(ada),"_linux.RData",sep=""))
}

quit()