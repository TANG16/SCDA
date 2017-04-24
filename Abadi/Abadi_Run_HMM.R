# Set working directory: ####
setwd("Abadi")

save_on = TRUE
# Load required packages and fix the random seed
library(rjags)
library(coda)
library(lattice)
set.seed(134522)

# MCMC details: ####

ada=100
iter=10000
th=1
cha=2

# Load data: ###
source("Abadi_Data_HMM.R")
# Set initial parameter values: ####
source("Abadi_StartingVals_HMM.R")

# Run the MCMC: ####
tstart = proc.time()
mod <- jags.model('Abadi_Bugs_HMM.R',data,inits,n.chains=cha,n.adapt=ada)
temp = proc.time()-tstart
time_init <- temp

tstart = proc.time()
output1 <- coda.samples(mod,params,n.iter=iter,thin=th)
temp = proc.time()-tstart
time_sample <- temp

if (save_on) {
  save(mod, output1, time_init, time_sample, file = "Abadi_HMM_iter10000.RData")
}

summary(output1)

mat1 = as.matrix(output1[1])
mat2 = as.matrix(output1[2])

# Retrieve the names ####
mat1_names <- colnames(mat1)


frame1 <- do.call(rbind.data.frame, output1)
frame2 <- do.call(rbind.data.frame, output2)

str(output1) # List of 2
str(output2) # List of 12


plot(output1[,2])
output1[j] [,"par_name"]
# it returns the j-th chain for the parameter of interest
plot(output1[1][,2])
plot(output1[2][,2])
plot(output1[1][,2])

summary(output1)

mat1 <- as.matrix(output1)