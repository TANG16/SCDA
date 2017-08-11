setwd("BKM/StateSpace_checks/")
rm(list=ls())

# Load required packages and fix the random seed
library(rjags)
library(coda)
library(lattice)
set.seed(134522)

save_on = FALSE

# MCMC details: ####

# # ada=1000
# # iter=3000
# # th=1
# # cha=2
# ada=100
# iter=2000
# th=1
# cha=3
ada=100
iter=10000
th=1
cha=2

# Read data
source("BKM_Data_scaled.R")
# Set parameters and inital values
source("BKM_StartingVals_scaled.R") 


# Compile the model: ####
tstart=proc.time()
mod <- jags.model('BKM_Bugs_scaled.R',data,inits,n.chains=cha,n.adapt=ada)
temp=proc.time()-tstart
time_init <- temp #0.14
if (save_on) {
  save(mod, time_init, file = "BKM_model_scaled.RData")
} 

# Compile the REVERSED model: ####
tstart=proc.time()
mod2 <- jags.model('BKM_Bugs_scaled_reversed.R',data,inits,n.chains=cha,n.adapt=ada)
temp=proc.time()-tstart
time_init2 <- temp #0.12
if (save_on) {
  save(mod2, time_init2, file = "BKM_model_scaled_reversed.RData")
} 



# Run the MCMC: #### 
tstart=proc.time()
output1 <- coda.samples(mod,params,n.iter=iter,thin=th)
temp=proc.time()-tstart
time_sample <- temp  # 8.86 

if (save_on) {
  save(output1, time_sample, file = paste("BKM_iter",toString(iter),"_ada",toString(ada),"_scaled.RData",sep=""))
}

# Run the MCMC for the reversed model: #### 
tstart=proc.time()
output2 <- coda.samples(mod2,params,n.iter=iter,thin=th)
temp=proc.time()-tstart
time_sample2 <- temp # 8.66 

if (save_on) {
  save(output2, time_sample2, file = paste("BKM_iter",toString(iter),"_ada",toString(ada),"_scaled_reversed.RData",sep=""))
}




# Collect the results ####
mat1 = as.matrix(output1[1])
mat2 = as.matrix(output1[2])
mat3 = as.matrix(output1[3])
summary(output1)

mat1_names <- colnames(mat1)
 


par(mfrow=c(2,1),oma=c(0,0,1.5,0))
plot(colMeans(mat1[,(1:T)]), type='l', xlab ="", ylab="", sub="N1")
plot(colMeans(mat1[,(T+1):(T+T)]), type='l', xlab ="", ylab="", sub="Na")
mtext("Posterior means", outer=TRUE, cex=1)


# Trace plots #### 

par(mfrow=c(3,3), oma = c(0, 0, 1.5, 0))
for (i in c(73:81)){
  plot(mat1[,i], type="l", xlab ="", ylab="", sub=mat1_names[i])
}


par(mfrow=c(3,3), oma = c(0, 0, 1.5, 0))
for (i in c(73:81)){
  acf(mat1[,i], main=mat1_names[i])
}

par(mfrow=c(3,3), oma = c(0, 0, 1.5, 0))
for (i in c(1:9)){
  plot(mat1[,4*(i-3)+9], type="l", xlab ="", ylab="", sub=mat1_names[4*(i-3)+9])
}

par(mfrow=c(3,3), oma = c(0, 0, 1.5, 0))
for (i in c(1:9)){
  plot(mat1[,36+4*(i-3)+9], type="l", xlab ="", ylab="", sub=mat1_names[36+4*(i-3)+9])
}


par(mfrow=c(3,3), oma = c(0, 0, 1.5, 0))
for (i in c(1:9)){
  acf(mat1[,36+4*(i-3)+9], main=mat1_names[36+4*(i-3)+9])
}


 