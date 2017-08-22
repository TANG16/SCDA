# Load required packages and fix the random seed
# setwd("BKM/Reversed")
# setwd("../Desktop/Reversed")
# install.packages(rjags)
# install.packages(coda)
# install.packages(lattice)
# rm(list=ls())

library(rjags)
library(coda)
library(lattice)
set.seed(134522)

save_on = TRUE
# MCMC details:  ####

ada=100
iter=10000
# ada=10000
# iter=10000
th=1
cha=1

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

##########################################

# Collect the results ####
mat_30_75 = as.matrix(output_HMM_30_75[1])
# mat2 = as.matrix(output1[2])
# summary(output1)
mat1_names <- colnames(mat_30_75)


par(mfrow=c(1,1),oma=c(0,0,1.5,0))
plot(colMeans(mat_30_75[,(1:36)]), type='l', xlab ="", ylab="", sub="N1")
mtext("Posterior means", outer=TRUE, cex=1)


par(mfrow=c(3,3), oma = c(0, 0, 1.5, 0))
for (i in c(37:45)){
  plot(mat_30_75[,i], type="l", xlab ="", ylab="", sub=mat1_names[i])
}


par(mfrow=c(3,3), oma = c(0, 0, 1.5, 0))
for (i in c(1:9)){
  plot(mat_30_75[,4*(i-3)+9], type="l", xlab ="", ylab="", sub=mat1_names[4*(i-3)+9])
}
##########################################

# Collect the results ####
mat_50_45 = as.matrix(output_HMM_50_45[1])
# mat2 = as.matrix(output1[2])
# summary(output1)
mat1_names <- colnames(mat_50_45)


par(mfrow=c(1,1),oma=c(0,0,1.5,0))
plot(colMeans(mat_50_45[,(1:36)]), type='l', xlab ="", ylab="", sub="N1")
mtext("Posterior means", outer=TRUE, cex=1)


par(mfrow=c(3,3), oma = c(0, 0, 1.5, 0))
for (i in c(37:45)){
  plot(mat_50_45[,i], type="l", xlab ="", ylab="", sub=mat1_names[i])
}



par(mfrow=c(3,3), oma = c(0, 0, 1.5, 0))
for (i in c(1:9)){
  plot(mat_50_45[,4*(i-3)+9], type="l", xlab ="", ylab="", sub=mat1_names[4*(i-3)+9])
}