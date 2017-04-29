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
iter=2000
th=1
cha=2

# Load data: ###
source("Abadi_Data.R")
# Set initial parameter values: ####
source("Abadi_StartingVals.R")

# Initialise the model: ####
tstart = proc.time()
mod <- jags.model('Abadi_Bugs.R',data,inits,n.chains=cha,n.adapt=ada)
temp = proc.time()-tstart
time_init <- temp

# Run the MCMC simulations: ####
tstart = proc.time()
output1 <- coda.samples(mod,params,n.iter=iter,thin=th)
temp = proc.time()-tstart
time_sample <- temp

if (save_on) {
  save(mod, output1, time_init, time_sample, file = paste("Results/Abadi_iter",toString(iter),".RData",sep=""))
}


summary(output1)
# frame1 <- do.call(rbind.data.frame, output1)
# str(output1) # List of 2
mat1 <- as.matrix(output1[[1]])
mat2 <- as.matrix(output1[[2]])

mat1_names <- colnames(mat1) 
mat1_names[155:161]
mat1_names[1:ti] # 1st years
mat1_names[(ti+(1:ti))] # survivors 
mat1_names[(2*ti+(1:ti))] # immigrants

# output2 <- jags.samples(mod,params,n.iter=iter)
# summary(output2)
# frame2 <- do.call(rbind.data.frame, output2)
# str(output2) # List of 12



par(mfrow=c(3,3),oma=c(0,0,1.5,0))
for (ii in (1:7)){
  plot(mat1[,154+ii],type='l',col='blue')
  lines(mat2[,154+ii],col='red')
}
mtext("Regression parameters", outer=TRUE, cex=1)


par(mfrow=c(3,3),oma=c(0,0,1.5,0))
for (ii in c(1,13,26)){
  plot(mat1[,ii],type='l', col='blue', xlab ="", ylab="", sub=paste("N1[",toString(ii),"]",sep=""))
  lines(mat2[,ii],col='red')
  
  plot(mat1[,ti+ii],type='l', col='blue', xlab ="", ylab="", sub=paste("NadSurv[",toString(ii),"]",sep=""))
  lines(mat2[,ti+ii],col='red')  

  plot(mat1[,2*ti+ii],type='l', col='blue', xlab ="", ylab="", sub=paste("Nadimm[",toString(ii),"]",sep=""))
  lines(mat2[,2*ti+ii],col='red')      
}
mtext("Populations", outer=TRUE, cex=1)



mat1_names[155:161]
mat1_names[1:ti] # 1st years
mat1_names[(ti+(1:ti))] # survivors 
mat1_names[(2*ti+(1:ti))] # immigrants

# CCF_1_SURV = rep(0,ti)
# CCF_SURV_IMM = rep(0,ti)
# CCF_1_IMM = rep(0,ti)
# for (ii in (1:ti)){
#   CCF_1_SURV[ii] = sum(ccf(mat1[,ii], mat1[,ti+ii], lag.max = 40, plot = FALSE)$acf)
#   CCF_SURV_IMM[ii] = sum(ccf(mat1[,ti+ii], mat1[,2*ti+ii],lag.max = 40, plot = FALSE)$acf)
#   CCF_1_IMM[ii] = sum(ccf(mat1[,ii], mat1[,2*ti+ii], lag.max = 40, plot = FALSE)$acf)
# }
# 
# par(mfrow=c(3,1),oma=c(0,0,1.5,0))
# barplot(CCF_1_SURV,xlab ="", ylab="", sub="N1 vs NadSurv")
# barplot(CCF_SURV_IMM,xlab ="", ylab="", sub="NadSurv vs Nadimm")
# barplot(CCF_1_IMM,xlab ="", ylab="", sub="N1 vs Nadimm")
# mtext("Sum of CCF (cross-correl of two univar series)", outer=TRUE, cex=1)
