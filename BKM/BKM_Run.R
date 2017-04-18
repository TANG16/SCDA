# Load required packages and fix the random seed
library(rjags)
library(coda)
library(lattice)
set.seed(134522)


# MCMC details: ####

ada=100
iter=1000
th=1
cha=2

# Read data
source("BKM_Data.R")
# Set parameters and inital values
source("BKM_StartingVals.R")

# Run the MCMC:
# tstart=proc.time()
mod <- jags.model('BKM_Bugs.R',data,inits,n.chains=cha,n.adapt=ada)
# temp=proc.time()-tstart
# tend <- temp

output1 <- coda.samples(mod,params,n.iter=iter,thin=th)
sd(output1[[1]][,1])
sd(output1[[1]][,"N1[1]"])
sd(output1[[1]][,100])
sd(output1[[1]][,"sigy"])
mean(output1[[1]][,"sigy"])
mean(output1[[2]][,"sigy"])
sd(output1[[2]][,"sigy"])

mean(output1[[1]][,"alpha1"])
mean(output1[[1]][,"beta1"])
mean(output1[[1]][,"alphaa"])
mean(output1[[1]][,"betaa"])
mean(output1[[1]][,"alphar"])
mean(output1[[1]][,"betar"])
mean(output1[[1]][,"alphal"])
mean(output1[[1]][,"betal"])

summary(output1)


mat1 = as.matrix(output1[1])
IF_ess <- nrow(mat1)/ apply(mat1, 2, effectiveSize)
ESS <- apply(mat1, 2, effectiveSize)

par(mfrow=c(1,2))
barplot(ESS)
hist(ESS)


par(mfrow=c(1,2))
barplot(IF)
hist(IF)


IF = matrix(data=NA,nrow=ncol(mat1),ncol=1)
for (ii in 1:ncol(mat1)){
  acf_curr = acf(mat1[,ii],lag=iter,plot=FALSE)
  acf_curr = acf_curr$acf
  IF_curr = 1 + 2*sum(acf_curr[-1])
  IF[ii] = IF_curr
}

