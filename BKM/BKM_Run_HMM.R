# Load required packages and fix the random seed
setwd("BKM")
library(rjags)
library(coda)
library(lattice)
set.seed(134522)


# MCMC details: ####

ada=100
iter=1000
th=1
cha=2

# Read data ###
source("BKM_Data_HMM.R")
# Set parameters and inital values
source("BKM_StartingVals_HMM.R")

alpha1 = 1
beta1 =-2
index = alpha1 + beta1*f
phi1 = exp(index)/(1+exp(index))

alpha1 = 0.54
beta1 = -0.19
index = alpha1 + beta1*f
phi12 = exp(index)/(1+exp(index))

alphar = -2
betar = -0.7      
index = alphar + betar*stdT
rho = exp(index)

par(mfrow=c(1,2))
plot(f,type='l',col='blue')
lines(phi1,col='red')
lines(phi12,col='magenta')
lines(rho,col='green')
lines(phi1*0)
plot(Na,type='l',col='blue')


# Run the MCMC: ###
# tstart=proc.time()
mod <- jags.model('BKM_Bugs_HMM.R',data,inits,n.chains=cha,n.adapt=ada)
# temp=proc.time()-tstart
# tend <- temp

output1 <- coda.samples(mod,params,n.iter=iter,thin=th)

save(output1, file = "BKM_HMM_try_iter1000.RData")
mat1 = as.matrix(output1[1])

# Retrieve the names ####
mat1_names <- colnames(mat1) 

mat1_names[1] # "G[1,3]"
mat1_names[35*100]      # "G[100,37]"
mat1_names[35*100+1]    # "Na[1]"
mat1_names[35*100+1+35] # "Na[36]"
mat1_names[35*100+1+36] # "P[1,3]"
mat1_names[35*100+1+36+35*100-1] #"P[100,37]"
mat1_names[7036]


# Check the transition probabilities ####
Gamma_last = matrix(mat1[1000,1:(35*100)], nrow = 100, ncol = 35, byrow = FALSE)
Gamma_last_10 = diag(Gamma_last[,12])
sum(Gamma_last_10) #1
sum_Gamma_last = colSums(Gamma_last)

par(mfrow=c(2,2))
plot(matrix(mat1[100,1:(35*100)], nrow = 100, ncol = 35, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(Gamma)", sub='MCMC iter: 100')
for (t in 4:35){
  lines(matrix(mat1[100,1:(35*100)], nrow = 100, ncol = 35, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[500,1:(35*100)], nrow = 100, ncol = 35, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 500')
for (t in 4:35){
  lines(matrix(mat1[500,1:(35*100)], nrow = 100, ncol = 35, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[900,1:(35*100)], nrow = 100, ncol = 35, byrow = FALSE)[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 900')
for (t in 4:35){
  lines(matrix(mat1[900,1:(35*100)], nrow = 100, ncol = 35, byrow = FALSE)[,t],type='l')
}
 
plot(Gamma_last[,3],type='l', xlab ="",ylab="diag(Gamma)", sub='MCMC iter: 1000')
for (t in 4:35){
  lines(Gamma_last[,t],type='l')
}

# Track plots ####
par(mfrow=c(3,3))
plot(mat1[,"sigy"], type="l", xlab ="", ylab="", sub="sigy")
plot(mat1[,"alphar"], type="l", xlab ="", ylab="", sub="alpha r")
plot(mat1[,"betar"], type="l", xlab ="", ylab="", sub="beta r")

plot(mat1[,"Na[10]"], type="l", xlab ="", ylab="", sub="Na[10]")
plot(mat1[,"Na[23]"], type="l", xlab ="", ylab="", sub="Na[23]")
plot(mat1[,"Na[36]"], type="l", xlab ="", ylab="", sub="Na[36]")
# transition probability to N1=20
plot(mat1[,"G[20,10]"], type="l", xlab ="", ylab="", sub="Gamma[20,10], k=20,t=10")
plot(mat1[,"G[20,20]"], type="l", xlab ="", ylab="", sub="Gamma[20,20], k=20,t=20")
plot(mat1[,"G[20,30]"], type="l", xlab ="", ylab="", sub="Gamma[20,30], k=20,t=30")

# Refer to particular variables ####
sd(output1[[1]][,1])
sd(output1[[1]][,"Na[10]"])
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



# ESS and IF ####
hist(mat1[,20])
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