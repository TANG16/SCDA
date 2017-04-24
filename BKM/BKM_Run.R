setwd("BKM")
rm(list=ls())

# Load required packages and fix the random seed
library(rjags)
library(coda)
library(lattice)
set.seed(134522)

save_on = FALSE
scaled_on = TRUE

# MCMC details: ####

# ada=1000
# iter=3000
# th=1
# cha=2
ada=1000
iter=3000
th=1
cha=2


# Read data
if (scaled_on){
  source("BKM_Data_scaled.R")
  # Set parameters and inital values
  source("BKM_StartingVals_scaled.R") 
}else {
  source("BKM_Data.R")
  # Set parameters and inital values
  source("BKM_StartingVals.R")
}


# Compile the model: ####
if (scaled_on){
  tstart=proc.time()
  mod <- jags.model('BKM_Bugs_scaled.R',data,inits,n.chains=cha,n.adapt=ada)
  temp=proc.time()-tstart
  time_init <- temp 
  if (save_on) {
    save(mod, time_init, file = "BKM_model_scaled.RData")
  } 
}else {
  tstart=proc.time()
  mod <- jags.model('BKM_Bugs.R',data,inits,n.chains=cha,n.adapt=ada)
  temp=proc.time()-tstart
  time_init <- temp # ada = 100 --> PC: 1.23; ada = 1000 --> PC: 6.29
  if (save_on) {
    save(mod, time_init, file = "BKM_model.RData")
  }
}


# Run the MCMC: #### 
tstart=proc.time()
output1 <- coda.samples(mod,params,n.iter=iter,thin=th)
temp=proc.time()-tstart
time_sample <- temp # PC org:   5.23

if (save_on) {
  if (scaled_on){
    save(output1, time_sample, file = paste("BKM_iter",toString(iter),"_ada",toString(ada),"_scaled.RData",sep=""))
}else {
    save(output1, time_sample, file = paste("BKM_iter",toString(iter),"_ada",toString(ada),".RData",sep=""))

  } 
}

# Collect the results ####
mat1 = as.matrix(output1[1])
mat2 = as.matrix(output1[2])
summary(output1)

# Refer to particular variables ####
sd(output1[[1]][,1])
mean(output1[[1]][,"Na[1]"])
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


# Trace plots ####
par(mfrow=c(2,3))
plot(mat1[,"sigy"], type="l", xlab ="", ylab="", sub="sigy")
plot(mat1[,"alphar"], type="l", xlab ="", ylab="", sub="alpha r")
plot(mat1[,"betar"], type="l", xlab ="", ylab="", sub="beta r")

plot(mat1[,"Na[10]"], type="l", xlab ="", ylab="", sub="Na[10]")
plot(mat1[,"Na[23]"], type="l", xlab ="", ylab="", sub="Na[23]")
plot(mat1[,"Na[36]"], type="l", xlab ="", ylab="", sub="Na[36]")


par(mfrow=c(2,3))
plot(mat2[,"sigy"], type="l", xlab ="", ylab="", sub="sigy")
plot(mat2[,"alphar"], type="l", xlab ="", ylab="", sub="alpha r")
plot(mat2[,"betar"], type="l", xlab ="", ylab="", sub="beta r")

plot(mat2[,"Na[10]"], type="l", xlab ="", ylab="", sub="Na[10]")
plot(mat2[,"Na[23]"], type="l", xlab ="", ylab="", sub="Na[23]")
plot(mat2[,"Na[36]"], type="l", xlab ="", ylab="", sub="Na[36]")


par(mfrow=c(1,1))
plot(output1[[2]][,"sigy"])

plot(mat1[,"sigy"], type="l", xlab ="", ylab="", sub="sigy")
lines(mat2[,"sigy"], type="l", col='blue', xlab ="", ylab="", sub="sigy")

plot(mat1[,"sigy"]/100, type="l", xlab ="", ylab="", sub="sigy")
lines(mat2[,"sigy"]/100, type="l", col='blue', xlab ="", ylab="", sub="sigy")



par(mfrow=c(2,1))
acf(output1[[1]][,"Na[3]"], main="Na[3], Chain 1")
acf(output1[[2]][,"Na[3]"], main="Na[3], Chain 2")

par(mfrow=c(2,1))
acf(output1[[1]][,"Na[13]"], main="Na[13], Chain 1")
acf(output1[[2]][,"Na[13]"], main="Na[13], Chain 2")

par(mfrow=c(2,1))
acf(output1[[1]][,"Na[23]"], main="Na[23], Chain 1")
acf(output1[[2]][,"Na[23]"], main="Na[23], Chain 2")


par(mfrow=c(2,1))
acf(output1[[1]][,"Na[33]"], main="Na[33], Chain 1")
acf(output1[[2]][,"Na[33]"], main="Na[33], Chain 2")


par(mfrow=c(4,2))
acf(output1[[1]][,"Na[3]"], main="Na[3], Chain 1")
acf(output1[[2]][,"Na[3]"], main="Na[3], Chain 2")

acf(output1[[1]][,"Na[13]"], main="Na[13], Chain 1")
acf(output1[[2]][,"Na[13]"], main="Na[13], Chain 2")

acf(output1[[1]][,"Na[23]"], main="Na[23], Chain 1")
acf(output1[[2]][,"Na[23]"], main="Na[23], Chain 2")

acf(output1[[1]][,"Na[33]"], main="Na[33], Chain 1")
acf(output1[[2]][,"Na[33]"], main="Na[33], Chain 2")



par(mfrow=c(4,2))
acf(output1[[1]][,"alphar"], main="alphar, Chain 1")
acf(output1[[2]][,"alphar"], main="alphar, Chain 2")

acf(output1[[1]][,"betar"], main="betar, Chain 1")
acf(output1[[2]][,"betar"], main="betar, Chain 2")

acf(output1[[1]][,"alphal"], main="alphal, Chain 1")
acf(output1[[2]][,"alphal"], main="alphal, Chain 2")

acf(output1[[1]][,"betal"], main="betal, Chain 1")
acf(output1[[2]][,"betal"], main="betal, Chain 2")


par(mfrow=c(4,2))
acf(output1[[1]][,"alpha1"], main="alpha1, Chain 1")
acf(output1[[2]][,"alpha1"], main="alpha1, Chain 2")

acf(output1[[1]][,"beta1"], main="beta1, Chain 1")
acf(output1[[2]][,"beta1"], main="beta1, Chain 2")

acf(output1[[1]][,"alphaa"], main="alphaa, Chain 1")
acf(output1[[2]][,"alphaa"], main="alphaa, Chain 2")

acf(output1[[1]][,"betaa"], main="betaa, Chain 1")
acf(output1[[2]][,"betaa"], main="betaa, Chain 2")

# IF and ESS ####
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