setwd("BKM")
rm(list=ls())

# Load required packages and fix the random seed
library(rjags)
library(coda)
library(lattice)
set.seed(134522)

save_on = FALSE
scaled_on = FALSE

# MCMC details: ####

# # ada=1000
# # iter=3000
# # th=1
# # cha=2
# ada=100
# iter=2000
# th=1
# cha=3
ada=0
iter=50000
th=1
cha=3

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
    save(mod, time_init, file = "Results/BKM_model_scaled.RData")
  } 
}else {
  tstart=proc.time()
  mod <- jags.model('BKM_Bugs.R',data,inits,n.chains=cha,n.adapt=ada)
  temp=proc.time()-tstart
  time_init <- temp # ada = 100 --> PC: 1.23; ada = 1000 --> PC: 6.29
  if (save_on) {
    save(mod, time_init, file = "Results/BKM_model.RData")
  }
}


# Run the MCMC: #### 
tstart=proc.time()
output1 <- coda.samples(mod,params,n.iter=iter,thin=th)
temp=proc.time()-tstart
time_sample <- temp # PC org:   5.23

if (save_on) {
  if (scaled_on){
    save(output1, time_sample, file = paste("Results/BKM_iter",toString(iter),"_ada",toString(ada),"_scaled.RData",sep=""))
}else {
    save(output1, time_sample, file = paste("Results/BKM_iter",toString(iter),"_ada",toString(ada),".RData",sep=""))
  } 
}

# Collect the results ####
mat1 = as.matrix(output1[1])
mat2 = as.matrix(output1[2])
mat3 = as.matrix(output1[3])
summary(output1)

mat1_names <- colnames(mat1)
colMeans(mat1[,73:80])
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



sigy <- round(c(mean(mat2[501:2000,"sigy"]),sd(mat2[501:2000,"sigy"])),4)

# alphas betas  ####
alphar <- round(c(mean(mat2[501:2000,"alphar"]),sd(mat2[501:2000,"alphar"])),4)
betar <- round(c(mean(mat2[501:2000,"betar"]),sd(mat2[501:2000,"betar"])),4) 
alphal <- round(c(mean(mat2[501:2000,"alphal"]),sd(mat2[501:2000,"alphal"])),4)
betal <- round(c(mean(mat2[501:2000,"betal"]),sd(mat2[501:2000,"betal"])),4) 
alpha1 <- round(c(mean(mat2[501:2000,"alpha1"]),sd(mat2[501:2000,"alpha1"])),4)
beta1 <- round(c(mean(mat2[501:2000,"beta1"]),sd(mat2[501:2000,"beta1"])),4) 
alphaa <- round(c(mean(mat2[501:2000,"alphaa"]),sd(mat2[501:2000,"alphaa"])),4)
betaa <- round(c(mean(mat2[501:2000,"betaa"]),sd(mat2[501:2000,"betaa"])),4) 
# Nas ####
Na3 <- round(c(mean(mat2[501:2000,"Na[3]"]),sd(mat2[501:2000,"Na[3]"])),4)
Na13 <- round(c(mean(mat2[501:2000,"Na[13]"]),sd(mat2[501:2000,"Na[13]"])),4)
Na23 <- round(c(mean(mat2[501:2000,"Na[23]"]),sd(mat2[501:2000,"Na[23]"])),4)
Na33 <- round(c(mean(mat2[501:2000,"Na[33]"]),sd(mat2[501:2000,"Na[33]"])),4)

if (scaled_on){
  DA_scaled_Resultsparams <- matrix(c(sigy,alphar,betar,alphal,betal,alpha1,beta1,alphaa,betaa,Na3,Na13,Na23,Na33),
                              ncol=2,byrow=T,dimnames=list(c("sigy","alphar","betar","alphal","betal","alpha1","beta1","alphaa","betaa", "Na3","Na13","Na23","Na33"),
                                                           c("Mean Scaled DA", "SD Scaled DA")))
  save(DA_scaled_Resultsparams, file="BKM_DA_scaled_Resultsparams.RData")
} else {
    DA_Resultsparams <- matrix(c(sigy,alphar,betar,alphal,betal,alpha1,beta1,alphaa,betaa,Na3,Na13,Na23,Na33),
                                      ncol=2,byrow=T,dimnames=list(c("sigy","alphar","betar","alphal","betal","alpha1","beta1","alphaa","betaa", "Na3","Na13","Na23","Na33"),
                                                                   c("Mean DA", "SD DA")))
    save(DA_Resultsparams, file="BKM_DA_Resultsparams.RData")
}



par(mfrow=c(2,1),oma=c(0,0,1.5,0))
plot(colMeans(mat1[,(1:T)]), type='l', xlab ="", ylab="", sub="N1")
plot(colMeans(mat1[,(T+1):(T+T)]), type='l', xlab ="", ylab="", sub="Na")
mtext("Posterior means", outer=TRUE, cex=1)


# Trace plots ####
par(mfrow=c(2,3), oma = c(0, 0, 1.5, 0))
plot(mat1[,"sigy"], type="l", xlab ="", ylab="", sub="sigy")
plot(mat1[,"alphar"], type="l", xlab ="", ylab="", sub="alpha r")
plot(mat1[,"betar"], type="l", xlab ="", ylab="", sub="beta r")

plot(mat1[,"Na[13]"], type="l", xlab ="", ylab="", sub="Na[13]")
plot(mat1[,"Na[23]"], type="l", xlab ="", ylab="", sub="Na[23]")
plot(mat1[,"Na[33]"], type="l", xlab ="", ylab="", sub="Na[33]")
mtext("Chain 1", outer=TRUE, cex=1)

# Chain 2
par(mfrow=c(2,3), oma = c(0, 0, 1.5, 0))
plot(mat2[,"sigy"], type="l", xlab ="", ylab="", sub="sigy")
plot(mat2[,"alphar"], type="l", xlab ="", ylab="", sub="alpha r")
plot(mat2[,"betar"], type="l", xlab ="", ylab="", sub="beta r")

plot(mat2[,"Na[13]"], type="l", xlab ="", ylab="", sub="Na[13]")
plot(mat2[,"Na[23]"], type="l", xlab ="", ylab="", sub="Na[23]")
plot(mat2[,"Na[33]"], type="l", xlab ="", ylab="", sub="Na[33]")
mtext("Chain 2", outer=TRUE, cex=1)

# Chain 3
par(mfrow=c(2,3), oma = c(0, 0, 1.5, 0))
plot(mat3[,"sigy"], type="l", xlab ="", ylab="", sub="sigy")
plot(mat3[,"alphar"], type="l", xlab ="", ylab="", sub="alpha r")
plot(mat3[,"betar"], type="l", xlab ="", ylab="", sub="beta r")

plot(mat3[,"Na[13]"], type="l", xlab ="", ylab="", sub="Na[13]")
plot(mat3[,"Na[23]"], type="l", xlab ="", ylab="", sub="Na[23]")
plot(mat3[,"Na[33]"], type="l", xlab ="", ylab="", sub="Na[33]")
mtext("Chain 3", outer=TRUE, cex=1)


par(mfrow=c(1,1))
plot(output1[[2]][,"sigy"])

plot(mat1[,"sigy"], type="l", xlab ="", ylab="", sub="sigy")
lines(mat2[,"sigy"], type="l", col='blue', xlab ="", ylab="", sub="sigy")

plot(mat1[,"sigy"]/100, type="l", xlab ="", ylab="", sub="sigy")
lines(mat2[,"sigy"]/100, type="l", col='blue', xlab ="", ylab="", sub="sigy")



par(mfrow=c(2,1))
acf(output1[[1]][501:2000,"Na[3]"], main="Na[3], Chain 1")
acf(output1[[2]][501:2000,"Na[3]"], main="Na[3], Chain 2")

par(mfrow=c(2,1))
acf(output1[[1]][501:2000,"Na[13]"], main="Na[13], Chain 1")
acf(output1[[2]][501:2000,"Na[13]"], main="Na[13], Chain 2")

par(mfrow=c(2,1))
acf(output1[[1]][501:2000,"Na[23]"], main="Na[23], Chain 1")
acf(output1[[2]][501:2000,"Na[23]"], main="Na[23], Chain 2")


par(mfrow=c(2,1))
acf(output1[[1]][501:2000,"Na[33]"], main="Na[33], Chain 1")
acf(output1[[2]][501:2000,"Na[33]"], main="Na[33], Chain 2")


par(mfrow=c(4,3))
acf(output1[[1]][501:2000,"Na[3]"], main="Na[3], Chain 1")
acf(output1[[2]][501:2000,"Na[3]"], main="Na[3], Chain 2")
acf(output1[[3]][501:2000,"Na[3]"], main="Na[3], Chain 3")

acf(output1[[1]][501:2000,"Na[13]"], main="Na[13], Chain 1")
acf(output1[[2]][501:2000,"Na[13]"], main="Na[13], Chain 2")
acf(output1[[3]][501:2000,"Na[13]"], main="Na[13], Chain 3")

acf(output1[[1]][501:2000,"Na[23]"], main="Na[23], Chain 1")
acf(output1[[2]][501:2000,"Na[23]"], main="Na[23], Chain 2")
acf(output1[[3]][501:2000,"Na[23]"], main="Na[23], Chain 3")

acf(output1[[1]][501:2000,"Na[33]"], main="Na[33], Chain 1")
acf(output1[[2]][501:2000,"Na[33]"], main="Na[33], Chain 2")
acf(output1[[3]][501:2000,"Na[33]"], main="Na[33], Chain 3")



par(mfrow=c(4,3))
acf(output1[[1]][501:2000,"alphar"], main="alphar, Chain 1")
acf(output1[[2]][501:2000,"alphar"], main="alphar, Chain 2")
acf(output1[[3]][501:2000,"alphar"], main="alphar, Chain 3")

acf(output1[[1]][501:2000,"betar"], main="betar, Chain 1")
acf(output1[[2]][501:2000,"betar"], main="betar, Chain 2")
acf(output1[[3]][501:2000,"betar"], main="betar, Chain 3")

acf(output1[[1]][501:2000,"alphal"], main="alphal, Chain 1")
acf(output1[[2]][501:2000,"alphal"], main="alphal, Chain 2")
acf(output1[[3]][501:2000,"alphal"], main="alphal, Chain 3")

acf(output1[[1]][501:2000,"betal"], main="betal, Chain 1")
acf(output1[[2]][501:2000,"betal"], main="betal, Chain 2")
acf(output1[[3]][501:2000,"alphal"], main="alphal, Chain 3")


par(mfrow=c(4,3))
acf(output1[[1]][501:2000,"alpha1"], main="alpha1, Chain 1")
acf(output1[[2]][501:2000,"alpha1"], main="alpha1, Chain 2")
acf(output1[[3]][501:2000,"alpha1"], main="alpha1, Chain 3")

acf(output1[[1]][501:2000,"beta1"], main="beta1, Chain 1")
acf(output1[[2]][501:2000,"beta1"], main="beta1, Chain 2")
acf(output1[[3]][501:2000,"beta1"], main="beta1, Chain 3")

acf(output1[[1]][501:2000,"alphaa"], main="alphaa, Chain 1")
acf(output1[[2]][501:2000,"alphaa"], main="alphaa, Chain 2")
acf(output1[[3]][501:2000,"alphaa"], main="alphaa, Chain 3")

acf(output1[[1]][501:2000,"betaa"], main="betaa, Chain 1")
acf(output1[[2]][501:2000,"betaa"], main="betaa, Chain 2")
acf(output1[[3]][501:2000,"betaa"], main="betaa, Chain 3")

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





par(mfrow=c(3,3), oma = c(0, 0, 1.5, 0))

plot(mat1[,"Na[1]"], type="l", xlab ="", ylab="", sub="Na[1]")
plot(mat1[,"Na[3]"], type="l", xlab ="", ylab="", sub="Na[3]")
plot(mat1[,"Na[9]"], type="l", xlab ="", ylab="", sub="Na[9]")

plot(mat1[,"Na[10]"], type="l", xlab ="", ylab="", sub="Na[10]")
plot(mat1[,"Na[13]"], type="l", xlab ="", ylab="", sub="Na[13]")
plot(mat1[,"Na[17]"], type="l", xlab ="", ylab="", sub="Na[17]")

plot(mat1[,"Na[23]"], type="l", xlab ="", ylab="", sub="Na[23]")
plot(mat1[,"Na[33]"], type="l", xlab ="", ylab="", sub="Na[33]")
plot(mat1[,"Na[36]"], type="l", xlab ="", ylab="", sub="Na[36]")




par(mfrow=c(3,3), oma = c(0, 0, 1.5, 0))

plot(mat1[,"N1[1]"], type="l", xlab ="", ylab="", sub="N1[1]")
plot(mat1[,"N1[3]"], type="l", xlab ="", ylab="", sub="N1[3]")
plot(mat1[,"N1[9]"], type="l", xlab ="", ylab="", sub="N1[9]")

plot(mat1[,"N1[10]"], type="l", xlab ="", ylab="", sub="N1[10]")
plot(mat1[,"N1[13]"], type="l", xlab ="", ylab="", sub="N1[13]")
plot(mat1[,"N1[17]"], type="l", xlab ="", ylab="", sub="N1[17]")

plot(mat1[,"N1[23]"], type="l", xlab ="", ylab="", sub="N1[23]")
plot(mat1[,"N1[33]"], type="l", xlab ="", ylab="", sub="N1[33]")
plot(mat1[,"N1[36]"], type="l", xlab ="", ylab="", sub="N1[36]")




par(mfrow=c(3,3), oma = c(0, 0, 1.5, 0))
for (i in c(73:81)){
  plot(mat1[,i], type="l", xlab ="", ylab="", sub=mat1_names[i])
}


par(mfrow=c(3,3), oma = c(0, 0, 1.5, 0))
for (i in c(1:9)){
  plot(mat1[,4*(i-3)+9], type="l", xlab ="", ylab="", sub=mat1_names[4*(i-3)+9])
}

par(mfrow=c(3,3), oma = c(0, 0, 1.5, 0))
for (i in c(1:9)){
  plot(mat1[,36+4*(i-3)+9], type="l", xlab ="", ylab="", sub=mat1_names[36+4*(i-3)+9])
}
 