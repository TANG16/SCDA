# Load required packages and fix the random seed
rm(list=ls())

setwd("BKM")
library(rjags)
library(coda)
library(lattice)
set.seed(134522)

save_on = FALSE
# MCMC details: ####

ada=100
iter=1000
th=1
cha=3

# Read data ###
source("BKM_Data_HMM.R")
# Set parameters and inital values
source("BKM_StartingVals_HMM.R")

if (1==0){
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
}

# Compile the model: ####
tstart=proc.time()
mod <- jags.model('BKM_Bugs_HMM.R',data,inits,n.chains=cha,n.adapt=ada)
temp=proc.time()-tstart
time_HMM_init <- temp # ada=100 PC: 357.39 ~ 6min --> 188.28 ~3 min with 0!
# ada = 1000 --> PC: 1798.68 
# ada = 100 --> laptop: 659
if (save_on) {
  save(mod, time_HMM_init, file = paste("BKM_HMM_model_ada",toString(ada),"_laptop.RData",sep=""))
}


# Run the MCMC: #### 
tstart=proc.time()
output1 <- coda.samples(mod,params,n.iter=iter,thin=th)
temp=proc.time()-tstart
time_HMM_sample <- temp # PC:   1843.02 ~31 min
# ada = 1000 --> 1741.71
# ada = 100 --> laptop: 6449 ~ 108 min
if (save_on) {
  save(output1, time_HMM_sample, file =  paste("BKM_HMM_try_iter",toString(iter),"_ada",toString(ada),"_laptop.RData",sep=""))
}


# Collect the results ####
mat1 = as.matrix(output1[1])
mat2 = as.matrix(output1[2])


# Retrieve the names ####
mat1_names <- colnames(mat1) 

mat1_names[1] # "G[1,3]"
mat1_names[35*100]      # "G[100,37]"
mat1_names[35*100+1]    # "Na[1]"
mat1_names[35*100+1+35] # "Na[36]"
mat1_names[35*100+1+36] # "P[1,3]"
mat1_names[35*100+1+36+35*100-1] #"P[100,37]"
mat1_names[7036]
(35*100+1+360):(35*100+1+36+35*100-1)


# PLOTS ####
source("BKM_HMM_plots.R")

# Check the transition probabilities ####
Gamma_last = matrix(mat1[1000,1:(35*100)], nrow = 100, ncol = 35, byrow = FALSE)
Gamma_last_10 = diag(Gamma_last[,12])
sum(Gamma_last_10) #1
sum_Gamma_last = colSums(Gamma_last)

###  plots and acfs
par(mfrow=c(1,1))
plot(output1[[1]][,"sigy"])
par(mfrow=c(1,1))
acf(output1[[2]][,"alphar"])
par(mfrow=c(1,1))
acf(output1[[2]][,"betar"])
par(mfrow=c(1,1))
acf(output1[[2]][,"Na[36]"])
par(mfrow=c(1,1))
acf(output1[[2]][,"Na[1]"])

par(mfrow=c(1,1))
acf(output1[[1]][,"Na[10]"])


# Refer to particular variables ####
sd(output1[[1]][,1])
sd(output1[[1]][,"Na[10]"])
sd(output1[[1]][,100])
sd(output1[[1]][,"sigy"])




### FULL DA posterior means and std####
# alpha1     0.5540 6.888e-02 
# alphaa     1.5678 6.338e-02 
# alphal    -4.5760 3.534e-02 
# alphar    -1.1723 6.579e-02 
# beta1     -0.1913 5.668e-02 
# betaa     -0.2465 3.786e-02 
# betal     -0.3650 3.981e-02 
# betar     -0.3360 3.254e-02 

# Na[3]    993.9475 26.52
# Na[13]  1789.2550 52.95
# Na[23]  1469.4150 49.46
# Na[33]   975.9720 59.51

### FULL DA SCALED adapt 1000 iter 3000

# Na[3]  104.7057, 5.74023, 
# Na[13] 170.9540, 6.73734, 
# Na[23] 158.9238,  7.29764, 
# Na[33]  94.5547,  5.72801, 
# 
# alpha1   0.5261,  0.06790, 
# alphaa   1.5357,  0.06691,
# alphal  -2.2752,  0.03572, 
# alphar  -1.1330,  0.09853, 
# beta1   -0.2034,  0.05903,
# betaa   -0.2501,  0.03805, 
# betal   -0.3564,  0.04183, 
# betar   -0.3061,  0.07159, 
# sigy   157.1837, 57.75430, 

# DA WITH HMM MEANS STD #####
# sigy #### 
sigy <- round(c(30755.8259, 8700, 
                157.1837, 57.75430, 
                mean(mat1[,"sigy"]),sd(mat1[,"sigy"]),
                mean(mat2[,"sigy"]),sd(mat2[,"sigy"]),
                mean(mat3[,"sigy"]),sd(mat3[,"sigy"])),4)

# alphas betas  ####
alphar <- round(c(-1.1723, 0.06579, 
                  -1.1330, 0.09853, 
                mean(mat1[,"alphar"]),sd(mat1[,"alphar"]),
                mean(mat2[,"alphar"]),sd(mat2[,"alphar"]),
                mean(mat3[,"alphar"]),sd(mat3[,"alphar"])),4)
betar <- round(c(-0.3360, 0.03254, 
                 -0.3061, 0.07159, 
                  mean(mat1[,"betar"]),sd(mat1[,"betar"]),
                  mean(mat2[,"betar"]),sd(mat2[,"betar"]),
                  mean(mat3[,"betar"]),sd(mat3[,"betar"])),4) 

alphal <- round(c(-4.5760, 0.3534, 
                  -2.2752, 0.03572, 
                  mean(mat1[,"alphal"]),sd(mat1[,"alphal"]),
                  mean(mat2[,"alphal"]),sd(mat2[,"alphal"]),
                  mean(mat3[,"alphal"]),sd(mat3[,"alphal"])),4)
betal <- round(c(-0.3650, 0.03981,
                 -0.3564,  0.04183,                 
                 mean(mat1[,"betal"]),sd(mat1[,"betal"]),
                 mean(mat2[,"betal"]),sd(mat2[,"betal"]),
                 mean(mat3[,"betal"]),sd(mat3[,"betal"])),4) 

alpha1 <- round(c(0.5540, 0.06888,
                  0.5261,  0.06790,                  
                  mean(mat1[,"alpha1"]),sd(mat1[,"alpha1"]),
                  mean(mat2[,"alpha1"]),sd(mat2[,"alpha1"]),
                  mean(mat3[,"alpha1"]),sd(mat3[,"alpha1"])),4)
beta1 <- round(c(-0.1913, 0.05668, 
                 -0.2034,  0.05903,
                 mean(mat1[,"beta1"]),sd(mat1[,"beta1"]),
                 mean(mat2[,"beta1"]),sd(mat2[,"beta1"]),
                 mean(mat3[,"beta1"]),sd(mat3[,"beta1"])),4) 

alphaa <- round(c(1.5678, 0.06338, 
                  1.5357,  0.06691,                  
                  mean(mat1[,"alphaa"]),sd(mat1[,"alphaa"]),
                  mean(mat2[,"alphaa"]),sd(mat2[,"alphaa"]),
                  mean(mat3[,"alphaa"]),sd(mat3[,"alphaa"])),4)
betaa <- round(c(-0.2465, 0.03786, 
                 -0.2501,  0.03805,
                 mean(mat1[,"betaa"]),sd(mat1[,"betaa"]),
                 mean(mat2[,"betaa"]),sd(mat2[,"betaa"]),
                 mean(mat3[,"betaa"]),sd(mat3[,"betaa"])),4) 
# Nas ####
Na3 <- round(c(993.9475, 26.52, 
               104.7057,  5.74023,
                  mean(mat1[,"Na[3]"]),sd(mat1[,"Na[3]"]),
                  mean(mat2[,"Na[3]"]),sd(mat2[,"Na[3]"]),
                  mean(mat3[,"Na[3]"]),sd(mat3[,"Na[3]"])),4)

Na13 <- round(c(1789.2550, 52.95, 
                170.9540, 6.73734,
                  mean(mat1[,"Na[13]"]),sd(mat1[,"Na[13]"]),
                  mean(mat2[,"Na[13]"]),sd(mat2[,"Na[13]"]),
                  mean(mat3[,"Na[13]"]),sd(mat3[,"Na[13]"])),4)

Na23 <- round(c(1469.4150, 49.46, 
                158.9238,  7.29764,
                  mean(mat1[,"Na[23]"]),sd(mat1[,"Na[23]"]),
                  mean(mat2[,"Na[23]"]),sd(mat2[,"Na[23]"]),
                  mean(mat3[,"Na[23]"]),sd(mat3[,"Na[23]"])),4)

Na33 <- round(c(975.9720, 59.51, 
                94.5547,  5.72801, 
                mean(mat1[,"Na[33]"]),sd(mat1[,"Na[33]"]),
                mean(mat2[,"Na[33]"]),sd(mat2[,"Na[33]"]),
                mean(mat3[,"Na[33]"]),sd(mat3[,"Na[33]"])),4)


HMM_Resultsparams <- matrix(c(sigy,alphar,betar,alphal,betal,alpha1,beta1,alphaa,betaa,Na3,Na13,Na23,Na33),
       ncol=8,byrow=T,dimnames=list(c("sigy","alphar","betar","alphal","betal","alpha1","beta1","alphaa","betaa",
                                      "Na3","Na13","Na23","Na33"),
                                    c("Mean Full DA", "SD Full DA","Mean Scaled DA", "SD Scaled DA", "Mean Chain1","SD Chain1","Mean Chain2","SD Chain2","Mean Chain3","SD Chain3")))



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