setwd('SV')
library(R.matlab)
library(coda)
# data <- readMat('SV_results.mat')
data_DA <- readMat('SV_DA_IF.mat')
M = 10000

H_DA_init <- as.mcmc(data_DA$H.DA.init)
IF_DA_init <- as.mcmc(data_DA$IF.DA[1:100])
ESS_DA_init <- M/IF_DA_init
  
  
# EFF = lapply(H_DA_init,effectiveSize)
# EFF = effectiveSize(H_DA_init[,1])
EFF <- rep(NaN,100)
for (i in c(1:100)){
  EFF[i] <- effectiveSize(H_DA_init[,i])
}
IFF <- M/EFF