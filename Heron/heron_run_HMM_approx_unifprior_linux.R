# setwd("Heron")
rm(list=ls())
library(rjags)
library(coda)
set.seed(1345221)


# ada=10000
# iter=500000
# th=1000
# cha=1
save_on = TRUE

ada= 1000
iter=1000000
th=1
cha=2 #2
# ada=1000
# iter=30000
# th=1
# cha=2

cat("Heron HMM\n")

cat(sprintf("ada = %i",ada),"\n")
cat(sprintf("iter = %i",iter),"\n")
cat(sprintf("thinning = %i",th),"\n")
cat(sprintf("chains = %i",cha),"\n")

source("heron_data_HMM_approx_unifprior.R")
source("heron_startingvals_HMM_approx_unifprior.R")

cat(sprintf("bin1 size = %i",bin_size1),"\n")
cat(sprintf("bin3 size = %i",bin_size3),"\n")

cat(sprintf("N bin min1 = %i",N_bin_min1),"\n")
cat(sprintf("N bin min3 = %i",N_bin_min3),"\n")

cat("Initialise the model:\n")
tstart = proc.time()
mod_HMM <- jags.model('heron_jags_HMM_approx_unifprior.R',data,inits,n.chains=cha,n.adapt=ada)
time_init_HMM = proc.time()-tstart

cat("Run the MCMC simulations:\n")
tstart = proc.time()
output_HMM <- coda.samples(mod_HMM,params,n.iter=iter,thin=th)
time_sample_HMM = proc.time()-tstart
 
# if (save_on) {
#   save(output1, time_sample_HMM, mod, time_init_HMM, file = paste("Results/Heron_HMM_approx_iter",toString(iter),"_ada",toString(ada),"_linux_bin1",toString(bin_size1),"_bin3",toString(bin_size3),"_unifprior.Rdata",sep=""))
# }




mat1_HMM = as.matrix(output1[[1]]) 
mat2names_HMM <- colnames(mat1_HMM)
mat_names_HMM <- colnames(mat1_HMM)

ESS_HMM = lapply(output_HMM,effectiveSize)
ESS1_HMM = as.matrix(ESS_HMM[[1]])
ESS2_HMM = as.matrix(ESS_HMM[[2]])

theta1_HMM = mat1_HMM[, 2*72 + c(1:4,7:10,5,6,11,12)]
theta2_HMM = mat2_HMM[, 2*72 + c(1:4,7:10,5,6,11,12)] 

X2_short1_HMM = mat1_HMM[, seq(1,72,6)]
X2_short2_HMM = mat2_HMM[, seq(1,72,6)]

mean_X2_1_HMM = colMeans(mat1_HMM[, c(1:72)])
mean_X2_2_HMM = colMeans(mat2_HMM[, c(1:72)])


X4_short1_HMM = mat1_HMM[, 72+seq(1,72,6)]
X4_short2_HMM = mat2_HMM[, 72+seq(1,72,6)]

mean_X4_1_HMM = colMeans(mat1_HMM[, c(73:144)])
mean_X4_2_HMM = colMeans(mat2_HMM[, c(73:144)])

# theta1_HMM = mat1_HMM[,2*48*72+ 2*72 + c(1:4,7:10,5,6,11,12)]
# theta2_HMM = mat2_HMM[,2*48*72 + 2*72 + c(1:4,7:10,5,6,11,12)] 
# 
# X2_short1_HMM = mat1_HMM[,2*48*72 + seq(1,72,6)]
# X2_short2_HMM = mat2_HMM[,2*48*72+ seq(1,72,6)]
# 
# mean_X2_1_HMM = colMeans(mat1_HMM[,2*48*72+ c(1:72)])
# mean_X2_2_HMM = colMeans(mat2_HMM[,2*48*72+ c(1:72)])
# 
# 
# X4_short1_HMM = mat1_HMM[,2*48*72 + 72+seq(1,72,6)]
# X4_short2_HMM = mat2_HMM[,2*48*72 + 72+seq(1,72,6)]
# 
# mean_X4_1_HMM = colMeans(mat1_HMM[,2*48*72 + c(73:144)])
# mean_X4_2_HMM = colMeans(mat2_HMM[,2*48*72 + c(73:144)])

if (save_on) {
  save(file = paste("Results/Heron_HMM_approx_iter",toString(iter),"_ada",toString(ada),"_bin1",toString(bin_size1),"_bin3",toString(bin_size3),"_unifprior_selected.Rdata",sep=""),
       mod_HMM, time_init_HMM, time_sample_HMM, mat_names_HMM,
       ESS1_HMM, ESS2_HMM, theta1_HMM, theta2_HMM, 
       X2_short1_HMM, X2_short2_HMM, X4_short1_HMM, X4_short2_HMM,
       mean_X2_1_HMM,mean_X2_2_HMM,mean_X4_1_HMM,mean_X4_2_HMM)
}

quit()
