# setwd("Heron/Adaptive")
# rm(list=ls())
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

cat("Heron HMM adaptive\n")

cat(sprintf("adapt = %i",ada),"\n")
cat(sprintf("iter = %i",iter),"\n")
cat(sprintf("thinning = %i",th),"\n")
cat(sprintf("chains = %i",cha),"\n")

# Read data ###
# source("heron_data_HMM_approx_unifprior_norm.R")
source("heron_data_HMM_approx_unifprior_norm_fun.R")
data_HMM_adapt <- data_fun(10, 5)
N_bin1 <- data_HMM_adapt$N_bin1
N_bin3 <- data_HMM_adapt$N_bin3



# source("heron_startingvals_HMM_approx_unifprior_norm.R")
source("heron_startingvals_HMM_approx_unifprior_norm_fun.R")

startingvals <- startingvals_fun(data_HMM_adapt)
inits <- startingvals$inits
params <- startingvals$params


##########################################
cat(paste("\n *** N_bin1 =", toString(N_bin1),", N_bin3 =", toString(N_bin3)," *** \n",sep=""))

cat("Initialise the model:\n")
tstart = proc.time()
mod_HMM_adapt <- jags.model('heron_jags_HMM_approx_unifprior_norm.R',data_HMM_adapt,
                       inits,n.chains=cha,n.adapt=ada)
time_init_HMM_adapt = proc.time()-tstart

cat("Run the MCMC simulations:\n")
tstart = proc.time()
output_HMM_adapt <- coda.samples(mod_HMM_adapt,params,n.iter=iter,thin=th)
time_sample_HMM_adapt = proc.time()-tstart

mat1_HMM_adapt = as.matrix(output_HMM_adapt[1]) 
mat2_HMM_adapt = as.matrix(output_HMM_adapt[2]) 
mat_names_HMM_adapt <- colnames(mat1_HMM_adapt)

ESS_HMM_adapt = lapply(output_HMM_adapt,effectiveSize)
ESS1_HMM_adapt = as.matrix(ESS_HMM_adapt[[1]])
ESS2_HMM_adapt = as.matrix(ESS_HMM_adapt[[2]])

theta1_HMM_adapt = mat1_HMM_adapt[,2*72 + c(1:4,7:10,5,6,11,12)]
theta2_HMM_adapt = mat2_HMM_adapt[,2*72 + c(1:4,7:10,5,6,11,12)] 

X2_short1_HMM_adapt = mat1_HMM_adapt[,seq(1,72,6)]
X2_short2_HMM_adapt = mat2_HMM_adapt[,seq(1,72,6)]

mean_X2_1_HMM_adapt = colMeans(mat1_HMM_adapt[,1:72])
mean_X2_2_HMM_adapt = colMeans(mat2_HMM_adapt[,1:72])


X4_short1_HMM_adapt = mat1_HMM_adapt[,72+seq(1,72,6)]
X4_short2_HMM_adapt = mat2_HMM_adapt[,72+seq(1,72,6)]

mean_X4_1_HMM_adapt = colMeans(mat1_HMM_adapt[,73:144])
mean_X4_2_HMM_adapt = colMeans(mat2_HMM_adapt[,73:144])

if (save_on) {
  save(file = paste("Heron_HMM_approx_iter",toString(iter),"_ada",toString(ada),"_N_q1",toString(N_bin1),"_N_q3",toString(N_bin3),"_unifprior_norm_selected.Rdata",sep=""),
       mod_HMM_adapt, time_init_HMM_adapt, time_sample_HMM_adapt, mat_names_HMM_adapt,
       ESS1_HMM_adapt, ESS2_HMM_adapt, theta1_HMM_adapt, theta2_HMM_adapt, 
       X2_short1_HMM_adapt, X2_short2_HMM_adapt, X4_short1_HMM_adapt, X4_short2_HMM_adapt,
       mean_X2_1_HMM_adapt,mean_X2_2_HMM_adapt,mean_X4_1_HMM_adapt,mean_X4_2_HMM_adapt)
}

quit()
