# setwd("Heron")
rm(list = ls())
library(rjags)
library(coda)
set.seed(1345221)

save_on = TRUE

# ada=10000
# iter=500000
# th=1000
# cha=1

ada=1000
iter=1000000
th=1
cha= 2


cat("Heron DA HMM init\n")

cat(sprintf("adapt = %i",ada),"\n")
cat(sprintf("iter = %i",iter),"\n")
cat(sprintf("thinning = %i",th),"\n")
cat(sprintf("chains = %i",cha),"\n")


source("heron_data.R")
# source("heron_startingvals.R")
load(file = "X2_X4_inits_from_HMM.RData")

cat("Initialise the model:\n")
tstart = proc.time()
mod_DA <- jags.model('heron_jags.R',data,inits,n.chains=cha,n.adapt=ada)
time_init_DA = proc.time()-tstart


cat("Run the MCMC simulations:\n")
tstart = proc.time()
output_DA <- coda.samples(mod_DA,params,n.iter=iter,thin=th)
time_sample_DA = proc.time()-tstart

# if (save_on) {
#   save(output_DA, time_sample_DA, mod_DA, time_init_DA,
#        file = paste("Results/Heron_DA_iter",toString(iter),"_ada",toString(ada),"_HMMinit_linux.RData",sep=""))
# }


# save selected output

mat1_DA = as.matrix(output_DA[1])
mat2_DA = as.matrix(output_DA[2])
mat_names_DA <- colnames(mat1_DA)

ESS_DA = lapply(output_DA,effectiveSize)
ESS1_DA = as.matrix(ESS_DA[[1]])
ESS2_DA = as.matrix(ESS_DA[[2]])

theta1_DA = mat1_DA[,2*72 + c(1:4,7:10,5,6,11,12)]
theta2_DA = mat2_DA[,2*72 + c(1:4,7:10,5,6,11,12)] 

X2_short1_DA = mat1_DA[,seq(1,72,6)]
X2_short2_DA = mat2_DA[,seq(1,72,6)]

mean_X2_1_DA = colMeans(mat1_DA[,1:72])
mean_X2_2_DA = colMeans(mat2_DA[,1:72])


X4_short1_DA = mat1_DA[,72+seq(1,72,6)]
X4_short2_DA = mat2_DA[,72+seq(1,72,6)]

mean_X4_1_DA = colMeans(mat1_DA[,73:144])
mean_X4_2_DA = colMeans(mat2_DA[,73:144])

if (save_on) {
  save(file = paste("Results/Heron_DA_iter",toString(iter),"_ada",toString(ada),"_HMMinit_selected.RData",sep=""),
       mod_DA, time_init_DA, time_sample_DA, mat_names_DA,
       ESS1_DA, ESS2_DA, theta1_DA, theta2_DA, 
       X2_short1_DA, X2_short2_DA, X4_short1_DA, X4_short2_DA,
       mean_X2_1_DA, mean_X2_2_DA,mean_X4_1_DA,mean_X4_2_DA)
}

# 
# if (save_on) {
#   save(output_DA, time_sample_DA, mod_DA, time_init_DA, 
#        file = paste("Results/Heron_DA_iter",toString(iter),"_ada",toString(ada),"_HMMinit_selected.RData",sep=""))
# }



#######################################################
# setwd("Heron")
rm(list = ls())
library(rjags)
library(coda)
set.seed(1345221)

save_on = TRUE

# ada=10000
# iter=500000
# th=1000
# cha=1

ada=1000
iter=1000000
th=1
cha= 2


cat("Heron DA original init\n")

cat(sprintf("adapt = %i",ada),"\n")
cat(sprintf("iter = %i",iter),"\n")
cat(sprintf("thinning = %i",th),"\n")
cat(sprintf("chains = %i",cha),"\n")


source("heron_data.R")
source("heron_startingvals.R")
# load(file = "X2_X4_inits_from_HMM.RData")

cat("Initialise the model:\n")
tstart = proc.time()
mod_DA <- jags.model('heron_jags.R',data,inits,n.chains=cha,n.adapt=ada)
time_init_DA = proc.time()-tstart


cat("Run the MCMC simulations:\n")
tstart = proc.time()
output_DA <- coda.samples(mod_DA,params,n.iter=iter,thin=th)
time_sample_DA = proc.time()-tstart

# if (save_on) {
#   save(output_DA, time_sample_DA, mod_DA, time_init_DA, 
#        file = paste("Results/Heron_DA_iter",toString(iter),"_ada",toString(ada),"_HMMinit_linux.RData",sep=""))
# }


# save selected output

mat1_DA = as.matrix(output_DA[1]) 
mat2_DA = as.matrix(output_DA[2]) 
mat_names_DA <- colnames(mat1_DA)

ESS_DA = lapply(output_DA,effectiveSize)
ESS1_DA = as.matrix(ESS_DA[[1]])
ESS2_DA = as.matrix(ESS_DA[[2]])

theta1_DA = mat1_DA[,2*72 + c(1:4,7:10,5,6,11,12)]
theta2_DA = mat2_DA[,2*72 + c(1:4,7:10,5,6,11,12)] 

X2_short1_DA = mat1_DA[,seq(1,72,6)]
X2_short2_DA = mat2_DA[,seq(1,72,6)]

mean_X2_1_DA = colMeans(mat1_DA[,1:72])
mean_X2_2_DA = colMeans(mat2_DA[,1:72])


X4_short1_DA = mat1_DA[,72+seq(1,72,6)]
X4_short2_DA = mat2_DA[,72+seq(1,72,6)]

mean_X4_1_DA = colMeans(mat1_DA[,73:144])
mean_X4_2_DA = colMeans(mat2_DA[,73:144])

if (save_on) {
  save(file = paste("Results/Heron_DA_iter",toString(iter),"_ada",toString(ada),"_selected.RData",sep=""),
       mod_DA, time_init_DA, time_sample_DA, mat_names_DA,
       ESS1_DA, ESS2_DA, theta1_DA, theta2_DA, 
       X2_short1_DA, X2_short2_DA, X4_short1_DA, X4_short2_DA,
       mean_X2_1_DA,mean_X2_2_DA,mean_X4_1_DA,mean_X4_2_DA)
}


quit()
