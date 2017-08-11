# setwd("Heron/Adaptive")
# setwd("Heron")
# rm(list=ls())
library(rjags)
library(coda)
set.seed(1345221)


# ada=10000
# iter=500000
# th=1000
# cha=1
save_on = FALSE # TRUE

ada= 1000
iter=30000
th=1
cha=2 #2

cat("Heron HMM\n")

cat(sprintf("adapt = %i",ada),"\n")
cat(sprintf("iter = %i",iter),"\n")
cat(sprintf("thinning = %i",th),"\n")
cat(sprintf("chains = %i",cha),"\n")

# source("heron_data_HMM_approx_unifprior_norm.R")
source("heron_data_HMM_approx_unifprior_norm_fun.R")
data_10_5 <- data_fun()
data_50_40 <- data_fun(50, 40)
data_100_70 <- data_fun(100,70)


# source("heron_startingvals_HMM_approx_unifprior_norm.R")
source("heron_startingvals_HMM_approx_unifprior_norm_fun.R")
 
startingvals <- startingvals_fun(data_10_5)
inits <- startingvals$inits
params <- startingvals$params


cat("Initialise the model:\n")
tstart = proc.time()
mod_HMM_10_5 <- jags.model('heron_jags_HMM_approx_unifprior_norm.R',data_10_5,inits,n.chains=cha,n.adapt=ada)
time_init_HMM_10_5 = proc.time()-tstart

if (save_on) {
  save(mod_10_5, time_init_HMM_10_5, file = paste("Results/Heron_HMM_approx_model_ada",toString(ada),"_linux_bin1",toString(bin_size1),"_bin3",toString(bin_size3),"_unifprior.RData",sep=""))
}

cat("Run the MCMC simulations:\n")
tstart = proc.time()
output_HMM_10_5 <- coda.samples(mod_HMM_10_5,params,n.iter=iter,thin=th)
time_sample_HMM_10_5 = proc.time()-tstart
 
if (save_on) {
  save(output_HMM_10_5, time_sample_HMM_10_5, mod_HMM_10_5, time_init_HMM_10_5, file = paste("Results/Heron_HMM_approx_iter",toString(iter),"_ada",toString(ada),"_linux_bin1",toString(bin_size1),"_bin3",toString(bin_size3),"_unifprior.RData",sep=""))
}

mat1_HMM_10_5 = as.matrix(output1[1]) 
mat2_HMM_10_5 = as.matrix(output1[2]) 
mat_names_HMM_adapt <- colnames(mat1_HMM_10_5)
 


par(mfrow=c(2,1),oma=c(0,0,1.5,0))
plot(colMeans(mat1_HMM_10_5[,(1:72)]), type='l', xlab ="", ylab="", sub="X2")
plot(colMeans(mat1_HMM_10_5[,72+(1:72)]), type='l', xlab ="", ylab="", sub="X4")  
mtext("Posterior means", outer=TRUE, cex=1)


par(mfrow=c(3,4),oma=c(0,0,1.5,0))
for (i in 1:12){
  plot(mat1_HMM_10_5[,2*72+i], type='l', xlab ="", ylab="", sub=mat_names_HMM_adapt[2*72+i])
}

par(mfrow=c(4,2),oma=c(0,0,1.5,0))
for (i in 1:4){
  plot(mat1_HMM_10_5[,20*(i-1)+20], type='l', xlab ="", ylab="", sub=mat_names_HMM_adapt[20*(i-1)+20])
  plot(mat1_HMM_10_5[,72 + 20*(i-1)+20], type='l', xlab ="", ylab="", sub=mat_names_HMM_adapt[72 + 20*(i-1)+20])
}


# more bins ####
cat("Initialise the model:\n")

tstart = proc.time()
mod_HMM_50_40 <- jags.model('heron_jags_HMM_approx_unifprior_norm.R',data_50_40,inits,n.chains=cha,n.adapt=ada)
time_init_HMM_50_40 = proc.time()-tstart

if (save_on) {
  save(mod_50_40, time_init_HMM_50_40, file = paste("Results/Heron_HMM_approx_model_ada",toString(ada),"_linux_bin1",toString(bin_size1),"_bin3",toString(bin_size3),"_unifprior.RData",sep=""))
}

cat("Run the MCMC simulations:\n")
tstart = proc.time()
output_HMM_50_40 <- coda.samples(mod_HMM_50_40,params,n.iter=iter,thin=th)
time_sample_HMM_50_40 = proc.time()-tstart

if (save_on) {
  save(output_HMM_50_40, time_sample_HMM_50_40, mod_HMM_50_40, time_init_HMM_50_40, file = paste("Results/Heron_HMM_approx_iter",toString(iter),"_ada",toString(ada),"_linux_bin1",toString(bin_size1),"_bin3",toString(bin_size3),"_unifprior.RData",sep=""))
}



mat1_HMM_10_5 = as.matrix(output_10_5[1]) 
mat2_HMM_10_5 = as.matrix(output_10_5[2]) 
mat_names_HMM_adapt <- colnames(mat1_HMM_10_5)

mat1_HMM_50_40 = as.matrix(output_50_40[1]) 
mat2_HMM_50_40 = as.matrix(output_50_40[2])  

mat1_HMM_100_70 = as.matrix(output_100_70[1]) 
mat2_HMM_100_70 = as.matrix(output_100_70[2])  

# PLOT MEANS X2, X4 #############################
plotname = "Figures/heron_HMMnorm_mean_X.png"
png(filename = plotname,width = 1121, height = 797, units = "px", 
    pointsize = 16, antialias = "cleartype")

par(mfrow=c(2,1),oma=c(0,0,1.5,0))
plot(colMeans(mat1_HMM_10_5[,(1:72)]), type='l', xlab ="", ylab="", sub="X2")
lines(colMeans(mat1_HMM_50_40[,(1:72)]), type='l', col="blue")
lines(colMeans(mat1_HMM_100_70[,(1:72)]), type='l', col="red")


plot(colMeans(mat1_HMM_10_5[,72+(1:72)]), type='l', xlab ="", ylab="", sub="X4")  
lines(colMeans(mat1_HMM_50_40[,72+(1:72)]), type='l', col="blue")
lines(colMeans(mat1_HMM_100_70[,72+(1:72)]), type='l', col="red")
mtext("Posterior means", outer=TRUE, cex=1)
legend(65,2750,legend = c("10,5","50,40","100,70"), lty=1, col=c("black","blue","red"))

dev.off()


# PLOT TRACE AND ACF X2, X4 #############################

# X2
plotname = "Figures/heron_HMMnorm_trace_X2.png"
png(filename = plotname,width = 1121, height = 797, units = "px", 
    pointsize = 16, antialias = "cleartype")

par(mfrow=c(3,4))
for (i in seq(1,72,6)){
  plot(mat1_HMM_10_5[,i], type="l", xlab = "", ylab="", sub = mat_names_HMM_adapt[i])
  lines(mat1_HMM_50_40[,i], type="l", col="blue")
  lines(mat1_HMM_100_70[,i], type="l", col="red")
}  

dev.off()


par(mfrow=c(3,4))
for (i in seq(1,72,6)){
  acf(mat1_HMM_10_5[,2*72+i] ,main = mat_names_HMM_adapt[2*72+i])
}

# X4
plotname = "Figures/heron_HMMnorm_trace_X4.png"
png(filename = plotname,width = 1121, height = 797, units = "px", 
    pointsize = 16, antialias = "cleartype")

par(mfrow=c(3,4))
for (i in seq(1,72,6)){
  plot(mat1_HMM_10_5[,72+i], type="l", xlab = "", ylab="", sub = mat_names_HMM_adapt[72+i])
  lines(mat1_HMM_50_40[,72+i], type="l", col="blue")
  lines(mat1_HMM_100_70[,72+i], type="l", col="red")
}  

dev.off()


par(mfrow=c(3,4))
for (i in seq(1,72,6)){
  acf(mat1_HMM_10_5[BurnIn:iter,2*72+i] ,main = mat_names_HMM_adapt[2*72+i])
}

# PLOT TRACE AND ACF PARAM #############################


plotname = "Figures/heron_HMMnorm_trace_param.png"
png(filename = plotname,width = 1121, height = 797, units = "px", 
    pointsize = 16, antialias = "cleartype")
par(mfrow=c(3,4))
for (i in c(1:4,7:10,5,6,11,12)){
  plot(mat1_HMM_10_5[,2*72+i], type="l", xlab = "", ylab="", sub = mat_names_HMM_adapt[2*72+i])
  lines(mat1_HMM_50_40[,2*72+i], type="l", col="blue")
  lines(mat1_HMM_100_70[,2*72+i], type="l", col="red")
}  

dev.off()


par(mfrow=c(3,4))
for (i in c(1:4,7:10,5,6,11,12)){
  acf(mat1_HMM_10_5[BurnIn:iter,2*72+i] ,main = mat_names_HMM_adapt[2*72+i])
}



## PRINT TIME ########################################
time_all = matrix(c(time_init_DA[3],time_init_HMM_10_5[3],time_init_HMM_50_40[3],time_init_HMM_100_70[3],
                    time_sample_DA[3],time_sample_HMM_10_5[3],time_sample_HMM_50_40[3],time_sample_HMM_100_70[3]),
                  byrow=FALSE,ncol= 2)
rownames(time_all) = c("DA","Q1=10, Q3=5","Q1=50, Q3=40","Q1=100, Q3=70") 

colnames(time_all) = c("Init","Sample")

myFile <- "time_all.txt"
write.table(round(time_all, 3), file=myFile, row.names=TRUE, col.names=TRUE)


# ESS ##################
ESS_10_5 = lapply(output_10_5,effectiveSize)
ESS_50_40 = lapply(output_50_40,effectiveSize)
ESS_100_70 = lapply(output_100_70,effectiveSize)


ESS_HMM <- matrix(unlist(ESS_10_5), ncol = 2, byrow = FALSE)
ESS_HMM <- cbind(ESS_HMM, matrix(unlist(ESS_50_40), ncol = 2, byrow = FALSE))
ESS_HMM <- cbind(ESS_HMM, matrix(unlist(ESS_100_70), ncol = 2, byrow = FALSE))
colnames(ESS_HMM) <- c("Q1=10, Q3=5","Q1=10, Q3=5",
                       "Q1=50, Q3=40","Q1=50, Q3=40",
                       "Q1=100, Q3=70","Q1=100, Q3=70")
rownames(ESS_HMM) <- mat_names_HMM_adapt
ESS_persec_HMM <- ESS_HMM/c(time_sample_HMM_10_5[3],time_sample_HMM_10_5[3],
                            time_sample_HMM_50_40[3],time_sample_HMM_50_40[3],
                            time_sample_HMM_100_70[3],time_sample_HMM_100_70[3])
ESS_HMM_X2 = ESS_HMM[1:72,]
ESS_HMM_X4 = ESS_HMM[73:144,]
ESS_HMM_param = ESS_HMM[145:156,]

ESS_persec_HMM_X2 = ESS_persec_HMM[1:72,]
ESS_persec_HMM_X4 = ESS_persec_HMM[73:144,]
ESS_persec_HMM_param = ESS_persec_HMM[145:156,]


save(ESS_HMM, ESS_persec_HMM,
     ESS_HMM_X2, ESS_HMM_X4, ESS_HMM_param,
     ESS_persec_HMM_X2, ESS_persec_HMM_X4, ESS_persec_HMM_param,
     file = paste("Results/Heron_HMMnorm_ESS.RData",sep=""))

