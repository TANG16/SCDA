setwd("Heron")
library(rjags)
library(coda)
set.seed(1345221)


# ada=10000
# iter=500000
# th=1000
# cha=1
save_on = TRUE

ada=100
iter=1000
th=1
cha=1


cat(sprintf("adapt = %i",ada),"\n")
cat(sprintf("iter = %i",iter),"\n")
cat(sprintf("thinning = %i",th),"\n")
cat(sprintf("chains = %i",cha),"\n")

source("heron_data_HMM.R")
source("heron_startingvals_HMM.R")

cat("Initialise the model:\n")
tstart = proc.time()
mod <- jags.model('heron_jags_HMM.R',data,inits,n.chains=cha,n.adapt=ada)
time_init_HMM = proc.time()-tstart

if (save_on) {
  save(mod, time_init_HMM, file = paste("Results/Heron_HMM_model_ada",toString(ada),".RData",sep=""))
}

cat("Run the MCMC simulations:\n")
tstart = proc.time()
output1 <- coda.samples(mod,params,n.iter=iter,thin=th)
time_sample_HMM = proc.time()-tstart
if (save_on) {
  save(output1, time_sample_HMM, file = paste("Results/Heron_HMM_iter",toString(iter),"_ada",toString(ada),".RData",sep=""))
}

summary(output1)
mat1 = as.matrix(output1[1])
mat1_names <- colnames(mat1)


tstart = proc.time()
output2 <- coda.samples(mod,params,n.iter=2*iter,thin=th)
time_sample2 = proc.time()-tstart
mat2 = as.matrix(output1[2])


# mat1_names[1:72] # "X2[1]" : "X2[72]"
# mat1_names[72+(1:72)] # "X4[1]" : "X4[72]"
# mat1_names[2*72+(1:4)] # "alpha[1]" "alpha[2]" "alpha[3]" "alpha[4]"
# mat1_names[2*72+5] # "alphal"
# mat1_names[2*72+6] # "alpharho"
# mat1_names[2*72+6+(1:4)] # "beta[1]" "beta[2]" "beta[3]" "beta[4]"
# mat1_names[2*72+11] # "betal"
# mat1_names[2*72+12] # "tauy"

mat1_names[1:80*72] # "G1[1,1]" : "G1[80,72]"
mat1_names[(80*72)+(1:(50*72))] # "G3[1,1]" : "G3[50,72]"
mat1_names[(80*72+50*72)+(1:(80*72))] # "P2[1,1]" : "P2[80,72]"
mat1_names[(80*72+50*72+80*72)+(1:(50*72))] # "P4[1,1]" : "P4[50,72]"
mat1_names[(80*72+50*72+80*72+50*72)+(1:(50*72))]  # "Q[1,1]" : "Q[50,72]"
mat1_names[(80*72+50*72+80*72+50*72+50*72)+(1:72)]  # "X2[1]" : "X2[72]"
mat1_names[(80*72+50*72+80*72+50*72+50*72+72)+(1:72)]  # "X4[1]" : "X4[72]"
mat1_names[(80*72+50*72+80*72+50*72+50*72+72+72)+(1:12)]  # "alpha[1]" : "tauy[50,72]"

# mat1_names[72+(1:72)] # "X4[1]" : "X4[72]"
# mat1_names[2*72+(1:4)] # "alpha[1]" "alpha[2]" "alpha[3]" "alpha[4]"
# mat1_names[2*72+5] # "alphal"
# mat1_names[2*72+6] # "alpharho"
# mat1_names[2*72+6+(1:4)] # "beta[1]" "beta[2]" "beta[3]" "beta[4]"
# mat1_names[2*72+11] # "betal"
# mat1_names[2*72+12] # "tauy"

# iter = 1000
par(mfrow=c(3,4))
for (i in c(1:4,7:10,5,6,11,12)){
  plot(mat1[,(80*72+50*72+80*72+50*72+50*72+72+72)+i] ,type="l", xlab = "", ylab="", sub = mat1_names[(80*72+50*72+80*72+50*72+50*72+72+72)+i])
}
par(mfrow=c(3,4))
for (i in c(1:4,7:10,5,6,11,12)){
  acf(mat1[,(80*72+50*72+80*72+50*72+50*72+72+72)+i] ,main = mat1_names[(80*72+50*72+80*72+50*72+50*72+72+72)+i])
}


# iter = 10000
# alphas
par(mfrow=c(3,4))
for (i in c(1:4,7:10,5,6,11,12)){
  plot(mat1[,2*72+i] ,type="l", xlab = "", ylab="", sub = mat1_names[2*72+i])
}
par(mfrow=c(3,4))
for (i in c(1:4,7:10,5,6,11,12)){
  acf(mat1[,2*72+i] ,main = mat1_names[2*72+i])
}

# X2
par(mfrow=c(3,4))
for (i in seq(1,72,6)){
  plot(mat1[,i] ,type="l", xlab = "", ylab="", sub = mat1_names[i])
}
# X4
par(mfrow=c(3,4))
for (i in seq(1,72,6)){
  plot(mat1[,72+i] ,type="l", xlab = "", ylab="", sub = mat1_names[72+i])
}



# X2
par(mfrow=c(3,4))
for (i in seq(1,72,6)){
  acf(mat1[,i], main = mat1_names[i])
}
# X4
par(mfrow=c(3,4))
for (i in seq(1,72,6)){
  acf(mat1[,72+i], main = mat1_names[72+i])
}


# iter = 5000
par(mfrow=c(2,1),oma=c(0,0,1.5,0))
plot(colMeans(mat1[,(80*72+50*72+80*72+50*72+50*72)+(1:72)]), type='l', xlab ="", ylab="", sub="X2")
plot(colMeans(mat1[,(80*72+50*72+80*72+50*72+50*72+72)+(1:72)]), type='l', xlab ="", ylab="", sub="X4")  
# for (i in 0:1){
#   plot(colMeans(mat1[,(80*72+50*72+80*72+50*72+50*72)+i*72+(1:72)]), type='l', xlab ="", ylab="", sub=paste("X",toString(2*i+2),sep=""))
# }
mtext("Posterior means", outer=TRUE, cex=1)



# iter = 10000
par(mfrow=c(2,1),oma=c(0,0,1.5,0))
plot(colMeans(mat2[,(1:72)]), type='l', xlab ="", ylab="", sub="X2")
plot(colMeans(mat2[,(72)+(1:72)]), type='l', xlab ="", ylab="", sub="X4")  
mtext("Posterior means", outer=TRUE, cex=1)


#ef=lapply(output,effectiveSize)
#ef[[1]][1]+ef[[2]][1]+ef[[3]][1]
#apply(rbind(ef[[1]]/tend[3],ef[[2]]/tend[3],ef[[3]]/tend[3]),2,mean) #samples per second
#summary(output)
#dens=density(c(output[[1]][,1]))
#mode=dens$x[which(dens$y==max(dens$y))]
#plot(output)
#bgr=gelman.diag(output,transform=TRUE)
#bgr$psrf
#bgr$mpsrf
#save.image("DA_snowshoehare_results.RData")



