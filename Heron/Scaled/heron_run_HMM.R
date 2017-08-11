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
mat2 = as.matrix(output1[2])
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
BurnIn = 5000;
# alphas
par(mfrow=c(3,4))
for (i in c(1:4,7:10,5,6,11,12)){
  plot(mat1[,2*72+i] ,type="l", xlab = "", ylab="", sub = mat1_names[2*72+i])
}
par(mfrow=c(3,4))
for (i in c(1:4,7:10,5,6,11,12)){
  acf(mat1[BurnIn:iter,2*72+i] ,main = mat1_names[2*72+i])
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
  acf(mat1[BurnIn:iter,i], main = mat1_names[i])
}
# X4
par(mfrow=c(3,4))
for (i in seq(1,72,6)){
  acf(mat1[BurnIn:iter,72+i], main = mat1_names[72+i])
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



######
range_G1 <- c(1:((N_bin1-N_bin_min1+1)*T))
range_G3 <- c(((N_bin1-N_bin_min1+1)*T) + 1:((N_bin3-N_bin_min3+1)*T))
range_P2 <- c(((N_bin1-N_bin_min1+1)*T) + (N_bin3-N_bin_min3+1)*T + 1:((N_bin3-N_bin_min3+1)*T))
range_P4 <- c(((N_bin1-N_bin_min1+1)*2*T) + (N_bin3-N_bin_min3+1)*T + 1:((N_bin3-N_bin_min3+1)*T))
range_Q <- c(((N_bin1-N_bin_min1+1)*2*T) + (N_bin3-N_bin_min3+1)*2*T + 1:((N_bin3-N_bin_min3+1)*T))

range_X2 <- c(((N_bin1-N_bin_min1+1)*2*T) + (N_bin3-N_bin_min3+1)*3*T + 1:T)
range_X4 <- c(((N_bin1-N_bin_min1+1)*2*T) + (N_bin3-N_bin_min3+1)*3*T + T + 1:T)
range_param <- c(((N_bin1-N_bin_min1+1)*2*T) + (N_bin3-N_bin_min3+1)*3*T + 2*T + 1:12)


mat1_names[range_G1] #"G1[1,1]" : "G1[50,72]"
mat1_names[range_G3] #"G3[1,1]" : "G3[50,72]"

mat1_names[range_P2] # "P2[1,1]" :"P4[50,72]"
mat1_names[range_P4] # "P4[1,1]" :"P4[50,72]"
mat1_names[range_Q]   # "Q[1,1]" :"Q[50,72]"
mat1_names[range_X2]  #  "X2[1]" : "X2[72]"
mat1_names[range_X4]  #  "X4[1]" : "X4[72]"
mat1_names[range_param] # "alpha[1]" "alpha[2]" "alpha[3]" "alpha[4]" "alphal"  "alpharho" "beta[1]"  "beta[2]"  "beta[3]"  "beta[4]"  "betal"    "tauy"



# Gamma1  ####
par(mfrow=c(2,2), oma = c(0, 0, 1.5, 0))
plot(matrix(mat1[1000,range_G1], nrow = (N_bin1-N_bin_min1+1), ncol = T, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(Gamma1)", sub='MCMC iter: 1000')
for (t in 2:T){
  lines(matrix(mat1[1000,range_G1], nrow = (N_bin1-N_bin_min1+1), ncol = T, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[5000,range_G1], nrow = (N_bin1-N_bin_min1+1), ncol = T, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(Gamma1)", sub='MCMC iter: 5000')
for (t in 2:T){
  lines(matrix(mat1[5000,range_G1], nrow = (N_bin1-N_bin_min1+1), ncol = T, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[7500,range_G1], nrow = (N_bin1-N_bin_min1+1), ncol = T, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(Gamma1)", sub='MCMC iter: 7500')
for (t in 2:T){
  lines(matrix(mat1[7500,range_G1], nrow = (N_bin1-N_bin_min1+1), ncol = T, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[10000,range_G1], nrow = (N_bin1-N_bin_min1+1), ncol = T, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(Gamma1)", sub='MCMC iter: 10000')
for (t in 2:T){
  lines(matrix(mat1[10000,range_G1], nrow = (N_bin1-N_bin_min1+1), ncol = T, byrow = FALSE)[,t],type='l')
}


# P4  ####
par(mfrow=c(2,2), oma = c(0, 0, 1.5, 0))
plot(matrix(mat1[1000,range_P2], nrow = (N_bin1-N_bin_min1+1), ncol = T, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(P2)", sub='MCMC iter: 1000')
for (t in 2:T){
  lines(matrix(mat1[1000,range_P2], nrow = (N_bin1-N_bin_min1+1), ncol = T, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[5000,range_P2], nrow = (N_bin1-N_bin_min1+1), ncol = T, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(P2)", sub='MCMC iter: 5000')
for (t in 2:T){
  lines(matrix(mat1[5000,range_P2], nrow = (N_bin1-N_bin_min1+1), ncol = T, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[7500,range_P2], nrow = (N_bin1-N_bin_min1+1), ncol = T, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(P2)", sub='MCMC iter: 7500')
for (t in 2:T){
  lines(matrix(mat1[7500,range_P2], nrow = (N_bin1-N_bin_min1+1), ncol = T, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[10000,range_P2], nrow = (N_bin1-N_bin_min1+1), ncol = T, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(P2)", sub='MCMC iter: 10000')
for (t in 2:T){
  lines(matrix(mat1[10000,range_P2], nrow = (N_bin1-N_bin_min1+1), ncol = T, byrow = FALSE)[,t],type='l')
}

mtext("P2", outer=TRUE, cex=1)



# Gamma3  ####
par(mfrow=c(2,2), oma = c(0, 0, 1.5, 0))
plot(matrix(mat1[1000,range_G3], nrow = (N_bin3-N_bin_min3+1), ncol = T, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(Gamma3)", sub='MCMC iter: 1000')
for (t in 2:T){
  lines(matrix(mat1[1000,range_G3], nrow = (N_bin3-N_bin_min3+1), ncol = T, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[5000,range_G3], nrow = (N_bin3-N_bin_min3+1), ncol = T, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(Gamma3)", sub='MCMC iter: 5000')
for (t in 2:T){
  lines(matrix(mat1[5000,range_G3], nrow = (N_bin3-N_bin_min3+1), ncol = T, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[7500,range_G3], nrow = (N_bin3-N_bin_min3+1), ncol = T, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(Gamma3)", sub='MCMC iter: 7500')
for (t in 2:T){
  lines(matrix(mat1[7500,range_G3], nrow = (N_bin3-N_bin_min3+1), ncol = T, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[10000,range_G3], nrow = (N_bin3-N_bin_min3+1), ncol = T, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(Gamma3)", sub='MCMC iter: 10000')
for (t in 2:T){
  lines(matrix(mat1[10000,range_G3], nrow = (N_bin3-N_bin_min3+1), ncol = T, byrow = FALSE)[,t],type='l')
}




# P4  ####
par(mfrow=c(2,2), oma = c(0, 0, 1.5, 0))
plot(matrix(mat1[1000,range_P4], nrow = (N_bin3-N_bin_min3+1), ncol = T, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(P4)", sub='MCMC iter: 1000')
for (t in 2:T){
  lines(matrix(mat1[1000,range_P4], nrow = (N_bin3-N_bin_min3+1), ncol = T, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[5000,range_P4], nrow = (N_bin3-N_bin_min3+1), ncol = T, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(P4)", sub='MCMC iter: 5000')
for (t in 2:T){
  lines(matrix(mat1[5000,range_P4], nrow = (N_bin3-N_bin_min3+1), ncol = T, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[7500,range_P4], nrow = (N_bin3-N_bin_min3+1), ncol = T, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(P4)", sub='MCMC iter: 7500')
for (t in 2:T){
  lines(matrix(mat1[7500,range_P4], nrow = (N_bin3-N_bin_min3+1), ncol = T, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[10000,range_P4], nrow = (N_bin3-N_bin_min3+1), ncol = T, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(P4)", sub='MCMC iter: 10000')
for (t in 2:T){
  lines(matrix(mat1[10000,range_P4], nrow = (N_bin3-N_bin_min3+1), ncol = T, byrow = FALSE)[,t],type='l')
}

mtext("P4", outer=TRUE, cex=1)



# Q  ####
par(mfrow=c(2,2), oma = c(0, 0, 1.5, 0))
plot(matrix(mat1[1000,range_Q], nrow = (N_bin3-N_bin_min3+1), ncol = T, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(Q)", sub='MCMC iter: 1000')
for (t in 2:T){
  lines(matrix(mat1[1000,range_Q], nrow = (N_bin3-N_bin_min3+1), ncol = T, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[5000,range_Q], nrow = (N_bin3-N_bin_min3+1), ncol = T, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(Q)", sub='MCMC iter: 5000')
for (t in 2:T){
  lines(matrix(mat1[5000,range_Q], nrow = (N_bin3-N_bin_min3+1), ncol = T, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[7500,range_Q], nrow = (N_bin3-N_bin_min3+1), ncol = T, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(Q)", sub='MCMC iter: 7500')
for (t in 2:T){
  lines(matrix(mat1[7500,range_Q], nrow = (N_bin3-N_bin_min3+1), ncol = T, byrow = FALSE)[,t],type='l')
}

plot(matrix(mat1[10000,range_Q], nrow = (N_bin3-N_bin_min3+1), ncol = T, byrow = FALSE)[,3],type='l', xlab ="", ylab="diag(Q)", sub='MCMC iter: 10000')
for (t in 2:T){
  lines(matrix(mat1[10000,range_Q], nrow = (N_bin3-N_bin_min3+1), ncol = T, byrow = FALSE)[,t],type='l')
}

mtext("Q", outer=TRUE, cex=1)


######

range_G1 <- c(1:((N_bin1-N_bin_min1+1)*T))
G1 = as.matrix(output1[20000:30000,range_G1])

par(mfrow=c(2,2), oma = c(0, 0, 1.5, 0))
for(ii in c(1,5001,7501,10001)){
  plot(matrix(G1[[1]][ii,], nrow = (N_bin1-N_bin_min1+1), ncol = T, byrow = FALSE)[,1],type='l', xlab ="", ylab="diag(Gamma1)", sub=paste('MCMC iter:',toString(ii+19999),seq=""))
  for (t in 2:T){
    lines(matrix(G1[[1]][ii,], nrow = (N_bin1-N_bin_min1+1), ncol = T, byrow = FALSE)[,t],type='l')
  }
}




range_P2 <- c((N_bin1-N_bin_min1+1)*T + (1:((N_bin1-N_bin_min1+1)*T)))
P2 = as.matrix(output1[20000:30000,range_P2])

par(mfrow=c(2,2), oma = c(0, 0, 1.5, 0))
for(ii in c(1,5001,7501,10001)){
  plot(matrix(P2[[1]][ii,], nrow = (N_bin1-N_bin_min1+1), ncol = T, byrow = FALSE)[,1],type='l', xlab ="", ylab="diag(P2)", sub=paste('MCMC iter:',toString(ii+19999),seq=""))
  for (t in 2:T){
    lines(matrix(P2[[1]][ii,], nrow = (N_bin1-N_bin_min1+1), ncol = T, byrow = FALSE)[,t],type='l')
  }
}
