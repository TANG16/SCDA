rm(list = ls())
library(rjags)
library(coda)
set.seed(1345221)

# ada=10000
# iter=500000
# th=1000
# cha=1

ada=100
iter=5000
th=1
cha=1

source("heron_data_scaled.R")
source("heron_startingvals_scaled.R")

tstart = proc.time()
mod <- jags.model('heron_jags_scaled.R',data,inits,n.chains=cha,n.adapt=ada)
time_init = proc.time()-tstart

tstart = proc.time()
output1 <- coda.samples(mod,params,n.iter=iter,thin=th)
time_sample = proc.time()-tstart

summary(output1)
mat1 = as.matrix(output1[1])
mat1_names <- colnames(mat1) 

tstart = proc.time()
output2 <- coda.samples(mod,params,n.iter=10*iter,thin=th)
time_sample2 = proc.time()-tstart
mat2 = as.matrix(output2[1])


mat1_names[1:72] # "X1[1]" : "X1[72]"
mat1_names[72+(1:72)] # "X2[1]" : "X2[72]"
mat1_names[2*72+(1:72)] # "X3[1]" : "X3[72]"
mat1_names[3*72+(1:72)] # "X4[1]" : "X4[72]"
mat1_names[4*72+(1:4)] # "alpha[1]" "alpha[2]" "alpha[3]" "alpha[4]"
mat1_names[4*72+5] # "alphal"
mat1_names[4*72+6] # "alpharho"
mat1_names[4*72+6+(1:4)] # "beta[1]" "beta[2]" "beta[3]" "beta[4]"
mat1_names[4*72+11] # "betal"
mat1_names[4*72+12] # "tauy"


# iter = 5000
par(mfrow=c(3,4))
for (i in c(1:4,7:10,5,6,11,12)){
  plot(mat1[,4*72+i] ,type="l", xlab = "", ylab="", sub = mat1_names[4*72+i])
}
par(mfrow=c(3,4))
for (i in c(1:4,7:10,5,6,11,12)){
  acf(mat1[,4*72+i] ,main = mat1_names[4*72+i])
}


# iter = 10000
par(mfrow=c(3,4))
for (i in c(1:4,7:10,5,6,11,12)){
  plot(mat2[,4*72+i] ,type="l", xlab = "", ylab="", sub = mat1_names[4*72+i])
}
par(mfrow=c(3,4))
for (i in c(1:4,7:10,5,6,11,12)){
  acf(mat2[,4*72+i] ,main = mat1_names[4*72+i])
}



par(mfrow=c(4,1),oma=c(0,0,1.5,0))
for (i in 0:3){
  plot(colMeans(mat2[,i*72+(1:72)]), type='l', xlab ="", ylab="", sub=paste("X",toString(i+1),sep=""))
}
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



