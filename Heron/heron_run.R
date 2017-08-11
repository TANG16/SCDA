# setwd("Heron")
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

source("heron_data.R")
source("heron_startingvals.R")

tstart = proc.time()
mod <- jags.model('heron_jags.R',data,inits,n.chains=cha,n.adapt=ada)
time_init = proc.time()-tstart

tstart = proc.time()
output_DA <- coda.samples(mod,params,n.iter=iter,thin=th)
time_sample = proc.time()-tstart

summary(output_DA)
mat1_DA = as.matrix(output_DA[1])
mat2_DA = as.matrix(output_DA[2])
mat_names_DA <- colnames(mat1_DA) 
 

# mat_names_DA[1:72] # "X1[1]" : "X1[72]"
# mat_names_DA[72+(1:72)] # "X2[1]" : "X2[72]"
# mat_names_DA[2*72+(1:72)] # "X3[1]" : "X3[72]"
# mat_names_DA[3*72+(1:72)] # "X4[1]" : "X4[72]"
# mat_names_DA[4*72+(1:4)] # "alpha[1]" "alpha[2]" "alpha[3]" "alpha[4]"
# mat_names_DA[4*72+5] # "alphal"
# mat_names_DA[4*72+6] # "alpharho"
# mat_names_DA[4*72+6+(1:4)] # "beta[1]" "beta[2]" "beta[3]" "beta[4]"
# mat_names_DA[4*72+11] # "betal"
# mat_names_DA[4*72+12] # "tauy"

mat_names_DA[(1:72)] # "X2[1]" : "X2[72]"
mat_names_DA[72+(1:72)] # "X4[1]" : "X4[72]"
mat_names_DA[2*72+(1:4)] # "alpha[1]" "alpha[2]" "alpha[3]" "alpha[4]"
mat_names_DA[2*72+5] # "alphal"
mat_names_DA[2*72+6] # "alpharho"
mat_names_DA[2*72+6+(1:4)] # "beta[1]" "beta[2]" "beta[3]" "beta[4]"
mat_names_DA[2*72+11] # "betal"
mat_names_DA[2*72+12] # "tauy"


# iter = 10000
BurnIn = 5000;
par(mfrow=c(3,4))
for (i in c(1:4,7:10,5,6,11,12)){
  plot(mat1_DA[,2*72+i] ,type="l", xlab = "", ylab="", sub = mat_names_DA[2*72+i])
}
par(mfrow=c(3,4))
for (i in c(1:4,7:10,5,6,11,12)){
  acf(mat1_DA[BurnIn:iter,4*72+i] ,main = mat_names_DA[2*72+i])
}




#######################################################################

mat1_DA = as.matrix(output_DA[1])
mat_names_DA <- colnames(mat1_DA) 
BurnIn = 5000;


plotname = "Figures/heron_DA_trace_param_HMMinit.png"
png(filename = plotname,width = 1121, height = 797, units = "px", 
    pointsize = 16, antialias = "cleartype")
par(mfrow=c(3,4))
for (i in c(1:4,7:10,5,6,11,12)){
  plot(mat1_DA[,2*72+i] ,type="l", xlab = "", ylab="", sub = mat_names_DA[2*72+i])
}
dev.off()



par(mfrow=c(3,4))
for (i in c(1:4,7:10,5,6,11,12)){
  acf(mat1_DA[BurnIn:iter,2*72+i] ,main = mat_names_DA[2*72+i])
}


# X States: trace plots and acfs#####
# X2
plotname = "Figures/heron_DA_trace_X2_HMMinit.png"
png(filename = plotname,width = 1121, height = 797, units = "px", 
    pointsize = 16, antialias = "cleartype")
par(mfrow=c(3,4))
for (i in seq(1,72,6)){
  plot(mat1_DA[,i] ,type="l", xlab = "", ylab="", sub = mat_names_DA[i])
}
dev.off()

# X4
plotname = "Figures/heron_DA_trace_X4_HMM_init.png"
png(filename = plotname,width = 1121, height = 797, units = "px", 
    pointsize = 16, antialias = "cleartype")
par(mfrow=c(3,4))
for (i in seq(1,72,6)){
  plot(mat1_DA[,72+i] ,type="l", xlab = "", ylab="", sub = mat_names_DA[72+i])
}
dev.off()


# X2
par(mfrow=c(3,4))
for (i in seq(1,72,6)){
  acf(mat1_DA[BurnIn:iter,i], main = mat_names_DA[i])
}
# X4
par(mfrow=c(3,4))
for (i in seq(1,72,6)){
  acf(mat1_DA[BurnIn:iter,72+i], main = mat_names_DA[72+i])
}


#####

plotname = "Figures/heron_DA_mean_X_HMMinit.png"
png(filename = plotname,width = 1121, height = 797, units = "px", 
    pointsize = 16, antialias = "cleartype")

par(mfrow=c(2,1),oma=c(0,0,1.5,0))
plot(colMeans(mat1_DA[,(1:72)]), type='l', xlab ="", ylab="", sub="X2")
plot(colMeans(mat1_DA[,72+(1:72)]), type='l', xlab ="", ylab="", sub="X4")
mtext("Posterior means", outer=TRUE, cex=1)
dev.off()


# par(mfrow=c(1,1),oma=c(0,0,1.5,0))
# for (i in 0:3){
#   plot(colMeans(mat1_DA[,i*72+(1:72)]), type='l', xlab ="", ylab="", sub=paste("X",toString(i+1),sep=""))
# }
# mtext("Posterior means", outer=TRUE, cex=1)


# ESS #####
ef=lapply(output_DA,effectiveSize)
ef1 = ef[[1]]
ef2 = ef[[2]]

EF_X2 <- matrix(c(ef1[2:72]),ncol=1)

EF_X2 <- matrix(c(ef1[2:72],ef2[2:72]),ncol=2)
barplot(t(EF_X2), beside=TRUE, col=c("blue", "red"), xlab="X2")


EF_X4 <- matrix(c(ef1[72+(2:72)],ef2[72+(2:72)]),ncol=2)
barplot(t(EF_X4), beside=TRUE, col=c("blue", "red"), xlab="X4")


EF_params <- matrix(c(ef1[2*72+(1:12)],ef2[2*72+(1:12)]),ncol=2)
barplot(t(EF_params), beside=TRUE, col=c("blue", "red"), names.arg=mat_names_DA[2*72+(1:12)])


EF_params <- matrix(c(ef1[2*72+c(1:4,6:10,12)],ef2[2*72+c(1:4,6:10,12)]),ncol=2)
barplot(t(EF_params), beside=TRUE, col=c("blue", "red"), names.arg=mat_names_DA[2*72+c(1:4,6:10,12)])


ESS_DA <- matrix(unlist(ef), ncol = 2, byrow = FALSE)
colnames(ESS_DA) <- c("DA","DA")
rownames(ESS_DA) <- mat_names_DA
ESS_persec_DA <- ESS_DA/c(time_sample_DA[3],time_sample_DA[3])


ESS_DA_X2 = ESS_DA[1:72,]
ESS_DA_X4 = ESS_DA[73:144,]
ESS_DA_param = ESS_DA[145:156,]

ESS_persec_DA_X2 = ESS_persec_DA[1:72,]
ESS_persec_DA_X4 = ESS_persec_DA[73:144,]
ESS_persec_DA_param = ESS_persec_DA[145:156,]


save(ESS_DA, ESS_persec_DA,
     ESS_DA_X2, ESS_DA_X4, ESS_DA_param,
     ESS_persec_DA_X2, ESS_persec_DA_X4, ESS_persec_DA_param,     
     file = paste("Results/Heron_DA_ESS.RData",sep=""))


######
# ef[[2]][1]+ef[[3]][1]
#apply(rbind(ef[[1]]/tend[3],ef[[2]]/tend[3],ef[[3]]/tend[3]),2,mean) #samples per second
#summary(output)
#dens=density(c(output[[1]][,1]))
#mode=dens$x[which(dens$y==max(dens$y))]
#plot(output)
#bgr=gelman.diag(output,transform=TRUE)
#bgr$psrf
#bgr$mpsrf
#save.image("DA_snowshoehare_results.RData")




##### GELMAN AND RUBIN'S CONVERGENCE DIAGNOSTIC
load("C:/Users/aba228/Dropbox/Research Visit/Ruth King/Codes/Heron/Results/Heron_DA_iter30000_ada1000_HMMinit_linux.RData")
x = output_DA
load("C:/Users/aba228/Dropbox/Research Visit/Ruth King/Codes/Heron/Results/Heron_DA_iter30000_ada1000_linux.RData")
x1 = mcmc.list(as.mcmc(x), as.mcmc(output_DA[1]))
BGR_diag <- gelman.diag(x1, confidence = 0.95, transform=FALSE, autoburnin=TRUE,
            multivariate=TRUE)

plotname = "Figures/heron_DA_BGR_diag.png"
png(filename = plotname,width = 1121, height = 797, units = "px", 
    pointsize = 16, antialias = "cleartype")

par(mfrow=c(3,1))
#potential scale reduction factor
plot(BGR_diag$psrf[1:72,1],type='l',xlab="",ylab="",main="X2: BGR Point Est")
plot(BGR_diag$psrf[73:144,1],type='l',xlab="",ylab="",main="X4: BGR Point Est")
barplot(BGR_diag$psrf[145:156,1],main="Params: BGR Point Est")
dev.off()

gelman.plot(x1)

# x: An mcmc.list object with more than one chain, and with starting values that are overdispersed with respect to the posterior distribution.
# gelman.plot(x1, bin.width = 10, max.bins = 50,
            confidence = 0.95, transform = FALSE, autoburnin=TRUE, auto.layout = TRUE)