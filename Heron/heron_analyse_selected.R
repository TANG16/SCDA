# DA HMM inits ####
load("Heron/Results/Heron_DA_iter30000_ada1000_HMMinit_linux.RData")

mat1_DA_hmminit_short = as.matrix(output_DA[1]) 
# mat2_DA_hmminit_short = as.matrix(output_DA[2]) 
mat_names_DA = colnames(mat1_DA_hmminit_short)


theta1_DA_hmminit_short = mat1_DA_hmminit_short[,2*72 + c(1:4,7:10,5,6,11,12)]
# theta2_DA_hmminit_short = mat2_DA_hmminit_short[,2*72 + c(1:4,7:10,5,6,11,12)] 

X2_short1_DA_hmminit_short = mat1_DA_hmminit_short[,seq(1,72,6)]
# X2_short2_DA_hmminit_short = mat2_DA_hmminit_short[,seq(1,72,6)]

mean_X2_1_DA_hmminit_short = colMeans(mat1_DA_hmminit_short[,1:72])
# mean_X2_2_DA_hmminit_short = colMeans(mat2_DA_hmminit_short[,1:72])


X4_short1_DA_hmminit_short = mat1_DA_hmminit_short[,72+seq(1,72,6)]
# X4_short2_DA_hmminit_short = mat2_DA_hmminit_short[,72+seq(1,72,6)]

mean_X4_1_DA_hmminit_short = colMeans(mat1_DA_hmminit_short[,73:144])
# mean_X4_2_DA_hmminit_short = colMeans(mat2_DA_hmminit_short[,73:144])


load("Heron/Results/Heron_DA_iter1e+06_ada1000_HMMinit_selected.RData")
mod_DA_hmminit = mod_DA

time_init_DA_hmminit = time_init_DA
time_sample_DA_hmminit = time_sample_DA

X2_short1_DA_hmminit = X2_short1_DA
X2_short2_DA_hmminit = X2_short2_DA
X4_short1_DA_hmminit = X4_short1_DA
X4_short2_DA_hmminit = X4_short2_DA

mean_X2_1_DA_hmminit = mean_X2_1_DA
mean_X2_2_DA_hmminit = mean_X2_2_DA
mean_X4_1_DA_hmminit = mean_X4_1_DA
mean_X4_2_DA_hmminit = mean_X4_2_DA

ESS1_DA_hmminit = ESS1_DA 
ESS2_DA_hmminit = ESS2_DA 

theta1_DA_hmminit = theta1_DA
theta2_DA_hmminit = theta2_DA



# DA default inits ####
load("Heron/Results/Heron_DA_iter30000_ada1000_linux.RData")

mat1_DA_short = as.matrix(output_DA[1]) 
mat2_DA_short = as.matrix(output_DA[2]) 
mat_names_DA = colnames(mat1_DA_short)

theta1_DA_short = mat1_DA_short[,2*72 + c(1:4,7:10,5,6,11,12)]
theta2_DA_short = mat2_DA_short[,2*72 + c(1:4,7:10,5,6,11,12)] 

X2_short1_DA_short = mat1_DA_short[,seq(1,72,6)]
X2_short2_DA_short = mat2_DA_short[,seq(1,72,6)]

mean_X2_1_DA_short = colMeans(mat1_DA_short[,1:72])
mean_X2_2_DA_short = colMeans(mat2_DA_short[,1:72])


X4_short1_DA_short = mat1_DA_short[,72+seq(1,72,6)]
X4_short2_DA_short = mat2_DA_short[,72+seq(1,72,6)]

mean_X4_1_DA_short = colMeans(mat1_DA_short[,73:144])
mean_X4_2_DA_short = colMeans(mat2_DA_short[,73:144])


load("Heron/Results/Heron_DA_iter1e+06_ada1000_selected.RData")


# HMM adapt ####
load("Heron/Adaptive/Heron_HMM_approx_iter30000_ada1000_linux_COMB_unifprior_norm.RData")

# mat1_HMM_adapt = as.matrix(output_10_5[1]) 
# mat2_HMM_adapt = as.matrix(output_10_5[2]) 

mat1_HMM_adapt_short = as.matrix(output_100_70[1]) 
mat2_HMM_adapt_short = as.matrix(output_100_70[2]) 
 
theta1_HMM_adapt_short = mat1_HMM_adapt_short[,2*72 + c(1:4,7:10,5,6,11,12)]
theta2_HMM_adapt_short = mat2_HMM_adapt_short[,2*72 + c(1:4,7:10,5,6,11,12)] 

X2_short1_HMM_adapt_short = mat1_HMM_adapt_short[,seq(1,72,6)]
X2_short2_HMM_adapt_short = mat2_HMM_adapt_short[,seq(1,72,6)]

mean_X2_1_HMM_adapt_short = colMeans(mat1_HMM_adapt_short[,1:72])
mean_X2_2_HMM_adapt_short = colMeans(mat2_HMM_adapt_short[,1:72])


X4_short1_HMM_adapt_short = mat1_HMM_adapt_short[,72+seq(1,72,6)]
X4_short2_HMM_adapt_short = mat2_HMM_adapt_short[,72+seq(1,72,6)]

mean_X4_1_HMM_adapt_short = colMeans(mat1_HMM_adapt_short[,73:144])
mean_X4_2_HMM_adapt_short = colMeans(mat2_HMM_adapt_short[,73:144])

load("Heron/Adaptive/Heron_HMM_approx_iter1e+06_ada1000_N_q110_N_q350_unifprior_norm_selected.Rdata")



# X States: trace plots and acfs#####
# X2
# plotname = "Figures/heron_DA_trace_X2_HMMinit.png"
# png(filename = plotname,width = 1121, height = 797, units = "px", 
#     pointsize = 16, antialias = "cleartype")
par(mfrow=c(3,4))
for (i in c(1:12)){
  plot(X2_short1_DA[,i], col="blue" ,type="l", xlab = "", ylab="", sub = mat_names_DA[6*i-5])
  lines(X2_short2_DA[,i], col="green" ,type="l")  
}
# dev.off()

# X4
# plotname = "Figures/heron_DA_trace_X4_HMM_init.png"
# png(filename = plotname,width = 1121, height = 797, units = "px", 
#     pointsize = 16, antialias = "cleartype")
par(mfrow=c(3,4))
for (i in c(1:12)){
  plot(X4_short1_DA[,i], col="blue", type="l", xlab = "", ylab="", sub = mat_names_DA[6*i+72-5])
  lines(X4_short2_DA[,i], col="green" ,type="l")
}
# dev.off()


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


# Posterior means: states
par(mfrow=c(2,1))
plot(mean_X2_1_DA,type="l",col="blue")
lines(mean_X2_2_DA,type="l",col="green")
lines(mean_X2_1_DA_hmminit,type="l",col="red")
lines(mean_X2_2_DA_hmminit,type="l",col="pink")
lines(mean_X2_1_HMM_adapt,type="l",col="black")
lines(mean_X2_2_HMM_adapt,type="l",col="grey")

plot(mean_X4_1_DA,type="l",col="blue")
lines(mean_X4_2_DA,type="l",col="green")
lines(mean_X4_1_DA_hmminit,type="l",col="red")
lines(mean_X4_2_DA_hmminit,type="l",col="pink")
lines(mean_X4_1_HMM_adapt,type="l",col="black")
lines(mean_X4_2_HMM_adapt,type="l",col="grey")
legend(0, 3000, legend=c("DA","DA hmminit","HMM adapt"),  lwd=1, col=c("blue","red","black") )


# Posterior means: states
par(mfrow=c(2,1))
plot(mean_X2_1_HMM_adapt,type="l",col="black")
lines(mean_X2_2_HMM_adapt,type="l",col="grey")
lines(mean_X2_1_DA,type="l",col="blue")
lines(mean_X2_2_DA,type="l",col="green")
lines(mean_X2_1_DA_hmminit,type="l",col="red")
lines(mean_X2_2_DA_hmminit,type="l",col="pink")

plot(mean_X4_1_HMM_adapt,type="l",col="black")
lines(mean_X4_2_HMM_adapt,type="l",col="grey")
lines(mean_X4_1_DA,type="l",col="blue")
lines(mean_X4_2_DA,type="l",col="green")
lines(mean_X4_1_DA_hmminit,type="l",col="red")
lines(mean_X4_2_DA_hmminit,type="l",col="pink")
legend(55, 2500, legend=c("DA","DA hmminit","HMM adapt"),  lwd=1, col=c("blue","red","black") )



# Posterior means: long simulation (1,000,000) vs short (30,000) ####
# HMM adapt
par(mfrow=c(2,1))
plot(mean_X2_1_HMM_adapt,type="l",col="black")
lines(mean_X2_2_HMM_adapt,type="l",col="grey")
lines(mean_X2_1_HMM_adapt_short,type="l",col="blue")
lines(mean_X2_2_HMM_adapt_short,type="l",col="green") 

plot(mean_X4_1_HMM_adapt,type="l",col="black")
lines(mean_X4_2_HMM_adapt,type="l",col="grey")
lines(mean_X4_1_HMM_adapt_short,type="l",col="blue")
lines(mean_X4_2_HMM_adapt_short,type="l",col="green") 


# DA with HMM inits
par(mfrow=c(2,1))
plot(mean_X2_1_DA_hmminit,type="l",col="black")
lines(mean_X2_2_DA_hmminit,type="l",col="grey")
lines(mean_X2_1_DA_hmminit_short,type="l",col="blue")
lines(mean_X2_2_DA_hmminit_short,type="l",col="green") 

plot(mean_X4_1_DA_hmminit,type="l",col="black")
lines(mean_X4_2_DA_hmminit,type="l",col="grey")
lines(mean_X4_1_DA_hmminit_short,type="l",col="blue")
lines(mean_X4_2_DA_hmminit_short,type="l",col="green") 




# DA with defualt inits
par(mfrow=c(2,1))
plot(mean_X2_1_DA,type="l",col="black")
lines(mean_X2_2_DA,type="l",col="grey")
lines(mean_X2_1_DA_short,type="l",col="blue")
lines(mean_X2_2_DA_short,type="l",col="green") 

plot(mean_X4_1_DA,type="l",col="black")
lines(mean_X4_2_DA,type="l",col="grey")
lines(mean_X4_1_DA_short,type="l",col="blue")
lines(mean_X4_2_DA_short,type="l",col="green") 




# Posterior means: states
par(mfrow=c(2,1))
plot(mean_X2_1_HMM,type="l",col="black")
lines(mean_X2_2_HMM,type="l",col="grey")
lines(mean_X2_1_DA,type="l",col="blue")
lines(mean_X2_2_DA,type="l",col="green")
lines(mean_X2_1_DA_hmminit,type="l",col="red")
lines(mean_X2_2_DA_hmminit,type="l",col="pink")

plot(mean_X4_1_HMM,type="l",col="black")
lines(mean_X4_2_HMM,type="l",col="grey")
lines(mean_X4_1_DA,type="l",col="blue")
lines(mean_X4_2_DA,type="l",col="green")
lines(mean_X4_1_DA_hmminit,type="l",col="red")
lines(mean_X4_2_DA_hmminit,type="l",col="pink")
legend(55, 2500, legend=c("DA","DA hmminit","HMM fixed"),  lwd=1, col=c("blue","red","black") )






### Old save
load("C:/Users/aga/Desktop/Dropbox_sie_nie_miesci/Heron/Heron_HMM_approx_iter30000_ada1000_linux_bin129_bin319_unifprior.RData")

mat_names_HMM <- colnames(output1[[1]])
mat_names_HMM[(11664-72-72+1):(11664-72)]
mat_names_HMM[(11664-72+1):11664]
mat_names_HMM[11665:11676]

X2 <- output1[[1]][,(11664-72-72+1):(11664-72)]
X4 <- output1[[1]][,(11664-72+1):11664]
theta <- output1[[1]][,11665:11676]

rm(output1)

mat_names_HMM = mat_names_HMM[(11664-72-72+1):11676]

X2 = as.matrix(X2)
X4 = as.matrix(X4)
theta = as.matrix(theta)
save.image("Heron/Results/Heron_HMM_approx_iter30000_ada1000_linux_bin129_bin319_unifprior.RData")


X2_short1_HMM = X2[, seq(1,72,6)]
mean_X2_1_HMM = colMeans(X2)

X4_short1_HMM = X4[, seq(1,72,6)]
mean_X4_1_HMM = colMeans(X4)


# Posterior means: states
par(mfrow=c(2,1))
plot(mean_X2_1_HMM,type="l",col="black") 

plot(mean_X4_1_HMM,type="l",col="black")

colMeans(theta)
# alpha[1]      alpha[2]      alpha[3]      alpha[4]       
# -0.18  0.40  0.80  1.58
# alphal      alpharho       beta[1] 
#  -2.03  0.20  0.01
# beta[2]       beta[3]       beta[4]         betal          tauy 
# -0.16   -0.20 -0.18 -0.50  4.20e-06 

# trace plots params
par(mfrow=c(3,4))
for (i in c(1:4,7:10,5,11,6,12)){
  plot( theta[,i],xlab="",ylab="",type="l",col="black",sub=mat_names_HMM[72+72+i])
  
}
  
