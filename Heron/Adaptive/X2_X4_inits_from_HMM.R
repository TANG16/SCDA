# setwd("Heron")
load("Adaptive/Heron_HMM_approx_iter30000_ada1000_linux_COMB_unifprior_norm.RData")


T = 72

par(mfrow=c(2,1),oma=c(0,0,1.5,0))
plot(colMeans(mat1_HMM_100_70[,(1:T)]), type='l', xlab ="", ylab="", sub="X2")
lines(colMeans(mat1_HMM_100_70[15000:30000,(1:T)]), type='l',col='red')

plot(colMeans(mat2_HMM_100_70[,T+(1:T)]), type='l', xlab ="", ylab="", sub="X4")  
lines(colMeans(mat1_HMM_100_70[15000:30000,T+(1:T)]), type='l',col='red')
mtext("Posterior means", outer=TRUE, cex=1)
legend(65,2750,legend = c("10,5","50,40","100,70"), lty=1, col=c("black","blue","red"))



mat1_HMM_100_70 = as.matrix(output_100_70[1]) 
mat2_HMM_100_70 = as.matrix(output_100_70[2])  
mat_names_HMM_adapt <- colnames(mat1_HMM_100_70)
mat_names_HMM_adapt[T+T+(1:11)]

X2_init_HMM = round(colMeans(mat1_HMM_100_70[15000:30000,(1:T)]))
X4_init_HMM = round(colMeans(mat1_HMM_100_70[15000:30000,T+(1:T)]))


alpha <- colMeans(mat1_HMM_100_70[15000:30000,2*T+(1:4)])
beta <- colMeans(mat1_HMM_100_70[15000:30000,2*T+(7:10)])
alpharho <- mean(mat1_HMM_100_70[15000:30000,2*T+6])
alphal <- mean(mat1_HMM_100_70[15000:30000,2*T+5])
betal <- mean(mat1_HMM_100_70[15000:30000,2*T+11])
tauy <- mean(mat1_HMM_100_70[15000:30000,2*T+12])

fdays <- c(12, 26, 2, 9, 6, 8, 3, 2, 17, 3, 4, 8, 38, 12, 30, 3, 2, 15, 18, 46, 6,  
           3, 6, 15, 6, 13, 17, 22, 22, 4, 11, 9, 2, 1, 22, 57, 12, 13, 13, 2, 15, 
           13, 20, 12, 3, 1, 5, 0, 6, 9, 11, 33, 7, 2, 27, 5, 7, 25, 30, 16, 3, 1,  
           1, 12, 12, 9, 6, 3, 18, 10, 0, 1)

# Normalise fdays:
f <- (fdays-mean(fdays))/sd(fdays)


phi1 = double(T)
phi2 = double(T)
phi3 = double(T)
phi4 = double(T)
rho = double(T)

for(t in 1:(T)){
  ind =  alpha[1] + beta[1]*f[t]
  phi1[t] <-exp(ind)/(1+exp(ind)) # corresponds to the year 1963
  ind =  alpha[2] + beta[2]*f[t]
  phi2[t] <-exp(ind)/(1+exp(ind)) # corresponds to the year 1963
  ind =  alpha[3] + beta[3]*f[t]
  phi3[t] <-exp(ind)/(1+exp(ind)) # corresponds to the year 1963
  ind =  alpha[4] + beta[4]*f[t]
  phi4[t] <-exp(ind)/(1+exp(ind)) # corresponds to the year 1963
  rho[t] <- exp(alpharho) # We assume here that t=1
}

par(mfrow=c(2,2))
plot(phi1,type='l',main='phi1')
plot(phi2,type='l',main='phi2')
plot(phi3,type='l',main='phi3')
plot(phi4,type='l',main='phi4')




# X1_init=rep(1400,T)
X1_init=rep(2000,T)
X1_init[2:T] <- rpois((T-1),X4_init_HMM[1:(T-1)]*rho[1:(T-1)]*phi1[1:(T-1)])
 
X3_init=rep(1000,T)
X3_init[2:T] <- rbinom((T-1),round(X2_init_HMM[1:(T-1)]),phi3[1:(T-1)])
 
par(mfrow=c(4,1))
plot(X1_init,type='l',col='red',main='X1_init')
plot(X2_init_HMM,type='l',col='blue',main='X2_init')
plot(X3_init,type='l',col='red',main='X3_init')
plot(X4_init_HMM,type='l',col='blue',main='X4_init')


X1_init_HMM = X1_init
X3_init_HMM = X3_init



# alpha = c(-0.55, -0.1, 0.5, 1.1)
# beta = c(0.02, -0.15, -0.15, -0.2)
# alphal = 1.65
# betal = -0.2
# alpharho = 1.1
# tauy = 0.1000


inits <- function()(list(alpha=alpha, beta=beta, alphal=alphal, betal=betal, alpharho=alpharho, tauy=tauy,
                         X1=X1_init_HMM, X2=X2_init_HMM, X3=X3_init_HMM, X4=X4_init_HMM))


params <- c('alpha','alpharho','alphal', 'beta', 'betal', 'tauy', 'X2', 'X4')

save(file="X2_X4_inits_from_HMM.RData",inits, params, 
     alpha, beta, alphal, betal, alpharho, tauy,
     X2_init_HMM, X4_init_HMM, X1_init_HMM, X3_init_HMM)
