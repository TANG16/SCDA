# setwd("SV")
# rm(list=ls())
library(rjags)
library(coda)


phi <- 0.98
sigma <- 0.2
sigma2 <- sigma^2
# beta = 0.05;
beta <- 0.5;
mu <- 2*log(beta);

a1 = mu
P1 <- (sigma^2)/(1-phi^2)


################# simultated data ####
T <- 1000
h_true <- rep(NaN,T)

h_true[1] <- a1 + sqrt(P1)*rnorm(1)
for (t in c(2:T)){
  h_true[t] = mu + phi*(h_true[t-1]-mu) + sigma*rnorm(1)
}
y <- exp(h_true/2)*rnorm(T)
plot(y,type='l')
lines(h_true, type='l', col='red')

################# full DA ####

# inits <- function()(list(mu = 0, phi = 0.97, sigma = 0.15, h = var(y)*rep(1,T)))
inits <- function()(list(mu = 0, phi_star = (0.97+1)/2, sigma2_star = 1/(0.15^2), h = log(var(y))*rep(1,T)))
params <- c('mu','sigma2','phi','h')
data <- list(y=y,T=T)

sv_model_string <- "
model{
  # define the observation process
  for (t in 1:T){
    y[t] ~ dnorm(0, 1.0/exp(h[t]))
  }
  
  # define the state process
  for (t in 2:T){
    h[t] ~ dnorm(mean_h[t],sigma2_star)
    mean_h[t] <- mu + phi*(h[t-1] - mu) 
  }

  h[1] ~ dnorm(mean_h[1],sigma2_star)
  mean_h[1] <- mu + phi*(h0 - mu) 

  h0 ~ dnorm(mean_h0,1/P1)
  P1 <- sigma2/(1-phi^2)
  mean_h0 <- mu
    
  # define the priors for parameters
  mu ~ dnorm(0,0.1)

  # (phi+1)/2 ~ betapdf((phi+1)/2, 20, 1.5)
  phi <- 2*phi_star - 1
  phi_star ~ dbeta(20, 1.5)
  # (phi+1)/2 ~ dbeta(20, 1.5)

  # 1/sigma2 ~ gampdf(1./s2, 5/2, 0.05/2);
  sigma2 <- 1/sigma2_star
#   sigma2_star ~ dgamma(5/2, 2/0.05) # shape and rate
  sigma2_star ~ dgamma(5/2, 0.05/2) # shape and rate
  # 1/sigma2  ~  dgamma(5/2, 2/0.05)
}
"

### Prior for phi
xxx <- seq(0,1,by=0.01)
bet <- dbeta((xxx+1)/2,20,1.5)
plot(xxx,bet,type="l",xlab="",ylab="",sub='phi prior')



### JAGS
tstart = proc.time()
sv_model <- jags.model(textConnection(sv_model_string), 
                       data=data, inits=inits, n.chains = 2, n.adapt = 1000)
time_init_DA = proc.time()-tstart

tstart = proc.time()
output_sv <- coda.samples(sv_model, params, n.iter = 10000, thin=1)
time_sample_DA = proc.time()-tstart

################# OUTPUT ANALYSIS ####
# print(summary(output_sv))
mean(mat1_DA[,T+1])
mean(mat1_DA[,T+2])
mean(mat1_DA[,T+3])
mat1_DA = as.matrix(output_sv_DA[1]) 
mat2_DA = as.matrix(output_sv_DA[2]) 
mat_names_DA <- colnames(mat1_DA)

par(mfrow = c(3,1))
for (i in 1:3){
  plot(mat1_DA[,T+i],type='l',col='blue', xlab="", ylab="", sub=mat_names_DA[T+i])
  lines(mat2_DA[,T+i],type='l',col='red')
}


mean(mat2_DA[,T+1])
mean(mat2_DA[,T+2])
mean(mat2_DA[,T+3])

mean(mat1_DA[5000:10000,T+1])
mean(mat1_DA[5000:10000,T+2])
mean(mat1_DA[5000:10000,T+3])

mean(mat2_DA[5000:10000,T+1])
mean(mat2_DA[5000:10000,T+2])
mean(mat2_DA[5000:10000,T+3])


colMeans(mat1_DA[,T+3])

par(mfrow = c(1,1))
plot(colMeans(mat1_DA[5000:10000,1:T]),type='l',col='blue', xlab="", ylab="", sub="mean h")
lines(colMeans(mat2_DA[5000:10000,1:T]),type='l',col='red')
lines(h_true,type="l",col="black")


ESS_DA = lapply(output_sv_DA,effectiveSize)
ESS_DA1 = as.matrix(ESS_DA[[1]])
ESS_DA2 = as.matrix(ESS_DA[[2]])

################# HMM integration ####
# integrate out the odd states, impute the even ones
# assume T is even so for the last observation we have an imputation

N_bin <- 30
bin_range <- 4
Up_h <- 10
inits_hmm <- function()(list(mu = 0, phi_star = (0.97+1)/2, sigma2_star = 1/0.15, h = log(var(y))*rep(1,T/2)))
# params_hmm <- c('mu','sigma2','phi','h','G_odd','G_even','Q_odd','Q_even')
# params_hmm <- c('mu','sigma2','phi','G_odd','G_even')
# params_hmm <- c('mu','sigma2','phi','Q_odd','Q_even')
params_hmm <- c('mu','sigma2','phi','h')

#   N_bin = 30
#   bin_range = 4
#   bins = seq(-bin_range,bin_range,length=N_bin+1);
#   bin_midpoint = (bins[1:N_bin] + bins[2:(N_bin+1)])/2;

bin = rep(NaN, N_bin+1) # bins boundaries
for (i in 0:N_bin){
  bin[i+1] <- -bin_range + i*2*bin_range/N_bin # the (i+1)'th bin's midpoint
}

bin_mid = rep(NaN, N_bin)
# Bins' midpoints
for (i in 0:(N_bin-1)){
  bin_mid[i+1] <- -bin_range + (2*i+1)*bin_range/N_bin # the (i+1)'th bin's midpoint
}

data_hmm <- list(y=y, T=T, N_bin=N_bin, bin_range=bin_range, Up_h = Up_h, bin_mid=bin_mid, bin = bin, zeros = rep(0,round(T/2)))

sv_model_hmm <- jags.model('sv_model_hmm.R', 
                           data=data_hmm, inits=inits_hmm, n.chains = 2, n.adapt = 1000)

output_sv_HMM <- coda.samples(sv_model_hmm, params_hmm, n.iter = 100000, thin=1)
# print(summary(output_sv_HMM))

# save(file="SV_HMM_G_mats.RData",output_sv_HMM)
save(file="SV_HMM_Q_mats.RData",output_sv_HMM)

mat1_HMM = as.matrix(output_sv_HMM[1]) 
mat2_HMM = as.matrix(output_sv_HMM[2]) 
mat_names_HMM <- colnames(mat1_HMM)



ESS_HMM = lapply(output_sv_HMM,effectiveSize)
ESS_HMM1 = as.matrix(ESS_HMM[[1]])
ESS_HMM2 = as.matrix(ESS_HMM[[2]])


ESS_DA1[T+1]
ESS_DA2[T+1]
ESS_HMM1[(T/2)+1]
ESS_HMM2[(T/2)+1]


ESS_DA1[T+2]
ESS_DA2[T+2]
ESS_HMM1[(T/2)+2]
ESS_HMM2[(T/2)+2]


ESS_DA1[T+3]
ESS_DA2[T+3]
ESS_HMM1[(T/2)+3]
ESS_HMM2[(T/2)+3]

par(mfrow = c(1,1))
plot(ESS_DA1[1:T],type='l',col='blue', xlab="", ylab="")
lines(ESS_DA2[1:T],type='l',col='green', xlab="", ylab="")
lines(seq(2,T,by=2),ESS_HMM1[1:(T/2)],type='l',col='magenta', xlab="", ylab="")
lines(seq(2,T,by=2),ESS_HMM2[1:(T/2)],type='l',col='red', xlab="", ylab="")

time_init_DA[3]
time_init_HMM[3]

time_sample_DA[3]
time_sample_HMM[3]

######## Params trace plots
# YL = matrix(c(-2,0, 0.96,1, 0,0.1), ncol = 2, byrow=TRUE)
par(mfrow = c(3,1))
for (i in 1:3){
#   plot(param[i] + 0*mat2_HMM[,T/2+i],type='l',col='red',ylim= YL[i,], xlab="", ylab="", sub=mat_names_HMM[T/2+i])
  plot(mat1_HMM[,T/2+i],type='l',col='blue', xlab="", ylab="", sub=mat_names_HMM[T/2+i])
#   lines(mat1_HMM[,T/2+i],type='l',col='blue')
  lines(mat2_HMM[,T/2+i],type='l',col='green')
#   lines(param[i] + 0*mat2_HMM[,T/2+i],type='l',col='red')
}

mean(mat1_HMM[,(T/2)+1])
mean(mat1_HMM[,(T/2)+2])
mean(mat1_HMM[,(T/2)+3])

mean(mat2_HMM[,(T/2)+1])
mean(mat2_HMM[,(T/2)+2])
mean(mat2_HMM[,(T/2)+3])

par(mfrow = c(1,1))
plot(h_true,type='l')
# lines(colMeans(mat1_DA[1:5000,1:T]),col='red')
lines(seq(2,T,by=2),colMeans(mat1_HMM[,1:(T/2)]),col='blue')
lines(seq(2,T,by=2),colMeans(mat2_HMM[,1:(T/2)]),col='red')


par(mfrow = c(3,1))
for (i in 1:3){
  plot(mat1_HMM[7000:10000,T/2+i],type='l',col='blue', xlab="", ylab="", sub=mat_names_HMM[T/2+i])
  lines(mat2_HMM[7000:10000,T/2+i],type='l',col='red')
}

# ################# HMM WITH FIXED PARAMS ####
# 
# inits_hmm_FIX <- function()(list(h = log(var(y))*rep(1,T/2)))
# params_hmm_FIX <- c( 'h')
# data_hmm_FIX <- list(y=y, T=T, N_bin=N_bin, bin_range=bin_range, Up_h=Up_h, bin=bin,
#                      bin_mid=bin_mid, mu=mu, phi=phi, sigma2 = sigma2)
# 
# sv_model_hmm_FIX <- jags.model('sv_model_hmm_FIX.R', 
#                            data=data_hmm_FIX, inits=inits_hmm_FIX, n.chains = 2, n.adapt = 1000)
# output_sv_HMM_FIX <- coda.samples(sv_model_hmm_FIX, params_hmm, n.iter = 10000, thin=1)
# # print(summary(output_sv_HMM))
# 
# mat1_HMM_FIX  = as.matrix(output_sv_HMM_FIX [1]) 
# mat2_HMM_FIX  = as.matrix(output_sv_HMM_FIX [2]) 
# mat_names_HMM_FIX  <- colnames(mat1_HMM_FIX )
# 
# plot(h_true,type='l')
# lines(seq(2,T,by=2),colMeans(mat1_HMM_FIX[,1:(T/2)]),col='green')
# lines(seq(2,T,by=2),colMeans(mat1_HMM_FIX[5000:10000,1:(T/2)]),col='red')

################# HMM WITH ADAPTIVE INTERVALS ####

N_q = 10
qu <- c(0:(N_q-1))/N_q
mid <- qu+qu[2]/2
Up_h <- 10

inits_hmm <- function()(list(mu = 0, phi_star = (0.97+1)/2, sigma2_star = 1/0.15, h = log(var(y))*rep(1,T/2)))
params_hmm <- c('mu','sigma2','phi','h')

data_hmm_adapt <- list(y=y, T=T, Up_h = Up_h,
                 mid = mid, qu= qu, N_q = N_q)

sv_model_hmm_adapt <- jags.model('sv_model_hmm_adapt.R', 
                           data=data_hmm_adapt, inits=inits_hmm, n.chains = 2, n.adapt = 1000)

output_sv_HMM_adapt <- coda.samples(sv_model_hmm_adapt, params_hmm, n.iter = 10000, thin=1)


mat1_HMM_adapt = as.matrix(output_sv_HMM_adapt[1]) 
mat2_HMM_adapt = as.matrix(output_sv_HMM_adapt[2]) 
mat_names_HMM_adapt <- colnames(mat1_HMM_adapt)



par(mfrow = c(3,1))
for (i in 1:3){
  plot(mat1_HMM_adapt[,T/2+i],type='l',col='blue', xlab="", ylab="", sub=mat_names_HMM_adapt[T/2+i])
  lines(mat2_HMM_adapt[,T/2+i],type='l',col='red')
}


mean(mat1_HMM_adapt[,(T/2)+1])
mean(mat1_HMM_adapt[,(T/2)+2])
mean(mat1_HMM_adapt[,(T/2)+3])

mean(mat2_HMM_adapt[,(T/2)+1])
mean(mat2_HMM_adapt[,(T/2)+2])
mean(mat2_HMM_adapt[,(T/2)+3])


par(mfrow = c(1,1))
plot(h_true,type='l')
# lines(colMeans(mat1_DA[1:5000,1:T]),col='red')
lines(seq(2,T,by=2),colMeans(mat1_HMM_adapt[5000:10000,1:(T/2)]),col='blue')
lines(seq(2,T,by=2),colMeans(mat2_HMM_adapt[,1:(T/2)]),col='red')


par(mfrow = c(3,3))
for (i in seq(100,900,by=100)){
  plot(mat1_HMM_adapt[,i/2],type='l',col='blue', xlab="", ylab="", sub=mat_names_HMM_adapt[i/2])
  lines(mat2_HMM_adapt[,i/2],type='l',col='red')
  lines(h_true[i]+0*(mat2_HMM_adapt[,i/2]),type='l',col='black')  
}
