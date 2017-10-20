# setwd("SV")
rm(list = ls())
library(rjags)
library(coda)
set.seed(1345221)

save_on = TRUE


ada=1000
iter=10000
th=1
cha= 2 #2


phi <- 0.98
sigma <- 0.2
sigma2 <- sigma^2
sigma2_init <- 0.15^2
nu <- 5;
# beta = 0.05;
beta <- 0.5;
mu <- 2*log(beta);
rho <- (nu-2)/nu

param <- c(mu, phi, sigma2, nu)

a1 = mu
P1 <- (sigma^2)/(1-phi^2)


################# simultated data ####
T <- 2000
h_true <- rep(NaN,T)

h_true[1] <- a1 + sqrt(P1)*rnorm(1)
for (t in c(2:T)){
  h_true[t] = mu + phi*(h_true[t-1]-mu) + sigma*rnorm(1)
}
y <- sqrt(rho)*exp(h_true/2)*rt(T,nu)
# plot(y,type='l')
# lines(h_true, type='l', col='red')



################# HMM integration ####
# integrate out the odd states, impute the even ones
# assume T is even so for the last observation we have an imputation

N_bin <- 30; #20 # 10 #30
bin_range <- 4
Up_h <- 10
inits_hmm <- function()(list(mu = 0, phi_star = (0.97+1)/2, nu = 8,  sigma2_star = 1/sigma2_init, h = log(var(y))*rep(1,T/2)))
# params_hmm <- c('mu','sigma2','phi','h','G_odd','G_even','Q_odd','Q_even')
# params_hmm <- c('mu','sigma2','phi','G_odd','G_even')
# params_hmm <- c('mu','sigma2','phi','Q_odd','Q_even')
params_hmm <- c('mu','sigma2','phi','nu','h')

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


tstart = proc.time()
svt_model_HMM <- jags.model('svt_model_hmm.R', 
                           data=data_hmm, inits=inits_hmm, n.chains = cha, n.adapt = ada)
time_init_HMM = proc.time()-tstart

tstart = proc.time()
output_svt_HMM <- coda.samples(svt_model_HMM, params_hmm, n.iter = iter, thin=th)
time_sample_HMM = proc.time()-tstart


# save selected output

mat1_HMM = as.matrix(output_svt_HMM[1]) 
mat2_HMM = as.matrix(output_svt_HMM[2]) 
mat_names_HMM <- colnames(mat1_HMM)

ESS_HMM = lapply(output_svt_HMM,effectiveSize)
ESS1_HMM = as.matrix(ESS_HMM[[1]])
ESS2_HMM = as.matrix(ESS_HMM[[2]])

theta1_HMM = mat1_HMM[,(T/2)+(1:4)] 
theta2_HMM = mat2_HMM[,(T/2)+(1:4)] 

H_short1_HMM = mat1_HMM[,seq(100,T/2,by=100)]
H_short2_HMM = mat2_HMM[,seq(100,T/2,by=100)]

mean_H1_HMM = colMeans(mat1_HMM[,1:(T/2)])
mean_H2_HMM = colMeans(mat2_HMM[,1:(T/2)])


# if (save_on) {
#   save(file=paste("SVt_HMM_Nbin",toString(N_bin),"_T",toString(T),".RData",sep=""),
#        y, h_true, param, output_svt_HMM,svt_model_HMM,time_init_HMM,time_sample_HMM)
# }

if (save_on) {
  save(file=paste("SVt_HMM_Nbin",toString(N_bin),"_T",toString(T),"_selected.RData",sep=""),
       y, h_true, param, 
       svt_model_HMM, time_init_HMM, time_sample_HMM, mat_names_HMM,
       ESS1_HMM, ESS2_HMM, theta1_HMM, theta2_HMM, 
       H_short1_HMM, H_short2_HMM, mean_H1_HMM, mean_H2_HMM)
}


quit()