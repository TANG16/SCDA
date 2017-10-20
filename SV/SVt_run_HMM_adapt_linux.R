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
# beta = 0.05;
beta <- 0.5;
mu <- 2*log(beta);
nu <- 5
rho = (nu-2)/nu
param <- c(mu, phi, sigma2,nu)

sigma2_init = 0.15 # working case
# sigma2_init = 0.1^2

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


################# HMM WITH ADAPTIVE INTERVALS ####


N_q = 30 #20
qu <- c(0:(N_q-1))/N_q
mid <- qu+qu[2]/2
Up_h <- 10

inits_hmm_adapt <- function()(list(mu = 0, phi_star = (0.97+1)/2, sigma2_star = 1/sigma2_init, nu = 8, h = log(var(y))*rep(1,T/2)))
# params_hmm <- c('mu','sigma2','phi','h','bin_mid','Q_odd')
params_hmm_adapt <- c('mu','sigma2','phi', 'nu','h')

data_hmm_adapt <- list(y=y, T=T, Up_h = Up_h,
                       mid = mid, qu= qu, N_q = N_q,
                       zeros = rep(0,round(T/2)))

# data_hmm_adapt <- list(y=y, T=T, Up_h = Up_h,
#                        mid = mid, qu= qu, N_q = N_q,
#                        zeros = rep(0,round(T/2)),
#                        bin_max = Up_h)

tstart = proc.time()
svt_model_HMM_adapt <- jags.model('svt_model_hmm_adapt.R', 
                                 data=data_hmm_adapt, inits=inits_hmm_adapt, n.chains = cha, n.adapt = ada)
time_init_HMM_adapt = proc.time()-tstart

tstart = proc.time()
output_svt_HMM_adapt <- coda.samples(svt_model_HMM_adapt, params_hmm_adapt, n.iter = iter, thin=th)
time_sample_HMM_adapt = proc.time()-tstart



# save selected output

mat1_HMM_adapt = as.matrix(output_svt_HMM_adapt[1]) 
mat2_HMM_adapt = as.matrix(output_svt_HMM_adapt[2]) 
mat_names_HMM_adapt <- colnames(mat1_HMM_adapt)

ESS_HMM_adapt = lapply(output_svt_HMM_adapt,effectiveSize)
ESS1_HMM_adapt = as.matrix(ESS_HMM_adapt[[1]])
ESS2_HMM_adapt = as.matrix(ESS_HMM_adapt[[2]])

theta1_HMM_adapt = mat1_HMM_adapt[,(T/2)+(1:4)] 
theta2_HMM_adapt = mat2_HMM_adapt[,(T/2)+(1:4)] 

H_short1_HMM_adapt = mat1_HMM_adapt[,seq(100,T/2,by=100)]
H_short2_HMM_adapt = mat2_HMM_adapt[,seq(100,T/2,by=100)]

mean_H1_HMM_adapt = colMeans(mat1_HMM_adapt[,1:(T/2)])
mean_H2_HMM_adapt = colMeans(mat2_HMM_adapt[,1:(T/2)])


if (save_on) {
  save(file=paste("SVt_HMM_adapt_Nq",toString(N_q),"_T",toString(T),"_selected.RData",sep=""),
       y, h_true, param, 
       svt_model_HMM_adapt, time_init_HMM_adapt, time_sample_HMM_adapt, mat_names_HMM_adapt,
       ESS1_HMM_adapt, ESS2_HMM_adapt, theta1_HMM_adapt, theta2_HMM_adapt, 
       H_short1_HMM_adapt, H_short2_HMM_adapt, mean_H1_HMM_adapt, mean_H2_HMM_adapt)
}
# if (save_on) {
#   save(file=paste("SVt_HMM_adapt_Nq",toString(N_q),"_T",toString(T),".RData",sep=""),
#        y, h_true, param, output_svt_HMM_adapt, svt_model_HMM_adapt, time_init_HMM_adapt, time_sample_HMM_adapt)
# }

quit()