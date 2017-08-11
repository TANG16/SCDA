# setwd("SV")
rm(list = ls())
library(rjags)
library(coda)
set.seed(1345221)

save_on = TRUE


ada=1000
iter=1000
th=1
cha= 2 #2


phi <- 0.98
sigma <- 0.2
sigma2 <- sigma^2
# beta = 0.05;
beta <- 0.5;
mu <- 2*log(beta);

param <- c(mu, phi, sigma2)

# sigma2_init = 0.15 # working case
sigma2_init = 0.1^2

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
# plot(y,type='l')
# lines(h_true, type='l', col='red')


################# HMM WITH ADAPTIVE INTERVALS ####


N_q = 10
qu <- c(0:(N_q-1))/N_q
mid <- qu+qu[2]/2
Up_h <- 10

inits_hmm <- function()(list(mu = 0, phi_star = (0.97+1)/2, sigma2_star = 1/sigma2_init, h = log(var(y))*rep(1,T/2)))
params_hmm <- c('mu','sigma2','phi','h','bin_mid','Q_odd')

data_hmm_adapt <- list(y=y, T=T, Up_h = Up_h,
                       mid = mid, qu= qu, N_q = N_q,
                       zeros = rep(0,round(T/2)))

tstart = proc.time()
sv_model_HMM_adapt <- jags.model('sv_model_hmm_adapt.R', 
                                 data=data_hmm_adapt, inits=inits_hmm, n.chains = cha, n.adapt = ada)
time_init_HMM_adapt = proc.time()-tstart

tstart = proc.time()
output_sv_HMM_adapt <- coda.samples(sv_model_HMM_adapt, params_hmm, n.iter = iter, thin=th)
time_sample_HMM_adapt = proc.time()-tstart



if (save_on) {
  save(file=paste("SV_HMM_adapt_Nq",toString(N_q),".RData",sep=""),
       y, h_true, param, output_sv_HMM_adapt, sv_model_HMM_adapt, time_init_HMM_adapt, time_sample_HMM_adapt)
}

quit()