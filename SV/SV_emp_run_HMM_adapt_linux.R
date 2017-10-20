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

################# READ DATA ####
 
y <- read.csv('Data/Perc_Rets_GSPC_IBM_AAPL_MSFT_JPM_GE.csv')[,4]
time = c(1998,(2017 + 8/12));       
TT = length(y);
T = 2000;
y = y[(TT-T+1):TT];

sigma2_init = 0.15 


 
# plot(y,type='l')
# lines(h_true, type='l', col='red')


################# HMM WITH ADAPTIVE INTERVALS ####


N_q = 30
qu <- c(0:(N_q-1))/N_q
mid <- qu+qu[2]/2
Up_h <- 10

inits_hmm_adapt <- function()(list(mu = 0, phi_star = (0.97+1)/2, sigma2_star = 1/sigma2_init, h = log(var(y))*rep(1,T/2)))
# params_hmm <- c('mu','sigma2','phi','h','bin_mid','Q_odd')
# params_hmm_adapt <- c('mu','sigma2','phi','h','BIN100')
params_hmm_adapt <- c('mu','sigma2','phi','h')

data_hmm_adapt <- list(y=y, T=T, Up_h = Up_h,
                       mid = mid, qu= qu, N_q = N_q,
                       zeros = rep(0,round(T/2)))

# data_hmm_adapt <- list(y=y, T=T, Up_h = Up_h,
#                        mid = mid, qu= qu, N_q = N_q,
#                        zeros = rep(0,round(T/2)),
#                        bin_max = Up_h)

tstart = proc.time()
sv_model_HMM_adapt <- jags.model('sv_model_hmm_adapt.R', 
                                 data=data_hmm_adapt, inits=inits_hmm_adapt, n.chains = cha, n.adapt = ada)
time_init_HMM_adapt = proc.time()-tstart

tstart = proc.time()
output_sv_HMM_adapt <- coda.samples(sv_model_HMM_adapt, params_hmm_adapt, n.iter = iter, thin=th)
time_sample_HMM_adapt = proc.time()-tstart



# save selected output

mat1_HMM_adapt = as.matrix(output_sv_HMM_adapt[1])
mat2_HMM_adapt = as.matrix(output_sv_HMM_adapt[2])
mat_names_HMM_adapt <- colnames(mat1_HMM_adapt)

# BIN100_1 =  mat1_HMM_adapt[,1:300]
# BIN100_2 =  mat2_HMM_adapt[,1:300]

ESS_HMM_adapt = lapply(output_sv_HMM_adapt,effectiveSize)
ESS1_HMM_adapt = as.matrix(ESS_HMM_adapt[[1]])
ESS2_HMM_adapt = as.matrix(ESS_HMM_adapt[[2]])

theta1_HMM_adapt = mat1_HMM_adapt[,(T/2)+(1:3)] 
theta2_HMM_adapt = mat2_HMM_adapt[,(T/2)+(1:3)] 

H_short1_HMM_adapt = mat1_HMM_adapt[,seq(100,T/2,by=100)]
H_short2_HMM_adapt = mat2_HMM_adapt[,seq(100,T/2,by=100)]

mean_H1_HMM_adapt = colMeans(mat1_HMM_adapt[,1:(T/2)])
mean_H2_HMM_adapt = colMeans(mat2_HMM_adapt[,1:(T/2)])



# if (save_on) {
#   save(file=paste("SV_emp_HMM_adapt_Nq",toString(N_q),"_T",toString(T),"_selected_mid_BIN100.RData",sep=""),
#        mat_names_HMM_adapt, BIN100_1, BIN100_2)
# }

if (save_on) {
  save(file=paste("SV_emp_HMM_adapt_Nq",toString(N_q),"_T",toString(T),"_selected_mid.RData",sep=""),
       y, sv_model_HMM_adapt, time_init_HMM_adapt, time_sample_HMM_adapt, mat_names_HMM_adapt,
       ESS1_HMM_adapt, ESS2_HMM_adapt, theta1_HMM_adapt, theta2_HMM_adapt, 
       H_short1_HMM_adapt, H_short2_HMM_adapt, mean_H1_HMM_adapt, mean_H2_HMM_adapt)
}

quit()


 