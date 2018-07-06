# setwd("SV/UCSV")
rm(list = ls())
library(rjags)
library(coda)
set.seed(1345221)

save_on = TRUE


ada=10000
iter=100000
th=1
cha= 2 #2

################# READ DATA ####
# data = xlsread('OECD_G7CPI.xls', 'D51:J266'); 
# country_id = 7; %1: Canada; 2: France; 3: Germany; 4: Italy; 5: Japan; 6: UK; 7: US

D = 7
if (D == 0)
{
  y <- read.csv("USCPI.csv")[,1]
}else
{
  y <- read.csv("OECD_G7CPI.csv")[,D]
  T <- length(y)
  y <- log(y[2:T]/y[1:(T-1)])*400;       
}
T <- length(y)

par(mfrow=c(1,1))
plot(y,type='l',xlab='',ylab='')

################# full DA ####

# inits <- function()(list(mu = 0, phi = 0.97, sigma = 0.15, h = var(y)*rep(1,T)))
inits <- function()(list(h0 = 0, g0 = 0, omega_h = sqrt(0.2), omega_g = sqrt(0.2),
                         h = log(var(y))*rep(1,T),  g = log(var(y))*rep(1,T)))
params <- c('h0','g0','omega_h','omega_g', 'tau','h', 'g')
data <- list(y=y,T=T,tau0=0,V_tau = 10, V_h=10, V_g=10,
             zeros = rep(0,T/2))



### JAGS #####
tstart = proc.time()
ucsv_model_HMM <- jags.model(ucsv_model_hmm, 
                       data=data, inits=inits, n.chains = cha, n.adapt = ada)
time_init_HMM = proc.time()-tstart

tstart = proc.time()
output_ucsv_HMM <- coda.samples(ucsv_model_HMM, params, n.iter = iter, thin=th)
time_sample_HMM = proc.time()-tstart

# save selected output

mat1_HMM = as.matrix(output_ucsv_HMM[1]) 
mat2_HMM = as.matrix(output_ucsv_HMM[2]) 
mat_names_HMM <- colnames(mat1_HMM)

ESS_HMM = lapply(output_ucsv_HMM,effectiveSize)
ESS1_HMM = as.matrix(ESS_HMM[[1]])
ESS2_HMM = as.matrix(ESS_HMM[[2]])

mat_names_HMM[1:T] # g[1] : g[258]
mat_names_HMM[T+1] # g0
mat_names_HMM[(T+2):(2*T+1)] # h[1] : h[258]
mat_names_HMM[(2*T+2)] # h0
mat_names_HMM[(2*T+2+1)]  # omega_g
mat_names_HMM[(2*T+2+2)]  # omega_h
mat_names_HMM[(2*T+2+1):(3*T+4)] # tau[1] : tau[258]

theta1_HMM = mat1_HMM[,c(T+1,2*T+2,2*T+2+1,2*T+2+2)] 
theta2_HMM = mat2_HMM[,c(T+1,2*T+2,2*T+2+1,2*T+2+2)]
colnames(theta1_HMM) = c('g0','h0','omega_g','omega_h')
colnames(theta2_HMM) = c('g0','h0','omega_g','omega_h')

tau1_sel_HMM = mat1_HMM[,seq((2*T+2+1),(3*T+4),by=50)]
tau2_sel_HMM = mat2_HMM[,seq((2*T+2+1),(3*T+4),by=50)]

mean_tau1_HMM = colMeans(mat1_HMM[,(2*T+2+1):(3*T+4)])
mean_tau2_HMM = colMeans(mat2_HMM[,(2*T+2+1):(3*T+4)])


g1_sel_HMM = mat1_HMM[,seq(1,T,by=50)]
g2_sel_HMM = mat2_HMM[,seq(1,T,by=50)]

mean_g1_HMM = colMeans(mat1_HMM[,1:T])
mean_g2_HMM = colMeans(mat2_HMM[,1:T])


h1_sel_HMM = mat1_HMM[,seq((T+2),(2*T+1),by=50)]
h2_sel_HMM = mat2_HMM[,seq((T+2),(2*T+1),by=50)]

mean_h1_HMM = colMeans(mat1_HMM[,(T+2):(2*T+1)])
mean_h2_HMM = colMeans(mat2_HMM[,(T+2):(2*T+1)])



if (save_on) {
  save(file=paste("UCSV_emp_HMM_HMMta",toString(D),".RData",sep=""),
       y, ucsv_model_HMM, time_init_HMM, time_sample_HMM, mat_names_HMM,
       ESS1_HMM, ESS2_HMM, theta1_HMM, theta2_HMM, 
       tau1_sel_HMM, tau2_sel_HMM, mean_tau1_HMM, mean_tau2_HMM,
       g1_sel_HMM, g2_sel_HMM, mean_g1_HMM, mean_g2_HMM,
       h1_sel_HMM, h2_sel_HMM, mean_h1_HMM, mean_h2_HMM)
}
  
# quit()