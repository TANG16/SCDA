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
  y <- read.csv("OECD_G7CPI.csv",header=FALSE)[,D]
  T <- length(y)
  y <- log(y[2:T]/y[1:(T-1)])*400;       
}
T <- length(y)

plot(y,type='l',xlab='',ylab='')

################# full DA ####

g_sgn = -1

# inits <- function()(list(mu = 0, phi = 0.97, sigma = 0.15, h = var(y)*rep(1,T)))
inits <- function()(list(h0 = 0, g0 = 0, omega_h = sqrt(0.2), omega_g = sqrt(0.2),
                         h = log(var(y))*rep(1,T),  g = g_sgn*log(var(y))*rep(1,T)))
params <- c('h0','g0','omega_h','omega_g', 'tau','h', 'g')
data <- list(y=y,T=T,tau0=0,V_tau = 10, V_h=10, V_g=10)

ucsv_model_string <- "
model{
# define the observation process for inflation
for (t in 1:T){
  y[t] ~ dnorm(tau[t], 1.0/exp(h0 + omega_h*h[t]))
}

# define the state process for the trend
for (t in 2:T){
  tau[t] ~ dnorm(tau[t-1], 1.0/exp(g0 + omega_g*g[t]))
 }

# define the state process for the volatilities
for (t in 2:T){
  h[t] ~ dnorm(h[t-1],1)
  g[t] ~ dnorm(g[t-1],1)
}

# initialise
tau[1] ~ dnorm(tau0, 1/(V_tau*exp(g0 + omega_g*g[1])))
h[1] ~ dnorm(0,1/V_h)
g[1] ~ dnorm(0,1/V_g)
 

# define the priors for parameters
omega_h ~ dnorm(0,1/(0.2))
omega_g ~ dnorm(0,1/(0.2))
h0 ~ dnorm(0,1/10)
g0 ~ dnorm(0,1/10) 
}
"


### JAGS #####
tstart = proc.time()
ucsv_model_DA <- jags.model(textConnection(ucsv_model_string), 
                       data=data, inits=inits, n.chains = cha, n.adapt = ada)
time_init_DA = proc.time()-tstart

tstart = proc.time()
output_ucsv_DA <- coda.samples(ucsv_model_DA, params, n.iter = iter, thin=th)
time_sample_DA = proc.time()-tstart

# save selected output

mat1_DA = as.matrix(output_ucsv_DA[1]) 
mat2_DA = as.matrix(output_ucsv_DA[2]) 
mat_names_DA <- colnames(mat1_DA)

ESS_DA = lapply(output_ucsv_DA,effectiveSize)
ESS1_DA = as.matrix(ESS_DA[[1]])
ESS2_DA = as.matrix(ESS_DA[[2]])

mat_names_DA[1:T] # g[1] : g[258]
mat_names_DA[T+1] # g0
mat_names_DA[(T+2):(2*T+1)] # h[1] : h[258]
mat_names_DA[(2*T+2)] # h0
mat_names_DA[(2*T+2+1)]  # omega_g
mat_names_DA[(2*T+2+2)]  # omega_h
mat_names_DA[(2*T+2+1):(3*T+4)] # tau[1] : tau[258]

theta1_DA = mat1_DA[,c(T+1,2*T+2,2*T+2+1,2*T+2+2)] 
theta2_DA = mat2_DA[,c(T+1,2*T+2,2*T+2+1,2*T+2+2)]
colnames(theta1_DA) = c('g0','h0','omega_g','omega_h')
colnames(theta2_DA) = c('g0','h0','omega_g','omega_h')

tau1_sel_DA = mat1_DA[,seq((2*T+2+1),(3*T+4),by=50)]
tau2_sel_DA = mat2_DA[,seq((2*T+2+1),(3*T+4),by=50)]

mean_tau1_DA = colMeans(mat1_DA[,(2*T+2+1):(3*T+4)])
mean_tau2_DA = colMeans(mat2_DA[,(2*T+2+1):(3*T+4)])


g1_sel_DA = mat1_DA[,seq(1,T,by=50)]
g2_sel_DA = mat2_DA[,seq(1,T,by=50)]

mean_g1_DA = colMeans(mat1_DA[,1:T])
mean_g2_DA = colMeans(mat2_DA[,1:T])


h1_sel_DA = mat1_DA[,seq((T+2),(2*T+1),by=50)]
h2_sel_DA = mat2_DA[,seq((T+2),(2*T+1),by=50)]

mean_h1_DA = colMeans(mat1_DA[,(T+2):(2*T+1)])
mean_h2_DA = colMeans(mat2_DA[,(T+2):(2*T+1)])



if (save_on) {
  save(file=paste("UCSV_emp_DA_data",toString(D),".RData",sep=""),
       y, ucsv_model_DA, time_init_DA, time_sample_DA, mat_names_DA,
       ESS1_DA, ESS2_DA, theta1_DA, theta2_DA, 
       tau1_sel_DA, tau2_sel_DA, mean_tau1_DA, mean_tau2_DA,
       g1_sel_DA, g2_sel_DA, mean_g1_DA, mean_g2_DA,
       h1_sel_DA, h2_sel_DA, mean_h1_DA, mean_h2_DA)
}
  
# quit()