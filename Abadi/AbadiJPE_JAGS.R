# Set working directory: ####
# setwd("~/Documents/ATI/Interim/BUGS/Abadi")
# setwd("C:/Users/aga/Dropbox/Research Visit/Ruth King/Codes")

#install.packages("R2jags", dep=TRUE)
#library(R2jags)

# Load required packages and fix the random seed
library(rjags)
library(coda)
library(lattice)
set.seed(134522)

# MCMC details: ####

ada=100
iter=10000
th=1
cha=2

# Load data: ###
source("AbadiJPE_data.R")

# Set initial parameter values: ####

# Simulate values for population sizes from the model:

#N1init <- rep(2,ti)
#Nadimminit <- rep(2,ti)
#NadSurvinit <- rep(10,ti)
#Ntotinit <- N1init+Nadimminit+NadSurvinit

#for (i in 2:ti){
#  N1init[i] <- rpois(1,0.1*Ntotinit[i-1])
#  Nadimminit[i] <- rpois(1,Ntotinit[i-1]*0.1)
#  NadSurvinit[i] <- rbinom(1,Ntotinit[i-1],0.75)
#  Ntotinit[i] <- N1init[i] + Nadimminit[i] + NadSurvinit[i]
#}

#N1init[1] <- NA
#Nadimminit[1] <- NA
#NadSurvinit[1] <- NA

inits <- function(){list(v=rnorm(7),bp=rnorm((ti-1)),fec=runif((ti-1),0,5),
                           N1=c(NA,round(runif(ti-1,1,50),0)),
                           NadSurv=c(NA,round(runif(ti-1,5,10),0)),
                           Nadimm=c(NA,round(runif(ti-1,5,10),0)),
                           N1prior=5,NadSurvprior=20,Nadimmprior=5)}

#inits <- function(){list(v=rep(0,7),bp=rep(0,ti-1),fec=rep(1,(ti-1)),
#                           N1=N1init,NadSurv=NadSurvinit,
#                           Nadimm=Nadimminit)}

#inits <- function(){list(v=rep(0,7),bp=rep(0,ti-1),fec=rep(1,(ti-1)),
#                        N1=c(NA,rep(2,ti-1)),NadSurv=c(NA,rep(10,ti-1)),
#                        Nadimm=c(NA,rep(2,ti-1)),N1prior=2,NadSurvprior=10,Nadimmprior=2)}

# Set the parameters: ####

#params <- c('phij','phia','phijM','phiaM','fec','im','lambda',
#                'p','pM','NadSurv','Nadimm','Ntot','MEPHJUF',
#                'MEPHADF','MEFE','MEIM_H','MEIM_L','MEPHJUM',
#                'MEPHADM','v','bp')

params <- c('N1','Ntot','MEPHJUF',
            'MEPHADF','MEFE','MEIM_H','MEIM_L','MEPHJUM',
            'MEPHADM','v','bp','fec')


# Run the MCMC: ####

tstart=proc.time()
mod <- jags.model('AbadiJPE.jag.R',data,inits,n.chains=cha,n.adapt=ada)
temp=proc.time()-tstart
tend <- temp

output1 <- coda.samples(mod,params,n.iter=iter,thin=th)
summary(output1)

output2 <- jags.samples(mod,params,n.iter=iter)
summary(output2)


frame1 <- do.call(rbind.data.frame, output1)
frame2 <- do.call(rbind.data.frame, output2)

str(output1) # List of 2
str(output2) # List of 12


plot(output1[,2])
output1[j] [,"par_name"]
# it returns the j-th chain for the parameter of interest
plot(output1[1][,2])
plot(output1[2][,2])
plot(output1[1][,2])

summary(output1)

mat1 <- as.matrix(output1)