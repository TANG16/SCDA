
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
