inits <- function(){list(v=rnorm(7),bp=rnorm((ti-1)),fec=runif((ti-1),0,5),
                         NadSurv=c(NA,round(runif(ti-1,5,10),0)),
                         NadImm=c(NA,round(runif(ti-1,5,10),0)),
                         NadSurvprior=20,NadImmprior=5)}

# params <- c('NadImm','NadSurv','v','bp','fec')
params <- c('NadImm','NadSurv','v')
# params3 <- c('Q')
# params4 <- c('P')
