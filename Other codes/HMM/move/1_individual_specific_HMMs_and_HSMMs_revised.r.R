## Revised version 04 February 2015
## The previous version was very inefficient due to unnecessary loops. 
## This code produces the same estimates as the previous version 
# in a substantially shorter amount of time (the speed increases by a factor of about 30).


## this code has two parts: 
# 1.) ESTIMATION OF (INDIVIDUAL-SPECIFIC) HMM and 
# 2.) ESTIMATION OF (INDIVIDUAL-SPECIFIC) HSMM
## observations ("obs") need to be given in an n x 2 matrix, with first column giving the step lengths and second column giving the turning angles

## for illustration purposes the data stored in the file "obs.txt" (in the supplementary material to the paper) may be used: 
## read.table(".../obs.txt",header=T)->obs
## this file comprises the observed step lengths and turning angles for one of the bison


setwd("Other codes/HMM")

library(CircStats)
library(boot)

read.table("4_obs.txt",header=T)->obs

######## 1.) fit hidden Markov movement model (Weibull & wrapped Cauchy) ########

## function that transforms each of the (possibly constrained) parameters to the real line
move.HMM.pn2pw <- function(a,b,kappa,gamma,co){
  ta <- log(a)
  tb <- log(b)
  tkappa <- logit(kappa)
  tgamma <- logit(gamma)
  tco <- c(log((co[1])/(2*pi-co[1])),log((pi+co[2])/(pi-co[2])))
  parvect <- c(ta,tb,tkappa,tgamma,tco)
  return(parvect)
}

## inverse transformation back to the natural parameter space
move.HMM.pw2pn <- function(parvect){
  epar <- exp(parvect)
  a <- epar[1:2]
  b <- epar[3:4]
  kappa <- inv.logit(parvect[5:6])
  Gamma <- matrix(rep(NA,4),nrow=2) 
  Gamma[1,1] <- inv.logit(parvect[7])
  Gamma[1,2] <- 1-inv.logit(parvect[7])
  Gamma[2,2] <- inv.logit(parvect[8])
  Gamma[2,1] <- 1-inv.logit(parvect[8])
  delta <- solve(t(diag(2)-Gamma+1),c(1,1))             
  co <- c(2*pi*inv.logit(parvect[9]),pi*(epar[10]-1)/(epar[10]+1))
  return(list(a=a,b=b,kappa=kappa,Gamma=Gamma,delta=delta,co=co))
}                 

## function that computes minus the log-likelihood
move.HMM.mllk <- function(parvect,obs){
  n <- dim(obs)[1]
  lpn <- move.HMM.pw2pn(parvect)
  allprobs <- matrix(rep(NA,2*n),nrow=n)
  ind.step <- which(!is.na(obs[,1]))
  ind.angle <- which(!is.na(obs[,2]))
  for (j in 1:2){
    step.prob <- angle.prob <- rep(1,n)
    angle.prob[ind.angle] <- dwrpcauchy(obs[ind.angle,2],mu=lpn$co[j],rho=lpn$kappa[j])
    step.prob[ind.step] <- dweibull(obs[ind.step,1],shape=lpn$a[j],scale=lpn$b[j])
    allprobs[,j] <- angle.prob*step.prob
  } # j index
  foo <- lpn$delta 
  lscale <- 0
  for (i in 1:n){
    foo <- foo%*%lpn$Gamma*allprobs[i,]  
    sumfoo <- sum(foo)
    lscale <- lscale+log(sumfoo)
    foo <- foo/sumfoo
  }
  mllk <- -lscale
  mllk
}

## function that runs the numerical maximization of the above likelihood function and returns the results
move.HMM.mle <- function(obs,a0,b0,kappa0,gamma0,co0){
  parvect <- move.HMM.pn2pw(a0,b0,kappa0,gamma0,co0)
  mod <- nlm(move.HMM.mllk,parvect,obs,print.level=2)
  mllk <- mod$minimum
  pn <- move.HMM.pw2pn(mod$estimate)
  list(a=pn$a,b=pn$b,kappa=pn$kappa,Gamma=pn$Gamma,delta=pn$delta,co=pn$co,mllk=mllk)
}

## initial parameter values to be used in the numerical maximization (several different combinations should be tried in order to ensure hitting the global maximum)
a0 <- c(0.8,0.9) # Weibull shape parameters
b0 <- c(0.1,0.5) # Weibull scale parameters
kappa0 <- c(0.2,0.3) # wrapped Cauchy concentration parameters
co0 <- c(pi,0) # wrapped Cauchy mean parameters
gamma0 <- c(0.8,0.8) # diagonal entries of the transition probability matrix

## run the numerical maximization
moveHMM <- move.HMM.mle(obs,a0,b0,kappa0,gamma0,co0)
moveHMM

 




 

######## 2.) fit hidden semi-Markov movement model (Weibull & wrapped Cauchy & negative binomial) ########

## function that derives the t.p.m. of the HMM that represents the HSMM (see Langrock and Zucchini, 2011) 
gen.Gamma <- function(m,gam){
  Gamma <- diag(m[1]+m[2])*0
  probs1 <- dnbinom(0:(m[1]-1),size=gam[1],prob=gam[3])
  probs2 <- dnbinom(0:(m[2]-1),size=gam[2],prob=gam[4])
  den1 <- 1-c(0,pnbinom(0:(m[1]-1),size=gam[1],prob=gam[3]))
  den2 <- 1-c(0,pnbinom(0:(m[2]-1),size=gam[2],prob=gam[4]))
  if (length(which(den1<1e-10))>0) {probs1[which(den1<1e-10)] <- 1; den1[which(den1<1e-10)] <- 1}
  if (length(which(den2<1e-10))>0) {probs2[which(den2<1e-10)] <- 1; den2[which(den2<1e-10)] <- 1}
  ## state aggregate 1
  for (i in 1:(m[1])){
    Gamma[i,m[1]+1] <- probs1[i]/den1[i]
    ifelse(i!=m[1],Gamma[i,i+1]<-1-Gamma[i,m[1]+1],Gamma[i,i]<-1-Gamma[i,m[1]+1])
  }
  ## state aggregate 2
  for (i in 1:(m[2])){
    Gamma[m[1]+i,1] <- probs2[i]/den2[i]
    ifelse(i!=m[2],Gamma[m[1]+i,m[1]+i+1]<-1-Gamma[m[1]+i,1],Gamma[m[1]+i,m[1]+i]<-1-Gamma[m[1]+i,1])
  }
 Gamma
}

## function that transforms each of the (possibly constrained) parameters to the real line
move.HSMM.pn2pw <- function(a,b,kappa,gam,co){
  ta <- log(a)
  tb <- log(b)
  tkappa <- logit(kappa)
  tgam <- c(log(gam[1:2]),logit(gam[3:4]))
  tco <- c(log((co[1])/(2*pi-co[1])),log((pi+co[2])/(pi-co[2])))
  parvect <- c(ta,tb,tkappa,tgam,tco)
  return(parvect)
}

## inverse transformation back to the natural parameter space
move.HSMM.pw2pn <- function(parvect){
  epar <- exp(parvect)
  a <- epar[1:2]
  b <- epar[3:4]
  kappa <- inv.logit(parvect[5:6])
  gam <- c(epar[7:8],inv.logit(parvect[9:10]))    
  co <- c(2*pi*inv.logit(parvect[11]),pi*(epar[12]-1)/(epar[12]+1))
  return(list(a=a,b=b,kappa=kappa,gam=gam,co=co))
}

## function that computes minus the log-likelihood of the HSMM
move.HSMM.mllk <- function(parvect,obs){
  n <- dim(obs)[1]
  lpn <- move.HSMM.pw2pn(parvect)
  allprobs <- matrix(rep(NA,(m[1]+m[2])*n),nrow=n)
  gamma <- gen.Gamma(m,lpn$gam)
  if (!stat) {delta <- rep(1/(sum(m)),sum(m))}  
  if ( stat) {delta <- solve(t(diag(sum(m))-gamma+1),rep(1,sum(m)))}
  ind.step <- which(!is.na(obs[,1]))
  ind.angle <- which(!is.na(obs[,2]))
  for (j in 1:2){
    step.prob <- angle.prob <- rep(1,n)
    angle.prob[ind.angle] <- dwrpcauchy(obs[ind.angle,2],mu=lpn$co[j],rho=lpn$kappa[j])
    step.prob[ind.step] <- dweibull(obs[ind.step,1],shape=lpn$a[j],scale=lpn$b[j])
    allprobs[,((j-1)*m[1]+1):((j-1)*m[1]+m[j])] <- rep(angle.prob*step.prob,m[j])
  } # j index
  foo <- delta  
  lscale <- 0
  for (i in 1:n){
    foo <- foo%*%gamma*allprobs[i,]  
    sumfoo <- sum(foo)
    lscale <- lscale+log(sumfoo)
    foo <- foo/sumfoo
  }
  mllk <- -lscale
  mllk
}

## function that runs the numerical maximization of the above likelihood function and returns the results
move.HSMM.mle <- function(obs,a0,b0,kappa0,gam0,co0){
  parvect <- move.HSMM.pn2pw(a0,b0,kappa0,gam0,co0)
  mod <- nlm(move.HSMM.mllk,parvect,obs,print.level=2)
  mllk <- mod$minimum
  pn <- move.HSMM.pw2pn(mod$estimate)
  gamma <- gen.Gamma(m,pn$gam)
  delta <- solve(t(diag(sum(m))-gamma+1),rep(1,sum(m)))
  delta <- c(sum(delta[1:m[1]]),1-sum(delta[1:m[1]]))
  list(a=pn$a,b=pn$b,kappa=pn$kappa,gam=pn$gam,delta=delta,co=pn$co,mllk=mllk)
}

## initial parameter values used in the numerical maximization (several different combinations should be tried in order to ensure hitting the global maximum)
a0 <- moveHMM$a  
b0 <- moveHMM$b    
kappa0 <- moveHMM$kappa    
co0 <- moveHMM$co  
gam0 <- c(1,1,moveHMM$Gamma[1,2],moveHMM$Gamma[2,1])
m <- c(30,30)

## first maximize not assuming stationarity (otherwise errors are likely; this affects only the initial distribution)
stat <- FALSE
moveHSMM <- move.HSMM.mle(obs,a0,b0,kappa0,gam0,co0) 

## then maximize assuming stationarity
stat <- TRUE
a0 <- moveHSMM$a
b0 <- moveHSMM$b 
kappa0 <- moveHSMM$kappa
gam0 <- moveHSMM$gam
co0 <- moveHSMM$co
moveHSMM <- move.HSMM.mle(obs,a0,b0,kappa0,gam0,co0)
moveHSMM
