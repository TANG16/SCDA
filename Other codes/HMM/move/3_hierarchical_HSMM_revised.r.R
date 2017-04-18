

## Revised version 04 February 2015
## The previous version was very inefficient due to unnecessary loops. 
## This code produces the same estimates as the previous version in a substantially shorter amount of time (the speed increases by a factor of about 10).



## ESTIMATION OF HIERARCHICAL HSMM 
## observations ("OBS") need to be given in an n x (2*p) matrix, with p being the number of individuals, 
## and the observed step lengths and turning angles for individual i given in columns 2*(i-1)+1 and 2*(i-1)+2, respectively


library(CircStats)
library(boot)
             
               
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
  ta <- c(a[1:2],log(a[3:4]))
  tb <- log(b)
  tkappa <- logit(kappa)
  tgam <- c(log(gam[1:2]),log(gam[3:4]/(1-gam[3:4])))
  tco <- c(log((co[1])/(2*pi-co[1])),log((pi+co[2])/(pi-co[2])))
  parvect <- c(ta,tb,tkappa,tgam,tco)
  return(parvect)
}
            
## inverse transformation back to the natural parameter space            
move.HSMM.pw2pn <- function(parvect){
  epar <- exp(parvect)
  a <- c(parvect[1:2],epar[3:4])
  b <- epar[5:6]
  kappa <- inv.logit(parvect[7:8])
  gam <- c(exp(parvect[9:10]),exp(parvect[11:12])/(exp(parvect[11:12])+1))    
  co <- c(2*pi*inv.logit(parvect[13]),pi*(exp(parvect[14])-1)/(exp(parvect[14])+1))
  return(list(a=a,b=b,kappa=kappa,gam=gam,co=co))
}

## function that computes minus the log-likelihood
move.HSMM.mllk <- function(parvect,OBS){
  lpn <- move.HSMM.pw2pn(parvect)
  gamma <- gen.Gamma(m,lpn$gam)
  delta <- solve(t(diag(sum(m))-gamma+1),rep(1,sum(m)))
  mllk.all <- 0
  K <- 25
  a1 <- seq(re1min,re1max,length=K)
  a2 <- seq(re2min,re2max,length=K)
  a.m1 <- (a1[-1]+a1[-K])*0.5   
  a.m2 <- (a2[-1]+a2[-K])*0.5   
  am <- rbind(exp(a.m1),exp(a.m2)) 
  for (ani in 1:n.ind){  
    obs <- OBS[,((ani-1)*2+1):((ani-1)*2+2)]
    n <- max(which(!is.na(obs[,1])))
    obs <- obs[1:n,]
    one.animal <- matrix(rep(NA,(K-1)^2*2),ncol=2)
    for (j1 in 1:(K-1)){                 
    for (j2 in 1:(K-1)){                
      allprobs <- matrix(rep(NA,sum(m)*n),nrow=n)
      ind.step <- which(!is.na(obs[,1]))
      ind.angle <- which(!is.na(obs[,2]))
      for (j in 1:2){
        step.prob <- angle.prob <- rep(1,n)
        angle.prob[ind.angle] <- dwrpcauchy(obs[ind.angle,2],mu=lpn$co[j],rho=lpn$kappa[j])
        if (j==1) scale.re<-am[1,j1] else scale.re<-am[2,j2]
        step.prob[ind.step] <- dweibull(obs[ind.step,1],shape=lpn$b[j],scale=scale.re)
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
    llk <- lscale
    stan.fac1 <- pnorm(max(a1),lpn$a[1],lpn$a[3])-pnorm(min(a1),lpn$a[1],lpn$a[3])
    stan.fac2 <- pnorm(max(a2),lpn$a[2],lpn$a[4])-pnorm(min(a2),lpn$a[2],lpn$a[4])
    r <- (j1-1)*(K-1)+1+(j2-1)
    one.animal[r,1] <- llk
    one.animal[r,2] <- stan.fac1^(-1)*stan.fac2^(-1)*
      (pnorm(a1[j1+1],lpn$a[1],lpn$a[3])-pnorm(a1[j1],lpn$a[1],lpn$a[3]))*
      (pnorm(a2[j2+1],lpn$a[2],lpn$a[4])-pnorm(a2[j2],lpn$a[2],lpn$a[4]))
    }}                                       
    ma <- max(one.animal[,1]) 
    for (goo in 1:dim(one.animal)[1]){
      ifelse(ma-one.animal[goo,1]<700,one.animal[goo,1]<-exp(one.animal[goo,1]-ma),one.animal[goo,1]<-0)
    } # to avoid underflow
    mllk.all <- mllk.all-(log(sum(one.animal[,1]*one.animal[,2]))+ma)  
  } 
  return(mllk.all)
}
 
## function that runs the numerical maximization of the above likelihood function and returns the results
move.HSMM.mle <- function(OBS,a0,b0,kappa0,gamma0,co0){
  parvect <- move.HSMM.pn2pw(a0,b0,kappa0,gamma0,co0)
  mod <- nlm(move.HSMM.mllk,parvect,OBS,print.level=2,iterlim=1000) 
  mllk <- mod$minimum
  pn <- move.HSMM.pw2pn(mod$estimate)
  list(a=pn$a,b=pn$b,kappa=pn$kappa,H=mod$hessian,gamma=pn$gam,co=pn$co,mllk=mllk)
}

## initial parameter values to be used in the numerical maximization (several different combinations should be tried in order to ensure hitting the global maximum)
a0 <- c(log(0.08),log(0.46),0.11,0.05) # associated with the random effects distributions 
b0 <- c(0.8,0.85) # Weibull shape parameters
kappa0 <- c(0.2,0.4) # wrapped Cauchy concentration parameters
co0 <- c(pi,0) # wrapped Cauchy mean parameters
gamma0 <- c(0.39,3.75,0.14,0.59)  # negative binomial state dwell-time distribution parameters

## choose range over which to integrate the random effects' distributions
re1min <- log(0.05)
re1max <- log(0.12)
re2min <- log(0.39)
re2max <- log(0.59)

## size of state aggregates
m <- c(20,10)

## run the numerical maximization
moveHSMMre <- move.HSMM.mle(OBS,a0,b0,kappa0,gamma0,co0)
moveHSMMre

