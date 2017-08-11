#  MODEL DESCRIPTION
#***********************************************************************
### var G[N_max,N_max,T], P[N_max,N_max,T];
var loglam[T-1], Na[T], G[(N_max+1),T], P[(N_max+1),T]; # time verying vectors enought in our case

# Model

model{ 
  # THE COVARIATES ####
  # Define the priors for the logistic regression parameters
  alpha1 ~ dnorm(0,0.01)
  alphaa ~ dnorm(0,0.01)
  alphar ~ dnorm(0,0.01)
  alphal ~ dnorm(0,0.01)
  beta1 ~ dnorm(0,0.01)
  betaa ~ dnorm(0,0.01)
  betar ~ dnorm(0,0.01)
  betal ~ dnorm(0,0.01)
  
  # Define the observation error prior
  sigy <- 1/tauy
  tauy ~ dgamma(0.001,0.001)
  
  # Define the logistic regression equations
  # for(t in 1:(T-1)){
  for(t in 1:T){
    logit(phi1[t]) <- alpha1 + beta1*f[t] # corresponds to the year 1963
    logit(phia[t]) <- alphaa + betaa*f[t]
    log(rho[t]) <- alphar + betar*stdT[t] # We assume here that t=1
    logit(lambda[t]) <- alphal + betal*stdT[t]
  }
  
  # THE STATE SPACE MODEL ####
  
  # Define the initial population priors
  for(t in 1:2){
    #   N1[t] ~ dpois(200)
    # Na[t] ~ dbin(0.5,2000) # --> 1000+-500
    Na[t] ~ dbin(0.5,200) # --> 100+-50
  }
  
  #####
  # Zero trick for loglik of HMM ####
  #####
  for (t in 3:T){
    # Use a discrete Uniform prior so the only influence on the posterior distr is the Upper limit
    Na[t] ~ dcat(Na_prior[]) # Na_prior = rep(1/(Up+1), Up+1); entered as data
  }
  
  C <- 1000000
 
  for(t in 3:T){ 
    zeros[t] ~ dpois(phi[t])
    phi[t] <- -loglik[t] + C
    
    loglam[t-1] <- (log(Na[t-1-1]) + log(rho[t-1-1]) + log(phi1[t-1-1])) 

    for (i in 0:(N_max-1)){  # from 0!!! 
      G[i+1,t] <- dpois(i,exp(loglam[t-1]))
      P[i+1,t] <- ifelse((i + Na[t-1] - Na[t])>0,
                         dbin(min(Na[t],i+Na[t-1]),phia[t-1],i+Na[t-1]),
                         0)
    } 
    
    G[(N_max+1),t] <- max(0, 1- sum(G[1:(N_max),t]))
    P[(N_max+1),t] <- ifelse((N_max + Na[t-1] - Na[t])>0,
                             dbin( min(Na[t],N_max+Na[t-1]),phia[t-1],N_max+Na[t-1]),
                             0)
  
    loglik[t] <- log(sum(G[,t] * P[,t])) # piecewise multiplication enough here
  }
  
  # Define the observation process for the census/index data
  for(t in 3:T){
    y[t] ~ dnorm(Na[t],tauy)
  }
  
  # THE RECOVERY MODEL ####
  # Define the recovery likelihood 
  for(t in 1:T1){
    m[t, 1:(T2+1)] ~ dmulti(p[t,], rel[t])
  }
  
  # Calculate the cell probabilities for the recovery table 
  
  for(t1 in 1 : (T1-1)){
    # Calculate the diagonal
    p[t1, t1] <- lambda[t1]*(1-phi1[t1])
    
    # Calculate value one above the diagonal
    p[t1, t1+1] <- lambda[t1+1]* phi1[t1]*(1-phia[t1+1])
    
    # Calculate remaining terms above diagonal
    for(t2 in (t1+2):T2){
      for(t in (t1+1):(t2-1)){
        lphi[t1, t2, t] <- log(phia[t])
      }
      # Probabilities in table
      p[t1,t2] <- lambda[t2]*phi1[t1]*(1-phia[t2])*exp(sum(lphi[t1,t2,(t1+1):(t2-1)]))
    }
    for(t2 in 1:(t1-1)){
      # Zero probabilities in lower triangle of table
      p[t1, t2] <- 0
    }
    # Probability of an animal never being seen again
    p[t1, T2+1] <- 1 - sum(p[t1,1:T2])	
  }
  
  # Final row
  p[T1,T1] <- lambda[T1]*(1-phi1[T1])
  
  for(t in 1:(T1-1)){
    p[T1,t] <- 0
  }
  p[T1,T1+1] <- 1 - p[T1,T1]
}