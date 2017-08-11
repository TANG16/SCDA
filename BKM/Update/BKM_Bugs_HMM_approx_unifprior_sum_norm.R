# var bin[N_bin+1], Na_prior[Up+1], loglam[T-1], Na[T], G[(N_bin+1),T], P[(N_bin+1),T]; 

model{ 
  
  # THE STATE SPACE MODEL ####
   # Define the observation process for the census/index data
  for(t in 3:T){
    y[t] ~ dnorm(Na[t],tauy)
  }
  
  # 0-1 trick for the states process, the obervation process can stay the same
  
  for(t in 3:T){ 
    zeros[t] ~ dpois(phi[t])
    phi[t] <- -loglik[t] + C
    loglam[t-1] <- (log(Na[t-2]) + log(rho[t-2]) + log(phi1[t-2])) 
    
    for (i in 0:(N_bin-1)){

      midbin[i+1,t] <- qnorm(mid[i+1],exp(loglam[t-1]),1/exp(loglam[t-1]))
      
# G is now simply a constant for all values as we are using the quantiles of the distribution 
# and a normal approximation (so they are exact).
 
      P[i+1,t] <- dbin(Na[t],phia[t-1],round(midbin[i+1,t]) + Na[t-1])
    }
    
    loglik[t] <- log(sum(P[,t])) 
    # no G matrix since all values simply constant and so goes into constant of likelihood
  }
  
  C <- 1000000 

  # Define the initial population priors
  for(t in 1:2){    
    Na[t]  <- round(Na_cont[t])
    Na_cont[t] ~ dunif(0.5, 2009+0.5)
  }
  
  for (t in 3:T){
    # Use a discrete Uniform prior so the only influence on the posterior distr is the Upper limit
    # Na[t] ~ dcat(Na_prior[]) # Na_prior = rep(1/(Up+1), Up+1); entered as data
    Na[t]  <- round(Na_cont[t]) 
    Na_cont[t] ~ dunif(0.5, Up+0.5)
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

  
  
  # Define the logistic regression equations
  # for(t in 1:(T-1)){
  for(t in 1:(T-1)){
    logit(phi1[t]) <- alpha1 + beta1*f[t] # corresponds to the year 1963
    logit(phia[t]) <- alphaa + betaa*f[t]
    log(rho[t]) <- alphar + betar*stdT[t] # We assume here that t=1
    logit(lambda[t]) <- alphal + betal*stdT[t]
  }
  
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
  
}