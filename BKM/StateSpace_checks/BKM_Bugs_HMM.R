model{ 
  # THE COVARIATES ####
  # Define the priors for the logistic regression parameters
  alpha1 ~ dnorm(0,0.01)
  alphaa ~ dnorm(0,0.01)
  alphar ~ dnorm(0,0.01)
  beta1 ~ dnorm(0,0.01)
  betaa ~ dnorm(0,0.01)
  betar ~ dnorm(0,0.01)

  # Define the observation error prior
  sigy <- 1/tauy
  tauy ~ dgamma(0.001,0.001)
  
  # Define the logistic regression equations
  for(t in 1:T){
    logit(phi1[t]) <- alpha1 + beta1*f[t] 
    logit(phia[t]) <- alphaa + betaa*f[t]
    log(rho[t]) <- alphar + betar*stdT[t] 
  }
  
  # Define the initial population priors
  for(t in 1:2){
    Na[t] ~ dbin(0.5,200) # --> 100+-50
  }
  
  for (t in 3:T){
    # Use a discrete Uniform prior so the only influence on the posterior distr is the Upper limit
    Na[t] ~ dcat(Na_prior[]) # Na_prior = rep(1/(Up+1), Up+1); entered as data
  }
  
  C <- 1000000

  for(t in 3:T){ 
    zeros[t] ~ dpois(phi[t])
    phi[t] <- -loglik[t] + C
    
    loglam[t-1] <- (log(Na[t-1-1]) + log(rho[t-1-1]) + log(phi1[t-1-1])) 

    for (i in 0:(N_max-1)){  
      G[i+1,t] <- exp(-exp(loglam[t-1]) + i*loglam[t-1] - logfact(i))
  
      P[i+1,t] <- ifelse((i + Na[t-1] - Na[t])>0,
                         exp(Na[t]*log(phia[t-1]) + (i + Na[t-1] - Na[t])*log(1-phia[t-1]) + logfact(i + Na[t-1]) - logfact(abs(i + Na[t-1] - Na[t])) - logfact(Na[t])),
                         0)
    } 
    
    G[(N_max+1),t] <- max(0,1- sum(G[1:(N_max),t]))
    P[(N_max+1),t] <- ifelse((N_max + Na[t-1] - Na[t])>0,
                             exp(Na[t]*log(phia[t-1]) + (N_max + Na[t-1] - Na[t])*log(1-phia[t-1]) + logfact(N_max + Na[t-1]) - logfact(abs(N_max + Na[t-1] - Na[t])) - logfact(Na[t])),
                             0)
    loglik[t] <- log(sum(G[,t] * P[,t])) 
  }
  
  # Define the observation process for the census/index data
  for(t in 3:T){
    y[t] ~ dnorm(Na[t],tauy)
  }
}