var bin[N_bin+1], Na_prior[Up+1], loglam[T-1], Na[T], G[(N_bin+1),T], P[(N_bin+1),T]; 

model{ 
  # Bins' midpoints
  for (i in 0:(N_bin)){
    bin[i+1] <- 0.5*(bin_size*(2*i+1)-1) # ith bin's midpoint
  }
  for (i in 0:(Up)){
    Na_prior[i+1] <- 1/(Up+1)
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
  
  # Define the logistic regression equations
  for(t in 1:T){
    logit(phi1[t]) <- alpha1 + beta1*f[t] 
    logit(phia[t]) <- alphaa + betaa*f[t]
    log(rho[t]) <- alphar + betar*stdT[t] 
    logit(lambda[t]) <- alphal + betal*stdT[t]
  }
  
  for (t in 3:T){
    # Use a discrete Uniform prior so the only influence on the posterior distr is the Upper limit
    Na[t] ~ dcat(Na_prior[]) # Na_prior = rep(1/(Up+1), Up+1); entered as data
  }
  
  # THE STATE SPACE MODEL ####
  # 0-1 trick for the states process
  C <- 1000000

  # Define the initial priors
  Na_init <- round(2000/sc)
  for(t in 1:2){
    Na[t] ~ dbin(0.5,Na_init) 
  }

  for(t in 3:T){ 
    zeros[t] ~ dpois(phi[t])
    phi[t] <- -loglik[t] + C
    
    loglam[t-1] <- (log(Na[t-1-1]) + log(rho[t-1-1]) + log(phi1[t-1-1])) 
   
    for (i in 0:(N_bin)){ 
      G[i+1,t] <- exp(-exp(loglam[t-1]) + bin[i+1]*loglam[t-1] - logfact(bin[i+1]))      
      P[i+1,t] <- ifelse((bin[i+1] + Na[t-1] - Na[t])>0,
                         exp(Na[t]*log(phia[t-1]) + (bin[i+1] + Na[t-1] - Na[t])*log(1-phia[t-1]) + logfact(bin[i+1] + Na[t-1]) - logfact(abs(bin[i+1] + Na[t-1] - Na[t])) - logfact(Na[t])),
                         0)
    } 
    sumG[t] <- sum(G[,t])
    loglik[t] <- log(sum((G[,t]/sumG[t]) * P[,t])) 
  }
  
  # Define the observation process
  for(t in 3:T){
    y[t] ~ dnorm(Na[t],tauy)
  }
  
 # MORE CODE HERE
}