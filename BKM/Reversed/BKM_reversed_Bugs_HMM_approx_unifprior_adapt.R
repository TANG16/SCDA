var bin[N_bin], N1_prior[Up+1], loglam[N_bin,T-1], N1[T], G[N_bin,N_bin,T], P[N_bin,T], Q[N_bin,T]; 


model{ 
  
  # THE STATE SPACE MODEL ####
   # Define the observation process for the census/index data
#   for(t in 3:T){
#     y[t] ~ dnorm(Na[t],tauy)    
#   }
  
  # 0 trick for the observationa and the state process  
  for(t in 3:T){ 
    zeros[t] ~ dpois(phi[t])
    phi[t] <- -loglik[t] + C
#     sumG[t] <- sum(G[,t])
#     loglik[t] <- log(sum((G[,t]/sumG[t]) * P[,t])) # piecewise multiplication enough here
    loglik[t] <- log(sum(G[,,t-1] %*% P[,t])) + log(sum(G[,,t] %*% Q[,t]))
  } 

  for(t in 3:T){   
    #     loglam[t-1] <- (log(Na[t-1-1]) + log(rho[t-1-1]) + log(phi1[t-1-1])) 
    for (i in 1:N_bin){ 
#       G[i,t] <- exp(-exp(loglam[t-1]) + bin[i]*loglam[t-1] - logfact(bin[i]))
      P[i,t] <- exp(-exp(loglam[i,t-1]) + N1[t]*loglam[i,t-1] - logfact(N1[t]))

#       P[i,t] <- ifelse((bin[i] + Na[t-1] - Na[t])>0,
#                          exp(Na[t]*log(phia[t-1]) 
#                              + (bin[i] + Na[t-1] - Na[t])*log(1-phia[t-1]) 
#                              + logfact(bin[i] + Na[t-1]) 
#                              - logfact(abs(bin[i] + Na[t-1] - Na[t])) 
#                              - logfact(Na[t])),
#                          0)
      for (j in 1:N_bin){
        G[j,i,t] <- ifelse((N1[t-1] + bin[j] - bin[i])>0,
                         exp(bin[i]*log(phia[t-1]) 
                             + (N1[t-1] + bin[j] - bin[i])*log(1-phia[t-1]) 
                             + logfact(N1[t-1] + bin[j]) 
                             - logfact(abs(N1[t-1] + bin[j] - bin[i])) 
                             - logfact(bin[i])),
                         0)
      }
      Q[i,t] = dnorm(y[t],bin[i],tauy)
    } 
  }
  
  for(t in 3:T){ 
    for (i in 1:N_bin){ 
      loglam[i,t-1] <- (log(bin[i]) + log(rho[t-1]) + log(phi1[t-1])) 
    }
  }
  

  for (i in 1:N_bin){
    for (j in 1:N_bin){
      G[j,i,2] <- 1/N_bin
    }
  }

  C <- 1000000 

  # Define the initial population priors
  for(t in 1:2){    
#     Na[t]  <- round(Na_cont[t])
#     Na_cont[t] ~ dunif(961+ 0.5, 1039+0.5)
    N1[t]  <- round(N1_cont[t])
    N1_cont[t] ~ dunif(961/5+ 0.5, 1039/5+0.5)    
    # Na[t] ~ dbin(0.5,200) 
  }
  
  for (t in 3:T){
    # Use a discrete Uniform prior so the only influence on the posterior distr is the Upper limit
    # Na[t] ~ dcat(Na_prior[]) # Na_prior = rep(1/(Up+1), Up+1); entered as data
#     Na[t]  <- round(Na_cont[t]) 
#     Na_cont[t] ~ dunif(0.5, Up+0.5)
    N1[t]  <- round(N1_cont[t]) 
    N1_cont[t] ~ dunif(0.5, Up+0.5)    
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
  for(t in 1:T){
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
  
  
  # Bins' midpoints ####
  for (i in 1:N_bin){
    bin[i] <- bin_size*(i-1) + (bin_size-1)/2  # ith bin's midpoint
  }
  
}