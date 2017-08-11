var  P2[N_bin1,T], P4[N_bin3,T], Q[N_bin3,T];


model{ 
# HMM ####
# Impute X2 and X4, integrate out X1 and X3
# G - transition matrix for the integrated states
# P - augmented observation matrix for the imputed states    
# Q - observation matrix
    

#     for (t in 1:T){
#       zeros[t] ~ dpois(PHI[t])
#       PHI[t] <- -loglik[t] + C
#       loglik[t] <- log(exp(log(sum(G1[,t] * P2[,t])) + log(sum(G3[,t] * P4[,t] * Q[,t]))))
#     }
    
  for (t in 1:T){
    zeros[t] ~ dpois(PHI[t])
    PHI[t] <- -loglik[t] + C
    loglik[t] <- log(exp(log(sum(P2[,t])) + log(sum(P4[,t] * Q[,t]))))
    # no G matrices since all values simply constant and so go into constant of likelihoo
  }


  for(t in 3:T){ 
  
    for (i in 1:N_bin1){
      # loglam[t] = log(rho[t] * phi[t] * X4[t])
      # X1[t] ~ Poisson(exp(loglam[t-1])) 
      midbin1[i,t] <- qnorm(mid1[i],exp(loglam1[t-2]),1/exp(loglam1[t-2])) 
      
      # G is now simply a constant for all values as we are using the quantiles of the distribution 
      # and a normal approximation (so they are exact).
      
      # X2[t] ~ Binomial(X1[t-1], phi2[t-1])
      P2[i,t] <- dbin(X2[t],phi2[t-1],round(midbin1[i,t]))
    }
    
    
    for (i in 1:N_bin3){
      # X3[t] ~ Binomial(X2[t-1], phi3[t-1])
      midbin3[i,t] <- qnorm(mid3[i],X2[t-2]*phi3[t-2],1/(X2[t-2]*phi3[t-2]*(1-phi3[t-2])))
      
      # G is now simply a constant for all values as we are using the quantiles of the distribution 
      # and a normal approximation (so they are exact).
      
      # X4[t] ~ Binomial(X3[t-1]+X4[t-1], phi4[t-1])
      P4[i,t] <- dbin(X4[t], phi4[t-1], X4[t-1] + round(midbin3[i,t]))
      Q[i,t] <- sqrt(tauy)*exp(-0.5*tauy*pow((y[t] - (X2[t] + round(midbin3[i,t]) + X4[t])),2))/sqrt(2*pi)
    }  
  } 
 
  # prior for first transition/augmented observations/observation probabilities
  for (t in 1:2){
    for (i in 0:(N_bin1-1)){
      P2[i+1,t] <- 1/Up2 #X2_prior[1]
    }
    for (i in 0:(N_bin3-1)){
      P4[i+1,t] <- 1/Up4 #X4_prior[1]
      Q[i+1,t] <- 1 #sqrt(tauy)*exp(-0.5*tauy*pow((y[t] - (X2[t] + bin3[i+1] + X4[t])),2))/sqrt(2*pi)
    }    
  }   
  

  for (t in 1:(T-1)){
    loglam1[t] <- log(X4[t]) + log(rho[t]) + log(phi1[t])
  }

  for (t in 1:T){
    # Use a discrete Uniform prior so the only influence on the posterior distr is the Upper limit
    # X2[t] ~ dcat(X2_prior[]) # X2_prior = rep(1/(Up2+1), Up2+1); entered as data
    # X4[t] ~ dcat(X4_prior[]) # X4_prior = rep(1/(Up4+1), Up4+1); entered as data
    X2[t]  <- round(X2_cont[t]) 
    X2_cont[t] ~ dunif(0.5, Up2+0.5)
    X4[t]  <- round(X4_cont[t]) 
    X4_cont[t] ~ dunif(0.5, Up4+0.5)
  }  
   
  C <- 1000000
  pi <- 3.14159265359

# Define the recovery likelihood ####

  for(t in 1:ring1){
    m[t,1:(ring1+1)] ~ dmulti(p[t, ], rel[t]) 
    year[t] <- t + 27
  }
    
# Calculate the cell probabilities for the recovery table 
  
  for(t in 1:ring1){
    # Calculate the diagonal
    p[t, t] = lambda[t] * (1-phi1[year[t]])
  }

  # Calculate value one above the diagonal
  for(t in 1:(ring1-1)){
    p[t, t+1] = lambda[t+1] * phi1[year[t]]*(1-phi2[year[t]+1])
  }

  # Calculate value two above the diagonal
  for (t in 1:(ring1-2)) {
    p[t,t+2] = lambda[t+2] * phi1[year[t]]*phi2[year[t]+1]*(1-phi3[year[t]+2])
  }
  
  # Calculate remaining terms above diagonal 
  for (t in 1:(ring1-3)) {
    p[t,t+3] = lambda[t+2] * phi1[year[t]]*phi2[year[t]+1]*phi3[year[t]+3]*(1-phi4[year[t]+3])
  }
  
  for(t1 in 1:(ring1-4)){
    for(t2 in (t1+4):ring1){
      # for (t in (t1+1):(t2-1)) {
      #    lphi[t1, t2, t] = log(phi4[year[t]])
      # }
      # Probabilities in table
      p[t1,t2] = lambda[t2]*phi1[year[t1]]*phi2[year[t1]+1]*phi3[year[t1]+2]*(1-phi4[year[t2]])*
        exp(sum(log(phi4[year[(t1+3):(t2-1)]])))
      # exp(sum(lphi[t1,t2,(t1+1):(t2-1)])) 
    }
  }
  
  for(t1 in 1:ring1) {
    for(t2 in 1:(t1-1)){
      # Zero probabilities in lower triangle of table 
      p[t1, t2] <- 0
    }
    # Probability of an animal never being seen again 
    p[t1, (ring1+1)] <- 1 - sum(p[t1, 1:ring1])
  }
    
    
#  Define the logistic regression equations: ####
  
  for (t in 1:(T-1)) {
    logit(phi1[t]) = alpha[1] + beta[1]*f[t]
    logit(phi2[t]) = alpha[2] + beta[2]*f[t]
    logit(phi3[t]) = alpha[3] + beta[3]*f[t]
    logit(phi4[t]) = alpha[4] + beta[4]*f[t]
    log(rho[t]) = alpharho 
  }
  # We assume here that t=1 corresponds to the year 1927 
  for (t in 1:(ring1)) {
    logit(lambda[t]) = alphal + betal*(time[t])
  }    
  
  
#  Define the priors for the logistic regression parameters alpha1 ~ dnorm(0,0.01) ####
  alpha[1] ~ dnorm(0,0.01)
  alpha[2] ~ dnorm(0,0.01)
  alpha[3] ~ dnorm(0,0.01)
  alpha[4] ~ dnorm(0,0.01)
  alphal ~ dnorm(0,0.01)
  alpharho ~ dnorm(0,0.01)
  beta[1] ~ dnorm(0,0.01) 
  beta[2] ~ dnorm(0,0.01) 
  beta[3] ~ dnorm(0,0.01) 
  beta[4] ~ dnorm(0,0.01)
  betal ~ dnorm(0,0.01)
  
  #  Define the observation error prior sigy <- 1/tauy - initially assume N.
  
  tauy ~ dgamma(0.001,0.001)
  sigmay <- 1/tauy  
  
}

