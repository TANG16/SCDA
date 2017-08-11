# Model

model{ 
  #  Define the priors for the logistic regression parameters alpha1 ~ dnorm(0,0.01)
  alpha[1] ~ dnorm(0,0.01)
  alpha[2] ~ dnorm(0,0.01)
  alpha[3] ~ dnorm(0,0.01)
  alpha[4] ~ dnorm(0,0.01)
  alphal ~ dnorm(0,0.01)
  beta[1] ~ dnorm(0,0.01) 
  beta[2] ~ dnorm(0,0.01) 
  beta[3] ~ dnorm(0,0.01) 
  beta[4] ~ dnorm(0,0.01)
  betal ~ dnorm(0,0.01)

#  Define the logistic regression equations:
    
   for (t in 1:(T-1)) {
      logit(phi1[t]) = alpha[1] + beta[1]*f[t]
      logit(phi2[t]) = alpha[2] + beta[2]*f[t]
      logit(phi3[t]) = alpha[3] + beta[3]*f[t]
      logit(phi4[t]) = alpha[4] + beta[4]*f[t]
   }
# We assume here that t=1 corresponds to the year 1927 
  for (t in 1:(ring1)) {
      logit(lambda[t]) = alphal + betal*(t)
  }

# Define the recovery likelihood 

  for(t in 1:ring1){
    m[t,1:(ring1+1)] ~ dmulti(p[t, ], rel[t]) 
    
    year[t] <- t + 28
  }
    
# Calculate the cell probabilities for the recovery table 
  
  for(t in 1:ring1){
    
# Calculate the diagonal
    
    p[t, t] = lambda[t] * (1-phi1[year[t]])
    
    }

# Calculate value one above the diagonal

  for(t in 1:(ring1-1)){
    p[t, t+1] = lambda[t+1] * phi1[year[t]]*(1-phi2[(year[t]+1)])
  }

# Calculate value two above the diagonal

  for (t in 1:(ring1-2)) {
    p[t,t+2] = lambda[t+2] * phi1[year[t]]*phi2[(year[t]+1)]*(1-phi3[(year[t]+2)])
  }
  
  # Calculate remaining terms above diagonal 
  for (t in 1:(ring1-3)) {
    p[t,t+3] = lambda[t+3] * phi1[year[t]]*phi2[(year[t]+1)]*phi3[(year[t]+2)]*(1-phi4[(year[t]+3)])
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
  
}
