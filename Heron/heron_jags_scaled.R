# Model

model{ 
#  Define the priors for the logistic regression parameters alpha1 ~ dnorm(0,0.01)
  
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

#  Define the logistic regression equations:
    
   for (t in 1:(T-1)) {

      logit(phi1[t]) = alpha[1] + beta[1]*f[t]
      logit(phi2[t]) = alpha[2] + beta[2]*f[t]
      logit(phi3[t]) = alpha[3] + beta[3]*f[t]
      logit(phi4[t]) = alpha[4] + beta[4]*f[t]
      log(rho[t]) = alpharho 
    
# We assume here that t=1 corresponds to the year 1927 
      
      logit(lambda[t]) = alphal + betal*(time[t+1])

      }

#  Define the initial population priors for (a in 1:4){
     
    # X1[1] ~ dpois(2500)
    # X2[1] ~ dpois(1500) 
    # X3[1] ~ dpois(1000)
    # X4[1] ~ dpois(4000) 
    X1[1] ~ dpois(25)
    X2[1] ~ dpois(15) 
    X3[1] ~ dpois(10)
    X4[1] ~ dpois(40) 

#  Define the system process for the count data 

    for(t in 2:T){
  
#    mean1[t] <- rho[t-1]*phi1[t-1]*X4[t-1] 
#    mean2[t] <- phi2[t-1]*X1[t-1]
#    mean3[t] <- phi3[t-1]*X2[t-1]
#    mean4[t] <- phi4[t-1]*(X3[t-1]+X4[t-1])
    
#    tau1[t] <- 1/(rho[t-1]*phi1[t-1]*X4[t-1])
#    tau2[t] <- 1/(X1[t-1]*phi2[t-1]*(1-phi2[t-1]))
#    tau3[t] <- 1/(X2[t-1]*phi3[t-1]*(1-phi3[t-1]))
#    tau4[t] <- 1/((X3[t-1]+X4[t-1])*phi4[t-1]*(1-phi4[t-1]))
    
      X1[t] ~ dpois(rho[t-1]*phi1[t-1]*X4[t-1])
      X2[t] ~ dbin(phi2[t-1],X1[t-1])
      X3[t] ~ dbin(phi3[t-1],X2[t-1])
      X4[t] ~ dbin(phi4[t-1],(X3[t-1]+X4[t-1]))
    }
    
# Define the observation process for the count data 

    for(t in 2:T){
      y[t] ~ dnorm((X2[t]+X3[t]+X4[t]),tauy) 
#      y[t] ~ dpois(X2[t]+X3[t]+X4[t]) 
    }

# Define the recovery likelihood 

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
  
  for(t1 in 1:(ring1-3)){
    for(t2 in (t1+3):ring1){
      for (t in (t1+1):(t2-1)) {
         lphi[t1, t2, t] = log(phi4[year[t]])
      }
  # Probabilities in table
  p[t1,t2] = lambda[t2]*phi1[year[t1]]*(1-phi4[year[t2]])*
              exp(sum(lphi[t1,t2,(t1+1):(t2-1)])) 
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
