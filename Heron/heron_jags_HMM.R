# Model
# var G1[(N_max1+1),T], P2[(N_max1+1),(N_max3+1),T], G3[(N_max3+1),T], P4[(N_max3+1),T], Q[(N_max3+1),T];
var G1[(N_max1+1),T], P2[(N_max1+1),T], G3[(N_max3+1),T], P4[(N_max3+1),T], Q[(N_max3+1),T];


model{ 
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

#  Define the logistic regression equations: ####
    
   for (t in 1:(T-1)) {

      logit(phi1[t]) = alpha[1] + beta[1]*f[t]
      logit(phi2[t]) = alpha[2] + beta[2]*f[t]
      logit(phi3[t]) = alpha[3] + beta[3]*f[t]
      logit(phi4[t]) = alpha[4] + beta[4]*f[t]
      log(rho[t]) = alpharho 
    
# We assume here that t=1 corresponds to the year 1927 
      
      logit(lambda[t]) = alphal + betal*(time[t+1])

      }

# HMM ####
# Impute X2 and X4, integrate out X1 and X3
# G - transition matrix for the integrated states
# P - augmented observation matrix for the imputed states    
# Q - observation matrix
    
    C <- 1000000
    pi <- 3.14159265359
    for (t in 1:T){
      zeros[t] ~ dpois(PHI[t])
      PHI[t] <- -loglik[t] + C
      loglik[t] <- log(exp(log(sum(G1[,t] * P2[,t])) + log(sum(G3[,t] * P4[,t] * Q[,t]))))
    }
    
    for (t in 1:(T-1)){
      loglam1[t] <- log(X4[t]) + log(rho[t]) + log(phi1[t])
    }
    
    for (t in 1:T){
      # Use a discrete Uniform prior so the only influence on the posterior distr is the Upper limit
      X2[t] ~ dcat(X2_prior[]) # X2_prior = rep(1/(Up2+1), Up2+1); entered as data
      X4[t] ~ dcat(X4_prior[]) # X4_prior = rep(1/(Up4+1), Up4+1); entered as data
    }
    
    # prior for first transition/augmented observations/observation probabilities
    for (t in 1:2){
      for (i in 0:N_max1){
        G1[i+1,t] <- 1/N_max1 # diffuse initialisation
        P2[i+1,t] <- X2_prior[1]
      }
      for (i in 0:N_max3){
        G3[i+1,t] <- 1/N_max3 # diffuse initialisation
        P4[i+1,t] <- X4_prior[1]
        Q[i+1,t] <- sqrt(tauy)*exp(-0.5*tauy*pow((y[t] - (X2[t] + i + X4[t])),2))/sqrt(2*pi)
      }    
    } 
    
    # for (t in 2:T){
    for (t in 3:T){
        for (i in 0:(N_max1-1)){  # X1 (depends only on [imputed] X4)
        G1[i+1,t] <- exp(-exp(loglam1[t-2]) + i*loglam1[t-2] - logfact(i)) # 2nd order
        } 
      G1[(N_max1+1),t] <- max(0,1- sum(G1[(1:N_max1),t]))
      
      for (i in 0:(N_max3)){  # y & X3 (depends only on [imputed] X2)
        G3[i+1,t] <- ifelse((X2[t-1]-i)>0,
                            exp(i*log(phi3[t-1]) + (X2[t-1]-i)*log(1-phi3[t-1]) + logfact(X2[t-1]) - logfact(abs(i - X2[t-1])) - logfact(i)),
                            0)   
        
        Q[i+1,t] <- sqrt(tauy)*exp(-0.5*tauy*pow((y[t] - (X2[t] + i + X4[t])),2))/sqrt(2*pi)
      } 
 
      for (i in 0:(N_max1)){  # X2 (depends on [integrated] X1)
        P2[i+1,t] <- ifelse((i - X2[t])>0,
                            exp(X2[t]*log(phi2[t-1]) + (i - X2[t])*log(1-phi2[t-1]) + logfact(i) - logfact(abs(i - X2[t])) - logfact(X2[t])),
                            0)
      }

      for (i in 0:(N_max3)){  # X4 (depends on [integrated] X3)
        P4[i+1,t] <- ifelse((i + X4[t-1] - X4[t])>0,
                            exp(X4[t]*log(phi4[t-1]) + (i + X4[t-1] - X4[t])*log(1-phi4[t-1]) + logfact(i + X4[t-1]) - logfact(abs(i + X4[t-1] - X4[t])) - logfact(X4[t])),
                            0)
      }
    }
    
# #  Define the initial population priors for (a in 1:4){####
#     X1[1] ~ dpois(2500)
#     X2[1] ~ dpois(1500) 
#     X3[1] ~ dpois(1000)
#     X4[1] ~ dpois(4000) 
# 
# #  Define the system process for the count data 
#     for(t in 2:T){
# #    mean1[t] <- rho[t-1]*phi1[t-1]*X4[t-1] 
# #    mean2[t] <- phi2[t-1]*X1[t-1]
# #    mean3[t] <- phi3[t-1]*X2[t-1]
# #    mean4[t] <- phi4[t-1]*(X3[t-1]+X4[t-1])
#     
# #    tau1[t] <- 1/(rho[t-1]*phi1[t-1]*X4[t-1])
# #    tau2[t] <- 1/(X1[t-1]*phi2[t-1]*(1-phi2[t-1]))
# #    tau3[t] <- 1/(X2[t-1]*phi3[t-1]*(1-phi3[t-1]))
# #    tau4[t] <- 1/((X3[t-1]+X4[t-1])*phi4[t-1]*(1-phi4[t-1]))
#     
#       X1[t] ~ dpois(rho[t-1]*phi1[t-1]*X4[t-1])
#       X2[t] ~ dbin(phi2[t-1],X1[t-1])
#       X3[t] ~ dbin(phi3[t-1],X2[t-1])
#       X4[t] ~ dbin(phi4[t-1],(X3[t-1]+X4[t-1]))
#     }
#     
# # Define the observation process for the count data 
#     for(t in 2:T){
#       y[t] ~ dnorm((X2[t]+X3[t]+X4[t]),tauy) 
# #      y[t] ~ dpois(X2[t]+X3[t]+X4[t]) 
#     }

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
