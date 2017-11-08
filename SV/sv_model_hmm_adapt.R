model{
#   for (t in 1:10){
#     for (i in 1:N_q){
#       BIN100[i,t] <- bin_mid[i,round(100*t)] 
#     }
#   }
  
  # THE 0-1 TRICK
  for (t in 1:round(T/2)){
    zeros[t] ~ dpois(PHI[t])
    PHI[t] <- -loglik[t] + C
#     loglik[t] <- log(exp(log(Q_even[t]) + log(sum(G_odd[,t] * Q_odd[,t] * G_even[,t]))))
    loglik[t] <- log(exp(log(Q_even[t]) + log(sum(Q_odd[,t]* G_even[,t]))))
#     loglik[t] <- log(exp(log(Q_even[t]) + log(sum(bin_size[,t]*Q_odd[,t]* G_even[,t]))))
    
  }

  C <- 1000000
  
#   # bin sizes
#   for (t in 2:round(T/2)){  
#     bin_size[1,t] = bin_max
#     bin_size[N_q,t] = bin_max 
#     for (i in 2:(N_q-1)){ 
#       bin_size[i,t] = 0.5*(bin_mid[i+1,t] - bin_mid[i,t]) + 0.5*(bin_mid[i,t] - bin_mid[i-1,t])
#     }
#   }
# 
#   bin_size[1,1] = 1/N_q
#   bin_size[N_q,1] = 1/N_q 
#   for (i in 2:(N_q-1)){ 
#     bin_size[i,1] = 1/N_q
#   }

  # EVEN IMPUTED, ODD INTEGRATED
  for (t in 1:round(T/2)){
    Q_even[t] <- dnorm(y[2*t],0, 1.0/exp(h[t])) #  observation density for even observations given even states
  }

  for (t in 2:round(T/2)){ # we have T observations and T/2 imputed states
    #########  integration of odd
    for (i in 1:N_q){ # bins are DEMEANED VOLATILITIES
      # G_odd is now simply a constant for all values as we are using the quantiles of the distribution       
      #       G_odd[i,t] <- dnorm(bin_mid[i], phi*(h[t-1] - mu), sigma2_star) # transition probability to h_odd[t] from h_even[t-1]
#       G_odd[i,t] <- pnorm(bin[i], phi*(h[t-1] - mu), sigma2_star) 
#       - pnorm(bin[i] + mu , mu + phi*(h[t-1] - mu), sigma2_star) # transition probability to h_odd[t] from h_even[t-1]
#       bin_mid[i,t] <- qnorm(mid[i],phi*(h[t-1] - mu), sigma2_star) 
      bin_mid[i,t] <- phi*(h[t-1] - mu) + sigma*mid[i] 
      
      Q_odd[i,t] <- dnorm(y[2*t-1],0, 1.0/exp(bin_mid[i,t]  + mu)) # observation density for odd observations given odd states
      G_even[i,t] <- dnorm(h[t], mu + phi*(bin_mid[i,t]), sigma2_star) # transition probability to h_even[t] from h_odd[t]    
    }
  }
  

  ######### INTEGRATION INITIALISATION
  # for the integration of h[1] we need conditioning on h[0] - hence outside the main loop over time
  # integration of odd
  for (i in 1:N_q){ # bins are DEMEANED VOLATILITIES
    #     G_odd[i,1] = dnorm(bin_mid[i] + mu , mu + phi*(h0- mu), sigma2_star) # transition probability to h_odd[t] from h_even[t-1]
#     G_odd[i,1] = pnorm(bin[i+1] + mu , mu + phi*(h0- mu), sigma2_star) 
#     - pnorm(bin[i] + mu , mu + phi*(h0- mu), sigma2_star) # transition probability to h_odd[t] from h_even[t-1]
#     bin_mid[i,1] <- qnorm(mid[i],phi*(h0 - mu), sigma2_star) 
    bin_mid[i,1] <- phi*(h0 - mu) + sigma*mid[i] 
    
    Q_odd[i,1] = dnorm(y[1],0, 1.0/exp(bin_mid[i,1]  + mu)) # observation density for odd observations given odd states
    G_even[i,1] = dnorm(h[1], mu + phi*(bin_mid[i,1]), sigma2_star) # transition probability to h_even[t] from h_odd[t]    
  }
  
  for (t in 1:round(T/2)){
    # Use a discrete Uniform prior so the only influence on the posterior distr is the Upper limit  
    h[t] ~ dunif(-Up_h,Up_h)
  }  
  h0 ~ dnorm(mean_h0,1/P1)
  P1 <- sigma2/(1-phi^2)
  mean_h0 <- mu
    
  ######### PARAMETER PRIORS
  mu ~ dnorm(0,0.1)
  
  # (phi+1)/2 ~ betapdf((phi+1)/2, 20, 1.5)
  phi <- 2*phi_star - 1
  phi_star ~ dbeta(20, 1.5)
  
  # 1/sigma2 ~ gampdf(1./s2, 5/2, 0.05/2);
  sigma <- sqrt(sigma2)
  sigma2 <- 1/sigma2_star # variance = 1/precision
  sigma2_star ~ dgamma(5/2, 0.05/2) # shape and rate
}
