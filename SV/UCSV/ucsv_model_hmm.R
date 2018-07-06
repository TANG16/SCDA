model{
  
  # THE 0-1 TRICK
  for (t in 1:round(T/2)){
    zeros[t] ~ dpois(PHI[t])
    PHI[t] <- -loglik[t] + C
    loglik[t] <- ( log(exp(log(Q_even[t]) + log(sum(G_odd[,t] * Q_odd[,t] * G_even[,t])))) +
      log(exp(log(R_odd[t]) + log(sum(D_evend[,t] * R_even[,t] * D_odd[,t])))) )
    
  }
  
  C <- 1000000
  
  # define the observation process for inflation
  # for (t in 1:T){
  #   y[t] ~ dnorm(tau[t], 1.0/exp(h0 + omega_h*h[t]))
  # }
  
  # INTEGRATE H
  # EVEN H's IMPUTED, ODD INTEGRATED
  for (t in 1:round(T/2)){
    Q_even[t] <- dnorm(y[2*t],tau[2*t], 1.0/exp(h0 + omega_h*h[t])) #  observation density for even observations given even states
  }
  
  for (t in 2:round(T/2)){ # we have T observations and T/2 imputed states
    #########  integration of odd
    for (i in 1:N_bin_h){ 
      G_odd[i,t] <- dnorm(bin_mid_h[i], h[t-1] , 1) # transition probability to h_odd[t] from h_even[t-1]
      Q_odd[i,t] <- dnorm(y[2*t-1],tau[2*t-1], 1.0/exp(h0 + omega_h*bin_mid_h[i])) # observation density for odd observations given odd states
      G_even[i,t] <- dnorm(h[t],bin_mid_h[i], 1) # transition probability to h_even[t] from h_odd[t]    
    }
  }  
  

  # define the state process for the trend
  # for (t in 2:T){
  #   tau[t] ~ dnorm(tau[t-1], 1.0/exp(g0 + omega_g*g[t]))
  # }
  # ODD G's IMPUTED, EVEN INTEGRATED
  
  for (t in 2:round(T/2)){ # tau_3|tau_2
    R_odd[t] <- dnorm(tau[2*t-1],tau[2*t-2], 1.0/exp(g0 + omega_g*g[t])) #  observation density for even observations given even states
  }
  
  for (t in 2:round(T/2)){ # we have T observations and T/2 imputed states
    #########  integration of odd
    for (i in 1:N_bin_g){ 
      D_even[i,t] <- dnorm(bin_mid_g[i], g[t-1] , 1) # transition probability to g_even[t] from g_odd[t-1]
      R_even[i,t] <- dnorm(tau[2*t-2],tau[2*t-3], 1.0/exp(g0 + omega_g*bin_mid_g[i])) # observation density for odd observations given odd states
      D_odd[i,t] <- dnorm(g[t],bin_mid_g[i], 1) # transition probability to g_odd[t] from g_even[t]    
    }
  }  
  
  # define the state process for the volatilities
  # for (t in 2:T){
  #   h[t] ~ dnorm(h[t-1],1)
  #   g[t] ~ dnorm(g[t-1],1)
  # }
  
   
  
  
  # initialise
  # tau[1] ~ dnorm(tau0, 1/(V_tau*exp(g0 + omega_g*g[1])))
  # h[1] ~ dnorm(0,1/V_h)
  # g[1] ~ dnorm(0,1/V_g)
  
  for (i in 1:N_bin){  
    G_odd[i,1] <- dnorm(bin_mid_h[i], 0 , 1/V_h) # transition probability to h_odd[t] from h_even[t-1]
    Q_odd[i,1] <- dnorm(y[1],tau[1], 1.0/exp(h0 + omega_h*bin_mid_h[i])) # observation density for odd observations given odd states
    G_even[i,1] <- dnorm(h[1],bin_mid_h[i], 1) # transition probability to h_even[t] from h_odd[t]    
  }
  
  
  
  
  
  
  # define the priors for parameters
  omega_h ~ dnorm(0,1/(0.2))
  omega_g ~ dnorm(0,1/(0.2))
  h0 ~ dnorm(0,1/10)
  g0 ~ dnorm(0,1/10) 
}
