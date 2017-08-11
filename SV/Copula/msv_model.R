# from https://sourceforge.net/p/mcmc-jags/discussion/610037/thread/56b1c15f/
# A multivariate stochatstic volatility model with copulae for the dependence in the latent and model volatilites. 
# This is a simplified model, 2 series of length 1653, and 2 copulae.
model{
  
  C <- 10000    
  
  # Latent Likelihoods
  for (t in 1:n) {
    yisigma2[t] <- exp(-thy[t]);
    y[t] ~ dnorm(0, yisigma2[t]);
    xisigma2[t] <- exp(-thx[t]);
    x[t] ~ dnorm(0, xisigma2[t]);
  }
  
  # Initial Likelihoods
  thmeany[1] <- muy + phiy*(thetay0 - muy);
  thy[1] ~ dnorm(thmeany[1], itauy2);
  thmeanx[1] <- mux + phix*(thetax0 - mux);
  thx[1] ~ dnorm(thmeanx[1], itaux2);
  
  # Likelihoods
  for (t in 2:n) {
    thmeany[t] <- muy + phiy*(thy[t-1] - muy)
    thy[t] ~ dnorm(thmeany[t], itauy2)
    thmeanx[t] <- mux + phix*(thx[t-1] - mux)
    thx[t] ~ dnorm(thmeanx[t], itaux2)
  }
  
  # Y's Priors
  thetay0 ~ dnorm(muy, itauy2)
  phistary ~ dbeta(1, 1)
  phiy <- 2*phistary - 1
  itauy2 ~ dgamma(.001, 0.001);                                                      
  tauy <- sqrt(1/itauy2);  
  muy ~ dnorm(0, 0.001)
  betay <- exp(muy/2)
  
  # X's Priors
  thetax0 ~ dnorm(mux, itaux2)
  phistarx ~ dbeta(1, 1)
  phix <- 2*phistarx - 1
  itaux2 ~ dgamma(.001, 0.001);                                                      
  taux <- sqrt(1/itaux2);  
  mux ~ dnorm(0, 0.001)
  betax <- exp(mux/2)
  
  # Copula Parameter Priors
  rho ~ dunif(-1, 1)
  nu <- 1/inv_nu
  inv_nu ~ dunif(0, .5)
  rho2 ~ dunif(-1, 1)
  
  # Latent Volatility Copula
  for (t in 1:n){
    
    pseudo_x[t] <- x[t]*pow(xisigma2[t],0.5)
    pseudo_y[t] <- y[t]*pow(yisigma2[t],0.5)
    NZD_AUD_t_cop[t] <- log( (1/(2*3.14159*sqrt((1-pow(rho,2)))))
                             *(1/(dt(pseudo_x[t],0,1,nu)*dt(pseudo_y[t],0,1,nu)))
                             *pow(1 + (pow(pseudo_x[t],2) + pow(pseudo_y[t],2) + (2*rho*pseudo_x[t]*pseudo_y[t])/(nu*(1-pow(rho,2)))),-(nu + 2)/2) )
    loglik[t] <- NZD_AUD_t_cop[t]
    copulat[t] <- -loglik[t] + C
    zeros[t] ~ dpois(copulat[t])
  }
  
  # Model Volatility Copula
  for (t in 1:n){
    pseudo_x2[t] <- (thx[t]-thmeanx[t])/taux
    pseudo_y2[t] <- (thy[t]-thmeany[t])/tauy
    NZD_AUD_Gauss[t] <- log( (1/(sqrt(1 - pow(rho2,2))))*exp((pow(rho2,2)*pow(pseudo_x2[t],2)*pow(pseudo_y2[t],2) + 2*rho2*(pseudo_x2[t])*(pseudo_y2[t]))/(2*(1-pow(rho2,2)))) )
    loglik2[t] <- NZD_AUD_Gauss[t]
    copula2[t] <- -loglik2[t] + C
    zeros2[t] ~ dpois(copula2[t])
  }
}