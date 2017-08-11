model{ 
  # Define the observation process for the census/index data
  for(t in 3:T){
    y[t] ~ dnorm(Na[t],tauy)
  }
  
	for(t in 3:T){  
	  N1[t] ~ dpois(po[t])
  	Na[t] ~ dbin(bin2[t],bin1[t])
  	
  	bin1[t] <- N1[t-1]+Na[t-1]
  	bin2[t] <- phia[t-1]
  	po[t] <- Na[t-1]*rho[t-1]*phi1[t-1]
	}
  
  for(t in 1:2){
    N1[t] ~ dpois(20)
    Na[t] ~ dbin(0.5,200)
  }
  
  # Define the logistic regression equations
  for(t in 1:(T-1)){
    logit(phi1[t]) <- alpha1 + beta1*f[t] # corresponds to the year 1963
    logit(phia[t]) <- alphaa + betaa*f[t]
    log(rho[t]) <- alphar + betar*stdT[t] # We assume here that t=1
  }
  
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

}