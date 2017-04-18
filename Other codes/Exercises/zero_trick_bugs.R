# Use the "zeros trick" to model a normal distribution with unknown mean mu and 
# unknown standard deviation sigma, including predicting a new observation.

model{
  
  C <- 10000
  # Likelihood
  for (i in 1:8)
  {
    # z[i] <- 0
    z[i] ~ dpois(phi[i])
    # phi <- -log(L[i] + C)
    # L[i] <- myfun(y[i],theta) 
    #  L_i is set to a function of y_i and theta proportional to the likelihood p(y_i|theta)
    phi[i] <- log(sigma) + 0.5*pow((y[i]-mu)/sigma,2)
    # phi[i] <- myfun(mu,sigma,y[i])
  }
  # Priors
  y[8] ~ dunif(-100000, 100000)
  mu ~ dunif(-100,100)
  sigma ~ dunif(0,100)
  
}