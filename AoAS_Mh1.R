# The semi-complete data likelihood approach using the Jeffreys' prior specification for N.
model{
  Pi <- 3.14159265359
  # Priors:
  alpha ~ dnorm(0.0,0.01)
  tau ~ dgamma(0.01,0.01)
  sigma <- 1/sqrt(tau)
  
  for (i in 1:n) {
    y[i] ~ dbin(p[i],T)
    logitp[i] ~ dnorm(alpha,tau)
    logit(p[i]) <- logitp[i]
  }
  
  # Calculate probability of not being observed using Gauss-Hermite quadrature
  # q = number of quadrature points
  # weights and nodes correspond to q quadrature points; entered as data
  
  for(i in 1:q){
    probi[i] <- 1/sqrt(Pi)*weights[i]*(1/(1+exp(sqrt(2)*sigma*nodes[i]+alpha)))^T
  }
  prob<- sum(probi[])
  
  # Prior for N: Jeffreys' prior - this is incorporated in the zero trick below
  # in specifying the likelihood term
  # However a prior distribution is needed to be specified on N
  
  # Use a discrete Uniform prior so the only influence on the posterior
  # distribution is the upper limit
  
  n00 ~ dcat(prior[]) # prior = rep(1/(M+1-n),M+1-n); entered as data
  n0 <- n00 - 1
  N <- n + n0
  
  # Use zero trick for model likelihood
  # Note loggam(N) instead of loggam(N+1) because of Jeffreys' prior for N
  
  logzeroprob <- loggam(N) - loggam(n0+1) - loggam(n+1) + n0*log(prob)
  lambda <- -logzeroprob + 100000
  dummy ~ dpois(lambda) # dummy = 0; entered as data
}