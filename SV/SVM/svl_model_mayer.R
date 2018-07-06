  model{
  theta0 ~ dnorm(mu,itau2)
  thetamean[1] <- mu + phi*(theta0-mu)
  theta[1] ~ dnorm(thetamean[1],itau2)
  for (t in 2:n) {
    thetamean[t] <- mu + phi*(theta[t-1]-mu)
    theta[t] ~ dnorm(thetamean[t],itau2)
  }
  for (t in 1:(n-1)) {
    Ymean[t] <- rho/tau*exp(0.5*theta[t])*(theta[t+1] - mu - phi*(theta[t]-mu))
    Yisigma2[t] <- 1/(exp(theta[t])*(1-rho*rho))
    Y[t] ~ dnorm(Ymean[t],Yisigma2[t])
  }
  Ymean[n] <- mu-phi*(theta[n]-mu)
  Yisigma2[n] <-  1/(exp(theta[n]))
  Y[n] ~ dnorm(Ymean[n],Yisigma2[n])
  
  mu ~ dnorm(0,0.1)
  phistar ~ dbeta(20,1.5)
  itau2 ~ dgamma(2.5,0.025)
  beta <- exp(mu/2)
  phi  <- 2*phistar-1
  tau  <- sqrt(1/itau2)
  rho ~ dnorm(0, 1)
}