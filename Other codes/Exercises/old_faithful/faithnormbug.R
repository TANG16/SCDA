# BUGS model for ``Old Faithful'' log intervals (normal mixture): faithnormbug.
model {
  
  for (i in 1:n)
  {
    c[i]~dcat(q[])
    y[i]~dnorm(mu[c[i]],tau[c[i]])
  }
  
  for (j in 1:2)
  {
    tau[j]~dgamma(4,0.04)
  }
  
  mumean[1]<-mu0-0.2
  mumean[2]<-mu0+0.2
  mu[1]~dnorm(mumean[1],3.3) I(,mu[2]) # This imposes the order constraint.
  mu[2]~dnorm(mumean[2],3.3) I(mu[1],)
  
  mu0~dnorm(4.0,p.mu)
  p.mu<-1/0.3
  
  pi~dbeta(3,3)
  q[1]<-pi
  q[2]<-1-pi
}
