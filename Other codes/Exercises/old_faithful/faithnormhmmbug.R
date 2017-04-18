model {
  
  P[1]<-pi[1]/(pi[1]+pi[2])
  P[2]<-pi[2]/(pi[1]+pi[2])
  c[1]~dcat(P[])     # This is for the initial state.
  
  for (i in 2:n)
  {
    c[i]~dcat(q[,c[i-1]])        
  }
  
  for (i in 1:n)
  {
    y[i]~dnorm(mu[c[i]],tau[c[i]])
  }
  
  for (j in 1:2)
  {
    tau[j]~dgamma(4,0.04)
  }
  
  mumean[1]<-mu0-0.2
  mumean[2]<-mu0+0.2
  for (j in 1:2)
  {
    mudash[j]~dnorm(mumean[j],3.3) 
  }
  mu[1:2]<-sort(mudash)  # This imposes the order constraint.
  
  mu0~dnorm(4.0,p.mu)
  p.mu<-1/0.3
  
  pi[1]~dbeta(1,1)
  pi[2]~dbeta(1,1)
  q[1,2]<-pi[1]
  q[2,2]<-1-q[1,2]
  q[2,1]<-pi[2]
  q[1,1]<-1-q[2,1]
  
}
