m <- 30 #50 100 200
gmax <- 4
phi <- 0.98
sigma <- 0.2
beta <- 0.05 
P1 <- (sigma^2)/(1-phi^2)

T <- 1000
h <- rep(NaN,1000)

h[1] <- sqrt(P1)*rnorm(1)
for (t in c(2:T)){
  h[t] = phi*h[t-1] + sigma*rnorm(1)
}
y <- beta*exp(h/2)*rnorm(T)
plot(y,type='l')

loglikelihood.SV0<-function(y,phi,sigma,beta,m,gmax){	
  K <- m+1	
  b <- seq(- gmax,gmax,length=K) #endpoints
  bs <- (b[-1]+b[-K])*0.5 # midpoints
  sey <- beta*exp(bs/2) # st.dev. of y in each midpoint
  Gamma <- matrix(0,m,m)	
  for (i in 1:m) Gamma[i,] <- diff(pnorm(b,phi*bs[i],sigma)) #time constant tpm 
  Gamma <- Gamma/apply(Gamma,1,sum)	 #scale the rows of Gamma to sum up to 1
  foo <- solve(t(diag(m)-Gamma+1),rep(1,m)) # compute delta=delta*Gamma, the stationary distribution of the MC 
  llk <- 0	
  for (t in 1:length(y)){	
    foo <- foo%*%Gamma*dnorm(y[t],0,sey)	
    sumfoo <- sum(foo)
    llk <- llk+log(sumfoo)	
    foo <- foo/sumfoo
    }	
  return(llk)
}	#16
