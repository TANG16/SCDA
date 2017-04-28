phi1 = double(T-1)
phia = double(T-1)
lambda = double(T-1)
rho = double(T-1)

invlogit = function(x){exp(x)/(1+exp(x))}

for(t in 1:(T-1)){
  ind =  alpha1 + beta1*f[t]
  phi1[t] <-exp(ind)/(1+exp(ind)) # corresponds to the year 1963
  
  ind =  alphaa + betaa*f[t]
  phia[t] <- invlogit(ind)
  
  # log(rho[t]) <- alphar + betar*t # We assume here that t=1
  rho[t] <- exp(alphar + betar*stdT[t]) # We assume here that t=1
  
  # logit(lambda[t]) <- alphal + betal*(t+1)
  ind <- alphal + betal*stdT[t]
  lambda[t] <- invlogit(ind)
}

dN1 = double(T)
dNa = double(T)

for(t in 1:2){
  dN1[t] = dpois(N1[t],200)
  dNa[t] = dbinom(Na[t],2000,0.5)
}

for(t in 3:T){
  bin1 <- N1[t-1]+Na[t-1]
  bin2 <- phia[t-1]
  
  po <- Na[t-1]*rho[t-1]*phi1[t-1]
  
  dN1[t] = dpois(N1[t],po)
  dNa[t] = dbinom(Na[t],bin1,bin2)
}



# neg bin prior for N1 Na t=1,2
# mean = pr/(1-p)
# var = mean/(1-p)

p1=1-1/5000
r1 = 200*(1-p1)/p1


pa=1-1/1000
ra = 1000*(1-pa)/pa

rnbinom(1, r1, 1-p1)
