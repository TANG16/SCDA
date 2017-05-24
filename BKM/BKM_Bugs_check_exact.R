T=36

# covariates
time <- seq(1,T,1)                
stdT <- (time-mean(time))/sd(time)
f=c(0.1922, 0.3082, 0.3082, -0.9676, 0.5401, 0.3082, 1.1995, 0.1921, -0.8526,
    -1.0835, -0.6196, -1.1995, -0.5037, -0.1557, 0.0762, 2.628, -0.3877, -0.968,
    1.9318, -0.6196, -0.3877, 1.700, 2.2797, 0.6561, -0.8516, -1.0835, -1.0835,
    0.1922, 0.1922, -0.1557, -0.5037, -0.8516, 0.8880, -0.0398, -1.1995, 0)


# alpha1 = 1
# alphaa = 2
# alphar = -2
# alphal = -4
# beta1 =-2
# betaa = 0.1
# betar = -0.7
# betal = -0.3

alpha1 = 1
alphaa = 2
alphar = -1.106
alphal = -4
beta1 = -0.19
betaa = 0.1
betar = -0.299 
betal = -0.3

phi1 = double(T)
phia = double(T)
lambda = double(T)
rho = double(T)

invlogit = function(x){exp(x)/(1+exp(x))}

for(t in 1:(T)){
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
