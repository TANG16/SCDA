# alpha=rep(0,4)
# beta=rep(0,4)
# alphal=0
# betal=0
# alpharho=0


alpha = c(-0.55, -0.1, 0.5, 1.1)
beta = c(0.02, -0.15, -0.15, -0.2)
alphal = 1.65
betal = -0.2
alpharho = 1.1
tauy = 1000


X2init=rep(1200,T)/sc
X4init=rep(2000,T)/sc

# X1[1] ~ dpois(2500)
# X2[1] ~ dpois(1500) 
# X3[1] ~ dpois(1000)
# X4[1] ~ dpois(4000) 


X2init= c(1500,y[-1]*2/7);
X4init= c(4000, y[-1]*4/7);


X2_cont<-X2init
X4_cont<-X4init
X2 <- round(X2_cont)
X4 <- round(X4_cont)

X3_try <- rbinom(T,X2,phi3)
plot(X3_try ,type='l')

phi1 = double(T)
phi2 = double(T)
phi3 = double(T)
phi4 = double(T)
phi1 = double(T)
lambda = double(T)
rho = double(T)

invlogit = function(x){exp(x)/(1+exp(x))}
logfact <- function(x){sum(log(seq_len(x)))}

for(t in 1:(T)){
  ind =  alpha[1] + beta[1]*f[t]
  phi1[t] <-exp(ind)/(1+exp(ind)) # corresponds to the year 1963
  ind =  alpha[2] + beta[2]*f[t]
  phi2[t] <-exp(ind)/(1+exp(ind)) # corresponds to the year 1963
  ind =  alpha[3] + beta[3]*f[t]
  phi3[t] <-exp(ind)/(1+exp(ind)) # corresponds to the year 1963
  ind =  alpha[4] + beta[4]*f[t]
  phi4[t] <-exp(ind)/(1+exp(ind)) # corresponds to the year 1963
  rho[t] <- exp(alpharho) # We assume here that t=1
}

X1 <-X1init
X2 <-X2init
X3 <-X3init
X4 <-X4init


for(t in 2:T){
  X1[t] <- rpois(1,rho[t-1]*phi1[t-1]*X4[t-1])
  X2[t] <- rbinom(1,X1[t-1],phi2[t-1])
  X3[t] <- rbinom(1,X2[t-1],phi3[t-1])
  X4[t] <- rbinom(1,(X3[t-1]+X4[t-1]),phi4[t-1])
}


#####
loglam1 <- rep(NaN,T)
zeros <- rep(NaN,T)
PHI <- rep(NaN,T)
loglik <- rep(NaN,T)


for (t in 1:(T)){
  loglam1[t] <- log(X4[t]) + log(rho[t]) + log(phi1[t])
}


#####
G1 <- matrix(0,(N_bin1+1),T)
P2 <- matrix(0,(N_bin1+1),T)
# for (t in 2:T){
for (t in 3:T){
  for (i in 0:(N_bin1)){  # X1 (depends only on [imputed] X4)
    G1[i+1,t] <- exp(-exp(loglam1[t-2]) + bin1[i+1] *loglam1[t-2] - logfact(bin1[i+1] )) # 2nd order
  } 
}

for (t in 3:T){
  for (i in 0:(N_bin1)){  # X2 (depends on [integrated] X1)
    P2[i+1,t] <- ifelse((bin1[i+1] - X2[t])>0,
                        exp(X2[t]*log(phi2[t-1]) + (bin1[i+1] - X2[t])*log(1-phi2[t-1]) + logfact(bin1[i+1]) - logfact(abs(bin1[i+1]  - X2[t])) - logfact(X2[t])),
                        0)
  }
}


plot(P2[,3],type='l')
for (t in c(4:T)){lines(P2[,t])}


plot(G1[,3],type='l')
for (t in c(4:T)){lines(G1[,t])}


for (t in 3:T){
  print(log(sum(G1[,t] * P2[,t] )))
}
#####

N_bin3 =  49#  # reduction to be based on Gamma plots 
bin_size3 = 19#  odd so the midpoints are integer

G3 <- matrix(0,(N_bin3+1),T)
P4 <- matrix(0,(N_bin3+1),T)
Q <- matrix(0,(N_bin3+1),T)

# Bins' midpoints
bin3 = rep(0,N_bin3+1)
for (i in 0:(N_bin3)){
  bin3[i+1] <- 0.5*(bin_size3*(2*i+1)-1) # ith bin's midpoint
}

for (t in 3:T){
  for (i in 0:(N_bin3)){  # y & X3 (depends only on [imputed] X2)
    G3[i+1,t] <- ifelse((X2[t-1]-bin3[i+1] )>0,
                        exp(bin3[i+1]*log(phi3[t-2]) + (X2[t-2]-bin3[i+1] )*log(1-phi3[t-2]) + logfact(X2[t-2]) - logfact(abs(bin3[i+1] - X2[t-2])) - logfact(bin3[i+1])),
                        0)   
    
    Q[i+1,t] <- sqrt(tauy)*exp(-0.5*tauy*((y[t] - (X2[t] + bin3[i+1]  + X4[t]))^2))/sqrt(2*pi)
  } 
}

for (t in 3:T){
  for (i in 0:(N_bin3)){  # X4 (depends on [integrated] X3)
    P4[i+1,t] <- ifelse((bin3[i+1]  + X4[t-1] - X4[t])>0,
                        exp(X4[t]*log(phi4[t-1]) + (bin3[i+1]  + X4[t-1] - X4[t])*log(1-phi4[t-1]) + logfact(bin3[i+1]  + X4[t-1]) - logfact(abs(bin3[i+1]  + X4[t-1] - X4[t])) - logfact(X4[t])),
                        0)
  }
}



plot(P4[,3],type='l')
for (t in c(4:T)){lines(P4[,t])}


plot(G3[,3],type='l')
for (t in c(4:T)){lines(G3[,t])}

plot(Q[,3],type='l')
for (t in c(4:T)){lines(Q[,t])}

for (t in 3:T){
print(log(sum(G3[,t] * P4[,t] * Q[,t])))
}
#####
C <- 1000000

for (t in 1:T){
  loglik[t] <- log(exp(log(sum(G1[,t] * P2[,t])) + log(sum(G3[,t] * P4[,t] * Q[,t]))))
  PHI[t] <- -loglik[t]# + C
  zeros[t] <- dpois(1,PHI[t])
}

