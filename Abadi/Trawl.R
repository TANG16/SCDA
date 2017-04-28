#Algorithm 6 (Simulation from the bivariate logarithmic series distribution).
library(VGAM)# ??? Load the VGAM package in R.
N = 10
p1 = 0.4
p2 = 0.8

Sim_BLSD <- function (N, p1, p2){
p1 <- p1/(1-p2) #??? Calculate the parameters of the modified LSD.
delta1 <- log(1-p2)/log(1-p1-p2)
L <- rlog(N, p1) #??? Simulate N i.i.d. Log(p1) rvs.
B <- rbinom(N, 1, 1- delta1) #??? Simulate N i.i.d. Bernoulli(1 - ??1) rvs.
C1 <- L * B #??? Generate N i.i.d. ModLog(~p1, ??1) rvs.
C2 <- numeric(N)
for (i in 1 : N) { 
  c1 <- C1[i]
  if (c1 == 0) {
    C2[i] <- rlog(1, p2) #??? Simulate a Log(p2) rv.
  }
  if (c1 > 0){
    C2[i] <- rnbinom(1, size = c1, prob = 1 - p2) #??? Simulate a NB(c1, p2) rv.
  }
}
C <- cbind(C1,C2) #??? Combine the component vectors to an N × 2 matrix.
return
}