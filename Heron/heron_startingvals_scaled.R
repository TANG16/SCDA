X1init=rep(1400,T)/sc
X1init[2:T] <- rpois((T-1),c(4000/sc,y[2:(T-1)])*0.5)
X2init=rep(1200,T)/sc
X2init[2:T] <- rbinom((T-1),X1init[1:(T-1)],0.9)
X3init=rep(1000,T)/sc
X3init[2:T] <- rbinom((T-1),X2init[1:(T-1)],0.9)
X4init=rep(900,T)/sc
X4init[2:T] <- rbinom((T-1),X3init[1:(T-1)],0.9)

# tauy = 1000
tauy = 1

inits <- function()(list(alpha=rep(0,4),beta=rep(0,4),alphal=0,betal=0,alpharho=0,tauy=tauy,X1=X1init,X2=X2init,X3=X3init,X4=X4init))

params <- c('alpha','alpharho','alphal', 'beta', 'betal', 'X1', 'X2', 'X3', 'X4','tauy')