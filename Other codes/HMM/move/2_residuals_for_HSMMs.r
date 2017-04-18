


## this code computes the pseudo-residuals for the individual-specific HSMMs
## observations ("obs") need to be given in an n x 2 matrix, with first column giving the step lengths and second column giving the turning angles

library(CircStats)


## function that computes (the logarithm of) the forward probability (which is required in the computation of the residuals below) 
moveHMM.lalpha <- function(x,m,a,b,co,kappa,gamma)
{
delta <- solve(t(diag(sum(m))-gamma+1),rep(1,sum(m)))
n <- length(x[,1])
lalpha <- matrix(NA,sum(m),n)
allprobs <- matrix(rep(NA,sum(m)*n),nrow=n)
for (h in 1:n){
if (is.na(x[h,1])) {
allprobs[h,] <- rep(1,m[1]+m[2])
}
if (!is.na(x[h,1])) {
for (j in 1:2){
angle.prob <- ifelse(is.na(x[h,2]),1,dwrpcauchy(x[h,2],mu=co[j],rho=kappa[j]))
allprobs[h,((j-1)*m[1]+1):((j-1)*m[1]+m[j])] <- rep(dweibull(x[h,1],shape=a[j],scale=b[j])*angle.prob,m[j])
}}} 
foo <- delta
lscale <- 0
for (i in 1:n)
{
foo <- foo%*%gamma*allprobs[i,]
sumfoo <-sum(foo)
lscale <- lscale+log(sumfoo)
foo <- foo/sumfoo
lalpha[,i] <- log(foo)+lscale
}
list(la=lalpha)
}

## function that computes k-th pseudo-residual
moveHMM.psres <- function(x,k,m,a,b,co,kappa,gamma)
{
la <- moveHMM.lalpha(x,m,a,b,co,kappa,gamma)$la
la <- la[,dim(la)[2]-1]
c <- max(la)   ### to avoid underflow
ela <- exp(la-c) ### (scalar cancels in Pro1 and Pro2 computation)
pweibullvec <- rep(NA,sum(m))
pweibullvec[1:m[1]] <- rep(pweibull(obs[k,1],shape=a[1],scale=b[1]),m[1])
pweibullvec[(m[1]+1):sum(m)] <- rep(pweibull(obs[k,1],shape=a[2],scale=b[2]),m[2])
pcauchyvec <- rep(NA,sum(m))
if (!is.na(obs[k,2])){
pcauchyvec[1:m[1]] <- rep(integrate(dwrpcauchy,lower=0,upper=obs[k,2],mu=co[1],rho=kappa[1])$value,m[1])
pcauchyvec[(m[1]+1):sum(m)] <- rep(integrate(dwrpcauchy,lower=0,upper=obs[k,2],mu=co[2],rho=kappa[2])$value,m[2])
}
Pro1 <- t(ela)%*%(gamma/sum(ela))%*%pweibullvec
if (!is.na(obs[k,2])) Pro2<-t(ela)%*%(gamma/sum(ela))%*%pcauchyvec
npsr1 <- qnorm(Pro1)                            
ifelse(!is.na(obs[k,2]),npsr2<-qnorm(Pro2),npsr2<-NA)
c(npsr1,npsr2)
}

## use fitted HSMM  
m <- c(30,30)
moveHSMM$a -> a
moveHSMM$b -> b
moveHSMM$kappa -> kappa
moveHSMM$gam -> gam0
moveHSMM$co -> co
gamma <- gen.Gamma(m,gam0) ## function gen.Gamma from "individual_specific_HMMs_and_HSMMs.r"

n.obs <- dim(obs)[1]
for (k in 2:n.obs){if (!is.na(obs[k,2])) {if (obs[k,2]<0) obs[k,2]<-obs[k,2]+2*pi }} # transform angles from [-pi,pi] to [0,2pi]
res <- matrix(rep(NA,2*(n.obs-1)),ncol=2)

## compute pseudo-residuals
for (k in 2:n.obs){
if (!is.na(obs[k,1])) res[k-1,] <- moveHMM.psres(obs[max(1,k-250):k,],k,m,a,b,co,kappa,gamma)
}

## Q-Q plots of pseudo-residuals
par(mfrow=c(2,1))
qqnorm(res[,1],main="Q-Q plot Weibull-res.",xlab="",ylab="",xlim=c(-3,3),ylim=c(-3,3))
abline(a=0,b=1)
qqnorm(res[,2],main="Q-Q plot cauchy-res.",xlab="",ylab="",xlim=c(-3,3),ylim=c(-3,3))
abline(a=0,b=1)
 





 

