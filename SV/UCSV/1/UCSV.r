# UCSV.r
source("SVRW.r");
source("gamrand.r");
##############<load a file>##############
file = as.matrix(read.csv("USCPI.csv"));#
file = matrix(file[,1]);                #
file = matrix(file[-258]);              #
d = 5.673853997; #The first line        #
USCPI = matrix(rbind(d,file),,1);       #
#########################################
y = USCPI;
T = length(y);
nloop = 11000;
burnin = 1000;
## prior
Vtau = 9; Vh = 9;
atau = 10; ltau = .25^2*(atau-1);
ah = 10; lh = .2^2*(ah-1);
## initialize the Markov chain
omega2tau = .25^2;
omega2h = .2^2;
h = as.numeric(log(var(y)*.8)) * matrix(rep(1,T));
####<Making a sparse matrix>####
Sp=diag(T);                    #
d=matrix(rep(0,T),1);          #
sparse=rbind(d,Sp);            #
sparse<-sparse[-(T+1),]        #
################################
H = diag(T) - sparse;
## initialize for storage
store_omega2tau = matrix(rep(0,(nloop-burnin)),,1);
store_omega2h = matrix(rep(0,(nloop-burnin)),,1);
store_tau = matrix(rep(0,(nloop-burnin)*T),,T);
store_h = matrix(rep(0,(nloop-burnin)*T),,T);
## compute a few things
newatau = (T-1)/2 + atau;
newah = (T-1)/2 + ah;
for (loop in 1:nloop){
  ## sample tau
  invOmegatau = diag(T)*c(1/Vtau, 1/omega2tau*rep(1,T-1));
  invSigy = diag(T)*c(exp(-h));
  Ktau = t(H) %*% invOmegatau %*% H + invSigy;
  Ctau = t(chol(Ktau));
  tauhat = solve(Ktau,(invSigy %*% y));
  tau = tauhat + solve(t(Ctau), matrix(rnorm(T),,1));
  ## sample h
  ystar = log((y-tau)^2 + .0001 );  
  result = SVRW(ystar,h,omega2h,Vh);
  h = result[[1]];
  ## sample omega2tau
  newltau = ltau + sum((tau[2:nrow(tau)]-tau[1:(nrow(tau)-1)])^2)/2;
  omega2tau = 1/gamrand(newatau, newltau);
  ## sample omega2h
  newlh = lh + sum((h[2:nrow(h)]-h[1:(nrow(h)-1)])^2)/2;
  omega2h = 1/gamrand(newah, newlh);
  if (loop>burnin){
    i = loop-burnin;
    store_tau[i,] = t(tau);
    store_h[i,] = t(h);
    store_omega2tau[i,] = omega2tau;
    store_omega2h[i,] = omega2h; 
  }
}
tauhat = matrix(rowMeans(t(store_tau)));
hhat = matrix(colMeans(store_h));

win.graph() #------------>For Windows users
#quartz()#---------------->For Mac users

plot(tauhat,type='l',col='blue')

win.graph() #------------>For Windows users
#quartz()#---------------->For Mac users
plot(hhat,type='l',col='blue')