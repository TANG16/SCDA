# UC_EM.m
##############<load a file>##############
file = as.matrix(read.csv("USCPI.csv"));#
file = matrix(file[,1]);                #
file = matrix(file[-258]);              #
d = 5.673853997; #The first line        #
USCPI = matrix(rbind(d,file),,1);       #
#########################################
y = USCPI;
T = length(y);
omega2_0 = 9; # initial condition
omega = .5^2; # fix omega
####<Making a sparse matrix>####
S=diag(T);                     #
d=matrix(rep(0,T),1);          #
sparse=rbind(d,S);             #
sparse<-sparse[-(T+1),]        #
################################
H = diag(T) - sparse;
invOmega = diag(T)*c(1/omega2_0,1/omega*rep(1,T-1));

HinvOmegaH = t(H) %*% invOmega %*% H;
sigma2t = var(y)[1]; # initial guess
err = 1;
while (err> 10^(-4)){
  # E-step
  Kt = HinvOmegaH + diag(T)/sigma2t;
  taut = solve(Kt,(y/sigma2t));
  # M-step
  lam = matrix(eigen(Kt)$values);
  newsigma2t = (sum(1/lam) + t(y-taut) %*% (y-taut))/T;
  # update
  err = abs(sigma2t-newsigma2t);
  sigma2t = newsigma2t[1];
}  
Kt = HinvOmegaH + diag(T)/sigma2t;
taut = solve(Kt,(y/sigma2t));
win.graph() #------------>For Windows users
#quartz()#---------------->For Mac users
plot(taut,type='l',col='blue',xlab='',ylab='');
lines(y,col='blue');