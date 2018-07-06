# TVPAR.r
source("SURform.r");
nloop = 11000;
burnin = 1000;
##############<load a file>##############
file = as.matrix(read.csv("USCPI.csv"));#
file = matrix(file[,1]);                #
file = matrix(file[-258]);              #
d = 5.673853997; #The first line        #
USCPI = matrix(rbind(d,file),,1);       #
#########################################
p = 2;   # no of lags
y0 = matrix(USCPI[1:p]);
y = matrix(USCPI[p+1:(nrow(USCPI)-p)]);
T = length(y);
q = p+1; # dim of states
Tq = T*q;
          ## prior
asigma2 = 5; lsigma2 = 1*(asigma2-1);
aomega2 = 5; lomega2 = (aomega2-1)*matrix(c(0.5^2, 0.1^2*rep(1,p)));
invOmega0 = matrix(rep(1,q))/5;
          ## initialize
omega2 = .1*matrix(rep(1,q));
sigma2 = 1;
store_omega2 = matrix(rep(0,(nloop-burnin)*q),,q);
store_sigma2 = matrix(rep(0,(nloop-burnin)));
store_beta = matrix(rep(0,Tq));
          ## construct/compute a few things
X = cbind(matrix(rep(1,T)),rbind(y0[nrow(y0)],matrix(y[1:(nrow(y)-1)])),rbind(y0,matrix(y[1:(nrow(y)-2)])));
#############################################################
bigX = SURform(X);
####<Making a sparse matrix>####
Sp=diag(Tq);                    #
d=matrix(rep(0,Tq*q),q);        #
sparse=rbind(d,Sp);            #
sparse<-sparse[(1:Tq),]        #
################################
H = diag(Tq) - sparse;
newaomega2 = aomega2 + (T - 1)/2;
newasigma2 = asigma2 + T/2 ;
for (loop in 1:nloop){
  ## sample beta
  invS = diag(Tq)*c(t(invOmega0),matrix(rep(1/t(omega2),T-1),1));
  ######K is da problem
  K = t(H) %*% invS %*% H + t(bigX) %*% bigX/sigma2;
  C = t(chol(K));
  betahat = solve(K,(t(bigX) %*% y/sigma2));
  beta = betahat + solve(t(C),matrix(rnorm(Tq)));
  ## sample omega2
  erromega2 = matrix(H %*% beta,q,T);
  newlomega2 = lomega2 + rowSums(erromega2[,2:ncol(erromega2)]^2)/2;
  omega2 = 1/rgamma(3,shape=newaomega2,scale=1/newlomega2);
  ## sample sigma2
  newlsigma2 = lsigma2 + sum((y-bigX %*% beta)^2)/2;
  sigma2 = 1/rgamma(1,shape=newasigma2,scale=1/newlsigma2);
  if (loop>burnin){
    i=loop-burnin;
    store_beta = store_beta + beta;
    store_omega2[i,] = t(omega2);
    store_sigma2[i,] = sigma2;
  }
}
betahat = store_beta/(nloop-burnin);
sigma2hat = mean(store_sigma2);
omega2hat = mean(store_omega2);

#win.graph() #------------>For Windows users
quartz()#---------------->For Mac users

plot(betahat[seq(1,3*T,3)],type='l',col='blue'); #plot beta_{0t}
