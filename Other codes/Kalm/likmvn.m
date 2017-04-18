function f=likmvn(x,mu,Sigma)

%LIKMVN  Multivariate normal log-likelihood
%   LIKMVN(X,MU,SIGMA) computes the minus log-likelihood
%   function of a multivariate normal with mean MU and 
%   variance SIGMA at observation X.


%initialisations
x=x(:); mu=mu(:);
p=length(x);

%log-likelihood
u=chol(inv(Sigma))*(x-mu);
logl=-p/2*log(2*pi)-1/2*log(det(Sigma))-1/2*u'*u;
f=-logl;
