function f=alikj(thet,X1,XA,Xp,atype)

%ALIKJ Approximate joint ring-recovery and census likelihood
%   ALIKJ(THET,X1,XA,XP,ATYPE) performs an approximate combined analysis
%   of ring-recovery and census data by calculating an approximation
%   to their joint likelihood.
%
%    THET is a (column) vector of parameters, in form
%       [beta1 betaa betalam betap sigma]'
%    X1 are the covariates, in columns, on which PHI1 may depend,
%       in form [ones|covariate1|covariate2|...]
%    XA are the covariates, in columns, on which PHIA may depend,
%       in form [ones|covariate1|covariate2|...]
%    XP are the covariates, in columns, on which P may depend,
%       in form [ones|covariate1|covariate2|...]
%    Use empty values ([]) for non-existent/applicable fields


global Mu Sigma Census ncensus

%initialisations
[rX1 cX1]=size(X1); cX1=max(1,cX1);
[rXa cXA]=size(XA); cXA=max(1,cXA);
[rXp cXp]=size(Xp); cXp=max(1,cXp);

%parameters
beta1=thet(1:cX1);
betaa=thet(cX1+1:cX1+cXA);
betap=thet(cX1+cXA+1:cX1+cXA+cXp);
H=thet(cX1+cXA+cXp+1);

%joint log-likelihood
f=likmvn([beta1;betaa],Mu,Sigma)+...
  likc([beta1;betaa;betap;H],X1,XA,Xp,atype);
