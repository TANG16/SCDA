function f=likj(thet,Xrr1,Xc1,Xrra,Xca,Xlam,Xp,atype)

%LIKJ Joint ring-recovery and census log-likelihood
%   LIKJ(THET,XRR1,XC1,XRRA,XCA,XLAM,XP,ATYPE) performs a combined analysis
%   of ring-recovery and census data by calculating their
%   joint likelihood.
%
%    THET is a (column) vector of parameters, in form
%       [beta1 betaa betalam betap sigma]'
%    XRR1 are the rr covariates, in columns, on which PHI1 may depend,
%       in form [year|covariate1|covariate2|...]
%    XC1 are the census covariates, in columns, on which PHI1 may depend,
%       in form [year|covariate1|covariate2|...]
%    XRRA are the rr covariates, in columns, on which PHIA may depend,
%       in form [year|covariate1|covariate2|...]
%    XCA are the census covariates, in columns, on which PHIA may depend,
%       in form [year|covariate1|covariate2|...]
%    XLAM are the covariates, in columns, on which LAMBDA may depend,
%       in form [year|covariate1|covariate2|...]
%    XP are the covariates, in columns, on which P may depend,
%       in form [year|covariate1|covariate2|...]
%    Use empty values ([]) for non-existent/applicable fields


global Ringed Mobs unrec nrows ncols Census ncensus

%initialisations
[rXrr1 cXrr1]=size(Xrr1); cXrr1=max(1,cXrr1);
[rXc1 cXc1]=size(Xc1); cXc1=max(1,cXc1);
[rXrra cXrra]=size(Xrra); cXrra=max(1,cXrra);
[rXca cXca]=size(Xca); cXca=max(1,cXca);
[rXlam cXlam]=size(Xlam); cXlam=max(1,cXlam);
[rXp cXp]=size(Xp); cXp=max(1,cXp);

if (cXrra == ncols-1 & cXca == ncensus),
  cXrra = cXca;
end

if ~(cXrr1 == cXc1 & cXrra == cXca), 
   error('number of covariates for ring-recovery and census survival must be equal') 
end


%parameters
beta1=thet(1:cXrr1);                          %note cXrr1=cXc1 etc
betaa=thet(cXrr1+1:cXrr1+cXrra);
betalam=thet(cXrr1+cXrra+1:cXrr1+cXrra+cXlam);
betap=thet(cXrr1+cXrra+cXlam+1:cXrr1+cXrra+cXlam+cXp);
H=thet(cXrr1+cXrra+cXlam+cXp+1);

%joint log-likelihood
if (cXrra == ncensus), %ols-1 & cXca == ncensus),
  f=likrr([beta1;betaa(2:cXrra);betalam],Xrr1,Xrra,Xlam,atype)+...
    likc([beta1;betaa;betap;H],Xc1,Xca,Xp,atype);
else
  f=likrr([beta1;betaa;betalam],Xrr1,Xrra,Xlam,atype)+...
    likc([beta1;betaa;betap;H],Xc1,Xca,Xp,atype);
end
