function [h,f,g]=ksc(thet,Xphi1,Xphia,Xp,atype)

%KSC  Kalman smoother for census 
%   KSC(THET,XPHI1,XPHIA,XP,ATYPE) filters the data Census(t) generated from
%
%      Census(t) = [ Z ] alpha(t) + epsilon(t),   Var[epsilon(t)] = sigma^2
%
%   alpha(t + 1) = [ T ] alpha(t) + ita(t),     Var[ita(t)] = Q
%
%   using the Kalman filter; normal errors are assumed in both
%   equations. KSVVV also produces the log-likelihood and 
%   the one-step ahead predictors.
%
%    THET is a (column) vector of parameters, in form
%       [beta1 betaa betap sigma]'
%    XPHI1 are the covariates, in columns, on which PHI1 may depend,
%       in form [ones|covariate1|covariate2|...]
%    XPHIA are the covariates, in columns, on which PHIA may depend,
%       in form [ones|covariate1|covariate2|...]
%    XP are the covariates, in columns, on which P may depend,
%       in form [ones|covariate1|covariate2|...]
%    Use empty values ([]) for non-existent/applicable fields

%    NB Only the variance Q(t) is time-variant.


global Census ncensus

if isempty(Xphi1), Xphi1 = ones(ncensus,1); end
if isempty(Xphia), Xphia = ones(ncensus,1); end
if isempty(Xp), Xp = ones(ncensus,1); end

[rX1 cX1] = size(Xphi1);
[rXa cXa] = size(Xphia);
[rXp cXp] = size(Xp);

if ~(rX1 == ncensus & rXa == ncensus & rXp == ncensus),
   error('wrong size for covariate matrices for census')
end

%parameters
beta1=thet(1:cX1);                       %application specific,
betaa=thet(cX1+1:cX1+cXa);               %here for census with Q=Q(t)
betap=thet(cX1+cXa+1:cX1+cXa+cXp);       %time-variant and T=T(X)
sigma=exp(thet(cX1+cXa+cXp+1));          %depending on covariates

%Kalman filtering
if strcmp(atype,'time'),
  Z=[0 1];
  a0=zeros(2,1); P0=1e6*eye(2);             %diffuse prior
  for i=1:ncensus
    phi1=ilogit(Xphi1(i,:)*beta1);          %ie logit(phi1) = Xphi1*beta1
    phia=ilogit(Xphia(i,:)*betaa);          %ie logit(phiA) = Xphia*betaA
    p=exp(Xp(i,:)*betap);                   %ie log(p) = Xp*betap
    T(:,:,i)=[0 p*phi1; phia phia];
    q11=p*phi1*(1-p*phi1)*a0(2);
    q22=phia*(1-phia)*(a0(1)+a0(2));
    Q=diag([q11,q22]);
    %predict
    a=T(:,:,i)*a0; P=T(:,:,i)*P0*T(:,:,i)'+Q; P=(P+P')/2;
    apred(:,i)=a; Ppred(:,:,i)=P;
    %prediction errors
    g(i)=Z*a;
    v(i)=Census(i)-Z*a; F(i)=Z*P*Z'+sigma^2;
    %update
    K=P*Z'*inv(F(i));
    a0=a+K*v(i); P0=P-K*Z*P; P0=(P0+P0')/2;
    acorr(:,i)=a0; Pcorr(:,:,i)=P0;
  end
elseif strcmp(atype,'age'),
  Z=[0 ones(1,cXa)];
  a0=zeros(cXa+1,1); P0=1e6*eye(cXa+1);     %diffuse prior
  for i=1:ncensus
    phi1=ilogit(Xphi1(i,:)*beta1);          %ie logit(phi1) = Xphi1*beta1
    phia=ilogit(betaa);                     %constant throughout
    p=exp(Xp(i,:)*betap);                   %ie log(p) = Xp*betap
    T(:,:,i)=diag(phia,-1); T(1,2:cXa+1,i)=p*phi1;
    T(cXa+1,cXa+1,i)=phia(cXa);
    Q=diag(T(:,:,i).*(1-T(:,:,i))*a0); %Q(1,1)=T(1,:,i)*a0;
    %predict
    a=T(:,:,i)*a0; P=T(:,:,i)*P0*T(:,:,i)'+Q; P=(P+P')/2;
    apred(:,i)=a; Ppred(:,:,i)=P;
    %prediction errors
    g(i)=Z*a;
    v(i)=Census(i)-Z*a; F(i)=Z*P*Z'+sigma^2;
    %update
    K=P*Z'*inv(F(i));
    a0=a+K*v(i); P0=P-K*Z*P; P0=(P0+P0')/2;
    acorr(:,i)=a0; Pcorr(:,:,i)=P0;
  end
end

%log-likelihood
l=-ncensus*log(2*pi)/2-sum(log(F))/2-sum(v.*v./F)/2;
f=-l;

%Kalman smoothing
asmooth=a0; h(:,ncensus)=asmooth;
for i=ncensus-1:-1:1
 Pstar=Pcorr(:,:,i)*T(:,:,i+1)'*inv(Ppred(:,:,i+1));
 asmooth=acorr(:,i)+Pstar*(asmooth-T(:,:,i+1)*acorr(:,i));
 h(:,i)=asmooth;
end
h=h';
