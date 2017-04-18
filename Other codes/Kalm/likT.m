function ell = likT(par,Xphi1,Xphia,Xlam,atype)
% calculates the negative loglikelihood for a Morgan-Freeman recovery model
% par = (column) vector of parameters in logistic transform space

global Ringed Mobs unrec nrows ncols const Toe ellmax;

if isempty(Xphi1),
   Xphi1 = ones(nrows,1);
end
if isempty(Xphia),
   Xphia = ones(ncols-1,1);
end
if isempty(Xlam),
   Xlam = ones(ncols,1);
end

[dum1 nphi1] = size(Xphi1);
[dum2 nphia] = size(Xphia);
[dum3 nlam] = size(Xlam);

if ~(dum1 == nrows & dum2 == ncols-1 & dum3 == ncols),
   error('wrong size for covariate matrices')
end

P = feval('probT',par,Xphi1,Xphia,Xlam,atype);
Q = ones(nrows,1) - sum(P')';
 
% ell = -sum( sum(Mobs.*log(P+Toe))) - unrec'*log(Q) - const;

ell = -sum( sum(Mobs(P>0).*log(P(P>0)+Toe(P>0)))) - unrec'*log(Q) - const;

