function P = probT(par,X1,Xa,Xlam,atype)         
% creates the matrix of survival probabilities
% par = parameter vector in logistic transform space

global Ringed Mobs unrec nrows ncols const Toe ellmax;

if isempty(X1),
   X1 = ones(nrows,1);
end
if isempty(Xa),
   Xa = ones(ncols-1,1);
end
if isempty(Xlam),
   Xlam = ones(ncols,1);
end

[dum1 n1] = size(X1);
[dum2 na] = size(Xa);
[dum3 nlam] = size(Xlam);

if ~(dum1 == nrows & dum2 == ncols-1 & dum3 == ncols),
  fprintf('wrong size for covariate matrices in probT\n')
end
if ~(n1+na+nlam == length(par)),
  fprintf('wrong length for parameter vector in probT\n')
end

b1 = par(1:n1);
ba = par(n1+1:n1+na);
blam = par(n1+na+1:n1+na+nlam);

phi1 = ilogit(X1*b1);
phia = ilogit(Xa*ba) ;
lam = ilogit(Xlam*blam);

Phi1 = [diag(1-phi1) + triu(phi1*ones(1,nrows),1), phi1*ones(1,ncols-nrows)];

Phia = ones(nrows,ncols);
e = [1; zeros(nrows-1,1)];

if strcmp(atype,'age'),
  Phia = toeplitz(e,[1;1-phia]).*toeplitz(e,cumprod([1;1;phia(1:ncols-2)]));
elseif strcmp(atype,'time'),
% this looks like it ought to be fast, but it isn't
%  for i = 1:nrows,
%    Phia(i,:) = [zeros(1,i-1),1,((1-phia(i:ncols-1)).*cumprod([1;phia(i:ncols-2)]))'];
%  end
  for i = 1:nrows,
    for j = i+1:ncols,
       Phia(i,j) = (1-phia(j-1))*prod(phia(i:j-2));
    end
  end;
else fprintf('phiatype must be  age or time\n')
end

Lambda = [triu(ones(nrows,1)*lam(1:nrows)'), ones(nrows,1)*lam(nrows+1:ncols)'];

P = Phi1.*Phia.*Lambda;
