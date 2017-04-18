function p = pchisq(x,n)
% pchisq(x,df)
% distribution function for chisq on df d.f.

a = int2str(n);
na = length(a);
digit = a(na);
even = (digit=='0' | digit=='2' | digit=='4' | digit=='6' | digit=='8');

if even,
  seq = 1:n/2;
%  p = sum( exp(-x/2)*(x/2).^(seq-1)./gamma(seq) );
  p = sum( exp( (-x/2) + log(x/2).*(seq-1) - gammaln(seq) ) );

elseif n==1,
  p = 1-erf(sqrt(x/2));

else
  seq = (3:2:n)/2;
%  p = 1-erf(sqrt(x/2)) + sum( exp(-x/2)*(x/2).^(seq-1)./gamma(seq) );
  p = 1-erf(sqrt(x/2)) + sum( exp( (-x/2)+log(x/2).*(seq-1)-gammaln(seq) ) );
end


