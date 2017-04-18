function s = scoreJT(x0,X1_null,Xa_null,Xlam_null,X1_alt,Xa_alt,Xlam_alt,atype)
% score statistic using 
% expected information matrix instead of observed

global nrows ncols;

if isempty(X1_null),   X1_null = ones(nrows,1);   end
if isempty(Xa_null),   Xa_null = ones(ncols-1,1); end
if isempty(Xlam_null), Xlam_null = ones(ncols,1); end
if isempty(X1_alt),    X1_alt = ones(nrows,1);    end
if isempty(Xa_alt),    Xa_alt = ones(ncols-1,1);  end
if isempty(Xlam_alt),  Xlam_alt = ones(ncols,1);  end

x = stretch(x0,X1_null,Xa_null,Xlam_null,X1_alt,Xa_alt,Xlam_alt);
g = grad('likT',x,X1_alt,Xa_alt,Xlam_alt,atype);    
J = infoT(x,X1_alt,Xa_alt,Xlam_alt,atype);
s = g'/J*g;
