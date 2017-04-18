function jay = infoT(x,X1,Xa,Xlam,atype)  

global Ringed Mobs unrec nrows ncols const Toe ellamx;

delta = 1e-6;
t = length(x);
jay = zeros(t);

Delta = delta*eye(t);

P = probT(x,X1,Xa,Xlam,atype) + tril(ones(nrows,ncols),-1);
R = triu(Ringed*ones(1,ncols));
Q = ones(nrows,1) - sum(triu(P)')';                            

for r=1:t,
  dr = Delta(:,r);
  eval(['dP',int2str(r),' = (probT(x+dr,X1,Xa,Xlam,atype) - probT(x-dr,X1,Xa,Xlam,atype))/(2*dr(r));'])
  eval(['dQ',int2str(r),' = -(sum(dP',int2str(r),'''','))','''',';'])
end

for r=1:t,
  for s = 1:t,
    eval(['jay(r,s) = sum(sum(R.*dP',int2str(r),'.*dP',int2str(s),'./P)) + sum(Ringed.*dQ',int2str(r),'.*dQ',int2str(s),'./Q);'])
  end
end

% for r=1:t,
%   for s = 1:t,
%     dr = Delta(:,r);
%     ds = Delta(:,s);
%     dPr = (probT(x+dr,X1,Xa,Xlam,atype) - probT(x-dr,X1,Xa,Xlam,atype))/(2*dr(r));
%     dPs = (probT(x+ds,X1,Xa,Xlam,atype) - probT(x-ds,X1,Xa,Xlam,atype))/(2*ds(s));
%     dQr = -(sum(dPr'))';
%     dQs = -(sum(dPs'))';
%     jay(r,s) = sum(sum(R(P>0).*dPr(P>0).*dPs(P>0)./P(P>0)))+ sum(Ringed.*dQr.*dQs./Q);
%   end
% end

jay = (jay+jay')/2;

