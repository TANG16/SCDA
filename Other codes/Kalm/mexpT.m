function A = mexpT(par,Xphi1,Xphia,Xlam,atype)
% expected numbers of recoveries

global Ringed Mobs unrec nrows ncols const Toe ellmax;

P = feval('probT',par,Xphi1,Xphia,Xlam,atype);

A = Ringed*ones(1,ncols).*P;
