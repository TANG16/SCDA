% UCSV.m
load 'USCPI.csv';
y = USCPI;
T = length(y);
nloop = 11000;
burnin = 1000;
%% prior
Vtau = 9; Vh = 9;
atau = 10; ltau = .25^2*(atau-1);
ah = 10; lh = .2^2*(ah-1);
%% initialize the Markov chain
omega2tau = .25^2;
omega2h = .2^2;
h = log(var(y)*.8)*ones(T,1);
H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T);
%% initialize for storage
store_omega2tau = zeros(nloop - burnin,1);
store_omega2h = zeros(nloop - burnin,1);
store_tau = zeros(nloop - burnin,T);
store_h = zeros(nloop - burnin,T);
%% compute a few things
newatau = (T-1)/2 + atau;
newah = (T-1)/2 + ah;
for loop = 1:nloop
    %% sample tau
    invOmegatau = sparse(1:T,1:T, ...
        [1/Vtau 1/omega2tau*ones(1,T-1)]);
    invSigy = sparse(1:T,1:T,exp(-h));
    Ktau = H'*invOmegatau*H + invSigy;
    Ctau = chol(Ktau,'lower');
    tauhat = Ktau\(invSigy*y);
    tau = tauhat + Ctau'\randn(T,1);
    %% sample h
    ystar = log((y-tau).^2 + .0001 );
    h = SVRW(ystar,h,omega2h,Vh);
    %% sample omega2tau
    newltau = ltau + sum((tau(2:end)-tau(1:end-1)).^2)/2;
    omega2tau = 1/gamrand(newatau, newltau);
    %% sample omega2h
    newlh = lh + sum((h(2:end)-h(1:end-1)).^2)/2;
    omega2h = 1/gamrand(newah, newlh);
    if loop>burnin
        i = loop-burnin;
        store_tau(i,:) = tau';
        store_h(i,:) = h';
        store_omega2tau(i,:) = omega2tau;
        store_omega2h(i,:) = omega2h;
    end
end
tauhat = mean(store_tau)';
hhat = mean(store_h);
figure(1)
plot(tauhat)
figure(2)
plot(hhat)