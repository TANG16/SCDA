% TVPAR.m
nloop = 11000;
burnin = 1000;
load 'USCPI.csv';
p = 2;   % no of lags
y0 = USCPI(1:p);
y = USCPI(p+1:end);
T = length(y);
q = p+1; % dim of states
Tq = T*q;
         %% prior
asigma2 = 5; lsigma2 = 1*(asigma2-1);
aomega2 = 5; lomega2 = (aomega2-1)*[0.5^2; 0.1^2*ones(p,1)];
invOmega0 = ones(q,1)/5;
         %% initialize
omega2 = .1*ones(q,1);
sigma2 = 1;
store_omega2 = zeros(nloop-burnin,q);
store_sigma2 = zeros(nloop-burnin,1);
store_beta = zeros(Tq,1);
         %% construct/compute a few things
X = [ones(T,1) [y0(end); y(1:end-1)] [y0; y(1:end-2)]];
bigX = SURform(X);
H = speye(Tq) - sparse(q+1:Tq,1:(T-1)*q, ...
ones(1,(T-1)*q),Tq,Tq);
newaomega2 = aomega2 + (T - 1)/2;
newasigma2 = asigma2 + T/2 ;
for loop = 1:nloop
        %% sample beta
    invS = sparse(1:Tq,1:Tq,[invOmega0' ...
    repmat(1./omega2',1,T-1)]);
    K = H'*invS*H + bigX'*bigX/sigma2;
    C = chol(K,'lower');
    betahat = K\(bigX'*y/sigma2);
    beta = betahat + C'\randn(Tq,1);
        %% sample omega2
    erromega2 = reshape(H*beta,q,T);
    newlomega2 = lomega2 + sum(erromega2(:,2:end).^2,2)/2;
    omega2 = 1./gamrnd(newaomega2, 1./newlomega2);
        %% sample sigma2
    newlsigma2 = lsigma2 + sum((y-bigX*beta).^2)/2;
    sigma2 = 1/gamrnd(newasigma2, 1./newlsigma2);
    if loop>burnin
        i=loop-burnin;
        store_beta = store_beta + beta;
        store_omega2(i,:) = omega2';
        store_sigma2(i,:) = sigma2;
    end
end
betahat = store_beta/(nloop-burnin);
sigma2hat = mean(store_sigma2);
omega2hat = mean(store_omega2);
plot(betahat(1:3:3*T)) %plot beta_{0t}
