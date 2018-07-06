% UC_EM.m

load 'USCPI.csv';
y = USCPI;
T = length(y);
omega2_0 = 9; % initial condition
omega = .5^2; % fix omega
H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T);
invOmega = sparse(1:T,1:T,[1/omega2_0 ...
1/omega*ones(1,T-1)],T,T);
HinvOmegaH = H'*invOmega*H;
sigma2t = var(y); % initial guess
err = 1;
while err> 10^(-4)
        % E-step
    Kt = HinvOmegaH + speye(T)/sigma2t;
    taut = Kt\(y/sigma2t);
        % M-step
    lam = eig(Kt);
    newsigma2t = (sum(1./lam) + (y-taut)'*(y-taut))/T;
        % update
    err = abs(sigma2t-newsigma2t);
    sigma2t = newsigma2t;
end
Kt = HinvOmegaH + speye(T)/sigma2t;
taut = Kt\(y/sigma2t);
hold on
plot(taut);
plot(y);
hold off
