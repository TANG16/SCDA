clear all
close all

% sc = 100;
sc = 1;

[y, f, m, rel, time, T, T_ring] = Heron_data(sc);

N_bin1 = 50;
N_bin3 = 30;

%Lower bin values:
qu1 =  (0:(N_bin1-1))/N_bin1;
mid1 = qu1+qu1(2)/2;
  
% Lower bin values:
qu3 = (0:(N_bin3-1))/N_bin3;
mid3 = qu3+qu3(2)/2;
  
logfact = @(x) sum(log(1:x));

% alpha = zeros(4,1);
% beta = zeros(4,1);
% alphal = 0;
% betal = 0;
% alpharho = 0;
% tauy = 1;%10;%1000;

alpha = [-0.55, -0.1, 0.5, 1.1]';
beta = [0.02, -0.15, -0.15, -0.2]';
alphal = 1.65;
betal = -0.2;
alpharho = 1.1;
tauy = 0.001000;

ind = alpha(1) + beta(1)*f(1:T-1);
phi1 = exp(ind)./(1+exp(ind));
ind = alpha(2) + beta(2)*f(1:T-1);
phi2 = exp(ind)./(1+exp(ind));
ind = alpha(3) + beta(3)*f(1:T-1);
phi3 = exp(ind)./(1+exp(ind));
ind = alpha(4) + beta(4)*f(1:T-1);
phi4 = exp(ind)./(1+exp(ind));

ind = alpharho*ones(1,T-1); 
rho = exp(ind);

ind = alphal + betal*(time(1:T_ring));
lambda =  exp(ind)./(1+exp(ind));

C = 1000000;
% pi = 3.14159265359;

dummy = zeros(1,T);
PHI = zeros(1,T);
loglik = zeros(1,T);
loglam1 = zeros(1,T-1); 


Up2 = 4000;
Up4 = 6000; 

X2 = round([1500, y(2:end)*2/7]);
X4 = round([2000, y(2:end)*4/7]);


P2 = zeros(N_bin1,T);
P4 = zeros(N_bin3,T);
Q = zeros(N_bin3,T);

% prior for first transition/augmented observations/observation probabilities
for t=1:2
    for ii = 0:N_bin1
        P2(ii+1,t) = 1/Up2;
    end

    for ii = 0:N_bin3
        P4(ii+1,t) = 1/Up4;
        Q(ii+1,t) = 1;
    end    
end


midbin1 = zeros(N_bin1, T);
midbin3 = zeros(N_bin3, T);


for t = 1:(T-1)
    loglam1(t) = log(X4(t)) + log(rho(t)) + log(phi1(t));
end

for t = 3:T
    for ii = 0:(N_bin1-1) % X2 (depends only on [=integrated] X1)
        midbin1(ii+1,t) = norminv(mid1(ii+1), exp(loglam1(t-2)), sqrt(exp(loglam1(t-2))));
        P2(ii+1,t) = binopdf(X2(t), round(midbin1(ii+1,t)), phi2(t-1)); %binopdf(X,N,P) 
        %   P2(ii+1,t) = exp(X2(t)*log(phi2(t-1)) + (ii - X2(t))*log(1-phi2(t-1)) + logfact(ii) - logfact(abs(ii - X2(t))) - logfact(X2(t)));
    end

    for ii = 0:(N_bin3-1)  
        % X4 (depends on [=integrated] X3)
        midbin3(ii+1,t) = norminv(mid3(ii+1), X2(t-2)*phi3(t-2), sqrt(X2(t-2)*phi3(t-2)*(1-phi3(t-2))));
        P4(ii+1,t) = binopdf(X4(t), X4(t-1) + round(midbin3(ii+1,t)), phi4(t-1));
        %   P4(ii+1,t) = exp(X4(t)*log(phi4(t-1)) + (ii + X4(t-1) - X4(t))*log(1-phi4(t-1)) + logfact(ii + X4(t-1)) - logfact(abs(ii + X4(t-1) - X4(t))) - logfact(X4(t)));       
        Q(ii+1,t) = exp(-tauy*((y(t) - (X2(t) + round(midbin3(ii+1,t)) + X4(t))).^2)/2)*sqrt(tauy)/sqrt(2*pi);
    end
end

figure(200)
set(gcf, 'Units', 'normalized', 'Position', [0.1 0 0.9 1]);
subplot(2,3,1)
hold all
for t = 1:T
    plot(G1(:,t))
end
title('G1')

subplot(2,3,2)
hold all
for t = 1:T
    plot(G3(:,t))
end
title('G3')

subplot(2,3,4)
hold all
for t = 1:T
    plot(P2(:,t))
end
title('P2')

subplot(2,3,5)
hold all
for t = 1:T
    plot(P4(:,t))
end
title('P4')

subplot(2,3,6)
hold all
for t = 1:T
    plot(Q(:,t))
end
title('Q')

% [maxP2 , inP2] = max(P2);
% figure(11)
% plot(inP2)
% plot(X1,'r')
% plot(X2,'g')
% 
% figure(12)
% plot(phi2)

sumsum = zeros(1,T);
for t = 1:T
    sumsum(t) = exp(log(sum(P2(:,t))) + log(sum(P4(:,t) .* Q(:,t))));
    loglik(t) = log(sumsum(t));
    PHI(t) = -loglik(t) + C;
    dummy(t) = poisspdf(0,PHI(t));
end

LL = sum(loglik) 
% tauy = 10 --> LL = -Inf

% figure(2)
% plot(sumsum)
% 
% 
% theta_init = [alpha', beta', alpharho, tauy];
% loglik_heron_HMM_statespace = @(xx) -Heron_HMM_loglik_statespace(xx, N_bin1, N_bin3, y, f, X2, X4, X2_prior, X4_prior)/T;
% options = optimoptions('fminunc');%,'MaxFunEvals',2000,'MaxIter',1500);
% [theta_mle_ss, ~, ~, ~, ~, ~] = fminunc(loglik_heron_HMM_statespace, theta_init,options);
% 

theta_init = [alpha', beta', alphal, betal];
loglik_heron_recovery = @(xx) -Heron_RR_loglik(xx, m, f)/numel(m);
options = optimoptions('fminunc');%,'MaxFunEvals',2000,'MaxIter',1500);
[theta_mle_rr, ~, ~, ~, ~, ~] = fminunc(loglik_heron_recovery, theta_init,options);
 
% theta_mle_ss =[1.2875    0.2482    0.6782    1.2535   -0.0046   -0.0935   -0.2305    0.0577    0.3901   -0.0504]
