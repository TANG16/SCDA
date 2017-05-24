clear all
close all

sc=10;
[y, T, time, stdT, f, m, T1, T2] = BKM_Data_HMM(sc);

Na = round( [1000, 1000, 1092.23, 1100.01, 1234.32, 1460.85, 1570.38, 1819.79,...
    1391.27, 1507.60, 1541.44, 1631.21, 1628.60, 1609.33, 1801.68, 1809.08, 1754.74,...
    1779.48, 1699.13, 1681.39, 1610.46, 1918.45, 1717.07, 1415.69, 1229.02, 1082.02,...
    1096.61, 1045.84, 1137.03, 981.1, 647.67, 992.65, 968.62, 926.83, 952.96, 865.64]/sc);
% N1 = 400*ones(1,T)/sc;
N1 = [400   400   400   400   400   400   400   400   400   400   400   400   400   400   400   400    40   400   400 ...
 40   400   400    40    40   400   400   400   400   400   400   400   400   400    40   400   400]/sc;
N = [N1;Na];

% alpha1 = 1;
% alphaa = 2;
% alphar = -2;
% alphal = -4;
% beta1 =-2;
% betaa = 0.1;
% betar = -0.7;
% betal = -0.3;
% sigy = 1;

alpha1 = 1;
alphaa = 2;
alphar = -2;
alphal = -4;
beta1 = -0.19;
betaa = 0.1;
betar = -1.106;
betal = -0.3;
sigy = 1;
 
theta = [alpha1, alphaa, alphar, alphal, beta1, betaa, betar, betal, sigy];
[phi1, phia, rho, lambda] = BKM_covariates(theta,f,stdT);  


Up = 2000;
Na_prior = ones(1,Up+1)/(Up+1) ;
dummy = zeros(1,T);
N_max= 100-1;
C = 1000000;
G = zeros(N_max+1, T);
P = zeros(N_max+1, T);
loglam = zeros(T-1,1);
loglik = zeros(T,1);

% logfact = @(xx) sum(log(1:1:xx));
logfact = @(xx) log(factorial(xx))

for t = 3:(T)
    loglam(t-1) = log(Na(t-1-1)) + log(rho(t-1-1)) + log(phi1(t-1-1));

    for ii = 0:(N_max-1)
        % logfact is the log of the factorial: log(x!)
        G(ii+1,t) = exp(-exp(loglam(t)) + ii*loglam(t) - logfact(ii)); 
        if ((ii + Na(t-1) - Na(t)) > 0)
            P(ii+1,t) = exp(Na(t)*log(phia(t-1)) + (ii + Na(t-1) - Na(t))*log(1-phia(t-1)) + ...
                logfact(ii + Na(t-1)) - logfact(ii + Na(t-1) - Na(t)) - logfact(Na(t-1))); 
        else
            P(ii+1,t) = 0;
        end
    end
    G(N_max+1,t) = max(0,1 - sum(G(1:(N_max),t)));
    if ((N_max + Na(t-1) - Na(t)) > 0)
        P(N_max+1,t) = exp(Na(t)*log(phia(t-1)) + (N_max + Na(t-1) - Na(t))*log(1-phia(t-1)) + ...
            logfact(N_max + Na(t-1)) - logfact(N_max + Na(t-1) - Na(t)) - logfact(Na(t-1))); 
    else
        P(ii+1,t) = 0;
    end    
    loglik(t) = log(sum(G(:,t) .* P(:,t))); % piecewise multiplication enough here
end

figure(1)
hold all
for ii=3:T
    plot(G(:,ii))
end
figure(2)
hold all
for ii=3:T
    plot(P(:,ii))
end
%% Recovery Data
loglik_BKM_recovery = @(xx) -BKM_Bugs_recovery(xx, m, f_r, stdT_r)/numel(m);
theta_init = [1, 2, -4, -2, 0.1, -0.3];
% theta = [alpha1, alphaa, alphal, beta1, betaa, betal]

[theta_mle_r, ~, ~, ~, ~, Sigma_mle_r] = fminunc(loglik_BKM_recovery, theta_init);
Sigma_mle_r = inv(numel(m)*Sigma_mle_r);
% theta_mle = [0.5297    1.5241   -4.5708   -0.1659   -0.2937   -0.3452]
std_mle = sqrt(diag(Sigma_mle_r));
% std_mle = [0.0686    0.0699    0.0352    0.0617    0.0428    0.0390]
[~, phi1_r, phia_r] = BKM_Bugs_recovery(theta_mle_r, m, f, stdT_r);



%% Index Data
theta_init = [1, 2, -2, -2, 0.1, -0.7, 1];
theta = theta_init;
f = f_ss;
stdT = stdT_ss;
% theta = [alpha1, alphaa, alphar, beta1, betaa, betar, sigy]; 
% stdT = stdT_ss;
% f = f_ss;

loglik_BKM_statespace = @(xx) -BKM_Bugs_statespace(xx, y, f_ss, stdT_ss)/T;
options = optimoptions('fminunc');%,'MaxFunEvals',2000,'MaxIter',1500);
[theta_mle_ss, ~, ~, ~, ~, ~] = fminunc(loglik_BKM_statespace, theta_init,options);
[theta_mle_ss, ~, ~, ~, ~, ~] = fminunc(loglik_BKM_statespace, theta_mle_ss);
[theta_mle_ss, ~, ~, ~, ~, Sigma_mle_ss] = fminunc(loglik_BKM_statespace, theta_mle_ss);

[~, N_smooth, phi1_ss, phia_ss, rho] = BKM_Bugs_statespace(theta_mle_ss, y, f_ss, stdT_ss);




%% Combined
theta_init = [1, 2, -2, -4, -2, 0.1, -0.7, -0.3, 1];
% theta = [alpha1, alphaa, alphar, alphals, beta1, betaa, betar, betal, sigy]; 
theta = theta_init;
loglik_BKM_combined = @(xx) -BKM_Bugs_combined(xx, y, f, m)/(length(y)+numel(m));
% [loglik_C, loglik_m, loglik_KF, N_smooth, PHI1, PHIA, RHO, LAMBDA] = BKM_Bugs_combined(theta, y, f, m)
[theta_mle_c, ~, ~, ~, ~, Sigma_mle_c] = fminunc(loglik_BKM_combined, theta_init);
% [0.542 1.558 -1.161 -4.581 -0.163 -0.243 -0.355 -0.362 29218.134]
Sigma_mle_c = inv((numel(m)+length(y))*Sigma_mle_c);
% theta_mle = [0.5297    1.5241   -4.5708   -0.1659   -0.2937   -0.3452]
std_mle_c = sqrt(diag(Sigma_mle_c)) ;
% [0.068 0.071 0.093 0.035 0.061 0.039 0.045 0.040 8013.534]


% alpha1     0.5258 0.072
% alphaa     1.4841 0.079
% alphal    -4.5917 0.035 
% alphar    -1.1017 0.094 
% beta1     -0.1889 0.054 
% betaa     -0.2360 0.038 
% betal     -0.3574 0.039 
% betar     -0.2871 0.037 
% sigy   28822.1313 8817

[loglik_C, loglik_m, loglik_KF, N_smooth, PHI1, PHIA, RHO, LAMBDA] = BKM_Bugs_combined(theta_mle_c, y, f, m);