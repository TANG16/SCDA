clear all
close all
plot_on = false;

sc = 1;
[y, T, time, stdT, f, m, T1, T2] = BKM_Data_HMM(sc);

Na = round( [1000, 1000, 1092.23, 1100.01, 1234.32, 1460.85, 1570.38, 1819.79,...
    1391.27, 1507.60, 1541.44, 1631.21, 1628.60, 1609.33, 1801.68, 1809.08, 1754.74,...
    1779.48, 1699.13, 1681.39, 1610.46, 1918.45, 1717.07, 1415.69, 1229.02, 1082.02,...
    1096.61, 1045.84, 1137.03, 981.1, 647.67, 992.65, 968.62, 926.83, 952.96, 865.64]/sc);

alpha1 = 1;
alphaa = 2;
alphar = -2; 
if (sc == 1)
    alphal = -4;
else
    alphal = -1;
end
beta1 =-2;
betaa = 0.1;
betar = -0.7;
betal = -0.3;
sigy = 1;

% alpha1 = 1;
% alphaa = 2;
% alphar = -2;
% alphal = -4;
% beta1 = -0.19;
% betaa = 0.1;
% betar = -1.106;
% betal = -0.3;
% sigy = 1;

params = {'alpha1', 'alphaa', 'alphar', 'alphal', ...
    'beta1', 'betaa', 'betar', 'betal',...
    'sigy'};

theta_init = [alpha1, alphaa, alphar, alphal, beta1, betaa, betar, betal, sigy];
[phi1, phia, rho, lambda] = BKM_covariates(theta_init,f,stdT);  

D = size(theta_init,2);
prior.N = [2000/sc 0.5];
prior.S = [0.001,0.001];
prior.T_mu = 0*ones(D-1,1);
prior.T_sigma2 = 100*ones(D-1,1);
priorN = prior.N;


%% Set the proposals
% for the parameters
update_T = 'NRW'; % 'NRW' or 'URW'

% step sizes 
% given step size delta, std for URW is delta/sqrt(3), for NRW 1*delta
% 0.5 of posterior st. dev. turns out to be: 
% [0.04 0.04 0.1 0.02 0.03 0.02 0.06 0.02]
% from JAGS [0.04 0.04 0.05 0.02 0.03 0.02 0.03 0.02]
if strcmp(update_T,'URW')
%     delta.T = sqrt(3)*[0.04 0.04 0.1 0.02 0.03 0.02 0.06 0.02];
    delta.T = sqrt(3)*[0.04 0.04 0.05 0.02 0.03 0.02 0.03 0.02];
elseif strcmp(update_T,'NRW')
%     delta.T = [0.04 0.04 0.07 0.02 0.03 0.02 0.07 0.02];
%     delta.T = 0.1*ones(D-1,1);
%     delta.T = [0.1 0.05 0.05 0.1 0.1 0.05 0.05 0.1]; %63
%     delta.T = [0.1 0.2 0.2 0.1 0.1 0.2 0.2 0.1]; %63
%     delta.T = [0.1 0.02 0.02 0.1 0.1 0.02 0.02 0.1]; %65
%     delta.T = [0.1 0.02 0.05 0.1 0.1 0.05 0.02 0.1]; %65
%  0.3799    0.5898    0.3390    0.3541  0.3674    0.2788    0.6205    0.4199
%     delta.T = [0.1 0.04 0.05 0.1 0.1 0.02 0.07 0.1]; %65
%  0.3746    0.3781    0.3262    0.3647  0.3539    0.5333    0.2643    0.4159
%     delta.T = [0.04 0.04 0.07 0.02 0.03 0.02 0.07 0.02];
%  0.6587    0.3899    0.2509    0.8101  0.7014    0.5417    0.2695    0.8402
%     delta.T = [0.07 0.04 0.07 0.1 0.1 0.02 0.07 0.1];
%  0.4779    0.3774    0.2531    0.3678 0.3609    0.5396    0.2664    0.4191
    delta.T = [0.1 0.04 0.05 0.1 0.1 0.035 0.05 0.12];
%  0.3834    0.3871    0.3395    0.3626 0.3776    0.3843    0.3550    0.3661
else
    delta.T = 0.1*ones(D-1,1);
end
deltaT = delta.T;


if (sc == 10)
    delta.N = 30/sc + 0.5; % std of the unif RW update 
elseif (sc == 1)
    % OLD uniform delta.T
    % with 4.5 the mean acceptence probability was 82% 
    % with 10.5 the mean acceptence probability was 24%    
    % with 6.5 the mean acceptence probability was 80%    
    % with 7.5 the mean acceptence probability was 78%
    % with 8.5 the mean acceptence probability was 76% 
    % with 9.5 the mean acceptence probability was 48% 
    % rate 21% ?] AND CRASH! 
%     delta.N = 8 + 0.5; % std of the unif RW update 
    % NEW 0.5 post st dev delta.T
%     delta.N = 8 + 0.5; % mean acceptance 66 mean 66 
    delta.N = 130 + 0.5; % [mean acceptance rate 75 mean acceptance prob 75]
%     delta.N = 8+ 0.5; % mean acceptance rate 66 mean probability 66 
%     delta.N = 10 + 0.5; % mean acceptance rate 81 mean probability 82 
else
    delta.N = 1.5;
end
deltaN = delta.N;


% No bins
N_max = 69 * 10/sc;
IND = 0:N_max;
% Bins
%  Bins' midpoints  
N_bin = 30;
bin_size = 29;
% bin = zeros(N_bin,1);
% for ii = 0:N_bin
%     bin(ii+1) = 0.5*(bin_size*(2*ii+1)-1); % ith bin's midpoint
% end
bin = 0.5*(bin_size*(2*(1:N_bin)+1)-1)';
logfact = @(xx) sum(log(1:1:xx));
logfact = arrayfun(logfact,(0:7000)') ;

% oldlikhood = BKM_calclikhood_HMM(Na, theta_init, y, m, f, stdT, prior.N, N_max, logfact);
oldlikhood = BKM_calclikhood_HMM_bin(Na, theta_init, y, m, f, stdT, prior.N, bin, logfact);

M = 10000;
BurnIn = 1000;
N = Na;
theta = theta_init;
% theta(9) = 30000;
NN = zeros(T,M);
sample = zeros(M,9);
accept = zeros(M,T+D-1);
mean_A = zeros(M,1);

tic
% profile on
for ii = -BurnIn:M
    % Update the parameters in the model using function "updateparam": 
    % Set parameter values and log(likelihood) value of current state to be the output from
    % the MH step:
    
    if (mod(ii,1000)==0)
        fprintf('MH iter = %i\n',ii); toc;
    end
    if (ii>0)
%     [N, theta, acc, a_sum] = BKM_update_HMM(N, theta, prior, delta, y, m, f, stdT, N_max, logfact);
        [N, theta, acc, a_sum] = BKM_update_HMM_bin(N, theta, prior, delta, y, m, f, stdT, bin, logfact);
        NN(:,ii) = N;
        sample(ii,:)= theta; 
        accept(ii,:) = acc; 
        mean_A(ii) = a_sum;
    end
end
time_sampl = toc; 
mean_A = mean_A/(T+D-1);
mean_accept = mean(accept);

% name = 'BKM_HMM_results_sigma2.mat';
% save(name,'sample','NN','accept','theta_init','prior','delta','time_sampl', 'accept', 'mean_A');

% with poissrnd and binopdf 4.97 sec per draw
% with explicit formulae 0.5055 sec per draw
% with explicit formula with logfact as array 0.0170 sec per draw
% profile off
% profile viewer

std(sample(:,1:8)) %   0.1682    0.1449    0.1461    0.4110    0.2104    0.0788    0.0654    0.0928
std(sample(:,9))
std(sample(M/2+1:M,1:8)) %    0.0677    0.0652    0.0791    0.0351    0.0626    0.0463    0.0287    0.0396
std(sample(M/2+1:M,9)) %  8.9622e+03

mean(sqrt(var(NN,1,2))) %    60.3553
mean(sqrt(var(NN(:,M/2+1:M),1,2))) % 42.9882
min(sqrt(var(NN(:,M/2+1:M),1,2))) % 19.8906 ==> 1/2 post st dev = 9.5


BurnIn = 0;
BurnIn = 2000; %5000;
