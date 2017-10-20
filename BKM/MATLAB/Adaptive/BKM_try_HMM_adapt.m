clear all
close all
plot_on = false;

sc = 1;
[y, T, time, stdT, f, m, T1, T2] = BKM_Data_HMM(sc);

Na = round( [1000, 1000, 1092.23, 1100.01, 1234.32, 1460.85, 1570.38, 1819.79,...
    1391.27, 1507.60, 1541.44, 1631.21, 1628.60, 1609.33, 1801.68, 1809.08, 1754.74,...
    1779.48, 1699.13, 1681.39, 1610.46, 1918.45, 1717.07, 1415.69, 1229.02, 1082.02,...
    1096.61, 1045.84, 1137.03, 981.1, 647.67, 992.65, 968.62, 926.83, 952.96, 865.64]/sc);

% Na = [1013 1000 1119 1181 1351 1435 1494 1534 1535 1685 1930 2091 2332 2398,...
%     2452 2490 2178 2172 2299 2168 2079 2126 1900 1689 1486 1504 1531 1590,...
%     1523 1408 1346 1281 1269 1154 1043 1024];
%  alpha1 = 0.5, alphaa = 2, alphar = -1, alphal = -4,
%                          beta1 =-0.2, betaa = -0.25, betar = -0.14, betal = -0.37)
alpha1 = 0.5;
alphaa = 2;
alphar = -1; 
alphal = -4; 
beta1 = -0.2;
betaa = -0.25;
betar = -0.14;
betal = -0.37;
sigy = 1; 

params = {'alpha1', 'alphaa', 'alphar', 'alphal', ...
    'beta1', 'betaa', 'betar', 'betal',...
    'sigy'};

theta_init = [alpha1, alphaa, alphar, alphal, beta1, betaa, betar, betal, sigy];
% theta_init = [0.553014371426413 1.39106960372498 -0.964995157020072,...
%     -4.63076743939019 -0.376963675487789 -0.217044463641838,...
%     -0.337437817454555 -0.326383747699795 137274.149642374];
[phi1, phia, rho, lambda] = BKM_covariates(theta_init,f,stdT);  

D = size(theta_init,2);
prior.N = [2000/sc 0.5];
prior.S = [0.001,0.001];
prior.T_mu = 0*ones(D-1,1);
prior.T_sigma2 = 100*ones(D-1,1);
priorN = prior.N;


%% Set the proposals
% for the parameters
update_T = 'NRW';  

% step sizes 
% given step size delta, std for URW is delta/sqrt(3), for NRW 1*delta
% 0.5 of posterior st. dev. turns out to be: 
% [0.04 0.04 0.1 0.02 0.03 0.02 0.06 0.02]
% from JAGS [0.04 0.04 0.05 0.02 0.03 0.02 0.03 0.02]
delta.T = [0.1 0.04 0.05 0.1 0.1 0.035 0.05 0.12];
% delta.T = [0.04 0.04 0.05 0.02 0.03 0.02 0.03 0.02];

deltaT = delta.T;
%     delta.N = 130 + 0.5;  
% delta.N = 80 + 0.5;  %mean(mean_accept(1:T)) = 0.3952
delta.N = 100 + 0.5;  %mean(mean_accept(1:T)) = 

deltaN = delta.N;


% No bins
N_max = 69 * 10/sc;
IND = 0:N_max;
% Bins
%  Bins' midpoints  
N_bin = 30; 
bin_size = 29;
bin = 0.5*(bin_size*(2*(1:N_bin)+1)-1);
% Quanatiles
N_q = 20; % 30 time_sampl = 186.4378; 10  time_sampl = 192.1283

qu = (0:(N_q-1))/N_q;
qu_mid = qu + qu(2)/2; 
mid = norminv(qu_mid);  
  
logfact_fun = @(xx) sum(log(1:1:xx));
logfact = arrayfun(logfact_fun,0:7000) ;

% oldlikhood = BKM_calclikhood_HMM(Na, theta_init, y, m, f, stdT, prior.N, N_max, logfact);
% oldlikhood = BKM_calclikhood_HMM_bin(Na, theta_init, y, m, f, stdT, prior.N, bin, logfact);
oldlikhood = BKM_calclikhood_HMM_adapt(Na, theta_init, y, m, f, stdT, prior.N, mid, logfact);

M = 10000;
N = Na;
theta = theta_init;
% theta(9) = 30000;
NN = zeros(T,M);
sample = zeros(M,9);
accept = zeros(M,T+D-1);
mean_A = zeros(M,1);

tic
% profile on
for ii = 1:M
    % Update the parameters in the model using function "updateparam": 
    % Set parameter values and log(likelihood) value of current state to be the output from
    % the MH step:
    
    if (mod(ii,1000)==0)
        fprintf('MH iter = %i\n',ii); toc;
    end
%     [N, theta, acc, a_sum] = BKM_update_HMM(N, theta, prior, delta, y, m, f, stdT, N_max, logfact);
%     [N, theta, acc, a_sum] = BKM_update_HMM_bin(N, theta, prior, delta, y, m, f, stdT, bin, logfact);
    [N, theta, acc, a_sum] = BKM_update_HMM_adapt(N, theta, prior, delta, y, m, f, stdT, mid, logfact);
    NN(:,ii) = N;
    sample(ii,:)= theta; 
    accept(ii,:) = acc; 
    mean_A(ii) = a_sum;
end
time_sampl = toc; 
mean_A = mean_A/(T+D-1);
mean_accept = mean(accept);

% name = 'BKM_HMM_results_sigma2.mat';
% save(name,'sample','NN','accept','theta_init','prior','delta','time_sampl', 'accept', 'mean_A');
 
std(sample(:,1:8))  
std(sample(:,9))
std(sample(M/2+1:M,1:8))  
std(sample(M/2+1:M,9))  

mean(sqrt(var(NN,1,2)))  
mean(sqrt(var(NN(:,M/2+1:M),1,2)))  
min(sqrt(var(NN(:,M/2+1:M),1,2)))  