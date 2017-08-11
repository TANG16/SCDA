clear all
close all


plot_on = false;

% theta = [1, 0.97, 0.15^2];
% beta = 0.05;
beta = 0.5;
mu = 2*log(beta);
theta = [mu, 0.98, 0.2^2];
theta_true = theta;

T = 2000;
[y,h_true] = generate_SV(theta,T);

if plot_on
    plot(y)
    plot(h_true)
end

mu = theta(1);
phi = theta(2);
sigma2 = theta(3);

theta_init = [0, 0.95, 0.1^2];
kernel = @(xx) -loglik_h(h_true, xx)/T;
theta_mle = fminunc(kernel,theta_init);
theta_mle = fminunc(kernel,theta_mle);

hyper.S = [5/2, 0.05/2];
hyper.M = 10;
hyper.P = [20, 3/2];

% P0 = sqrt(sigma2/(1-phi^2)); % unconditional st. dev. 0.6120
% bin_range = 5*P0;
bin_range = 4;
N_bin = 30; %50;
bins = linspace(-bin_range,bin_range,N_bin+1);
bin_midpoint = (bins(1:N_bin) + bins(2:N_bin+1))/2;
% bins are the demeaned volatilities: b=h-mu
stdev_y = exp((mu+bin_midpoint)/2);

delta.h = 0.1;
delta_h = delta.h;
M = 10000;
BurnIn = 10000;
h_init = log(var(y))*ones(1,T); 


%% FULL DA - FULL CONDITIONAL DENSITIES FROM KIM ET AL. (1998)
h = h_init;
theta = theta_init;
H_DA = zeros(M,T);
theta_DA = zeros(M,3);

accept_DA = zeros(M,1);
A_sum_DA = zeros(M,1);
A_phi_DA = zeros(M,1);

fprintf('Full DA, full conditional updates\n')
tic
for ii = -BurnIn:1:M 
    if (mod(ii,1000) == 1)
        toc;
    end
    [h, acc, A_s] = update_h(y, h, theta, delta.h);
    [theta, A_phi] = update_theta_cond(h, theta, hyper);
    if (ii > 0)
        theta_DA(ii,:) = theta;
        H_DA(ii,:) = h;
        accept_DA(ii,1) = acc;
        A_sum_DA(ii,1) = A_s;
        A_phi_DA(ii,1) = A_phi;        
    end
end
time_DA = toc;

accept_DA = accept_DA/T;
A_sum_DA = A_sum_DA/T;

%% FULL DA - RANDOM WALK UPDATES
% delta.t = [0.1, 0.01, 0.01]; mean(A_theta_DA_RW):  0.8422    0.3067    0.1098
% delta.t = [0.3, 0.01, 0.005]; mean(A_theta_DA_RW):  0.5300    0.3449    0.1554
delta.t = [0.4, 0.01, 0.002]; % mean(A_theta_DA_RW):     0.5390    0.3013    0.4061


h = h_init;
theta = theta_init;
H_DA_RW = zeros(M,T);
theta_DA_RW = zeros(M,3);

accept_DA_RW = zeros(M,1);
A_sum_DA_RW = zeros(M,1);
A_theta_DA_RW = zeros(M,3);

fprintf('Full DA, RW updates\n')
tic
for ii = -BurnIn:1:M 
    if (mod(ii,1000) == 1)
        toc;
    end
    [h, acc, A_s] = update_h(y, h, theta, delta.h);
    [theta, A_theta] = update_theta_RW(y, h, theta, hyper, delta.t);
    if (ii > 0)
        theta_DA_RW(ii,:) = theta;
        H_DA_RW(ii,:) = h;
        
        accept_DA_RW(ii,1) = acc;
        A_sum_DA_RW(ii,1) = A_s;
        A_theta_DA_RW(ii,:) = A_theta;        
    end
end
time_DA_RW = toc;

accept_DA_RW = accept_DA_RW/T;
A_sum_DA_RW = A_sum_DA_RW/T;

%% SEMI DA:
% integrate out the odd h(t)'s and impute the even ones
% delta.t = [0.1, 0.01, 0.01];
delta.t = [0.4, 0.01, 0.002]; 

h = h_init;
theta = theta_init;

H_HMM_v2 = zeros(M,T);
theta_HMM_v2 = zeros(M,3);

accept_HMM_v2 = zeros(M,1);
A_sum_HMM_v2 = zeros(M,1);
A_theta_HMM_v2 = zeros(M,3);

fprintf('Semi DA, even states imputed\n')
tic
for ii = 1:1:M 
    if (mod(ii,1000) == 1)
        toc;
    end
    [h, acc, A_s] = update_h_HMM_v2(y, h, theta, delta.h, bins, bin_midpoint);
    [theta, A_theta] = update_theta_HMM_RW(h, theta, bins, bin_midpoint, hyper, delta.t);
  
    if (ii > 0)
        theta_HMM_v2(ii,:) = theta;
        
        H_HMM_v2(ii,:) = h;
        accept_HMM_v2(ii,1) = acc;
        A_sum_HMM_v2(ii,1) = A_s;
        A_theta_HMM_v2(ii,:) = A_theta;        
    end
end
accept_HMM_v2 = accept_HMM_v2/(T/2);
A_sum_HMM_v2 = A_sum_HMM_v2/(T/2);
time_HMM = toc; 


%% SEMI DA shift: 
% if shift integrate out the odd h(t)'s and impute the even ones
% if not shift integrate out the even h(t)'s and impute the odd ones
delta.t = [0.1, 0.01, 0.01];

h = h_init;
theta = theta_init;

H_HMM_v2_shift = zeros(M,T);
theta_HMM_v2_shift = zeros(M,3);

accept_HMM_v2_shift = zeros(M,1);
A_sum_HMM_v2_shift = zeros(M,1);
A_theta_HMM_v2_shift = zeros(M,3);

fprintf('Semi DA, shift\n')
tic
for ii = 1:1:M 
    if (mod(ii,1000) == 1)
        toc;
    end
    shift = mod(ii,2);
    [h, acc, A_s] = update_h_HMM_v2_shift(y, h, theta, delta.h, bins, bin_midpoint, shift);
    [theta, A_theta] = update_theta_RW(y, h, theta, bins, bin_midpoint, hyper, delta.t, shift);
  
    if (ii > 0)
        theta_HMM_v2_shift(ii,:) = theta;
        
        H_HMM_v2_shift(ii,:) = h;
        accept_HMM_v2_shift(ii,1) = acc;
        A_sum_HMM_v2_shift(ii,1) = A_s;
        A_theta_HMM_v2_shift(ii,:) = A_theta;        
    end
end
accept_HMM_v2_shift = accept_HMM_v2_shift/(T/2);
A_sum_HMM_v2_shift = A_sum_HMM_v2_shift/(T/2);
time_HMM_shift = toc; 



%% SAVE
name = 'SV_results_param.mat';
save(name,'y','h_true','theta_true',...
    'H_DA','H_HMM_v2','H_HMM_v2_shift',...
    'accept_DA','accept_HMM_v2','accept_HMM_v2_shift',...
    'A_sum_DA','A_sum_HMM_v2','A_sum_HMM_v2_shift',...
    'time_DA','time_HMM','time_HMM_shift',...
    'IF_DA','IF_HMM','IF_HMM_shift');

