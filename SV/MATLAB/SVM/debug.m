clear all
close all
addpath(genpath('../../MATLAB'));


T = 3000;
ind_h_sel = 2:50:T;

beta = -0.2;
mu = -1;
theta_true = [mu, 0.98, 0.01, beta];

hyper.S = [5/2, 0.05/2];
hyper.M = 10;
hyper.P = [20, 3/2];
hyper.B = 10;

bin_range = 4;
N_bin = 50; 
bins = linspace(-bin_range,bin_range,N_bin+1);
bin_midpoint = (bins(1:N_bin) + bins(2:N_bin+1))/2;
% bins are the demeaned volatilities: b=h-mu
N_q = 20; %30;
mid = (0:(N_q-1))/N_q;
mid = mid + mid(2)/2;
mid_inv = norminv(mid); 

M = 10000;
BurnIn = 1000; 

delta_HMM.h = 0.3;
delta_HMM.t = [0.05, 0.008, 0.0008, 0.07]; %      

% delta_HMM_adapt.h = 0.4;
% delta_HMM_adapt.t = [0.25, 0.009, 0.0009, 0.08]; %      

ss = 1;


s = RandStream('mt19937ar','Seed',ss);
RandStream.setGlobalStream(s); 

[y,h_true,y0] = generate_SVM(theta_true,T);
h_init = log(var(y))*ones(1,T); 


h = h_true;
theta = theta_true;  

loglik = loglik_h_HMM(y, h, theta, bins, bin_midpoint);
sumloglik = sum(loglik_h_HMM(y, h, theta, bins, bin_midpoint)); 

[theta, A_theta] = update_theta_HMM_RW(y, h, theta, bins, bin_midpoint, hyper, delta_HMM.t);
[h, acc, A_h] = update_h_HMM(y, h, theta, delta_HMM.h, bins, bin_midpoint);


% adaot     
% loglik = loglik_h_HMM_adapt(y, h, theta_true, mid_inv); 
% 
% [theta, A_theta] = update_theta_HMM_RW_adapt(y, h, theta, mid_inv,  hyper, delta_HMM_adapt.t);
% [h, acc, A_h] = update_h_HMM_adapt(y, h, theta, delta_HMM_adapt.h, mid_inv);



