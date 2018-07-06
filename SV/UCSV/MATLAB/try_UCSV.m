% function try_UCSV(arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7)
clear all
close all
% 
% s = RandStream('mt19937ar','Seed',1);
% RandStream.setGlobalStream(s); 

plot_on = false;

hyper.h0 = [0, 10];
hyper.g0 = [0, 10];
hyper.omegah = 0.2;
hyper.omegag = 0.2;
hyper.Vtau = 10;
hyper.Vh = 10;
hyper.Vg = 10;

country_id = 7; %1: Canada; 2: France; 3: Germany; 4: Italy; 5: Japan; 6: UK; 7: US

if (country_id == 0)
    T = 250;
    theta_true = [0 ,0, sqrt(0.2), sqrt(0.2)];
    theta = theta_true;
    [y, tau_true, h_true, g_true] = ucsv_sim(theta_true,T, hyper);
else
    data = csvread('../OECD_G7CPI.csv'); %, 'D51:J266');
    y = 400*log(data(2:end,country_id)./data(1:end-1,country_id))';       
    T = length(y);
end


M = 100000; %1000000;% 10000;
BurnIn = 1000 + 3000 + 6000;


h_init = 0.5*log(var(y))*ones(1,T); 
g_init = 0.5*log(var(y))*ones(1,T); 
tau_init = mean(y)*ones(1,T); %y; 
theta_init = [0.5*log(var(y)), 0.5*log(var(y)), sqrt(0.2), sqrt(0.2)];
% theta_init = [0, 0 , sqrt(0.2), sqrt(0.2)];

h = h_init;
g = g_init;
tau = tau_init;
theta = theta_init;

% load('../Post_means.mat')
% h = store_h;
% g = store_g;
% tau = store_tau;
% theta = mean_theta;

h = load('../results csv/mean_h1_DA.csv')';
g = load('../results csv/mean_g1_DA.csv')';
tau = load('../results csv/mean_tau1_DA.csv')';
theta = load('../results csv/mean_theta1_DA.csv')';


delta.t = 0.75; %0.5;%0.55
delta.h = 2; %0.39; %1.0; %0.60;
delta.g = 0.5; %1.5;
% delta.theta =  [0.05, 0.05, 0.75, 0.75]; %   0.0618    0.0735    1.1462    1.3674
delta.theta =  [0.01, 0.01, 0.05, 0.05]; %   0.0618    0.0735    1.1462    1.3674

TAU_DA_RW = zeros(M,T);
H_DA_RW = zeros(M,T);
G_DA_RW = zeros(M,T);
theta_DA_RW = zeros(M,4);

accept_DA_RW = zeros(M,1);
A_TAU_DA_RW = zeros(M,1);
A_H_DA_RW = zeros(M,1);
A_G_DA_RW = zeros(M,1);
A_theta_DA_RW = zeros(M,4);

fprintf('Full DA, RW updates\n')
tic
% for ii = -BurnIn:1:M 
for ii = 1:1:M 
    if (mod(ii,10000) == 0)
        fprintf('Iter = %i,  ',ii)
        toc;
%         fprintf('\n')
    end
    [h, A_h] = update_h(y, tau, h, theta, delta.h, hyper);
    [g, A_g] = update_g(tau, g, theta, delta.g, hyper); 
    [tau, A_tau] = update_tau(y, tau, h, g, theta, delta.t, hyper);

    [theta, A_theta] = update_theta_RW(y, tau, h, g, theta, hyper, delta.theta);
    
    if (ii > 0)
        theta_DA_RW(ii,:) = theta;
        TAU_DA_RW(ii,:) = tau;
        H_DA_RW(ii,:) = h;
        G_DA_RW(ii,:) = g;        

        A_TAU_DA_RW(ii,1) = A_tau;
        A_H_DA_RW(ii,1) = A_h;
        A_G_DA_RW(ii,1) = A_g;
        A_theta_DA_RW(ii,:) = A_theta;        
    end
end
time_DA_RW = toc;
% 
A_TAU_DA_RW = A_TAU_DA_RW/T;
A_H_DA_RW = A_H_DA_RW/T;
A_G_DA_RW = A_G_DA_RW/T;
 
% mean_H_DA_RW = mean(H_DA_RW);
% 
% 
% ESS_H_DA_RW_40 = ESS(H_DA_RW,40);
% ESS_H_DA_RW_100 = ESS(H_DA_RW,100);
% ESS_H_DA_RW_1000 = ESS(H_DA_RW,1000);
% ESS_H_DA_RW_sig = ESS(H_DA_RW,0);
% 
% H_subset_DA_RW = H_DA_RW(:,(2:20)*100);        



figure(111)
subplot(3,1,1)
plot(store_tau);hold on; plot([NaN,tau],'r'); hold off
title('tau')
subplot(3,1,2)
plot(store_h);hold on; plot([NaN,h],'r'); hold off
title('h')
subplot(3,1,3)
plot(store_g);hold on; plot([NaN,g],'r'); hold off
title('g')
legend('Joshua Chan''s code','JAGS')