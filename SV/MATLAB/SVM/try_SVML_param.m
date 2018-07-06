% clear all
close all
addpath(genpath('../../MATLAB'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

plot_on = false;
save_on = false;

arg0 = 0;% 2; % 0 -- sim; 2 -- IBM; 100 -- TOPIX (Omori and Zucchini)
arg1 = 1; % DA
arg2 = 1; % fixed
arg3 = 1; % adaptive
path = '';

if (arg0 == 0)
    mu = -1.0;
    phi = 0.98;
    sigma2 = 0.10^2;
    beta = 0; %  -0.2; % NO SV IN THE MEAN
    rho = -0.35;
    theta_true = [mu, phi, sigma2, beta, rho];
%     theta_true = [-1.00, 0.98, 0.1^2, -0.20, -0.35];
    T = 2000;
    [y,h_true] = generate_SVML(theta_true,T);
    [y0,h_true0] = generate_SV(theta_true,T);
    theta_init = theta_true;
    h_init = h_true; 
    data_name = '_sim';
elseif (arg0 == 2)
    y = load('../other/Perc_Rets_GSPC_IBM_MSFT.csv');
    y = y(:,arg0); % arg0 = 1 = GSPC;  arg0 = 2 = IBM
    y = y';
    time = [2000,2018];       
    T = length(y); 
    DATA_NAMES = {'_GSPC','_IBM','_MSFT'};
    data_name = DATA_NAMES{arg0}; %'_GSPC';  
%     path = '../Results/Empirical/';
    fprintf('Data loaded: %s\n',data_name);
    theta_true = [];
    h_true = [];   
%     theta_init = [-1, 0.98, 0.01, -0.2, -0.3];
    theta_init = [0.4, 0.96, 0.08, 0.0045, -0.2];
    h_init = log(var(y))*ones(1,T); 
elseif (arg0 == 100)
    y = xlsread('topix-9802.xls','B:B');
    y=100*diff(log(y))';
    T = length(y); 
    data_name = '_TOPIX';

    fprintf('Data loaded: %s\n',data_name);
    theta_true = [];
    h_true = [];   
    theta_init = [0.37 0.95 0.018 0.0 -0.35];
    h_init = log(var(y))*ones(1,T);     
end

 

hyper.S = [5/2, 0.05/2];
hyper.M = 10;
hyper.P = [20, 3/2];
hyper.B = 10;
hyper.R = 100;

% bins are the demeaned volatilities: b=h-mu
if (arg2 > 0)
    bin_range = 4;
    N_bin = 50; 
    bins = linspace(-bin_range,bin_range,N_bin+1);
    bin_midpoint = (bins(1:N_bin) + bins(2:N_bin+1))/2;
end
if (arg3 > 0)
    N_q = 50;
    mid = (0:(N_q-1))/N_q;
    mid = mid + mid(2)/2;
    mid_inv = norminv(mid);
end

M = 10000;
BurnIn = 10000;
ind_h_sel = 2:50:T;

if (arg0 == 0) % simulation
    delta.h = 0.25;
%     delta.t = [0.32, 0.01, 0.0008, 0.00, 0.04]; % NO UPDATE FOR BETA
    delta.t = [0.0, 0.01, 0.0008, 0.00, 0.04]; % NO UPDATE FOR MU & NO UPDATE FOR BETA
%     delta.t = [0.32, 0.01, 0.0008, 0.065, 0.04]; %  
    
    if (arg2 > 0)
        delta_HMM.h = 0.3;
        delta_HMM.t = [0.05, 0.008, 0.0008, 0.07, 0.1]; %      
    end
    if (arg3 > 0)    
        delta_HMM_adapt.h = 0.4;
%         delta_HMM_adapt.t = [0.25, 0.009, 0.0009, 0.00, 0.1]; %  NO UPDATE FOR BETA
        delta_HMM_adapt.t = [0.0, 0.009, 0.0009, 0.00, 0.1]; %  NO UPDATE FOR MU & NO UPDATE FOR BETA
%         delta_HMM_adapt.t = [0.25, 0.009, 0.0009, 0.08, 0.1]; % 
    end
elseif (arg0 == 100) % TOPIX data (Omori and Z-Book)
    delta.h = 0.25;
    delta.t = [0.32, 0.01, 0.0008, 0, 0.04]; %  NO UPDATE FOR BETA 
    if (arg3 > 0)    
        delta_HMM_adapt.h = 0.4;
        delta_HMM_adapt.t = [0.15, 0.015, 0.002, 0, 0.06]; % 
    end
elseif (arg0 == 2) % IBM
    delta.h = 0.7;
    delta.t = [0.3, 0.01, 0.004, 0.026, 0.04]; %   
    if (arg3 > 0)    
        delta_HMM_adapt.h = 0.95;
%         delta_HMM_adapt.t = [0.27, 0.013, 0.0065, 0.03, 0.065]; %   0.3500    0.2878    0.3415    0.2900    0.3597
        delta_HMM_adapt.t = [0.27, 0.011, 0.0065, 0.03, 0.065]; 
    end    
end

%% FULL DA - RANDOM WALK UPDATES 
if arg1 
    h = h_init;
    theta = theta_init;

    H_subset_DA_RW = zeros(M,length(ind_h_sel));
    mean_H_DA_RW = zeros(1,T);
    s_H_DA_RW  = zeros(1,T);
    theta_DA_RW = zeros(M,5);

    accept_DA_RW = zeros(M,1);
    A_H_DA_RW = zeros(M,1);
    A_theta_DA_RW = zeros(M,5);


    fprintf('Full DA, RW updates\n')
    tic
    for ii = -BurnIn:1:M 
        if (mod(ii,1000) == 0)
            tt = toc;
            fprintf('Iter = %i, Elapsed time = %6.4f s \n', ii, tt);
        end
        [h, acc, A_h] = update_h_leverage(y, h, theta, delta.h);
        [theta, A_theta] = update_theta_leverage(y, h, theta, hyper, delta.t);
        if (ii > 0)
            theta_DA_RW(ii,:) = theta;
            H_subset_DA_RW(ii,:) = h(ind_h_sel);

[mean_H_DA_RW, s_H_DA_RW, var_H_DA_RW] = h_seq_stats(h, mean_H_DA_RW, s_H_DA_RW, ii);                
 

            accept_DA_RW(ii,1) = acc;
            A_H_DA_RW(ii,1) = A_h;
            A_theta_DA_RW(ii,:) = A_theta;        
        end
    end
    time_DA_RW = toc;
    accept_DA_RW = accept_DA_RW/T;
    A_H_DA_RW = A_H_DA_RW/T;
    
mean(A_theta_DA_RW)
mean(A_H_DA_RW)
mean(theta_DA_RW)

    ESS_theta_DA_RW = ESS(theta_DA_RW,0); 
    ESS_H_DA_RW = ESS(H_subset_DA_RW,0); 

    if save_on
        name = [path,'SVML_results_param_DA_RW',data_name,'.mat'];
        save(name,...
            'y','h_init','theta_init','delta','hyper',...
            'mean_H_DA_RW','var_H_DA_RW',...
            'theta_DA_RW',...
            'accept_DA_RW',...
            'A_H_DA_RW',...
            'A_theta_DA_RW',...
            'H_subset_DA_RW',...
            'ESS_theta_DA_RW',...
            'ESS_H_DA_RW',...       
            'time_DA_RW');
    end
    
end

%% SEMI DA efficient implementation:
if (arg2 > 0)        % integrate out the odd h(t)'s and impute the even ones
    if (arg2 > 1)
        N_bin = arg2;
        bins = linspace(-bin_range,bin_range,N_bin+1);
        bin_midpoint = (bins(1:N_bin) + bins(2:N_bin+1))/2;
    end        

    h = h_init;
    theta = theta_init; 

    H_subset_HMM = zeros(M,length(ind_h_sel));       
    mean_H_HMM = zeros(1,T);
    s_H_HMM  = zeros(1,T);

    theta_HMM  = zeros(M,5);

    accept_HMM = zeros(M,1);
    A_H_HMM  = zeros(M,1);
    A_theta_HMM  = zeros(M,5);

    loglik = loglik_h_leverage_HMM_eff(y, h, theta, bins, bin_midpoint); 
    % this is the loglik for the integrals ONLY 
    % so the denominator for the theta updates (up to the prior terms)
    % does NOT include the observations logliks for the imputed h's 

    fprintf('Semi DA, even states imputed \n')

    tic
    for ii = -BurnIn:1:M 
        if (mod(ii,1000) == 0)
            tt = toc;
            fprintf('Iter = %i, Elapsed time = %6.4f s \n', ii, tt);
        end
        [h, acc, A_h, loglik] = update_h_leverage_HMM_eff(y, h, theta, delta_HMM.h, bins, bin_midpoint, loglik);
        [theta, A_theta, loglik] = update_theta_leverage_HMM_RW_eff(y, h, theta, bins, bin_midpoint, hyper, delta_HMM.t, loglik);

        if (ii > 0)
            theta_HMM (ii,:) = theta;
            H_subset_HMM(ii,:) = h(ind_h_sel);

[mean_H_HMM, s_H_HMM, var_H_HMM] = h_seq_stats(h, mean_H_HMM, s_H_HMM, ii);                


            accept_HMM(ii,1) = acc;
            A_H_HMM(ii,1) = A_h;
            A_theta_HMM(ii,:) = A_theta;        
        end
    end
    time_HMM = toc;

    accept_HMM = accept_HMM/(T/2);
    A_H_HMM = A_H_HMM/(T/2); 

    ESS_theta_HMM = ESS(theta_HMM,0);
    ESS_H_HMM = ESS(H_subset_HMM,0);

    if save_on
        name = [path,'SVML_results_param_HMM_Nbin',num2str(N_bin),data_name,'.mat'];
        save(name,...
            'y','h_init','theta_init','delta_HMM','hyper',...
            'H_subset_HMM',...
            'mean_H_HMM','var_H_HMM' ,...
            'theta_HMM',...
            'accept_HMM',...
            'A_H_HMM',...
            'A_theta_HMM',... 
            'ESS_theta_HMM',...
            'ESS_H_HMM',...    
            'time_HMM');
    end
end
   
%% ADAPTIVE
if (arg3 > 0)        % integrate out the odd h(t)'s and impute the even ones  
    h = h_init;
    theta = theta_init; 

% h = mean_H_DA_RW;
% theta = [0.37 0.95 0.018, 0.0, -0.36];

    H_subset_HMM_adapt = zeros(M,length(ind_h_sel));       
    mean_H_HMM_adapt = zeros(1,T);
    s_H_HMM_adapt  = zeros(1,T);

    theta_HMM_adapt  = zeros(M,5);

    accept_HMM_adapt = zeros(M,1);
    A_H_HMM_adapt  = zeros(M,1);
    A_theta_HMM_adapt  = zeros(M,5);

    loglik = loglik_h_leverage_HMM_adapt_eff(y, h, theta, mid_inv); 
    % this is the loglik for the integrals ONLY 
    % so the denominator for the theta updates (up to the prior terms)
    % does NOT include the observations logliks for the imputed h's 

    fprintf('Semi DA Adapted, even states imputed \n')

    tic
    for ii = -BurnIn:1:M 
        if (mod(ii,1000) == 0)
            tt = toc;
            fprintf('Iter = %i, Elapsed time = %6.4f s \n', ii, tt);
        end
        [theta, A_theta, loglik] = update_theta_leverage_HMM_adapt_eff(y, h, theta, mid_inv, hyper, delta_HMM_adapt.t, loglik);
        [h, acc, A_h, loglik] = update_h_leverage_HMM_adapt_eff_v2(y, h, theta, delta_HMM_adapt.h, mid_inv, loglik);
        if (ii > 0)
            theta_HMM_adapt(ii,:) = theta;
            H_subset_HMM_adapt(ii,:) = h(ind_h_sel);

[mean_H_HMM_adapt, s_H_HMM_adapt, var_H_HMM_adapt] = h_seq_stats(h, mean_H_HMM_adapt, s_H_HMM_adapt, ii);                

            accept_HMM_adapt(ii,1) = acc;
            A_H_HMM_adapt(ii,1) = A_h;
            A_theta_HMM_adapt(ii,:) = A_theta;        
        end
    end
    time_HMM_adapt = toc;

    accept_HMM_adapt = accept_HMM_adapt/(T/2);
    A_H_HMM_adapt = A_H_HMM_adapt/(T/2); 

mean(A_theta_HMM_adapt)
mean(A_H_HMM_adapt)
mean(theta_HMM_adapt)


    ESS_theta_HMM_adapt = ESS(theta_HMM_adapt,0);
    ESS_H_HMM_adapt = ESS(H_subset_HMM_adapt,0);    
    
    if save_on
        name = [path,'SVML_results_param_HMM_adapt_Nq',num2str(N_q),data_name,'.mat'];
        save(name,...
            'y','h_init','theta_init','delta_HMM_adapt','hyper',...
            'H_subset_HMM_adapt',...
            'mean_H_HMM_adapt','var_H_HMM_adapt' ,...
            'theta_HMM_adapt',...
            'accept_HMM_adapt',...
            'A_H_HMM_adapt',...
            'A_theta_HMM_adapt',... 
            'ESS_theta_HMM_adapt',...
            'ESS_H_HMM_adapt',...    
            'time_HMM_adapt');      
    end
end
    