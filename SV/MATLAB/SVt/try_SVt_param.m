% clear all
close all
addpath(genpath('../../MATLAB'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

plot_on = false;

% arg0 = 2;
% arg1 = 0;
% arg2 = 1;
path = '';

if (arg0 == 0)
    beta = 0.5;
    mu = 2*log(beta);
    theta = [mu, 0.98, 0.2^2, 7];
    theta_true = theta;
    T = 2000;
    [y,h_true] = generate_SVt(theta,T);
    data_name = '_sim';
else
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
end


theta_init = [0, 0.95, 0.1^2, 7]; 


hyper.S = [5/2, 0.05/2];
hyper.M = 10;
hyper.P = [20, 3/2];
hyper.N = 0.01;

 
bin_range = 4;
N_bin = 30; 
bins = linspace(-bin_range,bin_range,N_bin+1);
bin_midpoint = (bins(1:N_bin) + bins(2:N_bin+1))/2;
% bins are the demeaned volatilities: b=h-mu
 

M = 10000;
BurnIn = 10000;
h_init = log(var(y))*ones(1,T); 
 
%% FULL DA - RANDOM WALK UPDATES 
if arg1 
    if (arg0 == 0)
        delta.h = 0.1;
        delta.t = [0.6, 0.005, 0.0007, 1.0]; % 
    elseif (arg0 == 2)
        delta.h = 0.5;
        delta.t = [0.6, 0.005, 0.0007, 1.0]; %     
    end
    h = h_init;
    theta = theta_init;
    ind_h_sel = 2:50:T;

    H_subset_DA_RW = zeros(M,length(ind_h_sel));
    mean_H_DA_RW = zeros(1,T);
    s_H_DA_RW  = zeros(1,T);
    theta_DA_RW = zeros(M,4);

    accept_DA_RW = zeros(M,1);
    A_H_DA_RW = zeros(M,1);
    A_theta_DA_RW = zeros(M,4);
 

    fprintf('Full DA, RW updates\n')
    tic
    for ii = -BurnIn:1:M 
        if (mod(ii,1000) == 0)
            toc;
        end
        [h, acc, A_h] = update_h(y, h, theta, delta.h);
        [theta, A_theta] = update_theta_RW(y, h, theta, hyper, delta.t);
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

    ESS_theta_DA_RW_sig = ESS(theta_DA_RW,0); 
    ESS_H_DA_RW_sig = ESS(H_subset_DA_RW,0); 

    name = [path,'SVt_results_param_DA_RW',data_name,'.mat'];
    save(name,...
        'y','h_true','theta_true','delta','hyper',...
        'mean_H_DA_RW','var_H_DA_RW',...
        'theta_DA_RW',...
        'accept_DA_RW',...
        'A_H_DA_RW',...
        'A_theta_DA_RW',...
        'H_subset_DA_RW',...
        'ESS_theta_DA_RW_sig',...
        'ESS_H_DA_RW_sig',...       
        'time_DA_RW');

    
end

%% SEMI DA efficient implementation:
if (arg2 > 0)        % integrate out the odd h(t)'s and impute the even ones
        if (arg2 > 1)
            N_bin = arg2;
            bins = linspace(-bin_range,bin_range,N_bin+1);
            bin_midpoint = (bins(1:N_bin) + bins(2:N_bin+1))/2;
        end        

        
if (arg0 == 0)
    delta.h = 0.2;
    delta.t = [0.05, 0.005, 0.001, 0.5]; %

    h = h_init;
%     theta_init = [0, 0.95, 0.1^2, 7]; 
    theta_init = [-1, 0.95, 0.5, 6]; 

    theta = theta_init;
    % theta(4) = theta_true(4);

%     h = h_true;
%     theta = theta_true;
elseif (arg0 == 2)
    delta.h = 0.1;
    delta.t = [0.1 0.005 0.001 0.1];
    h = h_init;
    theta_init = [-1.3803    0.9817    0.0301    6.3280];
    theta = theta_init;
end


ind_h_sel = 2:50:T;

        H_subset_HMM_eff = zeros(M,length(ind_h_sel));       
        mean_H_HMM_eff = zeros(1,T);
        s_H_HMM_eff  = zeros(1,T);
        
        theta_HMM_eff  = zeros(M,4);

        accept_HMM_eff = zeros(M,1);
        A_H_HMM_eff  = zeros(M,1);
        A_theta_HMM_eff  = zeros(M,4);

        loglik = loglik_h_HMM_eff(y, h, theta, bins, bin_midpoint); 
        % this is the loglik for the integrals ONLY 
        % so the denominator for the theta updates (up to the prior terms)
        % does NOT include the observations logliks for the imputed h's 
                                                                   
        fprintf('Semi DA, even states imputed efficient implementation\n')
        
        tic
        for ii = -BurnIn:1:M 
%         for ii = 1:M 
            if (mod(ii,1000) == 1)
                tt = toc;
                fprintf('Iter = %i, Elapsed time = %6.4f s \n', ii, tt);
            end
%             [h, acc, A_h, loglik] = update_h_HMM_eff(y, h, theta, delta.h, bins, bin_midpoint, loglik);
%             [theta, A_theta, loglik] = update_theta_HMM_RW_eff(y, h, theta, bins, bin_midpoint, hyper, delta.t, loglik);
            [h, acc, A_h] = update_h_HMM_v2(y, h, theta, delta.h, bins, bin_midpoint);
            [theta, A_theta] = update_theta_HMM_RW(y, h, theta, bins, bin_midpoint, hyper, delta.t);

            if (ii > 0)
                theta_HMM_eff (ii,:) = theta;
                H_subset_HMM_eff(ii,:) = h(ind_h_sel);

[mean_H_HMM_eff, s_H_HMM_eff, var_H_HMM_eff] = h_seq_stats(h, mean_H_HMM_eff, s_H_HMM_eff, ii);                


                accept_HMM_eff(ii,1) = acc;
                A_H_HMM_eff(ii,1) = A_h;
                A_theta_HMM_eff(ii,:) = A_theta;        
            end
        end
        time_HMM_eff = toc;

        accept_HMM_eff = accept_HMM_eff/(T/2);
        A_H_HMM_eff = A_H_HMM_eff/(T/2); 
        
        ESS_theta_HMM_eff_sig = ESS(theta_HMM_eff,0);
        ESS_H_HMM_eff_sig = ESS(H_subset_HMM_eff,0);
          
        name = [path,'SVt_results_param_HMM_eff_Nbin',num2str(N_bin),data_name,'.mat'];
        save(name,...
            'y','h_true','theta_true','delta','hyper',...
            'H_subset_HMM_eff',...
            'mean_H_HMM_eff','var_H_HMM_eff' ,...
            'theta_HMM_eff',...
            'accept_HMM_eff',...
            'A_H_HMM_eff',...
            'A_theta_HMM_eff',... 
            'ESS_theta_HMM_eff_sig',...
            'ESS_H_HMM_eff_sig',...    
            'time_HMM_eff');
end
 