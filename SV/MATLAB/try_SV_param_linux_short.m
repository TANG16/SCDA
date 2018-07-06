function try_SV_param_linux_short(arg0,~,arg2,~,~,arg5,~,arg7)
% arg0 = 2; arg2 = 2;
addpath(genpath('../MATLAB'));

% try_SV_param_linux_short(1,~,2,~,~,1,~,1)
% clear all
close all

% arg0 = 2;

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

plot_on = false;

if (arg0 == 0) % simulate   
 %%
    beta = 0.5;
    mu = 2*log(beta);
     
%     theta = [-1, 0.95, 0.15^2];
    theta = [mu, 0.98, 0.2^2];
    theta_true = theta;
    params = {'\mu','\phi','\sigma^2'};

    T = 2000;
    [y,h_true] = generate_SV(theta,T);

    if plot_on
        plot(y)
        plot(h_true)
    end

    mu = theta(1);
    phi = theta(2);
    sigma2 = theta(3);
    data_name = '';
    path = 'Results/Simulation/';
    
else % empirical
    if false
        %%
        y = load('../Data/Perc_Rets_GSPC_IBM_AAPL_MSFT_JPM_GE.csv');
        y = y(:,arg0); % arg0 = 4 = MSFT
        time = [1998,(2017 + 8/12)];       
        TT = length(y);
        T = 2000;
        y = y((TT-T+1):TT)';
        DATA_NAMES = {'_GSPC','_IBM','_AAPL','_MSFT','_JPM','_GE'};
        data_name = [DATA_NAMES{arg0},'_old']; %'_MSFT';  
%         data_name = DATA_NAMES{arg0}; %'_MSFT';  
        path = 'Results/Empirical/';
        fprintf('Data loaded: %s\n',data_name);
        theta_true = [];
        h_true = [];
    else
        y = load('other/Perc_Rets_GSPC_IBM_MSFT.csv');
        y = y(:,arg0); % arg0 = 1 = GSPC;  arg0 = 2 = IBM
        time = [2000,2018];       
        T = length(y);
%         T = 2000;
%         y = y((TT-T+1):TT)';
        DATA_NAMES = {'_GSPC','_IBM','_MSFT'};
        data_name = DATA_NAMES{arg0}; %'_GSPC';  
%         path = 'Results/Empirical/';
        path = 'Results/Empirical_new/';
        fprintf('Data loaded: %s\n',data_name);
        theta_true = [];
        h_true = [];        
    end
end

hyper.S = [5/2, 0.05/2];
hyper.M = 10;
hyper.P = [20, 3/2];

% P0 = sqrt(sigma2/(1-phi^2)); % unconditional st. dev. 0.6120
% bin_range = 5*P0;
bin_range = 4;
N_bin = 30; %10; %25; %20; %30; %10; %15; %20; %25; %30; %15; % 20;%30; %50;
bins = linspace(-bin_range,bin_range,N_bin+1);
bin_midpoint = (bins(1:N_bin) + bins(2:N_bin+1))/2;
mid = (0:(N_bin-1))/N_bin;
mid = mid + mid(2)/2;
mid_inv = norminv(mid);
% bins are the demeaned volatilities: b=h-mu
 
M = 5*10000;
BurnIn = 10000; % + 3000 + 6000;

% q_L = 2.5;
% q_U = 97.5;
% q_no = M*q_L/100;
% H_tail_L = zeros(q_no,T);
% H_tail_U = zeros(q_no,T);
fprintf('*** M = %i and BurnIn = %i *** \n', M, BurnIn);

h_init = log(var(y))*ones(1,T); 
 
%% FULL DA - RANDOM WALK UPDATES
    if  (arg2 == 2)
%         [theta, delta] = SV_initialise(arg0, 2);
if false
    if (arg0 == 1)
        delta.t = [0.5, 0.005, 0.0015]; % mean(A_theta_DA_RW_eff) =   0.3441    0.3903    0.3417
        delta.h = 0.4; % mean(A_H_DA_RW_eff) =     0.3346
    elseif (arg0 == 2)
        delta.t = [0.15, 0.03, 0.02]; % mean(A_theta_DA_RW_eff) = 0.3378    0.3583    0.3149
        delta.h = 1.1;  %  mean(A_H_DA_RW_eff) = 0.3807
    end
else
    if (arg0 == 1)
        delta.t = [0.5, 0.005, 0.0015]; % mean(A_theta_DA_RW_eff) =   0.3441    0.3903    0.3417
        delta.h = 0.4; % mean(A_H_DA_RW_eff) =     0.3346
    elseif (arg0 == 2)
        delta.t = [0.3, 0.01, 0.005]; % mean(A_theta_DA_RW_eff) = 0.3378    0.3583    0.3149
        delta.h = 0.6;  %  mean(A_H_DA_RW_eff) = 0.3807
    end
end
h = h_init;
theta = [0 0.96 0.08];
ind_h_sel = 2:50:T;

        H_subset_DA_RW_eff = zeros(M,length(ind_h_sel));
        mean_H_DA_RW_eff = zeros(1,T);
        s_H_DA_RW_eff  = zeros(1,T);
        theta_DA_RW_eff = zeros(M,3);

        accept_DA_RW_eff = zeros(M,1);
        A_H_DA_RW_eff = zeros(M,1);
        A_theta_DA_RW_eff = zeros(M,3);
 
        loglik = loglik_h_eff(h, theta);
 
        fprintf('Full DA, RW_eff updates\n')
        tic
        for ii = -BurnIn:1:M 
             if (mod(ii,1000) == 1)
                toc;
            end
            [h, acc, A_h, loglik] = update_h_eff(y, h, theta, delta.h, loglik);
            [theta, A_theta, loglik] = update_theta_RW_eff(h, theta, hyper, delta.t, loglik);
            if (ii > 0)
                theta_DA_RW_eff(ii,:) = theta;
                H_subset_DA_RW_eff(ii,:) = h(ind_h_sel);

[mean_H_DA_RW_eff, s_H_DA_RW_eff, var_H_DA_RW_eff] = h_seq_stats(h, mean_H_DA_RW_eff, s_H_DA_RW_eff, ii);                

                accept_DA_RW_eff(ii,1) = acc;
                A_H_DA_RW_eff(ii,1) = A_h;
                A_theta_DA_RW_eff(ii,:) = A_theta;        
            end
        end
        time_DA_RW_eff = toc;

        accept_DA_RW_eff = accept_DA_RW_eff/T;
        A_H_DA_RW_eff = A_H_DA_RW_eff/T;
 
        ESS_theta_DA_RW_eff_sig = ESS(theta_DA_RW_eff,0); 
        ESS_H_DA_RW_eff_sig = ESS(H_subset_DA_RW_eff,0); 

        name = [path,'SV_results_param_DA_RW_eff',data_name,'.mat'];
        save(name,...
            'y','h_true','theta_true','delta','hyper',...
            'mean_H_DA_RW_eff','var_H_DA_RW_eff',...
            'theta_DA_RW_eff',...
            'accept_DA_RW_eff',...
            'A_H_DA_RW_eff',...
            'A_theta_DA_RW_eff',...
            'H_subset_DA_RW_eff',...
            'ESS_theta_DA_RW_eff_sig',...
            'ESS_H_DA_RW_eff_sig',...       
            'time_DA_RW_eff');
    end
 
%% SEMI DA efficient implementation:
    if (arg5 > 0)        % integrate out the odd h(t)'s and impute the even ones
        if (arg5 > 1)
            N_bin = arg5;
            bins = linspace(-bin_range,bin_range,N_bin+1);
            bin_midpoint = (bins(1:N_bin) + bins(2:N_bin+1))/2;
        end        
%         [theta, delta] = SV_initialise(arg0, 5);
if false
    if (arg0 == 0)
        delta.t = [0.25, 0.02, 0.008];  %            0.4275    0.3147    0.4102
        delta.h = 0.8;        
    else
        delta.t = [0.17, 0.04, 0.03];   %     0.3420    0.2995    0.3219
        delta.h = 1.7;         
    end
else
    if (arg0 == 1)
        if (N_bin == 30)
            delta.t = [0.1 0.006 0.0019]; % mean(A_theta_HMM_eff) =   0.3880    0.3416    0.3912
            delta.h = 0.52; % mean(A_H_HMM_eff) = 0.3679    
        elseif (N_bin == 25)
            delta.t = [0.1 0.006 0.0019]; 
            delta.h = 0.52;     
        else
            delta.t = [0.1 0.006 0.0019];   
            delta.h = 0.52;           
        end
    elseif (arg0 == 2)
        delta.h = 1.0; % mean(A_H_HMM_eff)
        if (N_bin == 30)
            delta.t = [0.3, 0.01, 0.005]; % mean(A_theta_HMM_eff) = 0.3378    0.3583    0.3149
        elseif (N_bin == 25)
            delta.t = [0.3, 0.01, 0.007]; 
        elseif (N_bin == 20)
            delta.t = [0.25, 0.01, 0.007];  
        elseif (N_bin <= 15)
            delta.t = [0.05, 0.051, 0.001]; 
            delta.h = 0.5; % mean(A_H_HMM_eff)
        end
    end
end
h = h_init;
theta = [0 0.96 0.08];

ind_h_sel = 2:50:T;

        H_subset_HMM_eff = zeros(M,length(ind_h_sel));       
        mean_H_HMM_eff = zeros(1,T);
        s_H_HMM_eff  = zeros(1,T);
        
        theta_HMM_eff  = zeros(M,3);

        accept_HMM_eff = zeros(M,1);
        A_H_HMM_eff  = zeros(M,1);
        A_theta_HMM_eff  = zeros(M,3);

        loglik = loglik_h_HMM_eff(y, h, theta, bins, bin_midpoint); 
        % this is the loglik for the integrals ONLY 
        % so the denominator for the theta updates (up to the prior terms)
        % does NOT include the observations logliks for the imputed h's 
                                                                   
        fprintf('Semi DA, even states imputed efficient implementation\n')
        
        tic
        for ii = -BurnIn:1:M 
            if (mod(ii,1000) == 1)
                toc;
            end
            [h, acc, A_h, loglik] = update_h_HMM_eff(y, h, theta, delta.h, bins, bin_midpoint, loglik);
            [theta, A_theta, loglik] = update_theta_HMM_RW_eff(y, h, theta, bins, bin_midpoint, hyper, delta.t, loglik);

            if (ii > 0)
                theta_HMM_eff (ii,:) = theta;
                H_subset_HMM_eff(ii,:) = h(ind_h_sel);

[mean_H_HMM_eff, s_H_HMM_eff, var_H_HMM_eff] = h_seq_stats(h, mean_H_HMM_eff, s_H_HMM_eff, ii);                

% if (ii <= q_no)
%     H_tail_L(ii,:) = h;
%     H_tail_U(ii,:) = h;
% else
%     [H_tail_L, H_tail_U] = update_tail(h,H_tail_L,H_tail_U);
%  end
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
          
        name = [path,'SV_results_param_HMM_eff_Nbin',num2str(N_bin),data_name,'.mat'];
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
 
%% SEMI DA ADAPTIVE EFFICIENT
    if (arg7 > 0)
        if (arg7 > 1)
            N_bin = arg7;
            mid = (0:(N_bin-1))/N_bin;
            mid = mid + mid(2)/2;
            mid_inv = norminv(mid);
        end
        
if (arg0 == 1)
    delta.t = [0.4 0.006 0.002]; % mean(A_theta_HMM_adapt_eff) =  0.3389    0.3451    0.3778
    delta.h = 0.54; % mean(A_H_HMM_adapt_eff) =  0.3658
elseif (arg0 == 2)
    delta.t = [0.25, 0.01, 0.006]; % mean(A_theta_HMM_adapt_eff) =  0.3751    0.3663    0.3090
    delta.h = 1.0; % mean(A_H_HMM_adapt_eff) = 0.3253
end
    h = h_init; 
    theta = [0 0.96 0.08];
    
    ind_h_sel = 2:50:T;

    % [theta, delta] = SV_initialise(arg0, 6);
        H_subset_HMM_adapt_eff = zeros(M,length(ind_h_sel));       
        mean_H_HMM_adapt_eff = zeros(1,T);
        s_H_HMM_adapt_eff  = zeros(1,T);

        theta_HMM_adapt_eff = zeros(M,3);

        accept_HMM_adapt_eff = zeros(M,1);
        A_H_HMM_adapt_eff = zeros(M,1);
        A_theta_HMM_adapt_eff = zeros(M,3);

        loglik = loglik_h_HMM_adapt_eff(y, h, theta, mid_inv); 

        fprintf('Semi DA adaptive efficient implementation\n')
        tic
        for ii = -BurnIn:1:M  
            if (mod(ii,1000) == 1)
                toc;
            end
            [h, acc, A_h, loglik] = update_h_HMM_adapt_eff(y, h, theta, delta.h, mid_inv, loglik);
            [theta, A_theta, loglik] = update_theta_HMM_RW_adapt_eff(y, h, theta, mid_inv, hyper, delta.t, loglik);

            if (ii > 0)
                theta_HMM_adapt_eff(ii,:) = theta;
                H_subset_HMM_adapt_eff(ii,:) = h(ind_h_sel);

[mean_H_HMM_adapt_eff, s_H_HMM_adapt_eff, var_H_HMM_adapt_eff] = ...
    h_seq_stats(h, mean_H_HMM_adapt_eff, s_H_HMM_adapt_eff, ii);   

                accept_HMM_adapt_eff(ii,1) = acc;
                A_H_HMM_adapt_eff(ii,1) = A_h;
                A_theta_HMM_adapt_eff(ii,:) = A_theta;        
            end
        end
        time_HMM_adapt_eff = toc;
        accept_HMM_adapt_eff = accept_HMM_adapt_eff/(T/2);
        A_H_HMM_adapt_eff = A_H_HMM_adapt_eff/(T/2);
            
        ESS_theta_HMM_adapt_eff_sig = ESS(theta_HMM_adapt_eff,0);     
        ESS_H_HMM_adapt_eff_sig = ESS(H_subset_HMM_adapt_eff,0);
        
%         name = [path,'SV_results_param_HMM_adapt_eff_Nbin',num2str(N_bin),data_name,'.mat'];
        name = [path,'SV_results_param_HMM_adapt_eff_Nbin',num2str(N_bin),data_name,'.mat'];
        save(name,...
            'y','h_true','theta_true','delta','hyper',...
            'H_subset_HMM_adapt_eff',...
            'mean_H_HMM_adapt_eff','var_H_HMM_adapt_eff',...
            'theta_HMM_adapt_eff',...
            'accept_HMM_adapt_eff',...
            'A_H_HMM_adapt_eff',...
            'A_theta_HMM_adapt_eff',...
            'ESS_theta_HMM_adapt_eff_sig',...
            'ESS_H_HMM_adapt_eff_sig',...    
            'time_HMM_adapt_eff');   
    end    
end