% clear all
close all
% addpath(genpath('../../MATLAB'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

plot_on = false;

arg0 = 0;
arg1 = 1;
arg2 = 1;
arg3 = 1;
path = '';

if (arg0 == 0)
    beta = -0.2;
    mu = -1;
    theta_true = [mu, 0.98, 0.01, beta];
    T = 2000;
    [y,h_true,y0] = generate_SVM(theta_true,T);
    theta_init = theta_true;
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
    theta_init = [-1, 0.98, 0.01, -0.2]; 
end



hyper.S = [5/2, 0.05/2];
hyper.M = 10;
hyper.P = [20, 3/2];
hyper.B = 10;

 
bin_range = 4;
N_bin = 50; 
bins = linspace(-bin_range,bin_range,N_bin+1);
bin_midpoint = (bins(1:N_bin) + bins(2:N_bin+1))/2;
% bins are the demeaned volatilities: b=h-mu
 N_q = 20;
mid = (0:(N_q-1))/N_q;
mid = mid + mid(2)/2;
mid_inv = norminv(mid);


M = 10000;
BurnIn = -1;

if (arg0 == 0)
    h_init = h_true;
else
    h_init = log(var(y))*ones(1,T); 
end

ind_h_sel = 2:50:T;

if (arg0 == 0)
    delta.h = 0.25;
    delta.t = [0.35, 0.01, 0.0008, 0.065]; % 

    delta_HMM.h = 0.4;
    delta_HMM.t = [0.2, 0.015, 0.0015, 0.085]; %      

    delta_HMM_adapt.h = 0.4;
%     delta_HMM_adapt.t = [0.3, 0.015, 0.00159, 0.10]; %
    delta_HMM_adapt.t = [0.25, 0.009, 0.0009, 0.08]; %      
%     delta_HMM_adapt.t = [0.25, 0.0012, 0.0012, 0.10]; %      

elseif (arg0 == 2)
    delta.h = 0.59;
    delta.t = [0.32, 0.0095, 0.004, 0.025]; % 

    delta_HMM.h = 0.3;
    delta_HMM.t = [0.05, 0.008, 0.0008, 0.07]; %      

    delta_HMM_adapt.h = 0.9;
    delta_HMM_adapt.t = [0.25, 0.01, 0.006, 0.022]; %      
end

%% FULL DA - RANDOM WALK UPDATES 
if arg1 
    h = h_init;
    theta = theta_init;

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
            tt = toc;
            fprintf('Iter = %i, Elapsed time = %6.4f s \n', ii, tt);
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
    
mean(A_theta_DA_RW)
mean(A_H_DA_RW)
mean(theta_DA_RW)



    ESS_theta_DA_RW = ESS(theta_DA_RW,0); 
    ESS_H_DA_RW = ESS(H_subset_DA_RW,0); 

    name = [path,'SVM_results_param_DA_RW',data_name,'.mat'];
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

    theta_HMM  = zeros(M,4);

    accept_HMM = zeros(M,1);
    A_H_HMM  = zeros(M,1);
    A_theta_HMM  = zeros(M,4);

    loglik = loglik_h_HMM_eff(y, h, theta, bins, bin_midpoint); 
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
        [h, acc, A_h, loglik] = update_h_HMM_eff(y, h, theta, delta_HMM.h, bins, bin_midpoint, loglik);
        [theta, A_theta, loglik] = update_theta_HMM_RW_eff(y, h, theta, bins, bin_midpoint, hyper, delta_HMM.t, loglik);

        if (ii > 0)
            theta_HMM(ii,:) = theta;
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

mean(A_theta_HMM)
mean(A_H_HMM)
mean(theta_HMM)
    ESS_theta_HMM = ESS(theta_HMM,0);
    ESS_H_HMM = ESS(H_subset_HMM,0);

    name = [path,'SVM_results_param_HMM_Nbin',num2str(N_bin),data_name,'.mat'];
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
   
%% SCDA ADAPTIVE
if (arg3 > 0)        % integrate out the odd h(t)'s and impute the even ones  
    h = h_init;
    theta = theta_init; 

    H_subset_HMM_adapt = zeros(M,length(ind_h_sel));       
    mean_H_HMM_adapt = zeros(1,T);
    s_H_HMM_adapt  = zeros(1,T);

    theta_HMM_adapt  = zeros(M,4);

    accept_HMM_adapt = zeros(M,1);
    A_H_HMM_adapt  = zeros(M,1);
    A_theta_HMM_adapt  = zeros(M,4);

    loglik = loglik_h_HMM_adapt_eff(y, h, theta, mid_inv); 
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
        [theta, A_theta, loglik] = update_theta_HMM_RW_adapt_eff(y, h, theta, mid_inv, hyper, delta_HMM_adapt.t, loglik);
        [h, acc, A_h, loglik] = update_h_HMM_adapt_eff(y, h, theta, delta_HMM_adapt.h, mid_inv, loglik);
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
    name = [path,'SVM_results_param_HMM_adapt_Nq',num2str(N_q),data_name,'.mat'];
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
    