function try_SV_param_linux(arg1,arg2,arg3,arg4,arg5, arg6)
% clear all
close all

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

plot_on = false;

% theta = [1, 0.97, 0.15^2];
% beta = 0.05;
beta = 0.5;
mu = 2*log(beta);
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

theta_init = [0, 0.95, 0.1^2];
theta_init2 = [0, 0.95, 0.16^2];
theta_init3 = [0, 0.95, 0.4^2];

kernel = @(xx) -loglik_h(h_true, xx)/T;
theta_mle = fminunc(kernel,theta_init);

hyper.S = [5/2, 0.05/2];
hyper.M = 10;
hyper.P = [20, 3/2];

% P0 = sqrt(sigma2/(1-phi^2)); % unconditional st. dev. 0.6120
% bin_range = 5*P0;
bin_range = 4;
N_bin = 30; %15; % 20;%30; %50;
bins = linspace(-bin_range,bin_range,N_bin+1);
bin_midpoint = (bins(1:N_bin) + bins(2:N_bin+1))/2;
mid = (0:(N_bin-1))/N_bin;
mid = mid + mid(2)/2;
% bins are the demeaned volatilities: b=h-mu

delta.h = 0.1;
delta_h = delta.h;
% delta.t = [0.4, 0.01, 0.002]; % mean(A_theta_DA_RW):     0.5390    0.3013    0.4061
                              % mean(A_theta_HMM):    0.3913    0.3599    0.5981
% delta.t = [0.4, 0.01, 0.003]; % mean(A_theta_HMM):  0.3931    0.3646    0.4835
delta.t = [0.4, 0.01, 0.004]; % mean(A_theta_HMM): 0.3929    0.3685    0.4065

M = 10000;
BurnIn = 1000;
h_init = log(var(y))*ones(1,T); 


%% FULL DA - FULL CONDITIONAL DENSITIES FROM KIM ET AL. (1998)
    if arg1
         
% %         delta.h = 0.1;
% % %         mean(A_H_DA) = 0.6089
%         delta.h = 0.5;
%         mean(A_H_DA) = 0.0621
        delta.h = 0.2;       
%         mean(A_H_DA) = 0.3734
        
        h = h_init;
        theta = theta_init;
        H_DA = zeros(M,T);
        theta_DA = zeros(M,3);

        accept_DA = zeros(M,1);
        A_H_DA = zeros(M,1);
        A_phi_DA = zeros(M,1);

        fprintf('Full DA, full conditional updates\n')
        tic
        for ii = -BurnIn:1:M 
            if (mod(ii,1000) == 1)
                toc;
            end
            [h, acc, A_h] = update_h(y, h, theta, delta.h);
            [theta, A_phi] = update_theta_cond(h, theta, hyper);
            if (ii > 0)
                theta_DA(ii,:) = theta;
                H_DA(ii,:) = h;
                accept_DA(ii,1) = acc;
                A_H_DA(ii,1) = A_h;
                A_phi_DA(ii,1) = A_phi;        
            end
        end
        time_DA = toc;

        accept_DA = accept_DA/T;
        A_H_DA = A_H_DA/T;
        mean_H_DA = mean(H_DA);

        ESS_theta_DA_40 = ESS(theta_DA,40);
        ESS_theta_DA_100 = ESS(theta_DA,100);
        ESS_theta_DA_1000 = ESS(theta_DA,1000);
        ESS_theta_DA_sig = ESS(theta_DA,0);

        ESS_H_DA_40 = ESS(H_DA,40);
        ESS_H_DA_100 = ESS(H_DA,100);
        ESS_H_DA_1000 = ESS(H_DA,1000);
        ESS_H_DA_sig = ESS(H_DA,0);
        
        H_subset_DA = H_DA(:,(2:10)*100);        
        
        name = 'SV_results_param_DA.mat';
        save(name,...
            'y','h_true','theta_true','delta','hyper',...
            'mean_H_DA',...
            'theta_DA',...
            'accept_DA',...
            'A_H_DA',...
            'A_phi_DA',...
            'H_subset_DA',...
            'ESS_theta_DA_40','ESS_theta_DA_100','ESS_theta_DA_1000','ESS_theta_DA_sig',...
            'ESS_H_DA_40','ESS_H_DA_100','ESS_H_DA_1000','ESS_H_DA_sig',...       
            'time_DA');    
    end

%% FULL DA - RANDOM WALK UPDATES
    if (arg2==1)
        % delta.t = [0.1, 0.01, 0.01]; mean(A_theta_DA_RW):  0.8422    0.3067    0.1098
        % delta.t = [0.3, 0.01, 0.005]; mean(A_theta_DA_RW):  0.5300    0.3449    0.1554
        % delta.t = [0.4, 0.01, 0.002]; % mean(A_theta_DA_RW):     0.5390    0.3013    0.4061

% % %         delta.h = 0.1;
% % % %         mean(A_H_DA_RW) =  0.6529
% %         delta.h = 0.5;
% % %         mean(A_H_DA_RW) =  0.0623       
%         delta.h = 0.2;
% %         mean(A_H_DA_RW) = 0.2614  
        delta.h = 0.15;
        
        h = h_init;
        theta = theta_init2;
        H_DA_RW = zeros(M,T);
        theta_DA_RW = zeros(M,3);

        accept_DA_RW = zeros(M,1);
        A_H_DA_RW = zeros(M,1);
        A_theta_DA_RW = zeros(M,3);

        fprintf('Full DA, RW updates\n')
        tic
        for ii = -BurnIn:1:M 
            if (mod(ii,1000) == 1)
                toc;
            end
            [h, acc, A_h] = update_h(y, h, theta, delta.h);
            [theta, A_theta] = update_theta_RW(h, theta, hyper, delta.t);
            if (ii > 0)
                theta_DA_RW(ii,:) = theta;
                H_DA_RW(ii,:) = h;

                accept_DA_RW(ii,1) = acc;
                A_H_DA_RW(ii,1) = A_h;
                A_theta_DA_RW(ii,:) = A_theta;        
            end
        end
        time_DA_RW = toc;

        accept_DA_RW = accept_DA_RW/T;
        A_H_DA_RW = A_H_DA_RW/T;
        mean_H_DA_RW = mean(H_DA_RW);

        ESS_theta_DA_RW_40 = ESS(theta_DA_RW,40);
        ESS_theta_DA_RW_100 = ESS(theta_DA_RW,100);
        ESS_theta_DA_RW_1000 = ESS(theta_DA_RW,1000);
        ESS_theta_DA_RW_sig = ESS(theta_DA_RW,0);
        
        ESS_H_DA_RW_40 = ESS(H_DA_RW,40);
        ESS_H_DA_RW_100 = ESS(H_DA_RW,100);
        ESS_H_DA_RW_1000 = ESS(H_DA_RW,1000);
        ESS_H_DA_RW_sig = ESS(H_DA_RW,0);
        
        H_subset_DA_RW = H_DA_RW(:,(2:10)*100);        

        name = 'SV_results_param_DA_RW.mat';
        save(name,...
            'y','h_true','theta_true','delta','hyper',...
            'mean_H_DA_RW',...
            'theta_DA_RW',...
            'accept_DA_RW',...
            'A_H_DA_RW',...
            'A_theta_DA_RW',...
            'H_subset_DA_RW',...
            'ESS_theta_DA_RW_40','ESS_theta_DA_RW_100','ESS_theta_DA_RW_1000','ESS_theta_DA_RW_sig',...
            'ESS_H_DA_RW_40','ESS_H_DA_RW_100','ESS_H_DA_RW_1000','ESS_H_DA_RW_sig',...       
            'time_DA_RW');
        
    elseif (arg2 == 2)
        delta.h = 0.15;
 
        h = h_init;
        theta = theta_init2;
        H_DA_RW_eff = zeros(M,T);
        theta_DA_RW_eff = zeros(M,3);

        accept_DA_RW_eff = zeros(M,1);
        A_H_DA_RW_eff = zeros(M,1);
        A_theta_DA_RW_eff = zeros(M,3);

%         h = h_true;
%         theta = theta_true;
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
                H_DA_RW_eff(ii,:) = h;

                accept_DA_RW_eff(ii,1) = acc;
                A_H_DA_RW_eff(ii,1) = A_h;
                A_theta_DA_RW_eff(ii,:) = A_theta;        
            end
        end
        time_DA_RW_eff = toc;

        accept_DA_RW_eff = accept_DA_RW_eff/T;
        A_H_DA_RW_eff = A_H_DA_RW_eff/T;
        mean_H_DA_RW_eff = mean(H_DA_RW_eff);

        ESS_theta_DA_RW_eff_40 = ESS(theta_DA_RW_eff,40);
        ESS_theta_DA_RW_eff_100 = ESS(theta_DA_RW_eff,100);
        ESS_theta_DA_RW_eff_1000 = ESS(theta_DA_RW_eff,1000);
        ESS_theta_DA_RW_eff_sig = ESS(theta_DA_RW_eff,0);

        ESS_H_DA_RW_eff_40 = ESS(H_DA_RW_eff,40);
        ESS_H_DA_RW_eff_100 = ESS(H_DA_RW_eff,100);
        ESS_H_DA_RW_eff_1000 = ESS(H_DA_RW_eff,1000);
        ESS_H_DA_RW_eff_sig = ESS(H_DA_RW_eff,0);
        
        H_subset_DA_RW_eff = H_DA_RW_eff(:,(2:10)*100);        

        name = 'SV_results_param_DA_RW_eff.mat';
        save(name,...
            'y','h_true','theta_true','delta','hyper',...
            'mean_H_DA_RW_eff',...
            'theta_DA_RW_eff',...
            'accept_DA_RW_eff',...
            'A_H_DA_RW_eff',...
            'A_theta_DA_RW_eff',...
            'H_subset_DA_RW_eff',...
            'ESS_theta_DA_RW_eff_40','ESS_theta_DA_RW_eff_100','ESS_theta_DA_RW_eff_1000','ESS_theta_DA_RW_eff_sig',...
            'ESS_H_DA_RW_eff_40','ESS_H_DA_RW_eff_100','ESS_H_DA_RW_eff_1000','ESS_H_DA_RW_eff_sig',...       
            'time_DA_RW_eff');
    end

%% SEMI DA:
    if arg3
        % integrate out the odd h(t)'s and impute the even ones
    % % % % % % % % % % % % %     delta.t = [0.4000 0.0100 0.0040]
    % % % % % % % % % % % % %     delta.h = 0.1
    % % % % % % % Starting from the true values
    % % % % % % % % % % % %     delta.t = [0.01, 0.025, 0.001]; 
    % % % % % % % % % % % %     delta.h = 0.01;
    % % % % % % % % % % % % %     mean(A_theta_HMM) = 0.9546    0.1611    0.7370
    % % % % % % % % % % % % %     mean(A_H_HMM) = 0.9839
    % % % % % % % % % % %     delta.t = [0.04, 0.01, 0.001]; 
    % % % % % % % % % % %     delta.h = 0.04;
    % % % % % % % % % % % %     mean(A_theta_HMM) =  0.8167    0.3762    0.7429
    % % % % % % % % % % % %     mean(A_H_HMM) = 0.9364
    % % % % % % % % % %     delta.t = [0.1, 0.01, 0.003]; 
    % % % % % % % % % %     delta.h = 0.06;
    % % % % % % % % % % %     mean(A_theta_HMM) =  0.6321    0.3579    0.4505  
    % % % % % % % % % % %     mean(A_H_HMM) =  0.9011
    % % % % % % % % %     delta.t = [0.2, 0.01, 0.003]; 
    % % % % % % % % %     delta.h = 0.08;
    % % % % % % % % % %     mean(A_theta_HMM) = 0.4784    0.3520    0.4328
    % % % % % % % % % %     mean(A_H_HMM)= 0.8650
    % % % % % % % %     delta.t = [0.3, 0.01, 0.004]; 
    % % % % % % % %     delta.h = 0.08;  
    % % % % % % % % %     mean(A_theta_HMM) = 0.3718    0.3490    0.3493   
    % % % % % % % % %     mean(A_H_HMM) = 0.8630
    % % % % % % %     delta.t = [0.3, 0.01, 0.004]; 
    % % % % % % %     delta.h = 0.1;  
    % % % % % % % %     mean(A_theta_HMM) =  0.4001    0.3520    0.3603
    % % % % % % % %     mean(A_H_HMM) = 0.8348
    % % % % % %     delta.t = [0.3, 0.01, 0.004]; 
    % % % % % %     delta.h = 0.5;
    % % % % % % %     mean(A_H_HMM) = 0.4003
    % % % % % % %     mean(A_theta_HMM)= 0.3685    0.3528    0.3345
    % % % % % % Starting from the true value of theta with h init
    % % % % %     delta.t = [0.3, 0.01, 0.004]; 
    % % % % %     delta.h = 0.5;
    % % % % % %     mean(A_H_HMM) = 0.4053
    % % % % % %     mean(A_theta_HMM) = 0.3789    0.3586    0.3486
    % % % % % Starting from the true value of sigma2, with mu and phi init and with h init
    % % % %     delta.t = [0.3, 0.01, 0.004]; 
    % % % %     delta.h = 0.5;
    % % % % %     mean(A_H_HMM) = 0.4070
    % % % % %     mean(A_theta_HMM) = 0.3844    0.3618    0.3443
    % % % % Starting from a higher (0.25^2) than the true (0.2^2) one value of sigma2, with mu and phi init and with h init
    % % %     delta.t = [0.3, 0.01, 0.004]; 
    % % %     delta.h = 0.5;
    % % % %     mean(A_H_HMM) = 0.4034
    % % % %     mean(A_theta_HMM) =  0.3717    0.3486    0.3422
    % % % Starting from a bit lower (0.19^2) than the true (0.2^2) one value of sigma2, with mu and phi init and with h init
    % %     delta.t = [0.3, 0.01, 0.004]; 
    % %     delta.h = 0.5;
    % % %     mean(A_H_HMM) = 0.4001
    % % %     mean(A_theta_HMM) =  0.3616    0.3561    0.3363
    % % Starting from a lower (0.17^2) than the true (0.2^2) one value of sigma2, with mu and phi init and with h init
    %     delta.t = [0.3, 0.01, 0.004]; 
    %     delta.h = 0.5;
    % %     mean(A_H_HMM) = 0.4011
    % %     mean(A_theta_HMM) =   0.3666    0.3534    0.3342
    % Starting from a lower (0.16^2) than the true (0.2^2) one value of sigma2, with mu and phi init and with h init
        delta.t = [0.3, 0.01, 0.004];
        if (N_bin >= 30)
            delta.h = 0.5;
        else
            delta.h = 0.1;            
        end
    %     mean(A_H_HMM) = 0.3998
    %     mean(A_theta_HMM) = 0.3686    0.3516    0.3347

        h = h_init;
        theta = theta_init;


        H_HMM = zeros(M,T);
        theta_HMM = zeros(M,3);

        accept_HMM = zeros(M,1);
        A_H_HMM = zeros(M,1);
        A_theta_HMM = zeros(M,3);

        theta = theta_init2;

        fprintf('Semi DA, even states imputed\n')
        tic
        for ii = -BurnIn:1:M 
            if (mod(ii,1000) == 1)
                toc;
            end
            [h, acc, A_h] = update_h_HMM_v2(y, h, theta, delta.h, bins, bin_midpoint);
            [theta, A_theta] = update_theta_HMM_RW(y, h, theta, bins, bin_midpoint, hyper, delta.t);

            if (ii > 0)
                theta_HMM(ii,:) = theta;

                H_HMM(ii,:) = h;
                accept_HMM(ii,1) = acc;
                A_H_HMM(ii,1) = A_h;
                A_theta_HMM(ii,:) = A_theta;        
            end
        end
        time_HMM = toc;
        accept_HMM = accept_HMM/(T/2);
        A_H_HMM = A_H_HMM/(T/2);

        mean_H_HMM = mean(H_HMM);

        ESS_theta_HMM_40 = ESS(theta_HMM,40);
        ESS_theta_HMM_100 = ESS(theta_HMM,100);
        ESS_theta_HMM_1000 = ESS(theta_HMM,1000);
        ESS_theta_HMM_sig = ESS(theta_HMM,0);     

        ESS_H_HMM_40 = ESS(H_HMM,40);
        ESS_H_HMM_100 = ESS(H_HMM,100);
        ESS_H_HMM_1000 = ESS(H_HMM,1000);
        ESS_H_HMM_sig = ESS(H_HMM,0);

        H_subset_HMM = H_HMM(:,(2:10)*100);
        
        name = ['SV_results_param_HMM_Nbin',num2str(N_bin),'.mat'];
        save(name,...
            'y','h_true','theta_true','delta','hyper',...
            'H_subset_HMM',...
            'mean_H_HMM',...
            'theta_HMM',...
            'accept_HMM',...
            'A_H_HMM',...
            'A_theta_HMM',...
            'ESS_theta_HMM_40','ESS_theta_HMM_100','ESS_theta_HMM_1000','ESS_theta_HMM_sig',...
            'ESS_H_HMM_40','ESS_H_HMM_100','ESS_H_HMM_1000','ESS_H_HMM_sig',...    
            'time_HMM');
    end

%% SEMI DA shift: 
    if arg4
        % if shift integrate out the odd h(t)'s and impute the even ones
        % if not shift integrate out the even h(t)'s and impute the odd ones

%         delta.h = 0.1; 
%         delta.t = [0.4, 0.01, 0.004]; 
% %         mean(A_H_HMM_shift) = 0.8325
% %         mean(A_theta_HMM_shift) = 0.3770    0.3657    0.3690
        delta.t = [0.4, 0.01, 0.004]; 
        delta.h = 0.5;
        h = h_init;
        theta = theta_init;

        H_HMM_shift = zeros(M,T);
        theta_HMM_shift = zeros(M,3);

        accept_HMM_shift = zeros(M,1);
        A_H_HMM_shift = zeros(M,1);
        A_theta_HMM_shift = zeros(M,3);

        theta = theta_init3;

        fprintf('Semi DA, shift\n')
        tic
        for ii = -BurnIn:1:M 
    %     for ii = 1:1:M 
            if (mod(ii,1000) == 1)
                toc;
            end
            shift = mod(ii,2);
            [h, acc, A_h] = update_h_HMM_v2_shift(y, h, theta, delta.h, bins, bin_midpoint, shift);
            [theta, A_theta] = update_theta_HMM_RW_shift(y, h, theta, bins, bin_midpoint, hyper, delta.t, shift);

            if (ii > 0)
                theta_HMM_shift(ii,:) = theta;

                H_HMM_shift(ii,:) = h;
                accept_HMM_shift(ii,1) = acc;
                A_H_HMM_shift(ii,1) = A_h;
                A_theta_HMM_shift(ii,:) = A_theta;        
            end
        end
        time_HMM_shift = toc;

        accept_HMM_shift = accept_HMM_shift/(T/2);
        A_H_HMM_shift = A_H_HMM_shift/(T/2);
        mean_H_HMM_shift = mean(H_HMM_shift);

        ESS_theta_HMM_shift_40 = ESS(theta_HMM_shift,40);
        ESS_theta_HMM_shift_100 = ESS(theta_HMM_shift,100);
        ESS_theta_HMM_shift_1000 = ESS(theta_HMM_shift,1000);

        ESS_H_HMM_shift_40 = ESS(H_HMM_shift,40);
        ESS_H_HMM_shift_100 = ESS(H_HMM_shift,100);
        ESS_H_HMM_shift_1000 = ESS(H_HMM_shift,1000);

        H_subset_HMM_shift = H_HMM_shift(:,(2:10)*100);
        H_subset_HMM_shift2 = H_HMM_shift(:,(2:10)*100-1);

        name = 'SV_results_param_HMM_shift.mat';
        save(name,...
            'y','h_true','theta_true','delta','hyper',...
            'mean_H_HMM_shift',...
            'theta_HMM_shift',...
            'accept_HMM_shift',...
            'A_H_HMM_shift',...
            'A_theta_HMM_shift',...
            'time_HMM_shift',...
            'H_subset_HMM_shift','H_subset_HMM_shift2',...
            'ESS_theta_HMM_shift_40','ESS_theta_HMM_shift_100','ESS_theta_HMM_shift_1000',...
            'ESS_H_HMM_shift_40','ESS_H_HMM_shift_100','ESS_H_HMM_shift_1000',...
            'time_HMM_shift');
    end
 

%% SEMI DA efficient implementation:
    if ((nargin == 5) || ((nargin == 6) && (arg5 == 1)))
        % integrate out the odd h(t)'s and impute the even ones 
% %         delta.t = [0.3, 0.01, 0.004]; 
%         delta.t = [0.1, 0.005, 0.001]; 
%         delta.h = 0.5; 
        if (N_bin >= 30)
            delta.t = [0.3, 0.01, 0.004]; 
            delta.h = 0.5;
        else
            delta.t = [0.03, 0.01, 0.004]; 
            delta.h = 0.1;            
        end
        
        h = h_init;
        theta = theta_init;

        H_HMM_eff  = zeros(M,T);
        theta_HMM_eff  = zeros(M,3);

        accept_HMM_eff = zeros(M,1);
        A_H_HMM_eff  = zeros(M,1);
        A_theta_HMM_eff  = zeros(M,3);

        theta = theta_init2;    
        
        loglik = loglik_h_HMM_eff(y, h, theta, bins, bin_midpoint); 
        % this is the loglik for the integrals ONLY 
        % so the denominator for the theta updates (up to the prior terms)
        % does NOT include the observations logliks for the imputed h's 
                                                                   
        fprintf('Semi DA, even states imputed efficient implementation\n')
%         profile on
        tic
        for ii = -BurnIn:1:M 
%         for ii = 1:1:10 
            if (mod(ii,1000) == 1)
                toc;
            end
            [h, acc, A_h, loglik] = update_h_HMM_eff(y, h, theta, delta.h, bins, bin_midpoint, loglik);
            [theta, A_theta, loglik] = update_theta_HMM_RW_eff(y, h, theta, bins, bin_midpoint, hyper, delta.t, loglik);

            if (ii > 0)
                theta_HMM_eff (ii,:) = theta;

                H_HMM_eff (ii,:) = h;
                accept_HMM_eff (ii,1) = acc;
                A_H_HMM_eff (ii,1) = A_h;
                A_theta_HMM_eff(ii,:) = A_theta;        
            end
        end
        time_HMM_eff = toc;
%         profile off
%         profile viewer
        accept_HMM_eff = accept_HMM_eff/(T/2);
        A_H_HMM_eff = A_H_HMM_eff/(T/2);

        mean_H_HMM_eff = mean(H_HMM_eff);

        ESS_theta_HMM_eff_40 = ESS(theta_HMM_eff,40);
        ESS_theta_HMM_eff_100 = ESS(theta_HMM_eff,100);
        ESS_theta_HMM_eff_1000 = ESS(theta_HMM_eff,1000);
        ESS_theta_HMM_eff_sig = ESS(theta_HMM_eff,0);

        ESS_H_HMM_eff_40 = ESS(H_HMM_eff,40);
        ESS_H_HMM_eff_100 = ESS(H_HMM_eff,100);
        ESS_H_HMM_eff_1000 = ESS(H_HMM_eff,1000);
        ESS_H_HMM_eff_sig = ESS(H_HMM_eff,0);
        
        H_subset_HMM_eff = H_HMM_eff(:,(2:10)*100);
        
        name = ['SV_results_param_HMM_eff_Nbin',num2str(N_bin),'.mat'];
        save(name,...
            'y','h_true','theta_true','delta','hyper',...
            'H_subset_HMM_eff',...
            'mean_H_HMM_eff',...
            'theta_HMM_eff',...
            'accept_HMM_eff',...
            'A_H_HMM_eff',...
            'A_theta_HMM_eff',...
            'ESS_theta_HMM_eff_40','ESS_theta_HMM_eff_100','ESS_theta_HMM_eff_1000','ESS_theta_HMM_eff_sig',...
            'ESS_H_HMM_eff_40','ESS_H_HMM_eff_100','ESS_H_HMM_eff_1000','ESS_H_HMM_eff_sig',...    
            'time_HMM_eff');
    end
   
    %% SEMI DA ADAPTIVE
    if (nargin == 6)
        delta.t = [0.3, 0.01, 0.004];
        if (N_bin >= 30)
            delta.h = 0.5;
        else
            delta.h = 0.1;            
        end
    %     mean(A_H_HMM_adapt) = 0.3998
    %     mean(A_theta_HMM_adapt) = 0.3686    0.3516    0.3347

        h = h_init;
        theta = theta_init;


        H_HMM_adapt = zeros(M,T);
        theta_HMM_adapt = zeros(M,3);

        accept_HMM_adapt = zeros(M,1);
        A_H_HMM_adapt = zeros(M,1);
        A_theta_HMM_adapt = zeros(M,3);

        theta = theta_init2;

        fprintf('Semi DA adaptive\n')
        tic
        for ii = -BurnIn:1:M 
            if (mod(ii,1000) == 1)
                toc;
            end
            [h, acc, A_h] = update_h_HMM_adapt(y, h, theta, delta.h, mid);
            [theta, A_theta] = update_theta_HMM_RW_adapt(y, h, theta, mid, hyper, delta.t);

            if (ii > 0)
                theta_HMM_adapt(ii,:) = theta;

                H_HMM_adapt(ii,:) = h;
                accept_HMM_adapt(ii,1) = acc;
                A_H_HMM_adapt(ii,1) = A_h;
                A_theta_HMM_adapt(ii,:) = A_theta;        
            end
        end
        time_HMM_adapt = toc;
        accept_HMM_adapt = accept_HMM_adapt/(T/2);
        A_H_HMM_adapt = A_H_HMM_adapt/(T/2);

        mean_H_HMM_adapt = mean(H_HMM_adapt);

        ESS_theta_HMM_adapt_40 = ESS(theta_HMM_adapt,40);
        ESS_theta_HMM_adapt_100 = ESS(theta_HMM_adapt,100);
        ESS_theta_HMM_adapt_1000 = ESS(theta_HMM_adapt,1000);
        ESS_theta_HMM_adapt_sig = ESS(theta_HMM_adapt,0);     

        ESS_H_HMM_adapt_40 = ESS(H_HMM_adapt,40);
        ESS_H_HMM_adapt_100 = ESS(H_HMM_adapt,100);
        ESS_H_HMM_adapt_1000 = ESS(H_HMM_adapt,1000);
        ESS_H_HMM_adapt_sig = ESS(H_HMM_adapt,0);

        H_subset_HMM_adapt = H_HMM_adapt(:,(2:10)*100);
        
        name = ['SV_results_param_HMM_adapt_Nbin',num2str(N_bin),'.mat'];
        save(name,...
            'y','h_true','theta_true','delta','hyper',...
            'H_subset_HMM_adapt',...
            'mean_H_HMM_adapt',...
            'theta_HMM_adapt',...
            'accept_HMM_adapt',...
            'A_H_HMM_adapt',...
            'A_theta_HMM_adapt',...
            'ESS_theta_HMM_adapt_40','ESS_theta_HMM_adapt_100','ESS_theta_HMM_adapt_1000','ESS_theta_HMM_adapt_sig',...
            'ESS_H_HMM_adapt_40','ESS_H_HMM_adapt_100','ESS_H_HMM_adapt_1000','ESS_H_HMM_adapt_sig',...    
            'time_HMM_adapt');   
    
    end
    
end