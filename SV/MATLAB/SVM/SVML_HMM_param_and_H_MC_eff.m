% clear all
close all
addpath(genpath('../../MATLAB'));

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s); 

plot_on = false;

arg0 = 0; % 0 - sim; 2 -- IBM
arg1 = 0; %1; % DA
arg2 = 1; % fixed
arg3 = 0;  % adaptive
path = '';

T = 2000; %3000;
ind_h_sel = 2:50:T;

S = 36; %24; %1; %50; %20; %24; %100;
if (S > 1)
    if (arg1 > 0)
        THETA_mean_DA = zeros(S,5);
        THETA_sd_DA = zeros(S,5);
        TIME_DA = zeros(S,1);
        A_TH_DA = zeros(S,5);
        A_H_DA = zeros(S,1);
        ESS_TH_DA = zeros(S,5);
        ESS_H_DA = zeros(S,length(ind_h_sel));
    end
    if (arg2 > 0)
        THETA_mean_HMM = zeros(S,5);
        THETA_sd_HMM = zeros(S,5);
        TIME_HMM = zeros(S,1);
        A_TH_HMM = zeros(S,5);
        A_H_HMM = zeros(S,1);
        ESS_TH_HMM = zeros(S,5);
        ESS_H_HMM = zeros(S,length(ind_h_sel));
    end
    if (arg3 > 0)
        THETA_mean_HMM_adapt = zeros(S,5);
        THETA_sd_HMM_adapt = zeros(S,5);
        TIME_HMM_adapt = zeros(S,1);
        A_TH_HMM_adapt = zeros(S,5);
        A_H_HMM_adapt = zeros(S,1);
        ESS_TH_HMM_adapt = zeros(S,5);
        ESS_H_HMM_adapt = zeros(S,length(ind_h_sel));
    end
end


mu = -1.0;
phi = 0.98;
sigma2 = 0.10^2;
beta = -0.2; % 0 -- NO SV IN THE MEAN
rho = -0.2;%3;
theta_true = [mu, phi, sigma2, beta, rho];
%     theta_true = [-1.00, 0.98, 0.1^2, -0.20, -0.35];
 
    
hyper.S = [5/2, 0.05/2];
hyper.M = 10;
hyper.P = [20, 3/2];
hyper.B = 10;
hyper.R = 10;


% bins are the demeaned volatilities: b=h-mu
if (arg2 > 0)
    bin_range = 4;
    N_bin = 40; 
    bins = linspace(-bin_range,bin_range,N_bin+1);
    bin_midpoint = (bins(1:N_bin) + bins(2:N_bin+1))/2;
else
    N_bin = NaN;    
end
if (arg3 > 0)
    N_q = 15;
    mid = (0:(N_q-1))/N_q;
    mid = mid + mid(2)/2;
    mid_inv = norminv(mid);
else
    N_q = NaN;
end

fprintf('S=%i, T=%i, N_bin=%i, N_q=%i \n',S,T,N_bin,N_q);

M = 10000;
BurnIn = -1; %10000;
 
% M = 5*10000;
% BurnIn = 10000; 

theta_init = theta_true;

if (arg0 == 0) % simulation
    delta.h = 0.25;
% %     delta.t = [0.3, 0.01, 0.0008, 0.065, 0.04]; %    
%     delta.t = [0.25, 0.015, 0.00085, 0.06, 0.05]; %    
    delta.t = [0.25, 0.015, 0.00085, 0.04, 0.02]; %    
    if (arg2 > 0)
%         delta_HMM.h = 0.3; % 0.2599
        delta_HMM.h = 0.25;
% % %         delta_HMM.t = [0.05, 0.008, 0.0008, 0.07, 0.1]; %   
% %         delta_HMM.t = [0.05, 0.008, 0.0008, 0.0, 0.01]; %   
        delta_HMM.t = [0.07, 0.008, 0.001, 0.04, 0.02]; %   
%         delta_HMM.t = [0.01, 0.001, 0.0001, 0.02, 0.04]; %   
    else 
        delta_HMM = [];
    end
    if (arg3 > 0)    
        delta_HMM_adapt.h = 0.35;
        delta_HMM_adapt.t = [0.5, 0.009, 0.001, 0.08, 0.085];
    else
        delta_HMM_adapt = [];
    end
end

ss = 1;%1;
            
parfor ss = 1:S
    fprintf('***** SEED = %i ***** \n', ss)
    s = RandStream('mt19937ar','Seed',ss);
    RandStream.setGlobalStream(s); 

    [y,h_true] = generate_SVML(theta_true,T);
    h_init = h_true; 


    %% FULL DA - RANDOM WALK UPDATES 
    if arg1 
        h = h_true;
        theta = theta_true; 
        
        H_subset_DA_RW = zeros(M,length(ind_h_sel));
        mean_H_DA_RW = zeros(1,T);
        s_H_DA_RW  = zeros(1,T);
        theta_DA_RW = zeros(M,5);

        accept_DA_RW = zeros(M,1);
        A_h_DA_RW = zeros(M,1);
        A_theta_DA_RW = zeros(M,5);

        fprintf('Full DA, RW updates\n')
        tic
        for ii = -BurnIn:1:M 
            if (mod(ii,5000) == 0)
                toc;
            end
            [h, acc, A_h] = update_h_leverage(y, h, theta, delta.h);
            [theta, A_theta] = update_theta_leverage(y, h, theta, hyper, delta.t);    
            if (ii > 0)
                theta_DA_RW(ii,:) = theta;
                H_subset_DA_RW(ii,:) = h(ind_h_sel);

    [mean_H_DA_RW, s_H_DA_RW, var_H_DA_RW] = h_seq_stats(h, mean_H_DA_RW, s_H_DA_RW, ii);                

                accept_DA_RW(ii,1) = acc;
                A_h_DA_RW(ii,1) = A_h;
                A_theta_DA_RW(ii,:) = A_theta;        
            end
        end
        time_DA_RW = toc;
        THETA_mean_DA(ss,:) = mean(theta_DA_RW);
        THETA_sd_DA(ss,:) = std(theta_DA_RW);    
        TIME_DA(ss) = time_DA_RW;

        A_TH_DA(ss,:) = mean(A_theta_DA_RW);
        A_H_DA(ss,:) = mean(A_h_DA_RW)/T;
        ESS_TH_DA(ss,:) = ESS(theta_DA_RW,0);
        ESS_H_DA(ss,:) = ESS(H_subset_DA_RW,0);
    end
 
    %% SEMI DA efficient implementation:
    if (arg2 > 0)        % integrate out the odd h(t)'s and impute the even ones  

        h = h_true;
        theta = theta_true; 
        
        H_subset_HMM = zeros(M,length(ind_h_sel));       
        mean_H_HMM = zeros(1,T);
        s_H_HMM  = zeros(1,T);

        theta_HMM  = zeros(M,5);

        accept_HMM = zeros(M,1);
        A_h_HMM  = zeros(M,1);
        A_theta_HMM  = zeros(M,5);

        loglik = loglik_h_leverage_HMM_eff(y, h, theta, bins, bin_midpoint); 
        % this is the loglik for the integrals ONLY 
        % so the denominator for the theta updates (up to the prior terms)
        % does NOT include the observations logliks for the imputed h's 

        fprintf('Semi DA, even states imputed \n')

        tic
        for ii = -BurnIn:1:M 
            if (mod(ii,5000) == 0)
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
                A_h_HMM(ii,1) = A_h;
                A_theta_HMM(ii,:) = A_theta;        
            end
        end
        time_HMM = toc;
        TIME_HMM(ss) = time_HMM;
        THETA_mean_HMM(ss,:) = mean(theta_HMM);
        THETA_sd_HMM(ss,:) = std(theta_HMM);    

        A_TH_HMM(ss,:) = mean(A_theta_HMM);
        A_H_HMM(ss,:) = mean(A_h_HMM/(T/2));
        ESS_TH_HMM(ss,:) = ESS(theta_HMM,0);
        ESS_H_HMM(ss,:) = ESS(H_subset_HMM,0);
    end

   %% SEMI DA adaptive efficient implementation:
    if (arg3 > 0)        % integrate out the odd h(t)'s and impute the even ones  
        h = h_true;
        theta = theta_true; 

        H_subset_HMM_adapt = zeros(M,length(ind_h_sel));       
        mean_H_HMM_adapt = zeros(1,T);
        s_H_HMM_adapt  = zeros(1,T);

        theta_HMM_adapt  = zeros(M,5);

        accept_HMM_adapt = zeros(M,1);
        A_h_HMM_adapt  = zeros(M,1);
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
            [theta, A_theta, loglik] = update_theta_leverage_HMM_adapt_eff(y, h, theta, ...
                mid_inv, hyper, delta_HMM_adapt.t, loglik);
            [h, acc, A_h, loglik] = update_h_leverage_HMM_adapt_eff_v2(y, h, theta, ...
                delta_HMM_adapt.h, mid_inv, loglik);
%             [theta, A_theta] = update_theta_leverage_HMM_adap(y, h, theta, ...
%                 mid_inv, hyper, delta_HMM_adapt.t);
%             [h, acc, A_h] = update_h_leverage_HMM_adapt(y, h, theta, ...
%                 delta_HMM_adapt.h, mid_inv);            
            if (ii > 0)
                theta_HMM_adapt(ii,:) = theta;
                H_subset_HMM_adapt(ii,:) = h(ind_h_sel);

[mean_H_HMM_adapt, s_H_HMM_adapt, var_H_HMM_adapt] = h_seq_stats(h, mean_H_HMM_adapt, s_H_HMM_adapt, ii);                

                accept_HMM_adapt(ii,1) = acc;
                A_h_HMM_adapt(ii,1) = A_h;
                A_theta_HMM_adapt(ii,:) = A_theta;        
            end
        end
        time_HMM_adapt = toc;
        
% accept_HMM_adapt = accept_HMM_adapt/(T/2);
% A_h_HMM_adapt = A_h_HMM_adapt/(T/2); 
% 
% mean(A_theta_HMM_adapt)
% mean(A_h_HMM_adapt)
% mean(theta_HMM_adapt)


        TIME_HMM_adapt(ss) = time_HMM_adapt;
        THETA_mean_HMM_adapt(ss,:) = mean(theta_HMM_adapt);
        THETA_sd_HMM_adapt(ss,:) = std(theta_HMM_adapt);  
       
        A_TH_HMM_adapt(ss,:) = mean(A_theta_HMM_adapt);
        A_H_HMM_adapt(ss,:) = mean(A_h_HMM_adapt/(T/2));
        ESS_TH_HMM_adapt(ss,:) = ESS(theta_HMM_adapt,0);
        ESS_H_HMM_adapt(ss,:) = ESS(H_subset_HMM_adapt,0);
    end

end

name = ['SVML_MC_est_H_T',num2str(T),'_Nbin',num2str(N_bin),'_Nq',num2str(N_q),'.mat'];
save(name,'-regexp',...
    'theta_true','^delta','hyper','BurnIn',...
    'T','S','N_bin','N_q',...
    '^THETA_mean_',...
    '^THETA_sd',...
    '^TIME',...
    '^A_TH',...
    '^A_H_',...
    '^ESS_TH',...
    '^ESS_H')


if false
	params = {'\mu','\phi','\sigma^2','\beta','\rho'};
    ColPal = [      1         1    0
               0.9290    0.6940    0.1250
               0.3010    0.7450    0.9330
                    0    0.4470    0.7410 
               0.8500    0.3250    0.0980
               0.6350    0.0780    0.1840];
        
    ff = figure(11);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);    
    for ii = 1:5
        subplot(2,3,ii)
%         b = bar([THETA_mean_DA(:,ii)]);
        b = bar([THETA_mean_DA(:,ii),THETA_mean_HMM(:,ii),THETA_mean_HMM_adapt(:,ii)]);
%         b = bar([THETA_mean_DA(:,ii),THETA_mean_HMM(:,ii)]);
%         bar(THETA_mean_HMM(:,ii));
        set(b(1), 'FaceColor',ColPal(1,:))
        set(b(2), 'FaceColor',ColPal(3,:))
        set(b(3), 'FaceColor',ColPal(5,:))
        xlim([0 (S+1)])
        hold on
        plot((0:(S+1)),mean(THETA_mean_DA(:,ii)) + 0*(0:(S+1)),'Color',ColPal(2,:),'linewidth',2)
        plot((0:(S+1)),mean(THETA_mean_HMM(:,ii)) + 0*(0:(S+1)),'Color',ColPal(4,:),'linewidth',2)
        plot(mean(THETA_mean_HMM_adapt(:,ii)) + 0*(0:(S+1)),'Color',ColPal(6,:),'linewidth',2)
        plot((0:(S+1)),theta_true(ii) + 0*(0:(S+1)),'k','linewidth',2)
        hold off
        title(params{ii})
    end
    legend('DA',['HMM ',num2str(N_bin)],['HMM adapt ',num2str(N_q)],'mean DA','mean HMM','mean HMM adapt','true')
%     legend('DA',['HMM ',num2str(N_bin)], 'mean DA','mean HMM','true')
    name = ['MC_est_H_T',num2str(T),'_Nbin',num2str(N_bin),'_Nq',num2str(N_q),'.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print(ff,name,'-dpng','-r0')
end