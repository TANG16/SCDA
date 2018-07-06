clear all
close all
addpath(genpath('../../MATLAB'));

plot_on = false;

arg0 = 0;
arg1 = 1;
arg2 = 1;
arg3 = 1;
path = '';
data_name = '_sim';   

T = 3000;
ind_h_sel = 2:50:T;

S = 1;% 50;
if (S > 1)
    THETA_mean_DA = zeros(S,4);
    THETA_sd_DA = zeros(S,4);
    TIME_DA = zeros(S,1);
    A_TH_DA = zeros(S,4);
    A_H_DA = zeros(S,1);
    ESS_TH_DA = zeros(S,4);
    ESS_H_DA = zeros(S,length(ind_h_sel));

    THETA_mean_HMM = zeros(S,4);
    THETA_sd_HMM = zeros(S,4);
    TIME_HMM = zeros(S,1);
    A_TH_HMM = zeros(S,4);
    A_H_HMM = zeros(S,1);
    ESS_TH_HMM = zeros(S,4);
    ESS_H_HMM = zeros(S,length(ind_h_sel));

    THETA_mean_HMM_adapt = zeros(S,4);
    THETA_sd_HMM_adapt = zeros(S,4);
    TIME_HMM_adapt = zeros(S,1);
    A_TH_HMM_adapt = zeros(S,4);
    A_H_HMM_adapt = zeros(S,1);
    ESS_TH_HMM_adapt = zeros(S,4);
    ESS_H_HMM_adapt = zeros(S,length(ind_h_sel));
end

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
N_q = 20;
mid = (0:(N_q-1))/N_q;
mid = mid + mid(2)/2;
mid_inv = norminv(mid);



fprintf('S=%i, T=%i, N_bin=%i, N_q=%i \n',S,T,N_bin,N_q);

M = 10000;
BurnIn = 1000;
theta_init = [0.2, 0.95, 0.1^2, 0.5]; 

delta.h = 0.25;
delta.t = [0.35, 0.01, 0.0008, 0.065]; % 

delta_HMM.h = 0.3;
delta_HMM.t = [0.05, 0.008, 0.0008, 0.07]; %      

delta_HMM_adapt.h = 0.4;
delta_HMM_adapt.t = [0.25, 0.009, 0.0009, 0.08]; %      
 
ss = 1;
            
parfor ss = 1:S
    fprintf('***** SEED = %i ***** \n', ss)
    s = RandStream('mt19937ar','Seed',ss);
    RandStream.setGlobalStream(s); 

    [y,h_true,y0] = generate_SVM(theta_true,T);
    h_init = log(var(y))*ones(1,T); 


    %% FULL DA - RANDOM WALK UPDATES 
    if arg1 
%     %     if (arg0 == 0)
%             delta.h = 0.25;
%             delta.t = [0.35, 0.01, 0.0008, 0.065]; % 
%     %     elseif (arg0 == 2)
%     %         delta.h = 0.5;
%     %         delta.t = [0.1, 0.005, 0.0007, 0.1]; %     
%     %     end
%         h = h_init;
%         theta = theta_init;
        h = h_true;
        theta = theta_true; 
        
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
%         for ii = 1:1:M 
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
        THETA_mean_DA(ss,:) = mean(theta_DA_RW);
        THETA_sd_DA(ss,:) = std(theta_DA_RW);    
        TIME_DA(ss) = time_DA_RW;

        A_TH_DA(ss,:) = mean(A_theta_DA_RW);
        A_H_DA(ss,:) = mean(A_H_DA_RW)/T;
        ESS_TH_DA(ss,:) = ESS(theta_DA_RW,0);
        ESS_H_DA(ss,:) = ESS(H_subset_DA_RW,0);
    end

    %% SEMI DA efficient implementation:
    if (arg2 > 0)        % integrate out the odd h(t)'s and impute the even ones  
%         delta_HMM.h = 0.1;
%         delta_HMM.t = [0.05, 0.005, 0.0005, 0.01]; %

    %         delta.h = 0.25;
    %         delta.t = [0.35, 0.01, 0.0008, 0.065]; %

    %     h = h_init;
%         theta = theta_init;
    % theta(3) = theta_true(3);

        h = h_true;
        theta = theta_true; 
        
        H_subset_HMM = zeros(M,length(ind_h_sel));       
        mean_H_HMM = zeros(1,T);
        s_H_HMM  = zeros(1,T);

        theta_HMM  = zeros(M,4);

        accept_HMM = zeros(M,1);
        A_H_HMM  = zeros(M,1);
        A_theta_HMM  = zeros(M,4);

%         loglik = loglik_h_HMM_eff(y, h, theta, bins, bin_midpoint); 
        % this is the loglik for the integrals ONLY 
        % so the denominator for the theta updates (up to the prior terms)
        % does NOT include the observations logliks for the imputed h's 

        fprintf('Semi DA, even states imputed \n')

        tic
%         for ii = -BurnIn:1:M 
        for ii = 1:M 
            if (mod(ii,1000) == 1)
                tt = toc;
                fprintf('Iter = %i, Elapsed time = %6.4f s \n', ii, tt);
            end
%             [h, acc, A_h, loglik] = update_h_HMM_eff(y, h, theta, delta.h, bins, bin_midpoint, loglik);
%             [theta, A_theta, loglik] = update_theta_HMM_RW_eff(y, h, theta, bins, bin_midpoint, hyper, delta.t, loglik);
            [h, acc, A_h] = update_h_HMM(y, h, theta, delta_HMM.h, bins, bin_midpoint);
            [theta, A_theta] = update_theta_HMM_RW(y, h, theta, bins, bin_midpoint, hyper, delta_HMM.t);

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
        TIME_HMM(ss) = time_HMM;
        THETA_mean_HMM(ss,:) = mean(theta_HMM);
        THETA_sd_HMM(ss,:) = std(theta_HMM);    

        A_TH_HMM(ss,:) = mean(A_theta_HMM);
        A_H_HMM(ss,:) = mean(A_H_HMM/(T/2));
        ESS_TH_HMM(ss,:) = ESS(theta_HMM,0);
        ESS_H_HMM(ss,:) = ESS(H_subset_HMM,0);
    end

   %% SEMI DA efficient implementation:
    if (arg3 > 0)        % integrate out the odd h(t)'s and impute the even ones  
%         delta_HMM_adapt.h = 0.1;
%         delta_HMM_adapt.t = [0.05, 0.005, 0.0005, 0.01]; %

    %         delta.h = 0.25;
    %         delta.t = [0.35, 0.01, 0.0008, 0.065]; %

    %     h = h_init;
%         theta = theta_init;
    % theta(3) = theta_true(3);

        h = h_true;
        theta = theta_true; 
        
        H_subset_HMM_adapt = zeros(M,length(ind_h_sel));       
        mean_H_HMM_adapt = zeros(1,T);
        s_H_HMM_adapt  = zeros(1,T);

        theta_HMM_adapt  = zeros(M,4);

        accept_HMM_adapt = zeros(M,1);
        A_H_HMM_adapt  = zeros(M,1);
        A_theta_HMM_adapt  = zeros(M,4);

%         loglik = loglik_h_HMM_adapt_eff(y, h, theta, bins, bin_midpoint); 
        % this is the loglik for the integrals ONLY 
        % so the denominator for the theta updates (up to the prior terms)
        % does NOT include the observations logliks for the imputed h's 

        fprintf('Semi DA Adapted, even states imputed \n')

        tic
%         for ii = -BurnIn:1:M 
        for ii = 1:M 
            if (mod(ii,1000) == 1)
                tt = toc;
                fprintf('Iter = %i, Elapsed time = %6.4f s \n', ii, tt);
            end
            [theta, A_theta] = update_theta_HMM_RW_adapt(y, h, theta, mid_inv,  hyper, delta_HMM_adapt.t);
            [h, acc, A_h] = update_h_HMM_adapt(y, h, theta, delta_HMM_adapt.h, mid_inv);
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
        TIME_HMM_adapt(ss) = time_HMM_adapt;
        THETA_mean_HMM_adapt(ss,:) = mean(theta_HMM_adapt);
        THETA_sd_HMM_adapt(ss,:) = std(theta_HMM_adapt);  
       
        A_TH_HMM_adapt(ss,:) = mean(A_theta_HMM_adapt);
        A_H_HMM_adapt(ss,:) = mean(A_H_HMM_adapt/(T/2));
        ESS_TH_HMM_adapt(ss,:) = ESS(theta_HMM_adapt,0);
        ESS_H_HMM_adapt(ss,:) = ESS(H_subset_HMM_adapt,0);
    end

end

name = ['MC_est_H_T',num2str(T),'_Nbin',num2str(N_bin),'_Nq',num2str(N_q),'.mat'];
save(name,'theta_true','delta','delta_HMM','delta_HMM_adapt','hyper','BurnIn',...
    'T','S','N_bin','N_q',...
    'THETA_mean_DA','THETA_mean_HMM','THETA_mean_HMM_adapt',...
    'THETA_sd_DA','THETA_sd_HMM','THETA_sd_HMM_adapt',...
    'TIME_DA','TIME_HMM','TIME_HMM_adapt',...
    'A_TH_DA','A_TH_HMM','A_TH_HMM_adapt',...
    'A_H_DA','A_H_HMM','A_H_HMM_adapt',...
    'ESS_TH_DA','ESS_TH_HMM','ESS_TH_HMM_adapt',...
    'ESS_H_DA','ESS_H_HMM','ESS_H_HMM_adapt')


if false
    ColPal = [      1         1    0
               0.9290    0.6940    0.1250
               0.3010    0.7450    0.9330
                    0    0.4470    0.7410 
               0.8500    0.3250    0.0980
               0.6350    0.0780    0.1840];
        
	params = {'\mu','\phi','\sigma^2','\beta'};
    ff = figure(1);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    
    for ii = 1:4
        subplot(2,2,ii)
        b = bar([THETA_mean_DA(:,ii),THETA_mean_HMM(:,ii),THETA_mean_HMM_adapt(:,ii)])
        set(b(1), 'FaceColor',ColPal(1,:))
        set(b(2), 'FaceColor',ColPal(3,:))
        set(b(3), 'FaceColor',ColPal(5,:))
        xlim([0 (S+1)])
        hold on
        plot(mean(THETA_mean_DA(:,ii)) + 0*(0:(S+1)),'Color',ColPal(2,:))
        plot(mean(THETA_mean_HMM(:,ii)) + 0*(0:(S+1)),'Color',ColPal(4,:))
        plot(mean(THETA_mean_HMM_adapt(:,ii)) + 0*(0:(S+1)),'Color',ColPal(6,:))
        plot(theta_true(ii) + 0*(0:(S+1)),'k')
        hold off
        title(params{ii})
    end
    legend('DA','HMM','HMM adapt','mean DA','mean HMM','mean HMM adapt','true')
    name = ['MC_est_H_T',num2str(T),'_Nbin',num2str(N_bin),'_Nq',num2str(N_q),'.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print(ff,name,'-dpng','-r0')
end