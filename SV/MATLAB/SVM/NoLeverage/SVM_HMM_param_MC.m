clear all
close all
addpath(genpath('../../MATLAB'));

plot_on = false;

arg0 = 0;
arg1 = 1;
arg2 = 1;
path = '';
data_name = '_sim';   


S = 50;
THETA_mean_DA = zeros(S,4);
THETA_sd_DA = zeros(S,4);
THETA_mean_HMM = zeros(S,4);
THETA_sd_HMM = zeros(S,4);
TIME_DA = zeros(S,1);
TIME_HMM = zeros(S,1);

beta = -0.2;
mu = -1;
theta_true = [mu, 0.98, 0.01, beta];
T = 2000;
    
hyper.S = [5/2, 0.05/2];
hyper.M = 10;
hyper.P = [20, 3/2];
hyper.B = 10;

bin_range = 4;
N_bin = 20; 
bins = linspace(-bin_range,bin_range,N_bin+1);
bin_midpoint = (bins(1:N_bin) + bins(2:N_bin+1))/2;
% bins are the demeaned volatilities: b=h-mu
fprintf('S=%i, T=%i, N_bin=%i \n',S,T,N_bin);

M = 10000;
BurnIn = 10000;
theta_init = [0.2, 0.95, 0.1^2, 0.5]; 

delta.h = 0.25;
delta.t = [0.35, 0.01, 0.0008, 0.065]; % 

delta_HMM.h = 0.1;
delta_HMM.t = [0.05, 0.005, 0.0005, 0.01]; %      
        
ind_h_sel = 2:50:T;
            
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
            if (mod(ii,1000) == 0)
                toc;
            end
%             [h, acc, A_h] = update_h(y, h, theta, delta.h);
            [theta, A_theta] = update_theta_RW(y, h, theta, hyper, delta.t);
            if (ii > 0)
                theta_DA_RW(ii,:) = theta;
                H_subset_DA_RW(ii,:) = h(ind_h_sel);

    [mean_H_DA_RW, s_H_DA_RW, var_H_DA_RW] = h_seq_stats(h, mean_H_DA_RW, s_H_DA_RW, ii);                

%                 accept_DA_RW(ii,1) = acc;
%                 A_H_DA_RW(ii,1) = A_h;
                A_theta_DA_RW(ii,:) = A_theta;        
            end
        end
        time_DA_RW = toc;
        THETA_mean_DA(ss,:) = mean(theta_DA_RW);
        THETA_sd_DA(ss,:) = std(theta_DA_RW);    
        TIME_DA(ss) = time_DA_RW;

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
%             [h, acc, A_h] = update_h_HMM(y, h, theta, delta_HMM.h, bins, bin_midpoint);
            [theta, A_theta] = update_theta_HMM_RW(y, h, theta, bins, bin_midpoint, hyper, delta_HMM.t);

            if (ii > 0)
                theta_HMM (ii,:) = theta;
                H_subset_HMM(ii,:) = h(ind_h_sel);

[mean_H_HMM, s_H_HMM, var_H_HMM] = h_seq_stats(h, mean_H_HMM, s_H_HMM, ii);                

%                 accept_HMM(ii,1) = acc;
%                 A_H_HMM(ii,1) = A_h;
                A_theta_HMM(ii,:) = A_theta;        
            end
        end
        time_HMM = toc;
        TIME_HMM(ss) = time_HMM;
        THETA_mean_HMM(ss,:) = mean(theta_HMM);
        THETA_sd_HMM(ss,:) = std(theta_HMM);    
    end

end

name = ['MC_fixed_H_T',num2str(T),'_Nbin',num2str(N_bin),'.mat'];
save(name,'theta_true','delta','hyper',...
    'T','S','N_bin',...
    'THETA_mean_DA','THETA_mean_HMM',...
    'THETA_sd_DA','THETA_sd_HMM',...
    'TIME_DA','TIME_HMM')

if false
	params = {'\mu','\phi','\sigma^2','\beta'};
    ff = figure(1);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    
    for ii = 1:4
        subplot(2,2,ii)
        bar([THETA_mean_DA(:,ii),THETA_mean_HMM(:,ii)])
        xlim([0 (S+1)])
        hold on
        plot(mean(THETA_mean_DA(:,ii)) + 0*(0:(S+1)),'b')
        plot(mean(THETA_mean_HMM(:,ii)) + 0*(0:(S+1)),'y')
        plot(theta_true(ii) + 0*(0:(S+1)),'r')
        hold off
        title(params{ii})
    end
    legend('DA','HMM','mean DA','mean HMM','true')
    name = ['MC_fixed_H_T',num2str(T),'_Nbin',num2str(N_bin),'.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print(ff,name,'-dpng','-r0')
end