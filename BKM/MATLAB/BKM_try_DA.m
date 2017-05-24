clear all
close all

plot_on = false;

sc = 1;
[y, T, time, stdT, f, m, T1, T2] = BKM_Data_HMM(sc);

Na = round( [1000, 1000, 1092.23, 1100.01, 1234.32, 1460.85, 1570.38, 1819.79,...
    1391.27, 1507.60, 1541.44, 1631.21, 1628.60, 1609.33, 1801.68, 1809.08, 1754.74,...
    1779.48, 1699.13, 1681.39, 1610.46, 1918.45, 1717.07, 1415.69, 1229.02, 1082.02,...
    1096.61, 1045.84, 1137.03, 981.1, 647.67, 992.65, 968.62, 926.83, 952.96, 865.64]/sc);
% N1 = 400*ones(1,T)/sc;
N1 = [400   400   400   400   400   400   400   400   400   400   400   400   400   400   400   400    40   400   400 ...
 40   400   400    40    40   400   400   400   400   400   400   400   400   400    40   400   400]/sc;
N = [N1;Na];

alpha1 = 1;
alphaa = 2;
alphar = -2;
alphal = -4;
beta1 =-2;
betaa = 0.1;
betar = -0.7;
betal = -0.3;
sigy = 1;

params = {'alpha1', 'alphaa', 'alphar', 'alphal', ...
    'beta1', 'betaa', 'betar', 'betal',...
    'sigy'};

theta_init = [alpha1, alphaa, alphar, alphal, beta1, betaa, betar, betal, sigy];
theta = theta_init;
D = size(theta,2);
prior.N = [200 2000 0.5];
prior.S = [0.001,0.001];
prior.T_mu = 0*ones(D-1,1);
prior.T_sigma2 = 100*ones(D-1,1);

delta.T = 0.1*ones(D-1,1);
delta.N = [20.5, 100.5]; %0.5 added to have a correct dicrete uniform distribution after rounding
oldlikhood = BKM_calclikhood(N, theta, y, m, f, stdT, prior.N);

[phi1, phia, rho, lambda] = BKM_covariates(theta,f,stdT);  
 
M = 10000;%50000;


%%
% update_N = 'U'; % or 'SP'
% N = [N1;Na];
% theta = theta_init;
% NN4 = zeros(2,T,M);
% sample4 = zeros(M,9);
% accept4 = zeros(M,1);
% tic
% for ii = 1:M
%     % Update the parameters in the model using function "updateparam": 
%     % Set parameter values and log(likelihood) value of current state to be the output from
%     % the MH step:
%     
%     if (mod(ii,1000)==0)
%         fprintf('MH iter = %i\n',ii); toc;
%     end
%     [N, theta, A] = BKM_update_URW(N, theta, prior, delta, y, m, f, stdT,update_N);
%     NN4(:,:,ii) = N;
%     sample4(ii,:)= theta; 
%     accept4(ii) = A; 
% end
% time = toc;
% name = 'BKM_DA_results4_long.mat';
% save(name,'sample4','NN4','accept4','theta_init','prior','delta','time');

%% WITH NORMAL PROPOSALS FOR THE REGRESSION COEFFS
update_N = 'SP'; % 'U' or 'SP'
N = [N1;Na];
theta = theta_init;
NN5 = zeros(2,T,M);
sample5 = zeros(M,9);
accept5 = zeros(M,1);
tic
for ii = 1:M
    % Update the parameters in the model using function "updateparam": 
    % Set parameter values and log(likelihood) value of current state to be the output from
    % the MH step:
    
    if (mod(ii,1000)==0)
        fprintf('MH iter = %i\n',ii); toc;
    end
    [N, theta, A] = BKM_update_NRW(N, theta, prior, delta, y, m, f, stdT, update_N);
%     [N, theta, A] = BKM_update_URW(N, theta, prior, delta, y, m, f, stdT, update_N);
    NN5(:,:,ii) = N;
    sample5(ii,:)= theta; 
    accept5(ii) = A; 
end
time_sampl = toc;
name = 'BKM_DA_results_U.mat';
save(name,'sample5','NN5','accept5','theta_init','prior','delta','time_sampl');

%%
if plot_on
    figure(11)
    set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
    for ii = 1:9
        subplot(3,3,ii)
        plot(sample5(:,ii))
        title(params{ii})
    end

    figure(10)
    set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
    for ii = 1:9
        subplot(3,3,ii)
        autocorr(sample4(:,ii),40)
        title(params{ii})
    end   
    
    figure(110)
    set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
    for ii = 1:9
        subplot(3,3,ii)
        autocorr(sample5(:,ii),40)
        title(params{ii})
    end
    
    figure(2)
    set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
    for ii = 1:9
        subplot(3,3,ii)
        plot(squeeze(NN5(1,4*(ii-3)+9,:)))
        title(['N1(',num2str(4*(ii-3)+9),')'])
    end


    figure(3)
    set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
    for ii = 1:9
        subplot(3,3,ii)
        plot(squeeze(NN5(2,4*(ii-3)+9,:)))
        title(['Na(',num2str(4*(ii-3)+9),')'])
    end


    figure(4)
    subplot(2,1,1)
    hold on
    plot(mean(NN4(1,:,:),3))
%     plot(mean(NN2(1,:,:),3),'r')
    hold off

    subplot(2,1,2)
    hold on
    plot(mean(NN4(2,:,:),3))
%     plot(mean(NN2(2,:,:),3),'r')
    hold off
end