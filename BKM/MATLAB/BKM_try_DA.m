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
if (sc == 1)
    alphal = -4;
else
    alphal = -1;
end
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

[phi1, phia, rho, lambda] = BKM_covariates(theta,f,stdT);  


prior.N = [200/sc 2000/sc 0.5];
prior.S = [0.001,0.001];
prior.T_mu = 0*ones(D-1,1);
prior.T_sigma2 = 100*ones(D-1,1);

priorN = prior.N;

logfact = @(xx) sum(log(1:1:xx));
logfact = arrayfun(logfact,0:7000);

oldlikhood = BKM_calclikhood(N, theta, y, m, f, stdT, prior.N, logfact);


%% Set the proposals
% for the states
update_N = 'U'; % 'U' or 'SP'
% for the parameters
update_T = 'NRW'; % 'NRW' or 'URW'

% step sizes 
% given step size delta, std for URW is delta/sqrt(3), for NRW 1*delta
% 0.5 of posterior st. dev. turns out to be: 
% [0.04 0.04 0.1 0.02 0.03 0.02 0.06 0.02]

if strcmp(update_T,'URW')
%     delta.T = sqrt(3)*[0.04 0.04 0.1 0.02 0.03 0.02 0.06 0.02];
    delta.T = sqrt(3)*[0.04 0.04 0.05 0.02 0.03 0.02 0.03 0.02];
elseif strcmp(update_T,'NRW')
%     delta.T =  [0.04 0.04 0.07 0.02 0.03 0.02 0.07 0.02];
%     delta.T = [0.04 0.04 0.05 0.02 0.03 0.02 0.03 0.02];
    delta.T = [0.1 0.04 0.05 0.1 0.1 0.035 0.05 0.12];
else
    delta.T = 0.1*ones(D-1,1);
end
deltaT = delta.T;

if strcmp(update_N,'U')
% if (sc == 1)
%     delta.N = [7, 10] + 0.5; %0.5 added to have a correct dicrete uniform distribution after rounding
% else
    delta.N = [20/sc, 100/sc] + 0.5;
%     delta.N = [50/sc, 140/sc] + 0.5;
    deltaN = delta.N;
% end
end


%% MH Algorithm
M = 10000;%50000;
NN = zeros(2,T,M);
sample = zeros(M,9);
accept = zeros(M,T+T+D-1);
mean_A = zeros(M,1);
theta_init(9) = 30000; % fixed for debugging
theta = theta_init;


tic
for ii = 1:M
    % Update the parameters in the model using function "updateparam": 
    % Set parameter values and log(likelihood) value of current state to be the output from
    % the MH step:
    
    if (mod(ii,1000)==0)
        fprintf('MH iter = %i\n',ii); toc;
    end
    if strcmp(update_T,'NRW') 
        [N, theta, acc, a_sum] = BKM_update_NRW(N, theta, prior, delta, y, m, f, stdT, update_N, logfact);
    else
        [N, theta, acc, a_sum] = BKM_update_URW(N, theta, prior, delta, y, m, f, stdT, update_N, logfact);
    end
    NN(:,:,ii) = N;
    sample(ii,:)= theta; 
    accept(ii,:) = acc; 
    mean_A(ii) = a_sum;
end
time_sampl = toc;
% accept = accept/(2*T+D-1);
mean_A = mean_A/(2*T+D-1);
name = ['BKM_DA_results_',update_T,'_',update_N,'_2.mat'];
save(name,'sample','NN','accept','theta_init','prior','delta','time_sampl', 'accept', 'mean_A');


BurnIn = 0
%%
if plot_on
    figure(1)
    set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
    for ii = 1:9
        subplot(3,3,ii)
        hold on
        plot(sample(BurnIn+1:M,ii))
%         plot(theta_init(ii)+0*sample(BurnIn+1:M,ii),'r')        
        hold off
        title(params{ii})
    end

    figure(10)
    set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
    for ii = 1:9
        subplot(3,3,ii)
        autocorr(sample(:,ii),40)
        title(params{ii})
    end   
    
    figure(110)
    set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
    for ii = 1:9
        subplot(3,3,ii)
        bar(acf(sample(:,ii),40))
        title(params{ii})
    end
    
    figure(2)
    set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
    for ii = 1:9
        subplot(3,3,ii)
        hold on
        plot(squeeze(NN(1,4*(ii-3)+9,:)))
%         plot(N1(4*(ii-3)+9) + 0*squeeze(NN(1,4*(ii-3)+9,:)),'r')
        hold off
        title(['N1(',num2str(4*(ii-3)+9),')'])
    end

    figure(22)
    set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
    for ii = 1:9
        subplot(3,3,ii)
        autocorr(squeeze(NN(1,4*(ii-3)+9,:)),40)
        title(['N1(',num2str(4*(ii-3)+9),')'])
    end

    figure(3)
    set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
    for ii = 1:9
        subplot(3,3,ii)
        hold on
        plot(squeeze(NN(2,4*(ii-3)+9,:)))
%         plot(Na(4*(ii-3)+9) + 0*squeeze(NN(2,4*(ii-3)+9,:)),'r')
        hold off
        title(['Na(',num2str(4*(ii-3)+9),')'])
    end

    figure(33)
    set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
    for ii = 1:9
        subplot(3,3,ii)
        autocorr(squeeze(NN(2,4*(ii-3)+9,:)),40)
        title(['Na(',num2str(4*(ii-3)+9),')'])
    end


    figure(33)
    set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
    for ii = 1:9
        subplot(3,3,ii)
        bar(acf(squeeze(NN(2,4*(ii-3)+9,:))',40))
        title(['Na(',num2str(4*(ii-3)+9),')'])
    end    
    
    figure(4)   
    subplot(2,1,1)
    hold on
    plot(mean(NN(1,:,:),3))
    hold off

    subplot(2,1,2)
    hold on
    plot(mean(NN(2,:,:),3))
    hold off
    
    
    figure(5)
    bar(sum((accept>0))/M)

end