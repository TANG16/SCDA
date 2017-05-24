clear all
close all
plot_on = false;

sc=10;
[y, T, time, stdT, f, m, T1, T2] = BKM_Data_HMM(sc);

Na = round( [1000, 1000, 1092.23, 1100.01, 1234.32, 1460.85, 1570.38, 1819.79,...
    1391.27, 1507.60, 1541.44, 1631.21, 1628.60, 1609.33, 1801.68, 1809.08, 1754.74,...
    1779.48, 1699.13, 1681.39, 1610.46, 1918.45, 1717.07, 1415.69, 1229.02, 1082.02,...
    1096.61, 1045.84, 1137.03, 981.1, 647.67, 992.65, 968.62, 926.83, 952.96, 865.64]/sc);

alpha1 = 1;
alphaa = 2;
alphar = -2;
alphal = -4;
beta1 =-2;
betaa = 0.1;
betar = -0.7;
betal = -0.3;
sigy = 1;

% alpha1 = 1;
% alphaa = 2;
% alphar = -2;
% alphal = -4;
% beta1 = -0.19;
% betaa = 0.1;
% betar = -1.106;
% betal = -0.3;
% sigy = 1;

params = {'alpha1', 'alphaa', 'alphar', 'alphal', ...
    'beta1', 'betaa', 'betar', 'betal',...
    'sigy'};

theta_init = [alpha1, alphaa, alphar, alphal, beta1, betaa, betar, betal, sigy];
[phi1, phia, rho, lambda] = BKM_covariates(theta_init,f,stdT);  

D = size(theta_init,2);
prior.N = [2000/sc 0.5];
prior.S = [0.001,0.001];
prior.T_mu = 0*ones(D-1,1);
prior.T_sigma2 = 100*ones(D-1,1);

delta.T = 0.1*ones(D-1,1); % std of the normal RW update
delta.N = 30/sc + 0.5; % std of the unif RW update 

priorN = prior.N;
deltaT = delta.T;
deltaN = delta.N;

N_max= 69;
oldlikhood = BKM_calclikhood_HMM(Na, theta_init, y, m, f, stdT, prior.N, N_max);


M = 10000;
N = Na;
theta = theta_init;
NN = zeros(T,M);
sample = zeros(M,9);
accept = zeros(M,1);
% tic
profile on
for ii = 1:M
    % Update the parameters in the model using function "updateparam": 
    % Set parameter values and log(likelihood) value of current state to be the output from
    % the MH step:
    
    if (mod(ii,1000)==0)
        fprintf('MH iter = %i\n',ii); toc;
    end
    [N, theta, A] = BKM_update_HMM(N, theta, prior, delta, y, m, f, stdT, N_max);
    NN(:,ii) = N;
    sample(ii,:)= theta; 
    accept(ii) = A; 
end
% time_sampl = toc; 
% with poissrnd and binopdf 4.97 sec per draw
% with explicit formulae 0.5055 sec per draw
profile off
profile viewer
name = 'BKM_HMM_results.mat';
save(name,'sample','NN','accept','theta_init','prior','delta','time_sampl');



%%
if plot_on
    figure(1)
    set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
    for ii = 1:9
        subplot(3,3,ii)
        plot(sample(:,ii))
        title(params{ii})
    end

    figure(10)
    set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
    for ii = 1:9
        subplot(3,3,ii)
        bar(acf(sample(10001:end,ii),40))
        title(params{ii})
    end   
    
    figure(110)
    set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
    for ii = 1:9
        subplot(3,3,ii)
        autocorr(sample(10001:end,ii),40)
        title(params{ii})
    end
    

    figure(3)
    set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
    for ii = 1:9
        subplot(3,3,ii)
        plot(NN(4*(ii-3)+9,:))
        title(['Na(',num2str(4*(ii-3)+9),')'])
    end

    figure(33)
    set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
    for ii = 1:9
        subplot(3,3,ii)
        autocorr(NN(4*(ii-3)+9,:),40)
        title(['Na(',num2str(4*(ii-3)+9),')'])
    end

    figure(4)
    set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
    hold on
    plot(mean(NN,2))
    hold off
end