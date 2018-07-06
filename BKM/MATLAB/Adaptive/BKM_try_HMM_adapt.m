function Results = BKM_try_HMM_adapt(N_q, M, BurnIn, save_on)
%     clear all
%     close all
%     N_q = 30;
%     M = 10000;
%     BurnIn = 20000; %50000;
    
    fprintf('*** BKM_HMM_adapt M=%i, BurnIn=%i, N_q = %i ***\n', M, BurnIn,N_q);
    
    addpath(genpath('../'));

    plot_on = false;

    sc = 1;
    [y, T, time, stdT, f, m, T1, T2] = BKM_Data_HMM(sc);

    Na = round( [1000, 1000, 1092.23, 1100.01, 1234.32, 1460.85, 1570.38, 1819.79,...
        1391.27, 1507.60, 1541.44, 1631.21, 1628.60, 1609.33, 1801.68, 1809.08, 1754.74,...
        1779.48, 1699.13, 1681.39, 1610.46, 1918.45, 1717.07, 1415.69, 1229.02, 1082.02,...
        1096.61, 1045.84, 1137.03, 981.1, 647.67, 992.65, 968.62, 926.83, 952.96, 865.64]/sc);

    alpha1 = 0.5;
    alphaa = 2;
    alphar = -1; 
    alphal = -4; 
    beta1 = -0.2;
    betaa = -0.25;
    betar = -0.14;
    betal = -0.37;
    sigy = 1; 

    
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
    priorN = prior.N;


    % Quanatiles
%     N_q = 20; 
    qu = (0:(N_q-1))/N_q;
    qu_mid = qu + qu(2)/2; 
    mid = norminv(qu_mid);  

    logfact_fun = @(xx) sum(log(1:1:xx));
%     logfact = arrayfun(logfact_fun,0:7000) ;
    logfact = arrayfun(logfact_fun,0:10000) ;
    %% Set the proposals
    % for the parameters
    update_T = 'NRW';  

    % step sizes 
    % given step size delta, std for URW is delta/sqrt(3), for NRW 1*delta
    % 0.5 of posterior st. dev. turns out to be: 
    % [0.04 0.04 0.1 0.02 0.03 0.02 0.06 0.02]
    % from JAGS [0.04 0.04 0.05 0.02 0.03 0.02 0.03 0.02]
    delta.T = [0.1 0.04 0.05 0.1 0.1 0.035 0.05 0.12];
    % delta.T = [0.04 0.04 0.05 0.02 0.03 0.02 0.03 0.02];

    deltaT = delta.T;
    %     delta.N = 130 + 0.5;  
    % delta.N = 80 + 0.5;  %mean(mean_accept(1:T)) = 0.3952
    % delta.N = 90 + 0.5;  %mean(mean_accept(1:T)) = 0.3770
    if (N_q == 20)
        delta.N = 95 + 0.5;  %mean(mean_accept(1:T)) =  0.3756
    elseif (N_q == 10)
    %     delta.N = 65 + 0.5;  %mean(mean_accept(1:T)) = 0.4376
        delta.N = 75 + 0.5;  %mean(mean_accept(1:T)) = 0.4141
    elseif (N_q == 30)
        delta.N = 95 + 0.5;  %mean(mean_accept(1:T)) =  0.3672
    else
        delta.N = 15 + 0.5;
    end

    deltaN = delta.N;

    oldlikhood = BKM_calclikhood_HMM_adapt(Na, theta_init, y, m, f, stdT, prior.N, mid, logfact);

%     M = 10000;
%     BurnIn = 1000;
    N = Na;
    theta = theta_init;
    % theta(9) = 30000;
    NN = zeros(T,M);
    Theta = zeros(M,9);
    accept = zeros(M,T+D-1);
    mean_A = zeros(M,1);

    tic
    % profile on
    for ii = -BurnIn:M
        % Update the parameters in the model using function "updateparam": 
        % Set parameter values and log(likelihood) value of current state to be the output from
        % the MH step:

        if (mod(ii,1000)==0)
            fprintf('MH iter = %i\n',ii); toc;
            fprintf('Sigma2 = %6.4f \n',theta(:,end));
        end
%         [N, theta, acc, a_sum] = BKM_update_HMM_adapt(N, theta, prior, delta, y, m, f, stdT, mid, logfact);
        [N, theta, acc, a_sum] = BKM_update_HMM_adapt_v2(N, theta, prior, delta, y, m, f, stdT, mid, logfact);
        if (ii > 0)
            NN(:,ii) = N;
            Theta(ii,:)= theta; 
            accept(ii,:) = acc; 
            mean_A(ii) = a_sum;
        end
    end
    time_sampl = toc; 
    mean_A = mean_A/(T+D-1);
    mean_accept = mean(accept);
    
    Results.NN = NN;
    Results.Theta = Theta;
    Results.accept = accept;
    Results.mean_A = mean_A;
    Results.time_sampl = time_sampl;
    if save_on
% %         name = ['../Results/BurnIn_',num2str(BurnIn),'/BKM_adapt_Nq',num2str(N_q),'.mat'];
%         name = ['Results/BurnIn_',num2str(BurnIn),'/BKM_adapt_Nq',num2str(N_q),'.mat'];
        name = ['/home/aba228/Documents/BKM/BKM_adapt_Nq',num2str(N_q),'_v2.mat'];
        save(name,'delta','prior','theta_init','Na','NN','Theta',...
            'accept','mean_A','mean_accept','time_sampl');
    end
end