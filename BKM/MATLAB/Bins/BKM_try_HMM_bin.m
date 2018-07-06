function Results = BKM_try_HMM_bin(N_bin, M, BurnIn, save_on)
    % clear all
    % close all
    % N_bin = 100;
    % M = 20000;
    % BurnIn = 20000;
    
    fprintf('*** BKM_HMM_bin M=%i, BurnIn=%i, N_bin = %i ***\n', M, BurnIn,N_bin);

    plot_on = false;

    addpath(genpath('../'));

    sc = 1;
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

    %%  Bins' midpoints  
%     switch N_bin 
%         case 30 % min = 0, max = 869
%             bin_size = 29; 
%         case 20 % min = 0, max = 859
%             bin_size = 43; 
%         case 10 % min = 0, max = 869
%             bin_size = 87; 
%     end
    switch N_bin 
        case 30 % min = 0, max = 689
            bin_size = 23; 
        case 20 % min = 0, max = 659
            bin_size = 33; 
        case 10 % min = 0, max = 669
            bin_size = 67; 
        case 100
            bin_size = 7; 
    end    
%     N_max = 66 * 10;
%     bin_size = N_max/N_bin + 1; 
    bin = 0.5*(bin_size*(2*(0:(N_bin-1))+1)-1)';

    logfact = @(xx) sum(log(1:1:xx));
    logfact = arrayfun(logfact,(0:7000)') ;
%     logfact = arrayfun(logfact,(0:10000)') ;

    %% Set the proposals
    % for the parameters
    update_T = 'NRW'; % 'NRW' or 'URW'

    % step sizes 
    % given step size delta, std for URW is delta/sqrt(3), for NRW 1*delta
    % 0.5 of posterior st. dev. turns out to be: 
    % [0.04 0.04 0.1 0.02 0.03 0.02 0.06 0.02]
    % from JAGS [0.04 0.04 0.05 0.02 0.03 0.02 0.03 0.02]
    if strcmp(update_T,'URW')
    %     delta.T = sqrt(3)*[0.04 0.04 0.1 0.02 0.03 0.02 0.06 0.02];
        delta.T = sqrt(3)*[0.04 0.04 0.05 0.02 0.03 0.02 0.03 0.02];
    elseif strcmp(update_T,'NRW')
        delta.T = [0.1 0.04 0.05 0.1 0.1 0.035 0.05 0.12];
    else
        delta.T = 0.1*ones(D-1,1);
    end
    deltaT = delta.T;

    if (N_bin == 10)
        delta.N = 90 + 0.5; 
%         delta.N = 80 + 0.5; 
%         delta.N = 60 + 0.5; 
    elseif (N_bin == 20)
        delta.N = 100 + 0.5; 
%         delta.N = 70 + 0.5; 
    elseif (N_bin == 30)
          delta.N = 100 + 0.5; 
%           delta.N = 80 + 0.5; 
%           delta.N = 60 + 0.5; 
          %delta.N = 40 + 0.5; 
    elseif (N_bin == 100)
        delta.N = 10 + 0.5; 
    end

    oldlikhood = BKM_calclikhood_HMM_bin(Na, theta_init, y, m, f, stdT, prior.N, bin, logfact);

    % M = 10000;
    % BurnIn = 5000; %1000;
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
        end
%         [N, theta, acc, a_sum] = BKM_update_HMM_bin(N, theta, prior, delta, y, m, f, stdT, bin, logfact);
        [N, theta, acc, a_sum] = BKM_update_HMM_bin_v2(N, theta, prior, delta, y, m, f, stdT, bin, logfact);
        if (ii>0)
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
%         name = ['../Results/BurnIn_',num2str(BurnIn),'/BKM_bin_Nbin',num2str(N_bin),'.mat'];
%         name = ['Results/BurnIn_',num2str(BurnIn),'/BKM_bin_Nbin',num2str(N_bin),'.mat'];
        name = ['/home/aba228/Documents/BKM/BKM_bin_Nbin',num2str(N_bin),'_v2.mat'];
        save(name,'delta','prior','theta_init','Na','NN','Theta',...
            'accept','mean_A','mean_accept','time_sampl','bin_size');
    end
end