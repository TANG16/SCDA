function Results = BKM_try_HMM_bin_NEW(N_bin, M, BurnIn, save_on)
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
%        alpha1        alphaa        alphal        alphar         beta1         betaa         betal         betar          sigy 
%     0.5490783     1.5673245    -4.5771131    -1.1760594    -0.1907766    -0.2472439    -0.3636677    -0.3421766 30440.2276841 
    alpha1 = 0.5490783;%1;
    alphaa = 1.5673245 ; %2;
    alphar = -1.1760594; %-2;
    alphal = -4.5771131 ; %-4;
    beta1 = -0.1907766; %-2;
    betaa = -0.2472439 ; %0.1;
    betar = -0.3421766; %-0.7;
    betal = -0.3636677  ; %-0.3;
    sigy = 30440;%1;

    params = {'alpha1', 'alphaa', 'alphar', 'alphal', ...
        'beta1', 'betaa', 'betar', 'betal',...
        'sigy'};

    theta_init = [alpha1, alphaa, alphar, alphal, beta1, betaa, betar, betal, sigy];
    [phi1, phia, rho, lambda] = BKM_covariates(theta_init,f,stdT);  

    D = size(theta_init,2);
    prior.N = [200/sc 2000/sc 0.5];
    prior.S = [0.001,0.001];
    prior.T_mu = 0*ones(D-1,1);
    prior.T_sigma2 = 100*ones(D-1,1);

    %%  Bins' midpoints  
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
    logfact = arrayfun(logfact,(0:10000)') ;

    %% Set the proposals
    % step sizes 
    % given step size delta, std for URW is delta/sqrt(3), for NRW 1*delta
    % 0.5 of posterior st. dev. turns out to be: 
    % [0.04 0.04 0.1 0.02 0.03 0.02 0.06 0.02]
    % from JAGS [0.04 0.04 0.05 0.02 0.03 0.02 0.03 0.02]
    delta.T = [0.1 0.04 0.05 0.1 0.1 0.035 0.05 0.12];
    
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


    N = Na;
    theta = theta_init;
    NN = zeros(T,M);
    Theta = zeros(M,9);
    accept = zeros(M,T+D-1);

    tic
    % profile on
    for ii = -BurnIn:M
        % Update the parameters in the model using function "updateparam": 
        % Set parameter values and log(likelihood) value of current state to be the output from
        % the MH step:

        if (mod(ii,1000)==0)
            fprintf('MH iter = %i\n',ii); toc;
        end
        [N, theta, acc] = BKM_update_HMM_bin_NEW(N, theta, prior, delta, y, m, f, stdT, bin, logfact);
        if (ii>0)
            NN(:,ii) = N;
            Theta(ii,:)= theta; 
            accept(ii,:) = acc; 
        end
    end
    time_sampl = toc; 

    Results.NN = NN;
    Results.Theta = Theta;
    Results.accept = accept;
    Results.time_sampl = time_sampl;
    
    if save_on
%         name = ['../Results/BurnIn_',num2str(BurnIn),'/BKM_bin_Nbin',num2str(N_bin),'.mat'];
%         name = ['Results/BurnIn_',num2str(BurnIn),'/BKM_bin_Nbin',num2str(N_bin),'.mat'];
%         name = ['/home/aba228/Documents/BKM/BKM_bin_Nbin',num2str(N_bin),'_v2.mat'];
        name = ['Results/BKM_bin_Nbin',num2str(N_bin),'.mat'];
        save(name,'delta','prior','theta_init','Na','NN','Theta',...
            'accept','mean_A','time_sampl','bin_size');
    end
end