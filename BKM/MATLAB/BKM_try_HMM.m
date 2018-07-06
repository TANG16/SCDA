function Results = BKM_try_HMM(M, BurnIn, save_on)
% clear all; close all; M = 10000; BurnIn = 10000; save_on = false;
% A version without bins: each value is itegrated out 
% (or, equiv., bin_size = 1)
    fprintf('*** BKM_HMM_exact ***\n');

    % clear all
    % close all
    plot_on = false;

    sc = 1;
    [y, T, time, stdT, f, m, T1, T2] = BKM_Data_HMM(sc);

    Na_init = round( [1000, 1000, 1092.23, 1100.01, 1234.32, 1460.85, 1570.38, 1819.79,...
        1391.27, 1507.60, 1541.44, 1631.21, 1628.60, 1609.33, 1801.68, 1809.08, 1754.74,...
        1779.48, 1699.13, 1681.39, 1610.46, 1918.45, 1717.07, 1415.69, 1229.02, 1082.02,...
        1096.61, 1045.84, 1137.03, 981.1, 647.67, 992.65, 968.62, 926.83, 952.96, 865.64]/sc);


    alpha1 = 0.5;
    alphaa = 1.5;
    alphar = -1.2; 
    alphal = -4.6;
    beta1 = -0.16;
    betaa = -0.20;
    betar = -0.3;
    betal = -0.4;
    sigy = 3e+04;
%%     
%     alpha1 = 1;
%     alphaa = 2;
%     alphar = -2; 
%     if (sc == 1)
%         alphal = -4;
%     else
%         alphal = -1;
%     end
%     beta1 =-2;
%     betaa = 0.1;
%     betar = -0.7;
%     betal = -0.3;
%     sigy = 1;

    % alpha1 = 1;
    % alphaa = 2;
    % alphar = -2;
    % alphal = -4;
    % beta1 = -0.19;
    % betaa = 0.1;
    % betar = -1.106;
    % betal = -0.3;
    % sigy = 1;
%%
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


    %% Set the proposals
    % for the parameters
    update_T = 'NRW'; % 'NRW' or 'URW'

    % step sizes 
    % given step size delta, std for URW is delta/sqrt(3), for NRW 1*delta
    % 0.5 of posterior st. dev. turns out to be: 
    % [0.04 0.04 0.1 0.02 0.03 0.02 0.06 0.02]
    % from JAGS [0.04 0.04 0.05 0.02 0.03 0.02 0.03 0.02]
    if strcmp(update_T,'URW')
         delta.T = sqrt(3)*[0.04 0.04 0.05 0.02 0.03 0.02 0.03 0.02];
    elseif strcmp(update_T,'NRW')
        delta.T = [0.1 0.04 0.05 0.1 0.1 0.035 0.05 0.12];
    else
        delta.T = 0.1*ones(D-1,1);
    end
    deltaT = delta.T;


    if (sc == 10)
        delta.N = 30/sc + 0.5; % std of the unif RW update 
    elseif (sc == 1)        
        delta.N = 60+0.5; %;100+0.5; %130 + 0.5; 
    else
        delta.N = 1.5;
    end
    deltaN = delta.N;


    % No bins
%     N_max = 69 * 10/sc;
    N_max = 68 * 10/sc -1;    
    IND = 0:N_max;

    logfact = @(xx) sum(log(1:1:xx));
%     logfact = arrayfun(logfact,(0:7000)') ;
    logfact = arrayfun(logfact,(0:10000)') ;

    oldlikhood = BKM_calclikhood_HMM(Na_init, theta_init, y, m, f, stdT, prior.N, N_max, logfact);
 
%     M = 10000;
%     BurnIn = 0; %1000;
    N = Na_init;
%     N = mean(squeeze(Results_DA.NN(2,:,:)),2)'; % DA init
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
%         [N, theta, acc, a_sum] = BKM_update_HMM(N, theta, prior, delta, y, m, f, stdT, N_max, logfact);
        [N, theta, acc, a_sum] = BKM_update_HMM_debug(N, theta, prior, delta, y, m, f, stdT, N_max, logfact);
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
mean_accept_Na = mean_accept(1:T);
mean(mean_accept_Na)
mean_accept_theta = mean_accept((T+1):end);
    % name = 'BKM_HMM_results_sigma2.mat';
    % save(name,'sample','NN','accept','theta_init','prior','delta','time_sampl', 'accept', 'mean_A');

    % with poissrnd and binopdf 4.97 sec per draw
    % with explicit formulae 0.5055 sec per draw
    % with explicit formula with logfact as array 0.0170 sec per draw
    % profile off
    % profile viewer

    Results.NN = NN;
    Results.Theta = Theta;
    Results.accept = accept;
    Results.mean_A = mean_A;
    Results.time_sampl = time_sampl;

    if save_on
%         name = ['Results/BurnIn_',num2str(BurnIn),'/BKM_exact_Nmax',num2str(N_max),'.mat'];
        name = ['/home/aba228/Documents/BKM/BKM_exact_Nmax',num2str(N_max),'_v2.mat'];
        save(name,'delta','prior','theta_init','Na','NN','Theta',...
            'accept','mean_A','mean_accept','time_sampl');
    end
end

 if false
    NN_old = NN;
    Theta_old = Theta;
 
    figure(10)
    set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
    for ii = 1:9
        subplot(3,3,ii)
        hold on
%         plot(Results_DA.Theta(1:M,ii))
%         plot(Theta(1:M,ii))
        plot(Results_DA2.Theta(1:M,ii))       
%         plot(theta_init(ii)+0*sample(1:M,ii),'r')        
        hold off
        title(params{ii})
    end
    
    
    figure(30)
    set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
    for ii = 1:36%9
%         subplot(3,3,ii)
        subplot(6,6,ii)
        hold on
%         plot(squeeze(Results_DA.NN(2,ii,1:M)))
%         plot(NN(ii,1:M))
        plot(squeeze(Results_DA2.NN(2,ii,1:M)))
        hold off
        title(['Na(',num2str(ii),')'])
    end
    
    figure(40)
    set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
    hold on
%     plot(mean(squeeze(Results_DA.NN(2,:,:)),2))
%     plot(mean(NN,2))
    plot(mean(squeeze(Results_DA2.NN(2,:,:)),2))
    hold off   
 end