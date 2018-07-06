function [NNa, Theta] = BKM_try_DA_generate(M, BurnIn, save_on, sd)% clear all
    % close all
    % M=10000; BurnIn=1000;

    s = RandStream('mt19937ar','Seed',sd);
    RandStream.setGlobalStream(s); 

    fprintf('*** BKM_DA generate ***\n');

    plot_on = false;

    sc = 1;
    [~, T, time, stdT, f, m, T1, T2] = BKM_Data_HMM(sc);

    [N_true, y] = BKM_genarate;
 
    N = N_true;
    
    params = {'alpha1', 'alphaa', 'alphar', 'alphal', ...
        'beta1', 'betaa', 'betar', 'betal',...
        'sigy'};

    theta_init = [0.5237    1.4728   -1.0636   -4.5935 ...
        -0.1528   -0.2435   -0.3597   -0.3546 ...
           6.0304e+04];
       
    theta = theta_init;
       
    [phi1, phia, rho, ~] = BKM_covariates(theta,f,stdT);

    D = size(theta,2);

    prior.N = [200/sc 2000/sc 0.5];
    prior.S = [0.001,0.001];
    prior.T_mu = 0*ones(D-1,1);
    prior.T_sigma2 = 100*ones(D-1,1);

    priorN = prior.N;

    logfact = @(xx) sum(log(1:1:xx));
    logfact = arrayfun(logfact,0:7000);
%     logfact = arrayfun(logfact,0:100000);

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
        delta.T = sqrt(3)*[0.04 0.04 0.05 0.02 0.03 0.02 0.03 0.02];
    elseif strcmp(update_T,'NRW')
        delta.T = [0.1 0.04 0.05 0.1 0.1 0.035 0.05 0.12];
    %   0.2842    0.3067    0.2388    0.3783    0.2618    0.3188    0.2730    0.3647    
    else
        delta.T = 0.1*ones(D-1,1);
    end
    deltaT = delta.T;

    if strcmp(update_N,'U')
    % if (sc == 1)
    %     delta.N = [7, 10] + 0.5; %0.5 added to have a correct dicrete uniform distribution after rounding
    % else
        delta.N = [60/sc, 100/sc] + 0.5;
        deltaN = delta.N;
    % end
    end

    %% MH Algorithm
    % M = 10000;%50000;
    % BurnIn = 0; %1000;
%          N = [N1;Na];

    NN = zeros(2,T,M);
    Theta = zeros(M,9);
    accept = zeros(M,T+T+D-1);
    mean_A = zeros(M,1);
    % theta_init(9) = 30000; % fixed for debugging
    theta = theta_init;

    tic
    for ii = -BurnIn:M
        % Update the parameters in the model using function "updateparam": 
        % Set parameter values and log(likelihood) value of current state to be the output from
        % the MH step:

        if (mod(ii,1000)==0)
            fprintf('MH iter = %i\n',ii); toc;
%             fprintf('Sigma2 = %6.4f \n',theta(:,end));            
        end
%         if strcmp(update_T,'NRW') 
%             [N, theta, acc, a_sum] = BKM_update_NRW(N, theta, prior, delta, y, m, f, stdT, update_N, logfact);
            [N, theta, acc, a_sum] = BKM_update_NRW_v2(N, theta, prior, delta, y, m, f, stdT, update_N, logfact);
%             [N, theta, acc, a_sum] = BKM_update_NRW_debug(N, theta, prior, delta, y, m, f, stdT, update_N, logfact);
%         else
%             [N, theta, acc, a_sum] = BKM_update_URW(N, theta, prior, delta, y, m, f, stdT, update_N, logfact);
%         end
        if (ii > 0)
            NN(:,:,ii) = N;
            Theta(ii,:)= theta; 
            accept(ii,:) = acc; 
            mean_A(ii) = a_sum;
        end
    end
    time_sampl = toc;
    % accept = accept/(2*T+D-1);
    mean_A = mean_A/(2*T+D-1);
    
    NNa = squeeze(NN(2,:,:));
    
    Results.NN = NN;
    Results.Theta = Theta;
    Results.accept = accept;
    Results.mean_A = mean_A;
    Results.time_sampl = time_sampl;

    if save_on
%         name = ['Results/BurnIn_',num2str(BurnIn),'/BKM_DA_',update_T,'_',update_N,'.mat'];
        name = ['/home/aba228/Documents/BKM/BKM_DA_',update_T,'_',update_N,'_v2_long.mat'];
        save(name,'Theta','NN','accept','theta_init','prior','delta','time_sampl', ...
            'accept', 'mean_A');
    end

    if false
        figure(1)
        subplot(2,1,1)
        hold on
        plot(squeeze(mean(NN(1,:,:),3)))
        plot(N_true(1,:),'r')
        hold off

        subplot(2,1,2)
        hold on
        plot(squeeze(mean(NN(2,:,:),3)))
        plot(N_true(2,:),'r')
        hold off


        figure(2)
        for ii = 1:9
            subplot(3,3,ii)
            hold on
            plot(Theta(:,ii))
            plot(theta_init(ii) + 0*Theta(:,ii),'r')
            hold off
            hold off
        end


        figure(3)
        for ii = 1:9
            subplot(3,3,ii)
            hold on
            plot(squeeze(NN(2,3*ii,:)))
            plot(N_true(2,3*ii) + 0*squeeze(NN(2,3*ii,:)),'r')
            hold off
            hold off
        end
    end

end