function [NN_HMM, Theta_HMM] = BKM_try_HMM_generate(M, BurnIn, save_on, sd)
    % close all
    % M=10000; BurnIn=1000;

s = RandStream('mt19937ar','Seed',sd);
RandStream.setGlobalStream(s); 

    fprintf('*** BKM_HMM generate ***\n');

    plot_on = false;

    sc = 1;
    [~, T, time, stdT, f, m, T1, T2] = BKM_Data_HMM(sc);

    [N_true, y] = BKM_genarate;
  
    params = {'alpha1', 'alphaa', 'alphar', 'alphal', ...
        'beta1', 'betaa', 'betar', 'betal',...
        'sigy'};

    theta_init = [0.5237    1.4728   -1.0636   -4.5935 ...
        -0.1528   -0.2435   -0.3597   -0.3546 ...
           6.0304e+04];
       
    theta = theta_init;
       
    [phi1, phia, rho, ~] = BKM_covariates(theta,f,stdT);

    D = size(theta,2);

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
         delta_HMM.T = sqrt(3)*[0.04 0.04 0.05 0.02 0.03 0.02 0.03 0.02];
    elseif strcmp(update_T,'NRW')
        delta_HMM.T = [0.1 0.04 0.05 0.1 0.1 0.035 0.05 0.12];
    else
        delta_HMM.T = 0.1*ones(D-1,1);
    end
    deltaT = delta_HMM.T;


    if (sc == 10)
        delta_HMM.N = 30/sc + 0.5; % std of the unif RW update 
    elseif (sc == 1)        
        delta_HMM.N = 60+0.5; %;100+0.5; %130 + 0.5; 
    else
        delta_HMM.N = 1.5;
    end
    deltaN = delta_HMM.N;


    % No bins
%     N_max = 69 * 10/sc;
    N_max = 68 * 10/sc -1;    
    IND = 0:N_max;

    logfact = @(xx) sum(log(1:1:xx));
%     logfact = arrayfun(logfact,(0:7000)') ;
    logfact = arrayfun(logfact,(0:10000)') ;
 
    N = N_true(2,:);
%     N = mean(squeeze(Results_DA.NN(2,:,:)),2)'; % DA init
    theta = theta_init;
    % theta(9) = 30000;
    NN_HMM = zeros(T,M);
    Theta_HMM = zeros(M,9);
    accept_HMM = zeros(M,T+D-1);
    mean_A_HMM = zeros(M,1);

    oldlikhood = BKM_calclikhood_HMM(N, theta, y, m, f, stdT, prior.N, N_max, logfact);

    tic
    % profile on
    for ii = -BurnIn:M
        % Update the parameters in the model using function "updateparam": 
        % Set parameter values and log(likelihood) value of current state to be the output from
        % the MH step:

        if (mod(ii,1000)==0)
            fprintf('MH iter = %i\n',ii); toc;
        end
        [N, theta, acc, a_sum] = BKM_update_HMM_v2(N, theta, prior, delta_HMM, y, m, f, stdT, N_max, logfact);
%         [N, theta, acc, a_sum] = BKM_update_HMM_debug(N, theta, prior, delta_HMM, y, m, f, stdT, N_max, logfact);
        if (ii>0)
            NN_HMM(:,ii) = N;
            Theta_HMM(ii,:)= theta; 
            accept_HMM(ii,:) = acc; 
            mean_A_HMM(ii) = a_sum;
        end
    end
    time_sampl_HMM = toc; 
    mean_A_HMM = mean_A_HMM/(T+D-1);
    mean_accept = mean(accept_HMM);
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

    Results.NN = NN_HMM;
    Results.Theta = Theta_HMM;
    Results.accept = accept_HMM;
    Results.mean_A = mean_A_HMM;
    Results.time_sampl = time_sampl_HMM;

    if save_on
%         name = ['Results/BurnIn_',num2str(BurnIn),'/BKM_exact_Nmax',num2str(N_max),'.mat'];
        name = ['/home/aba228/Documents/BKM/BKM_exact_Nmax',num2str(N_max),'_v2.mat'];
        save(name,'delta','prior','theta_init','Na','NN','Theta',...
            'accept','mean_A','mean_accept','time_sampl');
    end

if false
    figure(11)
    hold on
    plot(squeeze(mean(NN(2,:,:),3)))
    plot(mean(NN_HMM(:,:),2))
    plot(N_true(2,:),'r')
    hold off


    figure(22)
    for ii = 1:9
        subplot(3,3,ii)
        hold on
        plot(Theta(:,ii))    
        plot(Theta_HMM(:,ii))
        plot(theta_init(ii) + 0*Theta_HMM(:,ii),'r')
        title(params(ii))
        hold off
    end


    figure(33)
    for ii = 1:9
        subplot(3,3,ii)
        hold on
        plot(squeeze(NN(2,3*ii,:)))    
        plot(NN_HMM(3*ii,:))
        plot(N_true(2,3*ii) + 0*NN_HMM(3*ii,:),'r')
        title(['Na(',num2str(3*ii),')'])
        hold off
        hold off
    end

end

end
