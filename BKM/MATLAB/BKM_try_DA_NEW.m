function Results = BKM_try_DA_NEW(M, BurnIn, save_on)% clear all
    % close all
    M=100000; BurnIn=10000;
    
    fprintf('*** BKM_DA ***\n');

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
 
%  RN1 = load('C:/Users/ab507t/Desktop/N1_DA_mean.mat');
%  N1 = round(RN1.meanN1');
 
    N = [N1;Na];

    
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
    
    alpha1 = 1;
    alphaa = 2;
    alphar =-2;
    alphal = -4;
    beta1 = -2;
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

    prior.N = [200/sc 2000/sc 0.5];
    prior.S = [0.001,0.001];
    prior.T_mu = 0*ones(D-1,1);
    prior.T_sigma2 = 100*ones(D-1,1);


    logfact = @(xx) sum(log(1:1:xx));
    logfact = arrayfun(logfact,0:10000);

    oldlikhood = BKM_calclikhood(N, theta, y, m, f, stdT, prior.N, logfact);

    %% Set the proposals
%   mean(mean(accept(:,1:36)))  mean(mean(accept(:,37:72))) mean(accept(:,73:80))
%     delta.T = [0.1 0.04 0.05 0.1 0.1 0.035 0.05 0.12]; %    0.2431    0.2575    0.2002    0.3044    0.2211    0.2574    0.2241    0.2996
%     delta.T = [0.07 0.025 0.02 0.1 0.07 0.025 0.03 0.10]; %    0.3260    0.3708    0.4165    0.3031    0.3005    0.3411    0.3465    0.3454
    delta.T = [0.07 0.03 0.035 0.1 0.07 0.025 0.03 0.10];  %    0.3256    0.3298    0.2730    0.3033    0.3002    0.3387    0.3415    0.3505
%     delta.N = [60/sc, 100/sc] + 0.5; % 0.27 0.16
%     delta.N = [55/sc, 70/sc] + 0.5; %0.29 0.2238
%       delta.N = [50/sc, 60/sc] + 0.5; % 0.3170 0.2586
      delta.N = [50/sc, 40/sc] + 0.5;  % 0.3149 0.3607 

    
    
%     M = M/2;    
    %% MH Algorithm
    NN = zeros(2,T,M);
    Theta = zeros(M,9);
    accept = zeros(M,T+T+D-1);
    theta = theta_init;

    tic
    for ii = -BurnIn:M
        % Update the parameters in the model using function "updateparam": 
        % Set parameter values and log(likelihood) value of current state to be the output from
        % the MH step:

        if (mod(ii,1000)==0)
            fprintf('MH iter = %i\n',ii); toc;
        end
        
        [N, theta, acc] = BKM_update_NRW_NEW(N, theta, prior, delta, y, m, f, stdT, logfact);

        if (ii > 0)
            NN(:,:,ii) = N;
            Theta(ii,:)= theta; 
            accept(ii,:) = acc; 
        end
    end
    time_sampl = toc;
    
    Na = squeeze(NN(2,:,:));
    Results.NN = NN;
    Results.Theta = Theta;
    Results.accept = mean(accept);
    Results.time_sampl = time_sampl;


% R_output = load('C:\Users\ab507t\Desktop\BKM_DA_output.mat');
% Na_mean_R = mean(R_output.matDA(:,1:36));

    if save_on
%         name = ['Results/BurnIn_',num2str(BurnIn),'/BKM_DA_',update_T,'_',update_N,'.mat'];
%         name = ['/home/aba228/Documents/BKM/BKM_DA_',update_T,'_',update_N,'_v2_long.mat'];
        name = ['Results/BKM_DA_M',num2str(M),'.mat'];
        save(name,'Theta','NN','accept','theta_init','prior','delta','time_sampl', ...
            'accept');
    end
end