clear all
close all


plot_on = false;

% theta = [1, 0.97, 0.15^2];
% beta = 0.05;
beta = 0.5;
mu = 2*log(beta);
theta = [mu, 0.98, 0.2^2];
theta_true = theta;
theta_init = [0, 0.95, 0.1^2];

T = 1000;
[y,h_true] = generate_SV(theta,T);

if plot_on
    plot(y)
    plot(h_true)
end

mu = theta(1);
phi = theta(2);
sigma2 = theta(3);

% P0 = sqrt(sigma2/(1-phi^2)); % unconditional st. dev. 0.6120
% bin_range = 5*P0;
bin_range = 4;
N_bin = 30; %50;
bins = linspace(-bin_range,bin_range,N_bin+1);
bin_midpoint = (bins(1:N_bin) + bins(2:N_bin+1))/2;
% bins are the demeaned volatilities: b=h-mu
stdev_y = exp((mu+bin_midpoint)/2);

delta.h = 0.1;
delta_h = delta.h;
M = 10000;
h_init = var(y)*ones(1,T); 


%% FULL DA
h = h_init;
H_DA = zeros(M,T);
accept_DA = zeros(M,1);
A_sum_DA = zeros(M,1);

tic
for ii = 1:M
    [h, acc, A_s] = update_h(y-mean(y),h, theta, delta.h);
    if (mod(ii,1000) == 1)
        toc;
    end
    H_DA(ii,:) = h;
    accept_DA(ii,1) = acc;
    A_sum_DA(ii,1) = A_s;
end
accept_DA = accept_DA/T;
A_sum_DA = A_sum_DA/T;
time_DA = toc;

% NOT SURE WHETHER THIS IS CORRECT
IF_DA = zeros(1,T);
for t = 1:T
    IF_DA(1,t) = 1 + 2*(sum(autocorr(H_DA(:,t),T))-1); 
end
ESS_DA = M./IF_DA;

% to check with the R function effectiveSize
H_DA_init = H_DA(:,1:100);
save('SV_DA_IF.mat','H_DA_init','IF_DA')
clear H_DA_init


%%  SEMI: 
% integrate out the odd h(t)'s and impute the even ones
h = h_init;
theta = theta_init;
H_HMM_v2 = zeros(M,T);
accept_HMM_v2 = zeros(M,1);
A_sum_HMM_v2 = zeros(M,1);

tic
for ii = 1:M
    [h, acc, A_s] = update_h_HMM_v2(y, h, theta, delta.h, bins, bin_midpoint);
    if (mod(ii,1000) == 1)
        toc;
    end
    H_HMM_v2(ii,:) = h;
    accept_HMM_v2(ii,1) = acc;
    A_sum_HMM_v2(ii,1) = A_s;
end
accept_HMM_v2 = accept_HMM_v2/(T/2);
A_sum_HMM_v2 = A_sum_HMM_v2/(T/2);
time_HMM = toc;

% NOT SURE WHETHER THIS IS CORRECT
IF_HMM = zeros(1,T/2);
for t = 2:2:T
    IF_HMM(1,t/2) = 1 + 2*(sum(autocorr(H_HMM_v2(1:M,t),T))-1);
end
ESS_HMM = M./IF_HMM;


%% SEMI DA shift: 
% if shift integrate out the odd h(t)'s and impute the even ones
% if not shift integrate out the even h(t)'s and impute the odd ones
h = h_init;

H_HMM_v2_shift = zeros(M,T);
accept_HMM_v2_shift = zeros(M,1);
A_sum_HMM_v2_shift = zeros(M,1);

tic
for ii = 1:M
    shift = mod(ii,2);
    [h, acc, A_s] = update_h_HMM_v2_shift(y, h, theta, delta.h, bins, bin_midpoint, shift);
    if (mod(ii,1000) == 1)
        toc;
    end
    H_HMM_v2_shift(ii,:) = h;
    accept_HMM_v2_shift(ii,1) = acc;
    A_sum_HMM_v2_shift(ii,1) = A_s;
end
accept_HMM_v2_shift = accept_HMM_v2_shift/(T/2);
A_sum_HMM_v2_shift = A_sum_HMM_v2_shift/(T/2);
time_HMM_shift = toc; 

IF_HMM_shift = zeros(1,T);

for t = 1:T
    if mod(t,2) == 0 % even
        IF_HMM_shift(1,t) = 1 + 2*(sum(autocorr(H_HMM_v2_shift(2:2:M,t),T))-1); 
    else % odd
        IF_HMM_shift(1,t) = 1 + 2*(sum(autocorr(H_HMM_v2_shift(1:2:M,t),T))-1);         
    end
end



if plot_on
    figure(1720)
    hold on
    scatter(1:T,h,'b')
    % plot(2:2:T,h(2:2:T),'b')
    plot(h_true,'k')
    % plot(H_DA(M/2,:),'g')
    plot(mean(H_HMM_v2_shift(5000:ii,:),1),'g')
    hold off


    figure(1000)
    for ii = 2:10
        subplot(3,3,ii-1) 
        hold on
        plot(H_DA(:,100*ii),'r')
        plot(H_HMM_v2(:,100*ii),'g')
    %     scatter(2:2:T,H_HMM_v2_shift(2:2:T,100*ii),'b')
        hold off
        title(['t = ',num2str(100*ii)])
    end
    suptitle('Even')


    figure(1001)
    for ii = 2:10
        subplot(3,3,ii-1) 
        autocorr(H_DA(:,100*ii))
        set(gca,'ylabel',[])    
        title(['t = ',num2str(100*ii)])
    end
    suptitle('Even DA')


    figure(1002)
    for ii = 2:10
        subplot(3,3,ii-1) 
        autocorr(H_HMM_v2(:,100*ii))
        set(gca,'ylabel',[])    
        title(['t = ',num2str(100*ii)])
    end
    suptitle('Even HMM no shift')

    figure(1003)
    for ii = 2:10
        subplot(3,3,ii-1) 
        autocorr(H_HMM_v2_shift(2:2:M,100*ii))
        set(gca,'ylabel',[])    
        title(['t = ',num2str(100*ii)])
    end
    suptitle('Even HMM with shift')


    figure(1004)
    for ii = 2:10
        subplot(3,3,ii-1) 
        autocorr(H_HMM_v2_shift(1:2:T,100*ii-1))
        set(gca,'ylabel',[])    
        title(['t = ',num2str(100*ii-1)])
    end
    suptitle('Odd HMM with shift')

end



%% SAVE
name = 'SV_results.mat';
save(name,'y','h_true','theta_true',...
    'H_DA','H_HMM_v2','H_HMM_v2_shift',...
    'accept_DA','accept_HMM_v2','accept_HMM_v2_shift',...
    'A_sum_DA','A_sum_HMM_v2','A_sum_HMM_v2_shift',...
    'time_DA','time_HMM','time_HMM_shift',...
    'IF_DA','IF_HMM','IF_HMM_shift');

