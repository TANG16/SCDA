close all

mean(A_theta_DA_RW)
mean(A_theta_HMM)
mean(A_theta_HMM_adapt)

mean(A_H_DA_RW)
mean(A_H_HMM)
mean(A_H_HMM_adapt)

theta_true
mean(theta_DA_RW)
mean(theta_HMM)
mean(theta_HMM_adapt)


params = {'\mu','\phi','\sigma^2','\beta','\rho'};
figure(1)
for ii = 1:5
    subplot(2,3,ii)
    hold on 
    plot(theta_DA_RW(:,ii))
%     plot(theta_true(ii) + 0*theta_DA_RW(:,ii))
    plot(theta_HMM(:,ii))
%     plot(theta_true(ii) + 0*theta_HMM(:,ii))

%     plot(theta_HMM_adapt(:,ii))
%     plot(theta_true(ii) + 0*theta_HMM_adapt(:,ii))    
    hold off
    title(params{ii})
end

mean(theta_DA_RW)
mean(theta_HMM_adapt)

figure(10)
for ii = 1:5
    subplot(2,3,ii)
    autocorr(theta_DA_RW(:,ii),200)
    title(params{ii})   
end
suptitle('DA')

figure(13)
for ii = 1:5
    subplot(2,3,ii)
    autocorr(theta_HMM(:,ii),200)
    title(params{ii})   
end
suptitle('HMM')

Lags = [20, 200, 200, 20, 200];
figure(14)
for ii = 1:5
    subplot(2,3,ii)
    hold on
    plot(0:Lags(ii),autocorr(theta_DA_RW(:,ii),Lags(ii)))
    plot(0:Lags(ii),autocorr(theta_HMM(:,ii),Lags(ii)))    
    plot(0:Lags(ii),autocorr(theta_HMM_adapt(:,ii),Lags(ii)))
    hold off
    title(params{ii})   
end
legend('DA',['HMM bin',num2str(N_bin)],['HMM adapt',num2str(N_q)])



figure(2)
hold on
plot(1:length(mean_H_DA_RW),mean_H_DA_RW)
plot(2:2:length(mean_H_HMM),mean_H_HMM(2:2:length(mean_H_HMM)))
plot(2:2:length(mean_H_HMM_adapt),mean_H_HMM_adapt(2:2:length(mean_H_HMM_adapt)))
% plot(1:length(h_true),h_true)
hold off

if (ind0 == 0)
    figure(3)
    k = floor(size(H_subset_DA_RW,2)/9);    
    for ii = 1:9
        subplot(3,3,ii)
        hold on
        plot(H_subset_DA_RW(:,ii*k))
%         plot(h_true(ind_h_sel(3*ii)) + 0*H_subset_DA_RW(:,ii*3))
        plot(H_subset_HMM(:,ii*k))
%         plot(h_true(ind_h_sel(3*ii)) + 0*H_subset_HMM(:,ii*3))
        plot(H_subset_HMM_adapt(:,ii*k))
%         plot(h_true(ind_h_sel(3*ii)) + 0*H_subset_HMM_adapt(:,ii*3))
        hold off
        title(['Na(',num2str(ii*k*25),')'])        
    end
else
    figure(33)
    for ii = 1:9
        subplot(3,3,ii)
        hold on
        plot(H_subset_DA_RW(:,ii*10))
        plot(H_subset_HMM(:,ii*10))
        plot(H_subset_HMM_adapt(:,ii*10))
        hold off
    end
end

if (ind0 == 0)
    figure(12)
    for ii = 1:9
        subplot(3,3,ii)
        autocorr(H_subset_DA_RW(:,ii*5),200)
    end
else
    figure(12)
%     k = floor(T/9);
    k = floor(size(H_subset_DA_RW,2)/9);
    for ii = 1:9
        subplot(3,3,ii)
        hold on
        plot(0:200,autocorr(H_subset_DA_RW(:,ii*k),200))
        plot(0:200,autocorr(H_subset_HMM(:,ii*k),200))
        plot(0:200,autocorr(H_subset_HMM_adapt(:,ii*k),200))
        title(['Na(',num2str(ii*k*25),')'])
        hold off
    end 
    
    figure(13)
    k = floor(size(H_subset_DA_RW,2)/25);
    for ii = 1:25
        subplot(5,5,ii)
        hold on
        plot(0:200,autocorr(H_subset_DA_RW(:,ii*k),200))
        plot(0:200,autocorr(H_subset_HMM(:,ii*k),200))
        plot(0:200,autocorr(H_subset_HMM_adapt(:,ii*k),200))
        title(['Na(',num2str(ii*k*25),')'])
        hold off
    end 
end