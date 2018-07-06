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
%     plot(theta_HMM(:,ii))
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

figure(14)
for ii = 1:5
    subplot(2,3,ii)
    autocorr(theta_HMM_adapt(:,ii),200)
    title(params{ii})   
end
suptitle('HMM adapt')

figure(2)
hold on
plot(1:length(mean_H_DA_RW),mean_H_DA_RW)
plot(2:2:length(mean_H_HMM),mean_H_HMM(2:2:length(mean_H_HMM)))
plot(2:2:length(mean_H_HMM_adapt),mean_H_HMM_adapt(2:2:length(mean_H_HMM_adapt)))
% plot(1:length(h_true),h_true)
hold off

if (ind0 == 0)
    figure(3)
    for ii = 1:9
        subplot(3,3,ii)
        hold on
        plot(H_subset_DA_RW(:,ii*3))
%         plot(h_true(ind_h_sel(3*ii)) + 0*H_subset_DA_RW(:,ii*3))
        plot(H_subset_HMM(:,ii*3))
%         plot(h_true(ind_h_sel(3*ii)) + 0*H_subset_HMM(:,ii*3))
%         plot(H_subset_HMM_adapt(:,ii*3))
%         plot(h_true(ind_h_sel(3*ii)) + 0*H_subset_HMM_adapt(:,ii*3))
        hold off
    end
else
    figure(33)
    for ii = 1:9
        subplot(3,3,ii)
        hold on
        plot(H_subset_DA_RW(:,ii*10))
        plot(H_subset_HMM(:,ii*10))
%         plot(H_subset_HMM_adapt(:,ii*10))
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
    for ii = 1:9
        subplot(3,3,ii)
        autocorr(H_subset_DA_RW(:,ii*10),200)
    end
    suptitle('DA')   
    
    figure(24)
    for ii = 1:9
        subplot(3,3,ii)
        autocorr(H_subset_HMM(:,ii*10),200)
    end    
    suptitle('HMM')
    
    figure(36)
    for ii = 1:9
        subplot(3,3,ii)
        autocorr(H_subset_HMM_adapt(:,ii*10),200)
    end    
    suptitle('HMM adapt')    
end