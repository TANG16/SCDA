close all

mean(theta_DA_RW)
mean(theta_HMM)
mean(theta_HMM_adapt)

figure(1)
for ii = 1:4
    subplot(2,2,ii)
    hold on 
    plot(theta_DA_RW(:,ii))
    plot(theta_HMM(:,ii))
%     plot(theta_true(ii) + 0*theta_HMM(:,ii))

%     plot(theta_HMM_adapt(:,ii))
%     plot(theta_true(ii) + 0*theta_HMM_adapt(:,ii))    
    plot(theta_true(ii) + 0*theta_DA_RW(:,ii))
    hold off
end

figure(10)
for ii = 1:4
    subplot(2,2,ii)
    autocorr(theta_DA_RW(:,ii),200)
end

figure(2)
hold on
plot(1:length(mean_H_DA_RW),mean_H_DA_RW)
plot(2:2:length(mean_H_HMM),mean_H_HMM(2:2:length(mean_H_HMM)))
plot(2:2:length(mean_H_HMM_adapt),mean_H_HMM_adapt(2:2:length(mean_H_HMM_adapt)))
plot(1:length(h_true),h_true)
hold off

if (ind0 == 0)
    figure(3)
    for ii = 1:9
        subplot(3,3,ii)
        hold on
    %     plot(H_subset_DA_RW(:,ii*5))
    %     plot(h_true(ind_h_sel(5*ii)) + 0*H_subset_DA_RW(:,ii*5))
    %     plot(H_subset_HMM(:,ii*5))
    %     plot(h_true(ind_h_sel(5*ii)) + 0*H_subset_HMM(:,ii*5))
    %     plot(H_subset_HMM_adapt(:,ii*5))
    %     plot(h_true(ind_h_sel(5*ii)) + 0*H_subset_HMM_adapt(:,ii*5))
        hold off
    end
else
    figure(33)
    for ii = 1:9
        subplot(3,3,ii)
        hold on
    %     plot(H_subset_DA_RW(:,ii*10))
    %     plot(H_subset_HMM(:,ii*10))
        plot(H_subset_HMM_adapt(:,ii*10))
        hold off
    end
end

figure(12)
for ii = 1:9
    subplot(3,3,ii)
    autocorr(H_subset_DA_RW(:,ii*10),200)
%     autocorr(H_subset_HMM_adapt(:,ii*10),200)
end