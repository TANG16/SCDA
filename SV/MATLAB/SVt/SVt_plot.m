close all

figure(1)
for ii = 1:4
    subplot(2,2,ii)
    hold on 
%     plot(theta_DA_RW(:,ii))
%     plot(theta_true(ii) + 0*theta_DA_RW(:,ii))
    plot(theta_HMM_eff(:,ii))
%     plot(theta_true(ii) + 0*theta_HMM_eff(:,ii))
    hold off
end

figure(2)
% plot(1:length(mean_H_DA_RW),mean_H_DA_RW)
plot(2:2:length(mean_H_HMM_eff),mean_H_HMM_eff(2:2:length(mean_H_HMM_eff)))
hold on
plot(1:length(h_true),h_true)
hold off

figure(3)
for ii = 1:9
    subplot(3,3,ii)
%     plot(H_subset_DA_RW(:,ii*10))
    plot(H_subset_HMM_eff(:,ii*5))
    hold on
%     plot(h_true(ind_h_sel(5*ii)) + 0*H_subset_HMM_eff(:,ii*5))
    hold off
end

ind_h_sel = 2:50:length(y);
ind_h_sel(10*(1:8))


figure(12)
for ii = 1:4
    subplot(2,2,ii)
    autocorr(theta_DA_RW(:,ii),100)
end