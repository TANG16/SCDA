load('Codes\BKM\MATLAB\Results\BurnIn_10000_1\BKM_exact_Nmax690.mat', 'Theta')
load('Codes\BKM\MATLAB\Results\BurnIn_10000_1\BKM_bin_Nbin30.mat', 'Theta')
load('Codes\BKM\MATLAB\Results\BurnIn_10000_1\BKM_DA_NRW_U.mat', 'Theta')
ACF_DA = autocorr(Theta_DA(:,2),100);
ACF_30 = autocorr(Theta_30(:,2),100);
ACF_exact = autocorr(Theta_exact(:,2),100);

figure(1)
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
hold on
plot(ACF_DA,'k','linewidth',2)
plot(ACF_exact,'b','linewidth',2)
plot(ACF_30,'r','linewidth',2)
hold off
xlim([0 100])
xlabel('Lag','Interpreter','latex')
ylabel('ACF','Interpreter','latex')
title('\alpha_{a}')
legend('DA','SCDA Exact', 'SCDA B30')