load('SV_MC_est_H_T2000_Nbin30_Nq30.mat')


params = {'\mu','\phi','\sigma^2'};
ColPal = [ 0.6000    0.6000    0.6000
           0.9900    0.1900    0.1900
           0.3010    0.7450    0.9330
           0.2500    0.2500    0.2500 
           0.6350    0.0780    0.1840
                0    0.4470    0.7410];
           
ff = figure(11);
set(gcf,'units','normalized','outerposition',[0 0 1 1]);    
for ii = 1:3
%     subplot(3,1,ii)
    ff = figure(ii);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);    

    b = bar([THETA_mean_DA(:,ii),THETA_mean_HMM(:,ii),THETA_mean_HMM_adapt(:,ii)]);
    set(b(1), 'FaceColor',ColPal(1,:))
    set(b(2), 'FaceColor',ColPal(2,:))
    set(b(3), 'FaceColor',ColPal(3,:))

    xlim([0 (S+1)])
    hold on
    plot(0:(S+1),mean(THETA_mean_DA(:,ii)) + 0*(0:(S+1)),'Color',ColPal(4,:),'linewidth',2)
    plot(0:(S+1),mean(THETA_mean_HMM(:,ii)) + 0*(0:(S+1)),'Color',ColPal(5,:),'linewidth',2)
    plot(0:(S+1),mean(THETA_mean_HMM_adapt(:,ii)) + 0*(0:(S+1)),'Color',ColPal(6,:),'linewidth',2)
    plot(0:(S+1),theta_true(ii) + 0*(0:(S+1)),'y','linewidth',2)
    hold off
    title(params{ii})
    legend('DA',['HMM ',num2str(N_bin)],['HMM adapt ',num2str(N_q)], 'mean DA','mean HMM','mean HMM adapt','true')
end
name = ['MC_est_H_T',num2str(T),'_Nbin',num2str(N_bin),'_Nq',num2str(N_q),'.png']; %eps
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-dpng','-r0')