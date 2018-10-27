clear all
ind_meth = ones(8,1);
BurnIn = 20000;
ext = '';
BKM_collect_results;


[mean(Theta_DA(:,1:8))
% mean(Theta_HMM_adapt_10(:,1:8))
% mean(Theta_HMM_adapt_20(:,1:8))
% mean(Theta_HMM_adapt_30(:,1:8))
% mean(Theta_HMM_bin_10(:,1:8))
% mean(Theta_HMM_bin_20(:,1:8))
mean(Theta_HMM_bin_30(:,1:8))]

[mean(Theta_DA(1:4000,9))
% mean(Theta_HMM_adapt_10(:,9))
% mean(Theta_HMM_adapt_20(:,9))
% mean(Theta_HMM_adapt_30(:,9))
% mean(Theta_HMM_bin_10(:,9))
% mean(Theta_HMM_bin_20(:,9))
mean(Theta_HMM_bin_30(:,9))]


%  colMeans(mat_DA[,37:45])
%        alpha1        alphaa        alphal        alphar         beta1         betaa         betal         betar          sigy 
%     0.5490783     1.5673245    -4.5771131    -1.1760594    -0.1907766    -0.2472439    -0.3636677    -0.3421766 30440.2276841 
% > colMeans(mat_HMM_15_55[,37:45])
%        alpha1        alphaa        alphal        alphar         beta1         betaa         betal         betar          sigy 
%     0.5466424     1.5613021    -4.5779503    -1.1743723    -0.1878457    -0.2363512    -0.3640446    -0.3080615 26900.7221108 
% > colMeans(mat_HMM_29_29[,37:45])
%        alpha1        alphaa        alphal        alphar         beta1         betaa         betal         betar          sigy 
%     0.5435316     1.5520413    -4.5794260    -1.1639493    -0.1883214    -0.2365793    -0.3633601    -0.3080321 27033.6907352 
% > colMeans(mat_HMM_5_169[,37:45])
%        alpha1        alphaa        alphal        alphar         beta1         betaa         betal         betar          sigy 
%     0.5046276     1.4406817    -4.5972554    -1.0777130    -0.2134497    -0.2077289    -0.3534646    -0.3463269 29265.0988178 



bar([mean(NN_DA,2)'
% mean(NN_HMM_adapt_10,2)'
% mean(NN_HMM_adapt_20,2)'
% mean(NN_HMM_adapt_30,2)'
% mean(NN_HMM_bin_10,2)'
% mean(NN_HMM_bin_20,2)'
mean(NN_HMM_bin_30,2)']')

GG = 10;
for ii = 1:GG
    subplot(3,GG,ii)
    [acf1,lags,bounds,~]= autocorr(NN_DA(3+3*ii,:),100);
    xlabel(['t=',num2str(3+3*ii)])
    subplot(3,GG,GG+ii)
    [acf2,lags,bounds,~]= autocorr(NN_HMM_bin_30(3+3*ii,:),100);
    xlabel(['t=',num2str(3+3*ii)])
    subplot(3,GG,2*GG+ii)
    hold on
    plot(acf1)
    plot(acf2)
    hold off
    xlabel(['t=',num2str(3+3*ii)])
end



