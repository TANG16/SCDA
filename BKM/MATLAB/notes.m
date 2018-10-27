
for jj = 1:36
    subplot(4,9,jj)
    hold on
%     plot(NN(ii,1:50000))
    plot(NN_HMM(jj,1:(ii-1)))
    plot(squeeze(NN(2,jj,1:(ii-1))))
%     legend('HMM','DA')
end

NN_HMM = NN;
load('Results\BKM_DA_M100000.mat', 'NN')

hold on
% plot(mean(NN_HMM(:,10000:(ii-1)),2))
plot(mean(NN_HMM,2))
% plot(mean(NN_HMM2(:,1:(ii-1)),2))
% plot(mean(squeeze(NN(2,:,10000:(ii-1))),2))
plot(mean(squeeze(NN(2,:,:)),2))

mean(squeeze(NN(2,:,:)),2)



mean(Theta(1:(ii-1),1:8))
%     0.3477    0.9889   -0.7809   -2.9938   -0.1027   -0.1387   -0.2738   -0.2356

%        alpha1        alphaa        alphal        alphar         beta1         betaa         betal         betar          sigy 
%     0.5490783     1.5673245    -4.5771131    -1.1760594    -0.1907766    -0.2472439    -0.3636677    -0.3421766 30440.2276841 

figure(2)
for ii = 1:36
    subplot(4,9,ii)
    hold on
%     plot(NN(ii,1:50000))
    plot(squeeze(NN(1,ii,:)))
%     legend('HMM','DA')
end


figure(3)
for jj = 1:9
    subplot(3,3,jj)
    hold on
    plot(Theta(1:(ii-1),jj))
end
