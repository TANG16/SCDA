figure(1)
for ii = 1:9
    subplot(3,3,ii)
    plot(Theta(:,ii))
end

figure(2)
for ii = 1:12
    subplot(3,4,ii)
    plot(NN(3*ii,:))
end

figure(3)
bar(mean_accept)


figure(11)
for ii = 1:9
    subplot(3,3,ii)
    plot(sample(:,ii))
end

figure(22)
for ii = 1:12
    subplot(3,4,ii)
    plot(squeeze(NN(2,3*ii,:)))
end