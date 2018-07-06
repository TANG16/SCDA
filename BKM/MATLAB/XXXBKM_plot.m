%%
params = {'alpha1', 'alphaa', 'alphar', 'alphal', ...
        'beta1', 'betaa', 'betar', 'betal',...
        'sigy'};
    
    
sample = Results_DA.Theta;
mean_Theta_DA = mean(Results_DA.Theta);

if plot_on
    figure(1)
    set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
    for ii = 1:9
        subplot(3,3,ii)
        hold on
        plot(sample(BurnIn+1:M,ii))
%         plot(theta_init(ii)+0*sample(BurnIn+1:M,ii),'r')        
        hold off
        title(params{ii})
    end

    figure(10)
    set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
    for ii = 1:9
        subplot(3,3,ii)
        autocorr(sample(BurnIn+1:M,ii),40)
        title(params{ii})   
    end   
    
    figure(100)
    set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
    for ii = 1:9
        subplot(3,3,ii)
        bar(acf(sample(BurnIn+1:M,ii),40))
        title(params{ii})
    end
    

    figure(3)
    set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
    for ii = 1:9
        subplot(3,3,ii)
        hold on
        plot(NN(4*(ii-3)+9,BurnIn+1:M))
%         plot(Na(4*(ii-3)+9) + 0*NN(4*(ii-3)+9,BurnIn+1:M),'r')
        hold off
        title(['Na(',num2str(4*(ii-3)+9),')'])
    end

    figure(33)
    set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
    for ii = 1:9
        subplot(3,3,ii)
        autocorr(NN(4*(ii-3)+9,BurnIn+1:M),40)
        title(['Na(',num2str(4*(ii-3)+9),')'])
    end
    
    figure(333)
    set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
    for ii = 1:9
        subplot(3,3,ii)
        bar(acf(NN(4*(ii-3)+9,BurnIn+1:M),40))
        title(['Na(',num2str(4*(ii-3)+9),')'])
    end
    
    figure(4)
    set(gcf,'units','normalized','outerposition',[0.0 0.0 1.0 1.0]);
    hold on
    plot(mean(NN,2))
    hold off
end