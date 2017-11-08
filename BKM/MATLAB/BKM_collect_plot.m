for BurnIn = [0, 1000, 5000,10000]
    clearvars -except BurnIn
    close all

    BKM_collect_results;

    figures_path = ['Figures/BurnIn_',num2str(BurnIn),'/'];

    ColPalette= [0 0 0
        255 121 75
        255 51 0
        153 0 0 
        62 154 222
        0 102 204
        51 51 102]/255;

    params = {'alpha1', 'alphaa', 'alphar', 'alphal', ...
        'beta1', 'betaa', 'betar', 'betal',...
        'sigy'};

    method = {'DA','Adapt10','Adapt20','Adapt30','Bin10','Bin20','Bin30'};
    
    TH_all = who('-regexp', '^Theta');
    NN_all = who('-regexp', '^NN');
    accept_all = who('-regexp','mean_accept');
    time_sampl_all = who('-regexp','time');

    time = zeros(7,1);
    for jj = 1:7
        time(jj) = eval(char(time_sampl_all{jj}));
    end

    ff = figure(1);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
    for ii = 1:9
        subplot(3,3,ii)
        hold on
        for jj = 1:7
            xxx = eval(char(TH_all{jj}));
            plot(xxx(:,ii),'color',ColPalette(jj,:))
        end
        hold off
        title(params{ii})
    end
    leg = legend('DA','Adapt10','Adapt20','Adapt30','Bin10','Bin20','Bin30',...
        'Orientation','Horizontal','Position',[-0.1766 -0.0129    0.3534    0.0250]);

    name_fig = [figures_path,'BKM_ALL_param_trace.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print(ff,name_fig,'-dpng','-r0')



    ff = figure(2);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
    for ii = 1:12
        subplot(3,4,ii) 
        hold on
        for jj = 1:7
            xxx = eval(char(NN_all{jj}));
            plot(xxx(3*ii,:),'color',ColPalette(jj,:))
        end
        hold off
        title(['Na(',num2str(3*ii),')'])
    end
    legend('DA','Adapt10','Adapt20','Adapt30','Bin10','Bin20','Bin30',...
        'Orientation','Horizontal','Position',[-0.1766 -0.0129    0.3534    0.0250]);


    name_fig = [figures_path,'BKM_ALL_Na_trace.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print(ff,name_fig,'-dpng','-r0')



    ff = figure(3);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
    for ii = 1:8
        subplot(3,3,ii)
        ACC = zeros(1,7);
        hold on
        for jj = 1:7
            xxx = eval(char(accept_all{jj}));
            ACC(jj) = xxx(36+ii);     
            h = bar(jj,ACC(jj));
            set(h,'FaceColor',ColPalette(jj,:))
        end
        set(gca,'XTickLabel',[])
        hold off
        title(['Mean accept ',params{ii}])
    end
    s9 = subplot(3,3,9);
    delete(s9)
    legend('DA','Adapt10','Adapt20','Adapt30','Bin10','Bin20','Bin30')

    name_fig = [figures_path,'BKM_ALL_param_mean_accept.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print(ff,name_fig,'-dpng','-r0')




    ff = figure(4);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
    for ii = 1:12
        subplot(3,4,ii)
        ACC = zeros(1,7);
        hold on
        for jj = 1:7
            xxx = eval(char(accept_all{jj}));
            ACC(jj) = xxx(3*ii);     
            h = bar(jj,ACC(jj));
            set(h,'FaceColor',ColPalette(jj,:))
        end
        set(gca,'XTickLabel',[])
        hold off
        title(['Mean accept Na(',num2str(3*ii),')'])
    end 
    legend('DA','Adapt10','Adapt20','Adapt30','Bin10','Bin20','Bin30')


    name_fig = [figures_path,'BKM_ALL_Na_mean_accept.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print(ff,name_fig,'-dpng','-r0') %depsc
    
    
    %% ACF param
    for ii = 6:9
        ff = figure(ii*10) ;
        set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
        for jj = 1:7
            subplot(2,4,jj)
            xxx = eval(char(TH_all{jj}));
            autocorr(xxx(:,ii),100)
            title(method{jj})
        end         
        suptitle([params{ii},' (BurnIn = ', num2str(BurnIn),')'])
        
        name_fig = [figures_path,'BKM_ALL_param_',params{ii},'_acf.png']; %eps
        set(gcf,'PaperPositionMode','auto');
        print(ff,name_fig,'-dpng','-r0') 
    end
  
    %% ACF Na
    for ii = 1:12
        ff = figure(ii*100) ; 
        set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
        for jj = 1:7
            subplot(2,4,jj)
            xxx = eval(char(NN_all{jj}));
            autocorr(xxx(3*ii,:),100)
            title(method{jj})
        end         
        suptitle(['Na(',num2str(3*ii),') (BurnIn = ', num2str(BurnIn),')'])

        name_fig = [figures_path,'BKM_ALL_Na',num2str(3*ii),'_acf.png']; %eps
        set(gcf,'PaperPositionMode','auto');
        print(ff,name_fig,'-dpng','-r0') 
    end
    
    
    %% ESS
    T = 36;
    ESS_N_40 = zeros(7,T);
    ESS_N_100 = zeros(7,T);
    ESS_N_1000 = zeros(7,T);
    ESS_N_sig = zeros(7,T);
        
    for jj = 1:7
        ESS_N_40(jj,:) = ESS(eval(char(NN_all{jj}))',40);
        ESS_N_100(jj,:) = ESS(eval(char(NN_all{jj}))',100);
        ESS_N_1000(jj,:) = ESS(eval(char(NN_all{jj}))',1000);
%         ESS_N_sig(jj,:) = ESS(eval(char(NN_all{jj}))',0);
    end

    D = 9;
    ESS_TH_40 = zeros(7,D);
    ESS_TH_100 = zeros(7,D);
    ESS_TH_1000 = zeros(7,D);
    ESS_TH_sig = zeros(7,D);
        
    for jj = 1:7
        ESS_TH_40(jj,:) = ESS(eval(char(TH_all{jj})),40);
        ESS_TH_100(jj,:) = ESS(eval(char(TH_all{jj})),100);
        ESS_TH_1000(jj,:) = ESS(eval(char(TH_all{jj})),1000);
%         ESS_TH_sig(jj,:) = ESS(eval(char(TH_all{jj})),0);
    end
    
    figure(987)
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
    subplot(2,2,1)    
    hold on
    for jj = 1:7
        plot(ESS_N_40(jj,:),'color',ColPalette(jj,:))
    end
    hold off
    title('ESS L=40')

    subplot(2,2,2)    
    hold on
    for jj = 1:7
        plot(ESS_N_100(jj,:),'color',ColPalette(jj,:))
    end
    hold off
    title('ESS L=100')
    
    subplot(2,2,3)    
    hold on
    for jj = 1:7
        plot(ESS_N_1000(jj,:),'color',ColPalette(jj,:))
    end
    hold off
    title('ESS L=1000')
    
    legend('DA','Adapt10','Adapt20','Adapt30','Bin10','Bin20','Bin30')
    
%     plot(ESS_N_sig)
%     title('ESS L=sig')
    
end
 