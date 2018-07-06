for BurnIn = [0, 1000, 5000,10000]
    clearvars -except BurnIn
    close all
BurnIn = 10000;

    ext = ''; %'_3';
    
%     ind_meth = logical( [1,1,0,0,1,0,0,1]');
    ind_meth = logical( [1,1,1,1,1,1,1,1]');
%     ind_meth = true(8,1);
    BKM_collect_results;

    figures_path = ''; %['Results/BurnIn_',num2str(BurnIn),ext,'/'];
%     figures_path = ['Figures/BurnIn_',num2str(BurnIn),'/'];

    ColPalette= [0 0 0
        255 121 75
        255 51 0
        153 0 0 
        62 154 222
        0 102 204
        51 51 102
        0 255 0]/255;

     ColPalette= ColPalette(ind_meth,:);
    
    params = {'alpha1', 'alphaa', 'alphar', 'alphal', ...
        'beta1', 'betaa', 'betar', 'betal',...
        'sigy'};

    method = {'DA','Adapt10','Adapt20','Adapt30','Bin10','Bin20','Bin30', 'Exact'};
    method = method(ind_meth);
    K = length(method);
    
    TH_all = who('-regexp', '^Theta');
    NN_all = who('-regexp', '^NN');
    accept_all = who('-regexp','mean_accept');
    time_sampl_all = who('-regexp','time');

    time = zeros(K,1);
    for jj = 1:K
        time(jj) = eval(char(time_sampl_all{jj}));
    end
    
    %% param trace
    ff = figure(1);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
    for ii = 1:9
        subplot(3,3,ii)
        hold on
        for jj = 1:K
            xxx = eval(char(TH_all{jj}));
            plot(xxx(:,ii),'color',ColPalette(jj,:))
        end
        hold off
        title(params{ii})
    end
    leg = legend(method,...
        'Orientation','Horizontal');

    name_fig = [figures_path,'BKM_ALL_param_trace.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print(ff,name_fig,'-dpng','-r0')


    %% Na trace
    ff = figure(2);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
    for ii = 1:12
        subplot(3,4,ii) 
        hold on
        for jj = 1:K
            xxx = eval(char(NN_all{jj}));
            plot(xxx(3*ii,:),'color',ColPalette(jj,:))
        end
        hold off
        title(['Na(',num2str(3*ii),')'])
    end
    legend(method,...
        'Orientation','Horizontal');


    name_fig = [figures_path,'BKM_ALL_Na_trace.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print(ff,name_fig,'-dpng','-r0')

    
    %% Na mean
    ff = figure(2022);
%     set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.6]);
    set(gcf,'defaulttextinterpreter','latex');
    hold on
    for jj = 1:K
        xxx = mean(eval(char(NN_all{jj})),2);
        plot(xxx,'color',ColPalette(jj,:))
    end
    hold off
    set(gca,'TickLabelInterpreter','latex');
    legend(method,'Location','southoutside',...
        'Orientation','Horizontal','Interpreter','latex');

    name_fig = [figures_path,'BKM_ALL_Na_mean.eps']; %png
    set(gcf,'PaperPositionMode','auto');
    print(ff,name_fig,'-depsc','-r0')

    %% param mean accept
    ff = figure(3);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
    for ii = 1:(9-1)
        subplot(3,3,ii)
        ACC = zeros(1,K);
        hold on
        for jj = 1:K
            xxx = eval(char(accept_all{jj}));
            ACC(jj) = xxx(36+ii);     
            h = bar(jj,ACC(jj));
            set(h,'FaceColor',ColPalette(jj,:))
        end
        set(gca,'XTickLabel',[])
        hold off
        title(['Mean accept ',params{ii}])
    end
    legend(method)

    name_fig = [figures_path,'BKM_ALL_param_mean_accept.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print(ff,name_fig,'-dpng','-r0')



    %% Na mean accept
    ff = figure(4);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
    for ii = 1:12
        subplot(3,4,ii)
        ACC = zeros(1,K);
        hold on
        for jj = 1:K
            xxx = eval(char(accept_all{jj}));
            ACC(jj) = xxx(3*ii);     
            h = bar(jj,ACC(jj));
            set(h,'FaceColor',ColPalette(jj,:))
        end
        set(gca,'XTickLabel',[])
        hold off
        title(['Mean accept Na(',num2str(3*ii),')'])
    end 
    legend(method)


    name_fig = [figures_path,'BKM_ALL_Na_mean_accept.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print(ff,name_fig,'-dpng','-r0') %depsc
    
    
    %% ACF param
    for ii = 1:9
        ff = figure(ii*10) ;
        set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
        for jj = 1:K
            subplot(2,4,jj)
            xxx = eval(char(TH_all{jj}));
            stem(0:100,autocorr(xxx(:,ii),100))
            title(method{jj})
        end         
        suptitle([params{ii},' (BurnIn = ', num2str(BurnIn),')'])
        
        name_fig = [figures_path,'BKM_ALL_param_',params{ii},'_acf.png']; %eps
        set(gcf,'PaperPositionMode','auto');
        print(ff,name_fig,'-dpng','-r0') 
    end
  
    ff = figure(10*10);    
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
     ij = 0;
    for ii = [1,2,3,5,6,7]%1:9
%         subplot(3,3,ii)
ij = ij + 1;
        subplot(2,3,ij)
        hold on
        for jj = [1,7,8]; %1:K
%             subplot(2,4,jj)
            xxx = eval(char(TH_all{jj}));
            plot(0:100,autocorr(xxx(:,ii),100), 'color',ColPalette(jj,:))         
        end
        hold off
        title(params{ii})
%         title(method{jj}) 
%         suptitle([params{ii},' (BurnIn = ', num2str(BurnIn),')'])
    end    
    name_fig = [figures_path,'BKM_ALL_param__acf.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print(ff,name_fig,'-dpng','-r0')     
    
    %% ACF Na
    for ii = 1:12
        ff = figure(ii*100) ; 
        set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
        for jj = 1:K
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

    ff = figure(10*100) ; 
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
    for ii = 1:12
        subplot(3,4,ii)
        hold on
        for jj = [1,7,8] %1:K
%             subplot(2,4,jj)
            xxx = eval(char(NN_all{jj}));
            plot(0:100,autocorr(xxx(3*ii,:),100),'color',ColPalette(jj,:))
        end    
        hold off
        title(['Na(',num2str(3*ii),')'])
%         suptitle(['Na(',num2str(3*ii),') (BurnIn = ', num2str(BurnIn),')'])

%         name_fig = [figures_path,'BKM_ALL_Na',num2str(3*ii),'_acf.png']; %eps
%         set(gcf,'PaperPositionMode','auto');
%         print(ff,name_fig,'-dpng','-r0') 
    end
        
    
    %% ESS
    T = 36;
    ESS_N_40 = zeros(K,T);
    ESS_N_100 = zeros(K,T);
    ESS_N_1000 = zeros(K,T);
    ESS_N_sig = zeros(K,T);
        
    for jj = 1:K
        ESS_N_40(jj,:) = ESS(eval(char(NN_all{jj}))',40);
        ESS_N_100(jj,:) = ESS(eval(char(NN_all{jj}))',100);
        ESS_N_1000(jj,:) = ESS(eval(char(NN_all{jj}))',1000);
        ESS_N_sig(jj,:) = ESS(eval(char(NN_all{jj}))',0);
    end

    D = 9;
    ESS_TH_40 = zeros(8,D);
    ESS_TH_100 = zeros(8,D);
    ESS_TH_1000 = zeros(8,D);
    ESS_TH_sig = zeros(8,D);
        
    for jj = 1:K
        ESS_TH_40(jj,:) = ESS(eval(char(TH_all{jj})),40);
        ESS_TH_100(jj,:) = ESS(eval(char(TH_all{jj})),100);
        ESS_TH_1000(jj,:) = ESS(eval(char(TH_all{jj})),1000);
        ESS_TH_sig(jj,:) = ESS(eval(char(TH_all{jj})),0);
    end
    
    ff = figure(987);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
    subplot(2,2,1)    
    hold on
    for jj = 1:K
        plot(3:T,ESS_N_40(jj,3:T),'color',ColPalette(jj,:))
    end
    hold off
    xlabel('Time')
    title('ESS L=40')

    subplot(2,2,2)    
    hold on
    for jj = 1:K
        plot(3:T,ESS_N_100(jj,3:T),'color',ColPalette(jj,:))
    end
    hold off
    xlabel('Time')
    title('ESS L=100')
    
    subplot(2,2,3)    
    hold on
    for jj = 1:K
        plot(3:T,ESS_N_1000(jj,3:T),'color',ColPalette(jj,:))
    end
    hold off
    xlabel('Time')
    title('ESS L=1000')
    
    subplot(2,2,4)    
    hold on
    for jj = 1:K
        plot(3:T,ESS_N_sig(jj,3:T),'color',ColPalette(jj,:))
    end
    hold off
	xlabel('Time')
    title('ESS L=sig')    
    
    legend(method)
    

    name_fig = [figures_path,'BKM_ALL_Na_ESS.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print(ff,name_fig,'-dpng','-r0')     
%     plot(ESS_N_sig)
%     title('ESS L=sig')
    

ff = figure(98765);
subplot(2,2,1)
bar(ESS_TH_40')
set(gca, 'XTick', 1:9,'XTickLabel',params)
title('ESS L=40')

subplot(2,2,2)
bar(ESS_TH_100')
set(gca, 'XTick', 1:9,'XTickLabel',params)
title('ESS L=100')

subplot(2,2,3)
bar(ESS_TH_1000')
set(gca, 'XTick', 1:9,'XTickLabel',params)
title('ESS L=1000')

subplot(2,2,4)
bar(ESS_TH_sig')
set(gca, 'XTick', 1:9,'XTickLabel',params)
title('ESS L=sig')

legend(method);

    name_fig = [figures_path,'BKM_ALL_param_ESS.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print(ff,name_fig,'-dpng','-r0')   
end
 