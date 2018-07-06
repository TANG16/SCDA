ColPalette= [0 0 0
    255 121 75
    255 51 0
    153 0 0 
    62 154 222
    0 102 204
    255 0 0
    0 0 255]/255;



params = {'$\alpha_1$', '$\alpha_a$', '$\alpha_{\rho}$', '$\alpha_{\lambda}$', ...
            '$\beta_1$', '$\beta_a$', '$\beta_{\rho}$', '$\beta_{\lambda}$',...
            '$\sigma^{2}_{y}$'};

ff = figure(1010);    
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.55 0.65]);
set(gcf,'defaulttextinterpreter','latex');  
ij = 0;
for ii = [2,7]; %[1,2,3,5,6,7]%1:9
%         subplot(3,3,ii)
ij = ij + 1;
    subplot(2,2,ij)
    hold on
    for jj = [1,7,8]; %1:K
%             subplot(2,4,jj)
        xxx = eval(char(TH_all{jj}));
        plot(0:100,autocorr(xxx(:,ii),100),...
            'color',ColPalette(jj,:),'linewidth',2)                   
    end
    hold off
    ylabel('ACF','fontsize',11)
    xlabel('Lag','fontsize',11)  
    title(params{ii},'fontsize',11)
%         title(method{jj}) 
%         suptitle([params{ii},' (BurnIn = ', num2str(BurnIn),')'])
end    

for ii = [5,11]
    ij = ij + 1;
    subplot(2,2,ij)
    hold on
    for jj = [1,7,8] %1:K
%             subplot(2,4,jj)
        xxx = eval(char(NN_all{jj}));
        plot(0:100,autocorr(xxx(3*ii,:),100),...
            'color',ColPalette(jj,:),'linewidth',2)
    end    
    hold off
    ylabel('ACF','fontsize',11)
    xlabel('Lag','fontsize',11)    
    title(['Na(',num2str(3*ii),')'],'fontsize',11)
%         suptitle(['Na(',num2str(3*ii),') (BurnIn = ', num2str(BurnIn),')'])

%         name_fig = [figures_path,'BKM_ALL_Na',num2str(3*ii),'_acf.png']; %eps
%         set(gcf,'PaperPositionMode','auto');
%         print(ff,name_fig,'-dpng','-r0') 
end

legend({'DA','SCDA Bin30','SCDA Exact'},...
    'interpreter','latex','fontsize',11)