% clear all
BurnIn = 10000;

ext = '_Selected';
ext = '';

name = ['Results/BurnIn_',num2str(BurnIn),ext,'/BKM_DA_NRW_U.mat'];
Results_DA = load(name);
% name = ['Results/BurnIn_',num2str(BurnIn),ext,'/BKM_exact_Nmax690.mat'];
name = ['Results/BurnIn_',num2str(BurnIn),ext,'/BKM_exact_Nmax679.mat'];
Results_Exact = load(name);

            
for ii = 1:3
    name = ['Results/BurnIn_',num2str(BurnIn),ext,'/BKM_adapt_Nq',num2str(10*ii),'.mat'];
    Results_Adapt{ii} = load(name);
    name = ['Results/BurnIn_',num2str(BurnIn),ext,'/BKM_bin_Nbin',num2str(10*ii),'.mat'];
    Results_Bin{ii} = load(name);
end


method = {'DA','Adapt10','Adapt20','Adapt30','Bin10','Bin20','Bin30', 'Exact'};
Time = [Results_DA.time_sampl;
    Results_Adapt{1}.time_sampl;
    Results_Adapt{2}.time_sampl;
    Results_Adapt{3}.time_sampl;
    Results_Bin{1}.time_sampl;
    Results_Bin{2}.time_sampl;
    Results_Bin{3}.time_sampl;
    Results_Exact.time_sampl];
% save(['Results/BurnIn_',num2str(BurnIn),ext,'/Time_all.mat'], 'Time'); 



Accept_Na = [mean(Results_DA.accept(:,37:72));
    mean(Results_Adapt{1}.accept(:,1:36));
    mean(Results_Adapt{2}.accept(:,1:36));
    mean(Results_Adapt{3}.accept(:,1:36));
    mean(Results_Bin{1}.accept(:,1:36));
    mean(Results_Bin{2}.accept(:,1:36));
    mean(Results_Bin{3}.accept(:,1:36));
    mean(Results_Exact.accept(:,1:36))];


Accept_Theta = [mean(Results_DA.accept(:,73:80));
    mean(Results_Adapt{1}.accept(:,37:44));
    mean(Results_Adapt{2}.accept(:,37:44));
    mean(Results_Adapt{3}.accept(:,37:44));
    mean(Results_Bin{1}.accept(:,37:44));
    mean(Results_Bin{2}.accept(:,37:44));
    mean(Results_Bin{3}.accept(:,37:44));
    mean(Results_Exact.accept(:,37:44))];
% save(['Results/BurnIn_',num2str(BurnIn),ext,'/Accept_all.mat'], 'Accept_Na','Accept_Theta'); 


ColPalette= [0 0 0
    255 121 75
    255 51 0
    153 0 0 
    62 154 222
    0 102 204
    51 51 102
    0 255 0]/255;

params = {'$\alpha_1$', '$\alpha_a$', '$\alpha_{\rho}$', '$\alpha_{\lambda}$', ...
            '$\beta_1$', '$\beta_a$', '$\beta_{\rho}$', '$\beta_{\lambda}$',...
            '$\sigma^{2}_{y}$'};
%% Params trace
ff = figure(1);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
set(gcf,'defaulttextinterpreter','latex');
for ii = 1:9
    subplot(3,3,ii)
    hold on
    xxx = Results_DA.Theta(:,ii);
    plot(xxx,'color',ColPalette(1,:))


    for jj = 1:3
        xxx = Results_Adapt{1,jj};
        xxx = xxx.Theta(:,ii);
        plot(xxx,'color',ColPalette(1+jj,:))
    end

    for jj = 1:3
        xxx = Results_Bin{1,jj};
        xxx = xxx.Theta(:,ii);
        plot(xxx,'color',ColPalette(4+jj,:))
    end 

    xxx = Results_Exact.Theta(:,ii);
    plot(xxx,'color',ColPalette(8,:))
    hold off
    
    set(gca,'TickLabelInterpreter','latex');        
    title(params{ii},'interpreter','latex','FontSize',11)
end
ll = legend(method);
set(ll,'interpreter','latex','orientation','horizontal','FontSize',11);

name = ['Figures/BurnIn_10000/BKM_Theta_trace_all.eps'];
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-depsc','-r0')


%% Na trace
ff = figure(2);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
set(gcf,'defaulttextinterpreter','latex');
for ii = 1:9
    subplot(3,3,ii)
    hold on
    xxx = squeeze(Results_DA.NN(2,4*ii,:));
    plot(xxx,'color',ColPalette(1,:))


    for jj = 1:3
        xxx = Results_Adapt{1,jj};
        xxx = xxx.NN(4*ii,:);
        plot(xxx,'color',ColPalette(1+jj,:))
    end

    for jj = 1:3
        xxx = Results_Bin{1,jj};
        xxx = xxx.NN(4*ii,:);
        plot(xxx,'color',ColPalette(4+jj,:))
    end 

    xxx = Results_Exact.NN(4*ii,:);
    plot(xxx,'color',ColPalette(8,:))

    hold off
    set(gca,'TickLabelInterpreter','latex');    
    title(['$N_{a,',num2str(4*ii),'}$'],'interpreter','latex','FontSize',11)
end
ll = legend(method);
set(ll,'interpreter','latex','orientation','horizontal','FontSize',11);

name = ['Figures/BurnIn_10000/BKM_N_trace_all.eps'];
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-depsc','-r0')


%% Hist
figure(30)
for ii = 1:9
    xxx = squeeze(Results_DA.NN(2,:,:));
    subplot(3,3,ii)
    histogram(xxx(4*ii,:))
end

% lo = y - e;
% hi = y + e;
% hp = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)], 'r');
% hold on;
% hl = line(x,y);
% set(hp, 'facecolor', [1 0.8 0.8], 'edgecolor', 'none');
% set(hl, 'color', 'r', 'marker', 'x');
% HPDI =  patch([x; x(end:-1:1); x(1)],...
%     [lo; hi(end:-1:1); lo(1)],...
%     'color',ColPalette(1,:));

%     HPDI_m = hpdi(xxx',0.95);


%% Na means
ff = figure(3);
    tt = (1:36)';
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.5 0.95]);
    set(gcf,'defaulttextinterpreter','latex');
%     %%
    subplot(3,1,1)
    title('$N_{a}$ posterior means','interpreter','latex','FontSize',11)
    hold on
    
    xxx = squeeze(Results_DA.NN(2,:,:));    
    plot(mean(xxx,2),'color',ColPalette(1,:),'linewidth',1)
    
    for jj = 1:3
        xxx = Results_Adapt{1,jj};                  
        plot(mean(xxx.NN,2),'color',ColPalette(1+jj,:),'linewidth',1)
    end

    for jj = 1:3
        xxx = Results_Bin{1,jj};                
        plot(mean(xxx.NN,2),'color',ColPalette(4+jj,:),'linewidth',1)
    end
    
    xxx = Results_Exact.NN;   
    plot(mean(xxx,2),'color',ColPalette(8,:),'linewidth',1) 
    
    hold off 
    set(gca,'TickLabelInterpreter','latex');    
    xlim([1,36])
    ll = legend(method);
    set(ll,'interpreter','latex','orientation','horizontal','FontSize',11);
    
%     %%
    subplot(3,1,2)
    xlabel('$N_{a}$ posterior means and 95\% CI, DA and adaptive bins','interpreter','latex','FontSize',11)     
    hold on
    
    xxx = squeeze(Results_DA.NN(2,:,:)); 
    xxx_L = prctile(xxx',2.5)';
    xxx_U = prctile(xxx',97.5)';   
    HPDI =  patch([tt; tt(end:-1:1); tt(1)],...
        [xxx_L; xxx_U(end:-1:1); xxx_L(1)],...
        'r');
    set(HPDI, 'facecolor', ColPalette(1,:), 'facealpha',0.2, 'edgecolor', 'none');
    plot(mean(xxx,2),'color',ColPalette(1,:),'linewidth',2)

    for jj = 1:3
        xxx = Results_Adapt{1,jj};        
            xxx_L = prctile(xxx.NN',2.5)';
            xxx_U = prctile(xxx.NN',97.5)';      
            HPDI =  patch([tt; tt(end:-1:1); tt(1)],...
                [xxx_L; xxx_U(end:-1:1); xxx_L(1)],...
                'r');
            set(HPDI, 'facecolor', ColPalette(1+jj,:), 'facealpha',0.2, 'edgecolor', 'none');        
        plot(mean(xxx.NN,2),'color',ColPalette(1+jj,:),'linewidth',2)
    end    
    hold off   
    set(gca,'TickLabelInterpreter','latex');    
    xlim([1,36])
    
%     %%
    subplot(3,1,3)
    xlabel('$N_{a}$ posterior means and 95\% CI, DA and fixed bins','interpreter','latex','FontSize',11) 
    
    hold on
    
    xxx = squeeze(Results_DA.NN(2,:,:)); 
    xxx_L = prctile(xxx',2.5)';
    xxx_U = prctile(xxx',97.5)';   
    HPDI =  patch([tt; tt(end:-1:1); tt(1)],...
        [xxx_L; xxx_U(end:-1:1); xxx_L(1)],...
        'r');
    set(HPDI, 'facecolor', ColPalette(1,:), 'facealpha',0.2, 'edgecolor', 'none');
    plot(mean(xxx,2),'color',ColPalette(1,:),'linewidth',2)

    for jj = 1:3
        xxx = Results_Bin{1,jj};        
            xxx_L = prctile(xxx.NN',2.5)';
            xxx_U = prctile(xxx.NN',97.5)';      
            HPDI =  patch([tt; tt(end:-1:1); tt(1)],...
                [xxx_L; xxx_U(end:-1:1); xxx_L(1)],...
                'r');
            set(HPDI, 'facecolor', ColPalette(4+jj,:), 'facealpha',0.2, 'edgecolor', 'none');        
        plot(mean(xxx.NN,2),'color',ColPalette(4+jj,:),'linewidth',2)
    end 

%     xxx = Results_Exact.NN;
%         xxx_L = prctile(xxx',2.5)';
%         xxx_U = prctile(xxx',97.5)';   
%         xxx = mean(xxx,2);
%         HPDI =  patch([tt; tt(end:-1:1); tt(1)],...
%             [xxx_L; xxx_U(end:-1:1); xxx_L(1)],...
%             'r');
%         set(HPDI, 'facecolor',ColPalette(8,:), 'facealpha',0.2, 'edgecolor', 'none');
%     plot(mean(xxx,2),'color',ColPalette(8,:),'linewidth',2)
     
    hold off        
    set(gca,'TickLabelInterpreter','latex');    
    xlim([1,36])   
 
name = ['Figures/BurnIn_10000/BKM_N_mean_CI_all.eps'];
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-depsc','-r0')



%% ACF Theta
par_ind = [1,2,3,5,6,7];
L = 200;
ff = figure(1000);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.55 0.65]);
%     set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.45]);
    set(gcf,'defaulttextinterpreter','latex');  
%     xxx1 = Results_DA.Theta;    
%     xxx2 = Results_Adapt{1,1};  xxx2 = xxx2.Theta;              
%     xxx3 = Results_Bin{1,2};  xxx3 = xxx3.Theta; 
    for ii = 1:6
        subplot(2,3,ii)  
        hold on
        [aaa,~,bbb] = autocorr(Results_DA.Theta(:,par_ind(ii)),L);
        plot(0:L,aaa,'color',ColPalette(1,:)) 
        
        for jj = 1:3
            [aaa] = autocorr(Results_Adapt{1,jj}.Theta(:,par_ind(ii)),L);
            plot(0:L,aaa,'color',ColPalette(1+jj,:))               
        end
        for jj = 1:3
            [aaa] = autocorr(Results_Bin{1,jj}.Theta(:,par_ind(ii)),L);
            plot(0:L,aaa,'color',ColPalette(4+jj,:))               
        end       
        [aaa] = autocorr(Results_Exact.Theta(:,par_ind(ii)),L);
        plot(0:L,aaa,'color',ColPalette(8,:)) 
        
%         plot([0 0; L L],[bbb([1 1]) bbb([2 2])],'-b'); 
        hold off
%         [aaa1,~,bbb1] = autocorr(xxx1(:,par_ind(ii)),L);
%         [aaa2,~,bbb2] = autocorr(xxx2(:,par_ind(ii)),L);
%         [aaa3,~,bbb3] = autocorr(xxx3(:,par_ind(ii)),L);        
%     	subplot(2,3,ii) 
%         hold on 
%         plot(0:L,aaa1,'color',ColPalette(1,:)) 
%         plot(0:L,aaa2,'color',ColPalette(2,:)) 
%         plot(0:L,aaa3,'color',ColPalette(6,:))         
%         hold off
        xlim([0,L])
        set(gca,'TickLabelInterpreter','latex');    
        xlabel('Lag','Interpreter','latex');
        ylabel('ACF','Interpreter','latex');
        title([params{par_ind(ii)}],'Interpreter','latex');
    end
%     ll = legend(method{[1,2,6]});
    ll = legend(method);
    set(ll,'interpreter','latex','orientation','horizontal','FontSize',11);

name = ['Figures/BurnIn_10000/BKM_Theta_acf_Comb.png'];
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-dpng','-r0') 

name = ['Figures/BurnIn_10000/BKM_Theta_acf_Comb.eps'];
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-depsc','-r0')




ff = figure(10);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.45]);
    set(gcf,'defaulttextinterpreter','latex');  
    xxx = Results_DA.Theta;    
    for ii = 1:6
    	subplot(2,3,ii) 
        [~,~,~,hh] = autocorr(xxx(:,par_ind(ii)),200);
        set(hh,'MarkerSize',2,'color',0.2*(1-ColPalette(1,:))) 
        set(gca,'TickLabelInterpreter','latex');    
        xlabel('Lag','Interpreter','latex');
        ylabel('ACF','Interpreter','latex');
        title(['DA: ',params{par_ind(ii)}])
    end
name = ['Figures/BurnIn_10000/BKM_Theta_acf_DA.png'];
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-dpng','-r0')


ff = figure(31);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.45]);
    set(gcf,'defaulttextinterpreter','latex');
    xxx = Results_Bin{1,1};  xxx = xxx.Theta;              
    for ii = 1:6
    	subplot(2,3,ii) 
        [~,~,~,hh] = autocorr(xxx(:,par_ind(ii)),200);
        set(hh,'MarkerSize',2,'color',ColPalette(5,:)) 
        set(gca,'TickLabelInterpreter','latex');    
        xlabel('Lag','Interpreter','latex');
        ylabel('ACF','Interpreter','latex');        
        title(['Bin10: ',params{par_ind(ii)}])       
    end
name = ['Figures/BurnIn_10000/BKM_Theta_acf_Bin10.png'];
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-dpng','-r0')


ff = figure(32);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.45]);
    set(gcf,'defaulttextinterpreter','latex');
    xxx = Results_Bin{1,2};  xxx = xxx.Theta;              
    for ii = 1:6
    	subplot(2,3,ii) 
        [~,~,~,hh] = autocorr(xxx(:,par_ind(ii)),200);
        set(hh,'MarkerSize',2,'color',ColPalette(6,:)) 
        set(gca,'TickLabelInterpreter','latex');    
        xlabel('Lag','Interpreter','latex');
        ylabel('ACF','Interpreter','latex');        
        title(['Bin20: ',params{par_ind(ii)}])       
    end
name = ['Figures/BurnIn_10000/BKM_Theta_acf_Bin20.png'];
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-dpng','-r0')

    
ff = figure(33);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.45]);
    set(gcf,'defaulttextinterpreter','latex');
    xxx = Results_Bin{1,3};  xxx = xxx.Theta;              
    for ii = 1:6
    	subplot(2,3,ii) 
        [~,~,~,hh] = autocorr(xxx(:,par_ind(ii)),200);
        set(hh,'MarkerSize',2,'color',ColPalette(7,:)) 
        set(gca,'TickLabelInterpreter','latex');    
        xlabel('Lag','Interpreter','latex');
        ylabel('ACF','Interpreter','latex');        
        title(['Bin30: ',params{par_ind(ii)}])       
    end    
name = ['Figures/BurnIn_10000/BKM_Theta_acf_Bin30.png'];
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-dpng','-r0')

    
    
ff = figure(21);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.45]);
    set(gcf,'defaulttextinterpreter','latex');
    xxx = Results_Adapt{1,1};  xxx = xxx.Theta;              
    for ii = 1:6
    	subplot(2,3,ii) 
        [~,~,~,hh] = autocorr(xxx(:,par_ind(ii)),200);
        set(hh,'MarkerSize',2,'color',ColPalette(2,:)) 
        set(gca,'TickLabelInterpreter','latex');    
        xlabel('Lag','Interpreter','latex');
        ylabel('ACF','Interpreter','latex');        
        title(['Adapt10: ',params{par_ind(ii)}])       
    end
name = ['Figures/BurnIn_10000/BKM_Theta_acf_Adapt10.png'];
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-dpng','-r0')

        
ff = figure(22);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.45]);
    set(gcf,'defaulttextinterpreter','latex');
    xxx = Results_Adapt{1,2};  xxx = xxx.Theta;              
    for ii = 1:6
    	subplot(2,3,ii) 
        [~,~,~,hh] = autocorr(xxx(:,par_ind(ii)),200);
        set(hh,'MarkerSize',2,'color',ColPalette(3,:)) 
        set(gca,'TickLabelInterpreter','latex');    
        xlabel('Lag','Interpreter','latex');
        ylabel('ACF','Interpreter','latex');        
        title(['Adapt20: ',params{par_ind(ii)}])       
    end
name = ['Figures/BurnIn_10000/BKM_Theta_acf_Adapt20.png'];
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-dpng','-r0')


ff = figure(23);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.45]);
    set(gcf,'defaulttextinterpreter','latex');
    xxx = Results_Adapt{1,3};  xxx = xxx.Theta;              
    for ii = 1:6
    	subplot(2,3,ii) 
        [~,~,~,hh] = autocorr(xxx(:,par_ind(ii)),200);
        set(hh,'MarkerSize',2,'color',ColPalette(4,:)) 
        set(gca,'TickLabelInterpreter','latex');    
        xlabel('Lag','Interpreter','latex');
        ylabel('ACF','Interpreter','latex');        
        title(['Adapt30: ',params{par_ind(ii)}])       
    end
name = ['Figures/BurnIn_10000/BKM_Theta_acf_Adapt30.png'];
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-dpng','-r0')            
    
% ACF NA    
L = 200;
ff = figure(2000);
%     set(gcf,'units','normalized','outerposition',[0.1 0.1 0.55 0.65]);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.85]);
     set(gcf,'defaulttextinterpreter','latex');  
     for ii = 1:9
        subplot(3,3,ii)  
        hold on
        xxx = squeeze(Results_DA.NN(2,:,:));    
        [aaa,~,bbb] = autocorr(xxx(4*ii,:),L);
        plot(0:L,aaa,'color',ColPalette(1,:)) 
  
        
        for jj = 1:3
            xxx = Results_Adapt{1,jj};  xxx = xxx.NN;              
            [aaa] = autocorr(xxx(4*ii,:),L);
            plot(0:L,aaa,'color',ColPalette(1+jj,:))               
        end
        for jj = 1:3
            xxx = Results_Bin{1,jj};  xxx = xxx.NN; 
            [aaa] = autocorr(xxx(4*ii,:),L);            
            plot(0:L,aaa,'color',ColPalette(4+jj,:))               
        end    
        
        xxx = Results_Exact.NN;
        [aaa] = autocorr(xxx(4*ii,:),L);
        plot(0:L,aaa,'color',ColPalette(8,:)) 
        
        hold off
        xlim([0,L])
        set(gca,'TickLabelInterpreter','latex');    
        xlabel('Lag','Interpreter','latex');
        ylabel('ACF','Interpreter','latex');
        title(['$N_{a,',num2str(ii*4),'}$'])
    end
    ll = legend(method);
    set(ll,'interpreter','latex','orientation','horizontal','FontSize',11);

name = ['Figures/BurnIn_10000/BKM_Na_acf_Comb.png'];
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-dpng','-r0')

name = ['Figures/BurnIn_10000/BKM_Na_acf_Comb.eps'];
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-depsc','-r0')



ff = figure(100);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.75]);
    set(gcf,'defaulttextinterpreter','latex');  
    xxx = squeeze(Results_DA.NN(2,:,:));    
    for ii = 1:9
        subplot(3,3,ii) 
        [~,~,~,hh] = autocorr(xxx(4*ii,:),200);
        set(hh,'MarkerSize',2,'color',0.2*(1-ColPalette(1,:))) 
        set(gca,'TickLabelInterpreter','latex');        
        set(gca,'TickLabelInterpreter','latex');    
        xlabel('Lag','Interpreter','latex');
        ylabel('ACF','Interpreter','latex');
        title(['DA: $N_{a,',num2str(ii*4),'}$'])
    end
name = ['Figures/BurnIn_10000/BKM_Na_acf_DA.png'];
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-dpng','-r0')      


ff = figure(231);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.75]);
    set(gcf,'defaulttextinterpreter','latex');
    xxx = Results_Bin{1,1};  xxx = xxx.NN;              
    for ii = 1:9
    	subplot(3,3,ii) 
        [~,~,~,hh] = autocorr(xxx(4*ii,:),200);
        set(hh,'MarkerSize',2,'color',ColPalette(5,:))         
        set(gca,'TickLabelInterpreter','latex');    
        xlabel('Lag','Interpreter','latex');
        ylabel('ACF','Interpreter','latex');        
        title(['Bin10: $N_{a,',num2str(ii*4),'}$'])       
    end
name = ['Figures/BurnIn_10000/BKM_Na_acf_Bin10.png'];
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-dpng','-r0')       
    

ff = figure(232);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.75]);
    set(gcf,'defaulttextinterpreter','latex');
    xxx = Results_Bin{1,2};  xxx = xxx.NN;              
    for ii = 1:9
    	subplot(3,3,ii) 
        [~,~,~,hh] = autocorr(xxx(4*ii,:),200);
        set(hh,'MarkerSize',2,'color',ColPalette(6,:))         
        set(gca,'TickLabelInterpreter','latex');    
        xlabel('Lag','Interpreter','latex');
        ylabel('ACF','Interpreter','latex');        
        title(['Bin20: $N_{a,',num2str(ii*4),'}$'])       
    end
name = ['Figures/BurnIn_10000/BKM_Na_acf_Bin20.png'];
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-dpng','-r0')       


ff = figure(233);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.75]);
    set(gcf,'defaulttextinterpreter','latex');
    xxx = Results_Bin{1,3};  xxx = xxx.NN;              
    for ii = 1:9
    	subplot(3,3,ii) 
        [~,~,~,hh] = autocorr(xxx(4*ii,:),200);
        set(hh,'MarkerSize',2,'color',ColPalette(7,:))         
        set(gca,'TickLabelInterpreter','latex');    
        xlabel('Lag','Interpreter','latex');
        ylabel('ACF','Interpreter','latex');        
        title(['Bin30: $N_{a,',num2str(ii*4),'}$'])       
    end
name = ['Figures/BurnIn_10000/BKM_Na_acf_Bin30.png'];
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-dpng','-r0')       
        


 
ff = figure(221);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.75]);
    set(gcf,'defaulttextinterpreter','latex');
    xxx = Results_Adapt{1,1};  xxx = xxx.NN;              
    for ii = 1:9
    	subplot(3,3,ii) 
        [~,~,~,hh] = autocorr(xxx(4*ii,:),200);
        set(hh,'MarkerSize',2,'color',ColPalette(2,:))         
        set(gca,'TickLabelInterpreter','latex');    
        xlabel('Lag','Interpreter','latex');
        ylabel('ACF','Interpreter','latex');        
        title(['Adapt10: $N_{a,',num2str(ii*4),'}$'])       
    end
name = ['Figures/BurnIn_10000/BKM_Na_acf_Adapt10.png'];
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-dpng','-r0')       

 
ff = figure(221);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.75]);
    set(gcf,'defaulttextinterpreter','latex');
    xxx = Results_Adapt{1,2};  xxx = xxx.NN;              
    for ii = 1:9
    	subplot(3,3,ii) 
        [~,~,~,hh] = autocorr(xxx(4*ii,:),200);
        set(hh,'MarkerSize',2,'color',ColPalette(3,:))         
        set(gca,'TickLabelInterpreter','latex');    
        xlabel('Lag','Interpreter','latex');
        ylabel('ACF','Interpreter','latex');        
        title(['Adapt20: $N_{a,',num2str(ii*4),'}$'])       
    end
name = ['Figures/BurnIn_10000/BKM_Na_acf_Adapt20.png'];
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-dpng','-r0')       

           
ff = figure(223);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.75]);
    set(gcf,'defaulttextinterpreter','latex');
    xxx = Results_Adapt{1,3};  xxx = xxx.NN;              
    for ii = 1:9
    	subplot(3,3,ii) 
        [~,~,~,hh] = autocorr(xxx(4*ii,:),200);
        set(hh,'MarkerSize',2,'color',ColPalette(4,:))         
        set(gca,'TickLabelInterpreter','latex');    
        xlabel('Lag','Interpreter','latex');
        ylabel('ACF','Interpreter','latex');        
        title(['Adapt30: $N_{a,',num2str(ii*4),'}$'])       
    end
name = ['Figures/BurnIn_10000/BKM_Na_acf_Adapt30.png'];
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-dpng','-r0')       
           
    
    
    
    
    xxx = Results_Exact.NN;   
    plot(mean(xxx,2),'color',ColPalette(8,:),'linewidth',1) 
    
    hold off 
    set(gca,'TickLabelInterpreter','latex');    
    xlim([1,36])
    ll = legend(method);
    set(ll,'interpreter','latex','orientation','horizontal','FontSize',11);
    
name = ['Figures/BurnIn_10000/BKM_N_acf_DA_B20.eps'];
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-depsc','-r0')

    
%%

NN_new = NN;
Theta_new = Theta;

load('Results/BurnIn_10000_1/BKM_DA_NRW_U.mat', 'NN','Theta')
NN_DA1 = squeeze(NN(2,:,:));
Theta_DA1 = Theta;
load('Results/BurnIn_10000_2/BKM_DA_NRW_U.mat', 'NN','Theta')
NN_DA2 = squeeze(NN(2,:,:));
Theta_DA2 = Theta;
load('Results/BurnIn_10000_3/BKM_DA_NRW_U.mat', 'NN','Theta')
NN_DA3 = squeeze(NN(2,:,:));
Theta_DA3 = Theta;

load('Results/BurnIn_10000_1/BKM_adapt_Nq30.mat', 'NN','Theta')
NN_A30_1 = NN;
Theta_A30_1 = Theta;
load('Results/BurnIn_10000_2/BKM_adapt_Nq30.mat', 'NN','Theta')
NN_A30_2 = NN;
Theta_A30_2 = Theta;
load('Results/BurnIn_10000_3/BKM_adapt_Nq30.mat', 'NN','Theta')
NN_A30_3 = NN;
Theta_A30_3 = Theta;

load('Results/BurnIn_10000_1/BKM_adapt_Nq20.mat', 'NN','Theta')
NN_A20_1 = NN;
Theta_A20_1 = Theta;
load('Results/BurnIn_10000_2/BKM_adapt_Nq20.mat', 'NN','Theta')
NN_A20_2 = NN;
Theta_A20_2 = Theta;
load('Results/BurnIn_10000_3/BKM_adapt_Nq20.mat', 'NN','Theta')
NN_A20_3 = NN;
Theta_A20_3 = Theta;

load('Results/BurnIn_10000_1/BKM_adapt_Nq10.mat', 'NN','Theta')
NN_A10_1 = NN;
Theta_A10_1 = Theta;
load('Results/BurnIn_10000_2/BKM_adapt_Nq10.mat', 'NN','Theta')
NN_A10_2 = NN;
Theta_A10_2 = Theta;
load('Results/BurnIn_10000_3/BKM_adapt_Nq10.mat', 'NN','Theta')
NN_A10_3 = NN;
Theta_A10_3 = Theta;


load('Results/BurnIn_10000_1/BKM_bin_Nbin30.mat', 'NN','Theta')
NN_B30_1 = NN;
Theta_B30_1 = Theta;
load('Results/BurnIn_10000_2/BKM_bin_Nbin30.mat', 'NN','Theta')
NN_B30_2 = NN;
Theta_B30_2 = Theta;
load('Results/BurnIn_10000_3/BKM_bin_Nbin30.mat', 'NN','Theta')
NN_B30_3 = NN;
Theta_B30_3 = Theta;


load('Results/BurnIn_10000_1/BKM_bin_Nbin20.mat', 'NN','Theta')
NN_B20_1 = NN;
Theta_B20_1 = Theta;
load('Results/BurnIn_10000_2/BKM_bin_Nbin20.mat', 'NN','Theta')
NN_B20_2 = NN;
Theta_B20_2 = Theta;
load('Results/BurnIn_10000_3/BKM_bin_Nbin20.mat', 'NN','Theta')
NN_B20_3 = NN;
Theta_B20_3 = Theta;


load('Results/BurnIn_10000_1/BKM_bin_Nbin10.mat', 'NN','Theta')
NN_B10_1 = NN;
Theta_B10_1 = Theta;
load('Results/BurnIn_10000_2/BKM_bin_Nbin10.mat', 'NN','Theta')
NN_B10_2 = NN;
Theta_B10_2 = Theta;
load('Results/BurnIn_10000_3/BKM_bin_Nbin10.mat', 'NN','Theta')
NN_B10_3 = NN;
Theta_B10_3 = Theta;


load('Results/BurnIn_10000_1/BKM_exact_Nmax690.mat', 'NN','Theta')
NN_Ex_1 = NN;
Theta_Ex_1 = Theta;
load('Results/BurnIn_10000_2/BKM_exact_Nmax690.mat', 'NN','Theta')
NN_Ex_2 = NN;
Theta_Ex_2 = Theta;
load('Results/BurnIn_10000_3/BKM_exact_Nmax690.mat', 'NN','Theta')
NN_Ex_3 = NN;
Theta_Ex_3 = Theta;



params = {'alpha1', 'alphaa', 'alphar', 'alphal', ...
            'beta1', 'betaa', 'betar', 'betal',...
            'sigy'};


ff = figure(100);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
for ii = 1:9
    subplot(3,3,ii)
    hold on
    plot(Theta_DA1(:,ii))
    plot(Theta_DA2(:,ii))
    plot(Theta_DA3(:,ii))
    hold off
    title(['DA ', params{ii}])
end  


ff = figure(113);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
for ii = 1:9
    subplot(3,3,ii)
    hold on    
%     plot(Theta_A30_1(:,ii))
    plot(Theta_A30_2(:,ii))
    plot(Theta_A30_3(:,ii))   

    plot(Theta_DA1(:,ii),'k')    
    
    hold off
    title(['A30 ', params{ii}])
end  


ff = figure(112);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
for ii = 1:9
    subplot(3,3,ii)
    hold on    
    plot(Theta_A20_1(:,ii))
%     plot(Theta_A20_2(:,ii))
    plot(Theta_A20_3(:,ii))    
    hold off
    title(['A20 ', params{ii}])
end  



ff = figure(111);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
for ii = 1:9
    subplot(3,3,ii)
    hold on    
    plot(Theta_A10_1(:,ii))
    plot(Theta_A10_2(:,ii))
    plot(Theta_A10_3(:,ii))    
    hold off
    title(['A10 ', params{ii}])
end  



ff = figure(123);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
for ii = 1:9
    subplot(3,3,ii)
    hold on    
    plot(Theta_B30_1(:,ii))
    plot(Theta_B30_2(:,ii))
    plot(Theta_B30_3(:,ii))    
    hold off
    title(['B30 ', params{ii}])
end  


ff = figure(122);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
for ii = 1:9
    subplot(3,3,ii)
    hold on    
    plot(Theta_B20_1(:,ii))
    plot(Theta_B20_2(:,ii))
    plot(Theta_B20_3(:,ii))    
    hold off
    title(['B20 ', params{ii}])
end  


ff = figure(121);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
for ii = 1:9
    subplot(3,3,ii)
    hold on    
    plot(Theta_B10_1(:,ii))
    plot(Theta_B10_2(:,ii))
    plot(Theta_B10_3(:,ii))    
    hold off
    title(['B10 ', params{ii}])
end



ff = figure(130);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
for ii = 1:9
    subplot(3,3,ii)
    hold on    
    plot(Theta_Ex_1(:,ii))
%     plot(Theta_Ex_2(:,ii))
    plot(Theta_Ex_3(:,ii))    
    hold off
    title(['Ex ', params{ii}])
end  

% Na


ff = figure(200);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
for ii = 1:9
    subplot(3,3,ii) 
    hold on
    plot(NN_DA1(4*ii,:))
    plot(NN_DA2(4*ii,:))
    plot(NN_DA3(4*ii,:)) 
    hold off
    title(['DA Na(',num2str(4*ii),')'])
end
 
ff = figure(213);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
for ii = 1:9
    subplot(3,3,ii) 
    hold on     
    plot(NN_A30_1(4*ii,:))
    plot(NN_A30_2(4*ii,:))
    plot(NN_A30_3(4*ii,:))         
    hold off
    title(['A30 Na(',num2str(4*ii),')'])
end
 

ff = figure(212);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
for ii = 1:9
    subplot(3,3,ii) 
    hold on     
    plot(NN_A20_1(4*ii,:))
%     plot(NN_A20_2(4*ii,:))
    plot(NN_A20_3(4*ii,:))         
    hold off
    title(['A20 Na(',num2str(4*ii),')'])
end
 
 

ff = figure(211);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
for ii = 1:9
    subplot(3,3,ii) 
    hold on     
    plot(NN_A10_1(4*ii,:))
    plot(NN_A10_2(4*ii,:))
    plot(NN_A10_3(4*ii,:))         
    hold off
    title(['A10 Na(',num2str(4*ii),')'])
end





ff = figure(223);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
for ii = 1:9
    subplot(3,3,ii) 
    hold on     
    plot(NN_B30_1(4*ii,:))
    plot(NN_B30_2(4*ii,:))
    plot(NN_B30_3(4*ii,:))         
    hold off
    title(['B30 Na(',num2str(4*ii),')'])
end
 

ff = figure(221);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
for ii = 1:9
    subplot(3,3,ii) 
    hold on     
    plot(NN_B20_1(4*ii,:))
    plot(NN_B20_2(4*ii,:))
    plot(NN_B20_3(4*ii,:))         
    hold off
    title(['B20 Na(',num2str(4*ii),')'])
end
 


ff = figure(221);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
for ii = 1:9
    subplot(3,3,ii) 
    hold on     
    plot(NN_B10_1(4*ii,:))
    plot(NN_B10_2(4*ii,:))
    plot(NN_B10_3(4*ii,:)) 
    
    plot(NN_DA1(4*ii,:),'k')   
    
    
    hold off
    title(['B10 Na(',num2str(4*ii),')'])
end
 


ff = figure(230);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
for ii = 1:9
    subplot(3,3,ii) 
    hold on     
    plot(NN_Ex_1(4*ii,:))
    plot(NN_Ex_2(4*ii,:))
    plot(NN_Ex_3(4*ii,:))         
    hold off
    title(['Ex Na(',num2str(4*ii),')'])
end
 


%% SINGLE PLOT
ff = figure(1);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
for ii = 1:9
    subplot(3,3,ii)
    hold on
    plot(Theta(:,ii))    
    hold off
    title(params{ii})
end  


ff = figure(2);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
for ii = 1:9
    subplot(3,3,ii) 
    hold on
    plot(NN(4*ii,:))
    hold off
    title(['Na(',num2str(4*ii),')'])
end
 

