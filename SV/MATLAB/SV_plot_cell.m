% clear all
arg0 = 2;
DATA_NAMES = {'_GSPC','_IBM','_MSFT'};
data_name = DATA_NAMES{arg0}; %'_GSPC';  
y = load('other/Perc_Rets_GSPC_IBM_MSFT.csv');
y = y(:,arg0); % arg0 = 1 = GSPC;  arg0 = 2 = IBM
T = length(y); 
time = [2000, 2018];
% path = 'Results/Empirical/';
path = 'Results/Empirical_new/';

name = [path,'SV_results_param_DA_RW_eff',data_name,'.mat'];
Results_DA = load(name);
 
method = {'DA','Adapt10','Adapt20','Adapt30','Bin10','Bin20','Bin30'};
            
for ii = 1:3
    name = [path,'SV_results_param_HMM_adapt_eff_Nbin',num2str(10*ii),data_name,'.mat'];
    Results_Adapt{ii} = load(name);
    
    name = [path,'SV_results_param_HMM_eff_Nbin',num2str(10*ii),data_name,'.mat'];
    try
        Results_Bin{ii} = load(name);
    catch
        Results_Bin{ii} = [];
        Results_Bin{ii}.time_HMM_eff = []; 
        Results_Bin{ii}.theta_HMM_eff = []; 
        Results_Bin{ii}.H_subset_HMM_eff = []; 
        Results_Bin{ii}.ESS_theta_HMM_eff_sig = []; 
        Results_Bin{ii}.ESS_H_subset_HMM_eff_sig = [];         
        method{4+ii} = [];  
    end
end

method = method(~cellfun('isempty',method));
 
%% PRINT TABLES
Time = [Results_DA.time_DA_RW_eff;
        Results_Adapt{1}.time_HMM_adapt_eff;
        Results_Adapt{2}.time_HMM_adapt_eff;
        Results_Adapt{3}.time_HMM_adapt_eff;
        Results_Bin{1}.time_HMM_eff;
        Results_Bin{2}.time_HMM_eff;
        Results_Bin{3}.time_HMM_eff];

Time(:,2) = Time/Time(1);

TH_mean = [mean(Results_DA.theta_DA_RW_eff);
        mean(Results_Adapt{1}.theta_HMM_adapt_eff);
        mean(Results_Adapt{2}.theta_HMM_adapt_eff);
        mean(Results_Adapt{3}.theta_HMM_adapt_eff);
%         mean(Results_Bin{1}.theta_HMM_eff);
        mean(Results_Bin{2}.theta_HMM_eff);
        mean(Results_Bin{3}.theta_HMM_eff)];


TH_std = [std(Results_DA.theta_DA_RW_eff);
        std(Results_Adapt{1}.theta_HMM_adapt_eff);
        std(Results_Adapt{2}.theta_HMM_adapt_eff);
        std(Results_Adapt{3}.theta_HMM_adapt_eff);
        std(Results_Bin{2}.theta_HMM_eff);
        std(Results_Bin{3}.theta_HMM_eff)];

TH_ESS = [Results_DA.ESS_theta_DA_RW_eff_sig;
        Results_Adapt{1}.ESS_theta_HMM_adapt_eff_sig;
        Results_Adapt{2}.ESS_theta_HMM_adapt_eff_sig;
        Results_Adapt{3}.ESS_theta_HMM_adapt_eff_sig;
%         Results_Bin{1}.ESS_theta_HMM_eff_sig;
        Results_Bin{2}.ESS_theta_HMM_eff_sig;
        Results_Bin{3}.ESS_theta_HMM_eff_sig];
 
H_mean = [mean(Results_DA.H_subset_DA_RW_eff(:,10*(1:9)));
        mean(Results_Adapt{1}.H_subset_HMM_adapt_eff(:,10*(1:9)));
        mean(Results_Adapt{2}.H_subset_HMM_adapt_eff(:,10*(1:9)));
        mean(Results_Adapt{3}.H_subset_HMM_adapt_eff(:,10*(1:9)));
        mean(Results_Bin{2}.H_subset_HMM_eff(:,10*(1:9)));
        mean(Results_Bin{3}.H_subset_HMM_eff(:,10*(1:9)))];


H_std = [std(Results_DA.H_subset_DA_RW_eff(:,10*(1:9)));
        std(Results_Adapt{1}.H_subset_HMM_adapt_eff(:,10*(1:9)));
        std(Results_Adapt{2}.H_subset_HMM_adapt_eff(:,10*(1:9)));
        std(Results_Adapt{3}.H_subset_HMM_adapt_eff(:,10*(1:9)));
        std(Results_Bin{2}.H_subset_HMM_eff(:,10*(1:9)));
        std(Results_Bin{3}.H_subset_HMM_eff(:,10*(1:9)))];
    
H_ESS = [Results_DA.ESS_H_DA_RW_eff_sig(10*(1:9));
        Results_Adapt{1}.ESS_H_HMM_adapt_eff_sig(10*(1:9));
        Results_Adapt{2}.ESS_H_HMM_adapt_eff_sig(10*(1:9));
        Results_Adapt{3}.ESS_H_HMM_adapt_eff_sig(10*(1:9));
        Results_Bin{2}.ESS_H_HMM_eff_sig(10*(1:9));
        Results_Bin{3}.ESS_H_HMM_eff_sig(10*(1:9))];    

 print_table_results_theta_and_h(path, data_name, method,...
    TH_mean, TH_std, TH_ESS,  H_mean, H_std, H_ESS, Time);

%% PLOT
    
ColPalette = [0 0 0
    255 121 75
    255 51 0
    153 0 0 
    62 154 222
    0 102 204
    51 51 102]/255;



%% Trace Plots params 
params = {'\mu','\phi','\sigma^2'};
ff = figure(16);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.8]);
for ii = 1:3
    subplot(3,1,ii)
    hold on    
    plot(Results_DA.theta_DA_RW_eff(:,ii),'color', ColPalette(1,:))
    for jj = 1:3
        plot(Results_Adapt{jj}.theta_HMM_adapt_eff(:,ii),'color',ColPalette(1+jj,:))     
    end
    for jj = 2:3
        plot(Results_Bin{jj}.theta_HMM_eff(:,ii),'color',ColPalette(4+jj,:))     
    end
    hold off
    title(params{ii}) 
end 
ll = legend(method);
set(ll,'interpreter','latex','orientation','horizontal','FontSize',11);

name = ['figures/Empirical/SV_trace_params_comb',data_name,'.png']; %eps
set(gcf,'PaperPositionMode','auto');
% print(ff,name,'-depsc','-r0')
print(ff,name,'-dpng','-r0')



%% Trace Plots H
ind_h_sel = 2:50:T;
ff = figure(2);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
for ii = 1:9
    subplot(3,3,ii)
    hold on    
    plot(Results_DA.H_subset_DA_RW_eff(:,10*ii),'color', ColPalette(1,:))
    for jj = 1:3
        plot(Results_Adapt{jj}.H_subset_HMM_adapt_eff(:,10*ii),'color',ColPalette(1+jj,:))     
    end
    for jj = 2:3
        plot(Results_Bin{jj}.H_subset_HMM_eff(:,10*ii),'color',ColPalette(4+jj,:))     
    end
    hold off 
    title(['$t = ',num2str(ind_h_sel(10*ii)),'$'],'interpreter','latex')   
end
ll = legend(method);
set(ll,'interpreter','latex','orientation','horizontal','FontSize',11);

name = ['figures/Empirical/SV_trace_H_comb',data_name,'.png']; %eps
set(gcf,'PaperPositionMode','auto');
% print(ff,name,'-depsc','-r0')
print(ff,name,'-dpng','-r0')


%% Mean H
ff = figure(3);
set(gcf,'units','normalized','outerposition',[0.25 0.25 0.45 0.45]);
hold on
xx = linspace(time(1,1),time(1,2),T);
plot(xx,y, 'color',[0.8 0.8 0.8]);
    plot(xx,Results_DA.mean_H_DA_RW_eff,'color', ColPalette(1,:))
    for jj = 1:3
        plot(xx(2:2:T),Results_Adapt{jj}.mean_H_HMM_adapt_eff(2:2:T),'color',ColPalette(1+jj,:))     
    end
    for jj = 2:3
        plot(xx(2:2:T), Results_Bin{jj}.mean_H_HMM_eff(2:2:T),'color',ColPalette(4+jj,:))     
    end
hold off 
set(gca,'XLim',time);
set(gca,'TickLabelInterpreter','latex');    
ll = legend(['y',method]);
set(ll,'interpreter','latex','orientation','horizontal','FontSize',11);

name = ['figures/Empirical/SV_y_mean_H_comb',data_name,'.png']; %eps
set(gcf,'PaperPositionMode','auto');
% print(ff,name,'-depsc','-r0')
print(ff,name,'-dpng','-r0')
 


ff = figure(34);
set(gcf,'units','normalized','outerposition',[0.25 0.25 0.45 0.45]);
hold on
xx = linspace(time(1,1),time(1,2),T);
% plot(xx,y, 'color',[0.8 0.8 0.8]);
    plot(xx,Results_DA.mean_H_DA_RW_eff,'color', ColPalette(1,:))
    for jj = 1:3
        plot(xx(2:2:T),Results_Adapt{jj}.mean_H_HMM_adapt_eff(2:2:T),'color',ColPalette(1+jj,:))     
    end
    for jj = 2:3
        plot(xx(2:2:T), Results_Bin{jj}.mean_H_HMM_eff(2:2:T),'color',ColPalette(4+jj,:))     
    end
hold off 
set(gca,'XLim',time);
set(gca,'TickLabelInterpreter','latex');    
ll = legend(method);
set(ll,'interpreter','latex','orientation','horizontal','FontSize',11);

name = ['figures/Empirical/SV_mean_H_comb',data_name,'.png']; %eps
set(gcf,'PaperPositionMode','auto');
% print(ff,name,'-depsc','-r0')
print(ff,name,'-dpng','-r0')





%% ACF params
params = {'$\mu$','$\phi$','$\sigma^2$'};
L = [20, 200, 200];
ff = figure(1000);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.4]); 
set(gcf,'defaulttextinterpreter','latex');  

for ii = 1:3
    subplot(1,3,ii)
    hold on    
    [aaa,~,bbb] = autocorr(Results_DA.theta_DA_RW_eff(:,ii),L(ii));
    aaa_min = min([aaa; bbb(2)]);

    plot(0:L(ii),aaa,'color', ColPalette(1,:))
    for jj = 1:3
        [aaa] = autocorr(Results_Adapt{jj}.theta_HMM_adapt_eff(:,ii),L(ii));
        aaa_min = min([aaa_min;aaa]);       
        plot(0:L(ii),aaa,'color',ColPalette(1+jj,:))               
    end
    for jj = 2:3
        [aaa] = autocorr(Results_Bin{jj}.theta_HMM_eff(:,ii),L(ii));
        aaa_min = min([aaa_min;aaa]);       
        plot(0:L(ii),aaa,'color',ColPalette(4+jj,:))               
    end
    plot(0:L(ii),bbb(1)+ 0*(0:L(ii)),'--b')
    plot(0:L(ii),bbb(2)+ 0*(0:L(ii)),'--b')
    plot(0:L(ii),0*(0:L(ii)),'-k');
   
    hold off
    title(params{ii}) 
    
    xlim([0,L(ii)])
    ylim([aaa_min,1])
    
    set(gca,'TickLabelInterpreter','latex');    
    xlabel('Lag','Interpreter','latex');
    ylabel('ACF','Interpreter','latex');
    

    pos=get(gca,'position');  % get current axes position vector
    dh=0.1;                   % guess for height adjustment value
    pos(2)=pos(2)+dh;pos(4)=pos(4)-dh;  % raise bottom, reduce height
    set(gca,'position',pos)   % resize with that modification         
end 
ll = legend(method);
set(ll,'interpreter','latex','orientation','horizontal','FontSize',11);

name = ['figures/Empirical/Theta_acf_Comb',data_name,'.png']; %eps
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-dpng','-r0') 





%% ACF Plots H
L = 200;
ind_h_sel = 2:50:T;
ff = figure(27);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.7]);
for ii = 1:9
    subplot(3,3,ii)
    hold on    
    [aaa, ~, bbb] = autocorr(Results_DA.H_subset_DA_RW_eff(:,10*ii),L);
    plot(0:L,aaa,'color', ColPalette(1,:))
    aaa_min = min([aaa; bbb(2)]);
    for jj = 1:3
        aaa = autocorr(Results_Adapt{jj}.H_subset_HMM_adapt_eff(:,10*ii),L);
        plot(0:L,aaa,'color',ColPalette(1+jj,:)) 
        aaa_min = min([aaa_min;aaa]);
    end
    for jj = 2:3
        aaa = autocorr(Results_Bin{jj}.H_subset_HMM_eff(:,10*ii),L);
        plot(0:L,aaa,'color',ColPalette(4+jj,:))  
        aaa_min = min([aaa_min;aaa]);        
    end
    
    plot(0:L,bbb(1)+ 0*(0:L),'--b')
    plot(0:L,bbb(2)+ 0*(0:L),'--b')
    plot(0:L,0*(0:L),'-k');

    hold off 

    title(['$t = ',num2str(ind_h_sel(10*ii)),'$'],'interpreter','latex') 
    xlim([0,L])
    ylim([aaa_min,1])
    set(gca,'TickLabelInterpreter','latex');    
    xlabel('Lag','Interpreter','latex');
    ylabel('ACF','Interpreter','latex');  
end
ll = legend(method);
set(ll,'interpreter','latex','orientation','horizontal','FontSize',11);

name = ['figures/Empirical/H_acf_Comb',data_name,'.png']; %eps
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-dpng','-r0') 


%% ESS H

ff = figure(16);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.5]);
hold on    
plot(Results_DA.ESS_H_DA_RW_eff_sig,'color', ColPalette(1,:))
for jj = 1:3
    plot(Results_Adapt{jj}.ESS_H_HMM_adapt_eff_sig,'color',ColPalette(1+jj,:))     
end
for jj = 2:3
    plot(Results_Bin{jj}.ESS_H_HMM_eff_sig,'color',ColPalette(4+jj,:))     
end
hold off
 
ll = legend(method);
set(ll,'interpreter','latex','orientation','horizontal','FontSize',11);


figure(23)
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.5]);
XXX = zeros(6,91);
hold on    
XXX(1,:) = Results_DA.ESS_H_DA_RW_eff_sig;
for jj = 1:3
    XXX(1+jj,:) = Results_Adapt{jj}.ESS_H_HMM_adapt_eff_sig;
end
for jj = 2:3
    XXX(3+jj,:) = Results_Bin{jj}.ESS_H_HMM_eff_sig;
end
hold off



ff = figure(23);
xxx = 2:50:length(y);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.4]);
for ii=1:6
    subplot(1,6,ii)
    hold on
    bar(xxx,XXX(ii,:),'facecolor',ColPalette(ii,:));
    plot(xxx,mean(XXX(ii,:)) + 0*xxx,'r','linewidth',2)
    hold off
    ylim([0 ,600])
    xlim([xxx(1), xxx(91)])
    xlabel(method{ii},'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex');    
end


name = ['figures/Empirical/H_ESS_Comb',data_name,'.png']; %eps
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-dpng','-r0') 

