T = 2000;
M = 10000;
params = {'\mu','\phi','\sigma^2'};


%% Volatility
ff = figure(1);
set(gcf,'units','normalized','outerposition',[0.25 0.25 0.45 0.45]);
% ax1 = axes('Position',[0.25 0.25 0.3 0.3],'Visible','on');
hold on
plot(h_true,'r')
% plot(h,'b')
% plot(mean(H_HMM_shift((5000+1):M,:),1),'m')
% plot(2:2:T,mean(H_HMM_adapt(:,2:2:T),1),'b')
% plot(2:2:T,mean_H_HMM_adapt_eff(:,2:2:T),'b')
plot(mean_H_DA_RW_eff)
% plot(mean_H_DA_RW,'c')
plot(mean_H_DA_RW)
plot(2:2:T,mean_H_HMM(:,2:2:T))
plot(2:2:T,mean_H_HMM_adapt(:,2:2:T))
% plot(2:2:T,mean_H_HMM_adapt_eff(:,2:2:T),'r')
% plot(mean_H_HMM_shift,'m')
hold off
plotTickLatex2D;%('FontSize',12);
leg = legend('true','DA','DA RW','HMM','HMM shift');
set(leg,'Interpreter','latex');%,'FontSize',12,'position',[0.8 0.42 0.18 0.16])

if save_on
    name = ['figures/SV_volatility.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print_fail = 1;
    while print_fail 
        try                 
%             print(ff,name,'-depsc','-r0')
            print(ff,name,'-dpng','-r0')
            print_fail = 0;
        catch
            print_fail = 1;
        end
    end
end

ff = figure(11);
set(gcf,'units','normalized','outerposition',[0.25 0.25 0.45 0.45]);
% ax1 = axes('Position',[0.25 0.25 0.3 0.3],'Visible','on');
hold on
plot(h_true,'r')
% plot(h,'b')
% plot(mean_H_DA,'g')
% plot(mean_H_HMM_shift,'m')
plot(2:2:T,mean_H_HMM(2:2:T),'b')
hold off
plotTickLatex2D;%('FontSize',12);
leg = legend('true','DA','HMM shift','HMM');
set(leg,'Interpreter','latex');%,'FontSize',12,'position',[0.8 0.42 0.18 0.16])


%% Trace Plots
ff = figure(2);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
for t = 1:9
    subplot(3,3,t)
    hold on
%     plot(H_subset_DA(:,t)*0 + h_true((t+1)*100),'r')
%     plot(H_subset_DA(:,t),'b')
    plot(H_subset_DA_RW_eff(:,10*t))
%     plot(H_subset_HMM_eff(:,t),'r')
%     plot(H_subset_HMM_adapt(:,t),'g')
    plot(H_subset_HMM_adapt_eff(:,10*t))
    plot(H_subset_HMM_eff(:,10*t))
    hold off    
%     XL = xlabel('MCMC iteration');
% %     pos = get(XL,'pos'); % Read position [x y z]
% %     YL = get(gca,'ylim');
% %     pos(2) = 1.01*YL(1); 
% %     set(XL,'pos',pos) % Move label up 
%     title(['Time = ',num2str(ind_h_sel(10*t))])   
% %     plotTickLatex2D;%('FontSize',12);
end
leg = legend('DA','HMM adapt', 'HMM fixed');
set(leg,'Interpreter','latex');%,'FontSize',12,'position',[0.8 0.42 0.18 0.16])

if save_on
    name = ['figures/SV_trace.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print_fail = 1;
    while print_fail 
        try                 
%             print(ff,name,'-depsc','-r0')
            print(ff,name,'-dpng','-r0')
            print_fail = 0;
        catch
            print_fail = 1;
        end
    end
end


%% ACF
if ~exist('autocorr','builtin')
    autocorr = @(xx,ll) acf(xx,ll);
end


ff = figure(3);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
for t = 1:9
    subplot(3,3,t)
    autocorr(H_DA_RW_eff(:,10*t),200)
%     autocorr(H_subset_DA(:,t),200)
%     set(gca,'ylabel',[])    
    title(['Time = ',num2str(ind_h_sel(10*t))])
%     plotTickLatex2D;%('FontSize',12);
end
suptitle('Even DA')
 
if save_on
    name = ['figures/SV_DA_acf.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print_fail = 1;
    while print_fail 
        try                 
%             print(ff,name,'-depsc','-r0')
            print(ff,name,'-dpng','-r0')
            print_fail = 0;
        catch
            print_fail = 1;
        end
    end
end


ff = figure(33);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
for t = 1:9
    subplot(3,3,t)
    autocorr(H_subset_DA_RW_eff(:,t),200)
%     set(gca,'ylabel',[])    
    title(['Time = ',num2str(100*(t+1))])
%     plotTickLatex2D;%('FontSize',12);
end
suptitle('Even DA RW')
 
if save_on
    name = ['figures/SV_DA_RW_acf.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print_fail = 1;
    while print_fail 
        try                 
%             print(ff,name,'-depsc','-r0')
            print(ff,name,'-dpng','-r0')
            print_fail = 0;
        catch
            print_fail = 1;
        end
    end
end


ff = figure(4);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
for t = 1:9
    subplot(3,3,t)
    autocorr(H_subset_HMM_eff(:,t),200)
%     set(gca,'ylabel',[])    
    title(['Time = ',num2str(100*(t+1))])
%     plotTickLatex2D;%('FontSize',12);
end
suptitle('Even HMM fixed')
 
if save_on
    name = ['figures/SV_HMM_acf.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print_fail = 1;
    while print_fail 
        try                 
%             print(ff,name,'-depsc','-r0')
            print(ff,name,'-dpng','-r0')
            print_fail = 0;
        catch
            print_fail = 1;
        end
    end
end


ff = figure(5);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
for t = 1:9
    subplot(3,3,t)
    autocorr(H_subset_HMM_adapt(:,t),200)
%     set(gca,'ylabel',[])    
    title(['Time = ',num2str(100*(t+1))])
%     plotTickLatex2D;%('FontSize',12);
end
suptitle('Even HMM adapt')
 
if save_on
    name = ['figures/SV_HMM_adapt_acf.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print_fail = 1;
    while print_fail 
        try                 
%             print(ff,name,'-depsc','-r0')
            print(ff,name,'-dpng','-r0')
            print_fail = 0;
        catch
            print_fail = 1;
        end
    end
end


ff = figure(6);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
for t = 1:9
    subplot(3,3,t)
    autocorr(H_subset_HMM_shift2(1:2:M,t),40)
    set(gca,'ylabel',[])    
    title(['Time = ',num2str(100*(t+1)-1)])
%     plotTickLatex2D;%('FontSize',12);
end
suptitle('Odd HMM shift')

if save_on
    name = ['figures/SV_HMM_shift_odd_acf.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print_fail = 1;
    while print_fail 
        try                 
%             print(ff,name,'-depsc','-r0')
            print(ff,name,'-dpng','-r0')
            print_fail = 0;
        catch
            print_fail = 1;
        end
    end
end


ff = figure(7);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
for t = 1:9
    subplot(3,3,t)
    autocorr(H_subset_HMM_adapt(:,t),200)
%     set(gca,'ylabel',[])    
    title(['Time = ',num2str(100*(t+1))])
%     plotTickLatex2D;%('FontSize',12);
end
suptitle('Even HMM adaptive')



if save_on
    name = ['figures/SV_HMM_adapt.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print_fail = 1;
    while print_fail 
        try                 
%             print(ff,name,'-depsc','-r0')
            print(ff,name,'-dpng','-r0')
            print_fail = 0;
        catch
            print_fail = 1;
        end
    end
end



%% Trace Plots params 
params = {'\mu','\phi','\sigma^2'};
ff = figure(16);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.8]);
for ii = 1:3
    subplot(3,1,ii)
    hold on    
%     plot(theta_DA_RW(:,ii))
    plot(theta_DA_RW_eff(:,ii))
%     plot(theta_HMM(:,ii))       
%     plot(theta_HMM_eff(:,ii))      
%     plot(theta_HMM_adapt(:,ii))      
    plot(theta_HMM_adapt_eff(:,ii))      
%     plot(theta_true(ii) + 0*theta_HMM_adapt_eff(:,ii),'r')
%     plot(theta_true(ii) + 0*theta_DA_RW_eff(:,ii),'r')
%     plot(mean(theta_DA_RW(:,ii)) + 0*theta_DA_RW(:,ii),'r')
%     plot(mean(theta_DA_RW_eff(:,ii)) + 0*theta_DA_RW_eff(:,ii),'r')
    hold off
%     set(gca,'ylabel',[])    
    title(params{ii})
%     plotTickLatex2D;%('FontSize',12);
 end
suptitle('Trace plots params')
% leg = legend('DA','HMM fixed','HMM adapt','true');
% leg = legend('DA','DA RW','HMM','HMM shift','true');
set(leg,'Interpreter','latex');%,'FontSize',12,'position',[0.8 0.42 0.18 0.16])

if save_on
    name = ['figures/SV_trace_params.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print_fail = 1;
    while print_fail 
        try                 
%             print(ff,name,'-depsc','-r0')
            print(ff,name,'-dpng','-r0')
            print_fail = 0;
        catch
            print_fail = 1;
        end
    end
end


%% ACF params 
ff = figure(16);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.8]);
for ii = 1:3
    subplot(3,1,ii)
    autocorr(theta_DA_RW_eff(:,ii),40)
%     autocorr(theta_HMM_eff(:,ii),40)
%     set(gca,'ylabel',[])
%     title(params{ii})
 end
suptitle('ACF params DA')

if save_on
    name = ['figures/SV_DA_acf_param.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print_fail = 1;
    while print_fail 
        try                 
%             print(ff,name,'-depsc','-r0')
            print(ff,name,'-dpng','-r0')
            print_fail = 0;
        catch
            print_fail = 1;
        end
    end
end



ff = figure(18);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.8]);
for ii = 1:3
    subplot(3,1,ii)
    autocorr(theta_DA_RW(:,ii),200)
%     autocorr(theta_DA_RW_eff(:,ii),40)
%     set(gca,'ylabel',[])
%     title(params{ii})
 end
suptitle('ACF params DA RW')

if save_on
    name = ['figures/SV_DA_RW_acf_param.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print_fail = 1;
    while print_fail 
        try                 
%             print(ff,name,'-depsc','-r0')
            print(ff,name,'-dpng','-r0')
            print_fail = 0;
        catch
            print_fail = 1;
        end
    end
end



ff = figure(19);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.8]);
for ii = 1:3
    subplot(3,1,ii)
    autocorr(theta_HMM_eff(:,ii),40)
%     set(gca,'ylabel',[])
%     title(params{ii})
 end
suptitle('ACF params HMM')

if save_on
    name = ['figures/SV_HMM_acf_param.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print_fail = 1;
    while print_fail 
        try                 
%             print(ff,name,'-depsc','-r0')
            print(ff,name,'-dpng','-r0')
            print_fail = 0;
        catch
            print_fail = 1;
        end
    end
end



ff = figure(20);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.8]);
for ii = 1:3
    subplot(3,1,ii)
    autocorr(theta_HMM_adapt(:,ii),40)
    set(gca,'ylabel',[])
    title(params{ii})
 end
suptitle('ACF params HMM adaptive')

if save_on
    name = ['figures/SV_HMM_adapt_acf_param.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print_fail = 1;
    while print_fail 
        try                 
%             print(ff,name,'-depsc','-r0')
            print(ff,name,'-dpng','-r0')
            print_fail = 0;
        catch
            print_fail = 1;
        end
    end
end


%% ESS params

ESS_theta_40 = [ESS_theta_DA_40;
                ESS_theta_DA_RW_eff_40;
                ESS_theta_HMM_eff_40;
                ESS_theta_HMM_shift_40];

ESS_theta_100 = [ESS_theta_DA_100;
                ESS_theta_DA_RW_eff_100;
                ESS_theta_HMM_eff_100;
                ESS_theta_HMM_shift_100];         

ESS_theta_1000 = [ESS_theta_DA_1000;
                ESS_theta_DA_RW_eff_1000;
                ESS_theta_HMM_eff_1000;
                ESS_theta_HMM_shift_1000];
            

ESS_theta_sig = [ESS_theta_DA_sig;
                ESS_theta_DA_RW_eff_sig;
                ESS_theta_HMM_eff_sig];

ff = figure(99);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.8]);
bbb = bar(ESS_theta_sig');
bbb(1).FaceColor  = 'g';
bbb(2).FaceColor  = 'c';
bbb(3).FaceColor  = 'b'; 
set(gca, 'XTickLabel', params)
title('ESS significant')

leg = legend('DA','DA RW','HMM');
set(leg,'Interpreter','latex');%,'FontSize',12,'position',[0.8 0.42 0.18 0.16])
 
if save_on
    name = ['figures/SV_ESS_param_sig.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print_fail = 1;
    while print_fail 
        try                 
%             print(ff,name,'-depsc','-r0')
            print(ff,name,'-dpng','-r0')
            print_fail = 0;
        catch
            print_fail = 1;
        end
    end
end
           
            
ff = figure(100);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.8]);
subplot(3,1,1)
bbb = bar(ESS_theta_40');
bbb(1).FaceColor  = 'g';
bbb(2).FaceColor  = 'c';
bbb(3).FaceColor  = 'b';
bbb(4).FaceColor  = 'm';
set(gca, 'XTickLabel', params)
title('ESS 40')


subplot(3,1,2)
bbb = bar(ESS_theta_100');
bbb(1).FaceColor  = 'g';
bbb(2).FaceColor  = 'c';
bbb(3).FaceColor  = 'b';
bbb(4).FaceColor  = 'm';
set(gca, 'XTickLabel', params)
title('ESS 100')


subplot(3,1,3)
bbb = bar(ESS_theta_1000');
bbb(1).FaceColor  = 'g';
bbb(2).FaceColor  = 'c';
bbb(3).FaceColor  = 'b';
bbb(4).FaceColor  = 'm';
set(gca, 'XTickLabel', params)
title('ESS 1000')

leg = legend('DA','DA RW','HMM','HMM shift');
set(leg,'Interpreter','latex');%,'FontSize',12,'position',[0.8 0.42 0.18 0.16])

if save_on
    name = ['figures/SV_ESS_param.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print_fail = 1;
    while print_fail 
        try                 
%             print(ff,name,'-depsc','-r0')
            print(ff,name,'-dpng','-r0')
            print_fail = 0;
        catch
            print_fail = 1;
        end
    end
end


%% ESS H

ESS_H_40 = [ESS_H_DA_40;
                ESS_H_DA_RW_eff_40;
                ESS_H_HMM_eff_40;
                ESS_H_HMM_shift_40];

ESS_H_100 = [ESS_H_DA_100;
                ESS_H_DA_RW_eff_100;
                ESS_H_HMM_eff_100;
                ESS_H_HMM_shift_100];         

ESS_H_1000 = [ESS_H_DA_1000;
                ESS_H_DA_RW_eff_1000;
                ESS_H_HMM_eff_1000;
                ESS_H_HMM_shift_1000];
            

ESS_H_sig = [ESS_H_DA_sig;
                ESS_H_DA_RW_eff_sig;
                ESS_H_HMM_eff_sig];       
            
ff = figure(101);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.8]);
subplot(3,1,1) 
hold on
plot(ESS_H_40(1,2:2:T),'g')
plot(ESS_H_40(2,2:2:T),'c')
plot(ESS_H_40(3,2:2:T),'b')
hold off
title('ESS 40')


subplot(3,1,2) 
hold on
plot(ESS_H_100(1,2:2:T),'g')
plot(ESS_H_100(2,2:2:T),'c')
plot(ESS_H_100(3,2:2:T),'b')
hold off
title('ESS 100')


subplot(3,1,3) 
hold on
plot(ESS_H_sig(1,2:2:T),'g')
plot(ESS_H_sig(2,2:2:T),'c')
plot(ESS_H_sig(3,2:2:T),'b')
hold off
title('ESS significant')

leg = legend('DA','DA RW','HMM');
set(leg,'Interpreter','latex');%,'FontSize',12,'position',[0.8 0.42 0.18 0.16])


 
if save_on
    name = ['figures/SV_ESS_h.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print_fail = 1;
    while print_fail 
        try                 
%             print(ff,name,'-depsc','-r0')
            print(ff,name,'-dpng','-r0')
            print_fail = 0;
        catch
            print_fail = 1;
        end
    end
end

mean(ESS_H_DA_sig)
mean(ESS_H_DA_RW_eff_sig)
mean(ESS_H_HMM_eff_sig(2:2:T))

ESS_theta_DA_sig
ESS_theta_DA_RW_eff_sig
ESS_theta_HMM_eff_sig



%% DIFFERENT NO. OF BINS

N_bins = [10,15,20,25,30];
params = {'\mu','\phi','\sigma^2'};


ff = figure(1600);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
for jj = 1:5
    N_bin = N_bins(jj);
    load(['Results/Simulation/SV_results_param_HMM_adapt_eff_Nbin',num2str(N_bin),'.mat'])
    for ii = 1:3
        subplot(3,5,(ii-1)*5+jj)
        hold on         
        plot(theta_HMM_adapt_eff(:,ii),'b')      
        plot(theta_true(ii) + 0*theta_HMM_adapt_eff(:,ii),'r')
        hold off
    %     set(gca,'ylabel',[])    
        title([params{ii},' ',', Nbin = ',num2str(N_bin)])
    %     plotTickLatex2D;%('FontSize',12);
    end
end
suptitle('Trace plots params')
name = ['figures/SV_HMM_adapt_param_trace_bins_comp.png']; %eps
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-dpng','-r0')

ff = figure(2600);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
 for jj = 1:5
    N_bin = N_bins(jj);
    load(['Results/Simulation/SV_results_param_HMM_adapt_eff_Nbin',num2str(N_bin),'.mat'])
    for ii = 1:5
        subplot(5,5,(ii-1)*5+jj)
        hold on         
        plot(H_subset_HMM_adapt_eff(:,2*ii-1),'b')      
        plot(H_subset_HMM_adapt_eff(:,2*ii-1)*0 + h_true((2*ii-1+1)*100),'r')
        hold off
    %     set(gca,'ylabel ',[])
        title(['Time = ',num2str(100*(2*ii-1)),', Nbin = ',num2str(N_bin)])
    %     plotTickLatex2D;%('FontSize',12);
    end
end
suptitle('Trace plots h')
name = ['figures/SV_HMM_adapt_h_trace_bins_comp.png']; %eps
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-dpng','-r0')


ff = figure(3600);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
 for jj = 1:5
    N_bin = N_bins(jj);
    load(['Results/Simulation/SV_results_param_HMM_adapt_eff_Nbin',num2str(N_bin),'.mat'])
    for ii = 1:5
        subplot(5,5,(ii-1)*5+jj)
        hold on         
        autocorr(H_subset_HMM_adapt_eff(:,2*ii-1),200)      
        hold off
%         set(gca,'ylabel ',[])
%         set(gca,'xlabel ',[])  
        ylabel('')
        xlabel('')
        title(['Time = ',num2str(100*(2*ii-1)),', Nbin = ',num2str(N_bin)])
    %     plotTickLatex2D;%('FontSize',12);
    end
end
suptitle('Acf plots h')
name = ['figures/SV_HMM_adapt_h_acf_bins_comp.png']; %eps
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-dpng','-r0')

ff = figure(4600);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
for jj = 1:5
    N_bin = N_bins(jj);
    load(['Results/Simulation/SV_results_param_HMM_adapt_eff_Nbin',num2str(N_bin),'.mat'])
    for ii = 1:3
        subplot(3,5,(ii-1)*5+jj)
        hold on         
        autocorr(theta_HMM_adapt_eff(:,ii),200)      
         hold off
    %     set(gca,'ylabel',[])    
        title([params{ii},' ',', Nbin = ',num2str(N_bin)])
    %     plotTickLatex2D;%('FontSize',12);
    end
end
suptitle('Acf plots params')
name = ['path_figures/SV_HMM_adapt_param_acf_bins_comp.png']; %eps
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-dpng','-r0')



%% ACF METHOD COMP
% EXT = {'DA','DA_RW_eff','HMM_eff_Nbin30','HMM_adapt_eff_Nbin10'};
EXT = {'DA','DA_RW_eff','HMM_eff_Nbin30','HMM_adapt_eff_Nbin30'};
params = {'\mu','\phi','\sigma^2'};
DATA = {'', '_GSPC' , '_MSFT', '_IBM'};

for data = DATA
    if ~isempty(char(data))
        path_results = 'Results/Empirical';
        path_figures = 'figures/Empirical';
    else 
        path_results = 'Results/Simulation';
        path_figures = 'figures/Simulation';
    end

    ff = figure(7200);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
    for jj = 1:4
        ext = EXT{jj};
        clear -regexp ^theta
        load([path_results,'/SV_results_param_',ext,char(data),'.mat'],'-regexp','^theta_')
        clear 'theta_true'
        aaa = who('-regexp','^theta');
        theta = eval(char(aaa{1}));

        for ii = 1:3
            subplot(3,4,(ii-1)*4+jj)
            hold on         
            autocorr(theta(:,ii),200)      
             hold off
        %     set(gca,'ylabel',[])    
            title([params{ii},' ',ext])
        %     plotTickLatex2D;%('FontSize',12);
        end
    end
    suptitle(['Acf plots params, ',char(data)])
    name = [path_figures,'/SV_ALL_param_acf',char(data),'.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print(ff,name,'-dpng','-r0')


    ff = figure(7300);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.9 0.9]);
    for jj = 1:4
        ext = EXT{jj};
        clear -regexp ^H_
        load([path_results,'/SV_results_param_',ext,char(data),'.mat'],'-regexp','^H_subset','^time')
        aaa = who('-regexp','^H_');
        H_subset = eval(char(aaa{1}));

        for ii = 1:5
            subplot(5,4,(ii-1)*4+jj)
            hold on         
            autocorr(H_subset(:,2*2*ii-1),200);    
            hold off
    %         set(gca,'ylabel ',[])
    %         set(gca,'xlabel ',[])  
            ylabel('')
            xlabel('')
            title(['Time = ',num2str(100*(2*ii-1)),', ',ext])
        %     plotTickLatex2D;%('FontSize',12);
        end
    end
    suptitle(['Acf plots h','',char(data)])
    name = [path_figures,'/SV_ALL_h_acf',char(data),'.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print(ff,name,'-dpng','-r0')


    %% Trace Plots
    col = [0 0 1;
           0 1 0;
           1 0 1;
           1 0 0];
    T = 2000;
    % data = '_MSFT';

    ff = figure(2);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);    
    for t = 1:9
        subplot(3,3,t)
        hold on
        for jj = 1:4
            ext = EXT{jj};
            clear -regexp ^H_
            load([path_results,'/SV_results_param_',ext,char(data),'.mat'],'-regexp','^H_subset','^h_true','^time')
            aaa = who('-regexp','^H_');
            H_subset = eval(char(aaa{1}));   
            plot(H_subset(:,2*t),'Color',col(jj,:))
        end
        if isempty(char(data))
            plot(H_subset(:,t)*0 + h_true((t+1)*100),'k')
        end
        hold off
    %     XL = xlabel('MCMC iteration');
    % %     pos = get(XL,'pos'); % Read position [x y z]
    % %     YL = get(gca,'ylim');
    % %     pos(2) = 1.01*YL(1); 
    % %     set(XL,'pos',pos) % Move label up 
        title(['Time = ',num2str(100*(t+1))])   
    % %     plotTickLatex2D;%('FontSize',12);
    end
    leg = legend('DA','DA RW','HMM fix', 'HMM adapt');
    suptitle(['Trace plots h','',char(data)])
    set(leg,'Interpreter','latex');%,'
    name = [path_figures,'/SV_ALL_h_trace',char(data),'.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print(ff,name,'-dpng','-r0')


    params = {'\mu','\phi','\sigma^2'};
    ff = figure(16);
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.8]); 
    for t = 1:3
        subplot(3,1,t)
        hold on
        for jj = 1:4
            ext = EXT{jj};
            clear -regexp ^theta_
            load([path_results,'/SV_results_param_',ext,char(data),'.mat'],'-regexp','^theta_','^time')
            TH_true = theta_true;
            clear theta_true
            aaa = who('-regexp','^theta_');
            theta = eval(char(aaa{1}));   
            plot(theta(:,t),'Color',col(jj,:))
        end
        if isempty(char(data))
            plot(theta(:,t)*0 + TH_true(t),'k')
        end
        hold off
    %     XL = xlabel('MCMC iteration');
    % %     pos = get(XL,'pos'); % Read position [x y z]
    % %     YL = get(gca,'ylim');
    % %     pos(2) = 1.01*YL(1); 
    % %     set(XL,'pos',pos) % Move label up 
        title(params{t})   
    % %     plotTickLatex2D;%('FontSize',12);
    end
    leg = legend('DA','DA RW','HMM fix', 'HMM adapt');
    suptitle(['Trace plots param','',char(data)])
    set(leg,'Interpreter','latex');%,'
    name = [path_figures,'/SV_ALL_param_trace',char(data),'.png']; %eps
    set(gcf,'PaperPositionMode','auto');
    print(ff,name,'-dpng','-r0')
end



%% Na means
ColPalette= [0 0 0
    255 121 75
    255 51 0
    153 0 0 
    62 154 222
    0 102 204
    51 51 102
    0 255 0]/255;

ff = figure(3);
    tt = (1:T)';
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.6 0.5]);
    set(gcf,'defaulttextinterpreter','latex');

%     xlabel('$\mathbf{h}$ posterior means and 95\% CI','interpreter','latex','FontSize',11)     
    xlabel('$\mathbf{h}$ posterior means and $\pm$std','interpreter','latex','FontSize',11)     
    hold on
    
    xxx = mean_H_DA_RW_eff'; 
    xxx_L = xxx - sqrt(var_H_DA_RW_eff)';
    xxx_U = xxx + sqrt(var_H_DA_RW_eff)';
%     xxx_L = prctile(xxx',2.5)';
%     xxx_U = prctile(xxx',97.5)';   
    HPDI =  patch([tt; tt(end:-1:1); tt(1)],...
        [xxx_L; xxx_U(end:-1:1); xxx_L(1)],...
        'r');
    set(HPDI, 'facecolor', ColPalette(1,:), 'facealpha',0.2, 'edgecolor', 'none');
    plot(xxx,'color',ColPalette(1,:),'linewidth',1)

 
    hold off  
    
    set(gca,'TickLabelInterpreter','latex');    
    xlim([1,36])
    
%     %%
 
name = ['Figures/BurnIn_10000/BKM_N_mean_CI_all.eps'];
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-depsc','-r0')



bar([ESS_H_DA_RW_eff_sig;
    ESS_H_HMM_eff_sig;
    ESS_H_HMM_adapt_eff_sig])