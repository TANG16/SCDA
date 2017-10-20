T = 2000;
M = 10000;
params = {'\mu','\phi','\sigma^2'};


%% Volatility
ff = figure(1);
set(gcf,'units','normalized','outerposition',[0.25 0.25 0.45 0.45]);
% ax1 = axes('Position',[0.25 0.25 0.3 0.3],'Visible','on');
hold on
% plot(h_true,'r')
% plot(h,'b')
% plot(mean(H_DA((5000+1):M,:),1),'g')
% plot(mean(H_HMM_shift((5000+1):M,:),1),'m')
% plot(2:2:T,mean(H_HMM((5000+1):M,2:2:T),1),'b')
% plot(2:2:T,mean_H_HMM_adapt(:,2:2:T),'b')
plot(mean_H_DA,'g')
% plot(mean_H_DA_RW,'c')
plot(mean_H_DA_RW_eff,'c')
% plot(2:2:T,mean_H_HMM(:,2:2:T),'b')
plot(2:2:T,mean_H_HMM_eff(:,2:2:T),'r')
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
plot(mean_H_DA,'g')
plot(mean_H_HMM_shift,'m')
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
%     plot(H_subset_HMM_adapt(:,t)*0 + h_true((t+1)*100),'r')
    plot(H_subset_DA(:,t),'b')
%     plot(H_subset_DA_RW_eff(:,t),'c')
%     plot(H_subset_HMM_eff(:,t),'r')
    plot(H_subset_HMM_adapt(:,t),'m')
    hold off    
%     XL = xlabel('MCMC iteration');
% %     pos = get(XL,'pos'); % Read position [x y z]
% %     YL = get(gca,'ylim');
% %     pos(2) = 1.01*YL(1); 
% %     set(XL,'pos',pos) % Move label up 
%     title(['Time = ',num2str(100*(t+1))])   
% %     plotTickLatex2D;%('FontSize',12);
end
leg = legend('true','DA','DA RW','HMM');
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
    autocorr(H_subset_DA(:,t),200)
%     set(gca,'ylabel',[])    
    title(['Time = ',num2str(100*(t+1))])
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
suptitle('Even HMM no shift')
 
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
    autocorr(H_subset_HMM_shift(2:2:M,t),40)
    set(gca,'ylabel',[])    
    title(['Time = ',num2str(100*(t+1))])
%     plotTickLatex2D;%('FontSize',12);
end
suptitle('Even HMM shift')
 
if save_on
    name = ['figures/SV_HMM_shift_even_acf.png']; %eps
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
ff = figure(16);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.8]);
for ii = 1:3
    subplot(3,1,ii)
    hold on    
    plot(theta_DA(:,ii),'g')
%     plot(theta_DA_RW_eff(:,ii),'c')
%     plot(theta_HMM_eff(:,ii),'r')   
    plot(theta_HMM_adapt(:,ii),'m')      
%     plot(theta_true(ii) + 0*theta_DA_RW_eff(:,ii),'r')
    hold off
%     set(gca,'ylabel',[])    
%     title(params{ii})
%     plotTickLatex2D;%('FontSize',12);
 end
suptitle('Trace plots params')
leg = legend('DA','DA RW','HMM','HMM shift','true');
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
ff = figure(17);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.8]);
for ii = 1:3
    subplot(3,1,ii)
    autocorr(theta_DA(:,ii),40)
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
    autocorr(theta_DA_RW_eff(:,ii),40)
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
    autocorr(theta_HMM_shift(:,ii),40)
    set(gca,'ylabel',[])
    title(params{ii})
 end
suptitle('ACF params HMM shift')

if save_on
    name = ['figures/SV_HMM_shift_acf_param.png']; %eps
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
