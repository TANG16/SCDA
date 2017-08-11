%% Volatility
ff = figure(1);
set(gcf,'units','normalized','outerposition',[0.25 0.25 0.45 0.45]);
% ax1 = axes('Position',[0.25 0.25 0.3 0.3],'Visible','on');
hold on
plot(h_true,'r')
% plot(h,'b')
% plot(mean(H_DA((5000+1):M,:),1),'g')
% plot(mean(H_HMM_shift((5000+1):M,:),1),'m')
% plot(2:2:T,mean(H_HMM((5000+1):M,2:2:T),1),'b')
% plot(2:2:T,mean_H_HMM(:,2:2:T),'b')
plot(mean_H_HMM_shift,'m')

hold off
plotTickLatex2D;%('FontSize',12);
leg = legend('true','DA','HMM shift','HMM');
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
for t = 2:10
    subplot(3,3,t-1)
    hold on
    plot(H_HMM(1:M,t*100)*0 + h_true(t),'r')
%     plot(H_DA(1:M,t*100),'g')
    plot(H_HMM(1:M,t*100),'b')
%     plot(2:2:M,H_HMM_shift(2:2:M,t*100),'m')
    hold off    
    XL = xlabel('MCMC iteration');
%     pos = get(XL,'pos'); % Read position [x y z]
%     YL = get(gca,'ylim');
%     pos(2) = 1.01*YL(1); 
%     set(XL,'pos',pos) % Move label up 
    title(['Time = ',num2str(100*t)])   
%     plotTickLatex2D;%('FontSize',12);
end
leg = legend('true','DA','HMM','HMM shift');
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
ff = figure(3);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
for t = 2:10
    subplot(3,3,t-1)
    autocorr(H_DA(1:M,t*100),40)
    set(gca,'ylabel',[])    
    title(['Time = ',num2str(100*t)])
%     plotTickLatex2D;%('FontSize',12);
end
suptitle('Even DA no shift')
 
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



ff = figure(4);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
for t = 2:10
    subplot(3,3,t-1)
    autocorr(H_HMM(1:M,t*100),40)
    set(gca,'ylabel',[])    
    title(['Time = ',num2str(100*t)])
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
for t = 2:10
    subplot(3,3,t-1)
    autocorr(H_HMM_shift(2:2:M,t*100),40)
    set(gca,'ylabel',[])    
    title(['Time = ',num2str(100*t)])
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
for t = 2:10
    subplot(3,3,t-1)
    autocorr(H_HMM_shift(1:2:M,t*100-1),40)
    set(gca,'ylabel',[])    
    title(['Time = ',num2str(100*t-1)])
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

%% Trace Plots params 
ff = figure(16);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.8]);
for ii = 1:3
    subplot(3,1,ii)
    hold on
    plot(theta_HMM_shift(:,ii))
    plot(theta_true(ii) + 0*theta_HMM_shift(:,ii),'r')
    hold off
    set(gca,'ylabel',[])    
    title(params{ii})
 end
suptitle('Trace plots params')

%% ACF params 
ff = figure(17);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.8]);
for ii = 1:3
    subplot(3,1,ii)
    autocorr(theta_HMM_shift(:,ii),40)
    set(gca,'ylabel',[])
    set(gca,'ylabel',[])    
    title(params{ii})
 end
suptitle('ACF params')
