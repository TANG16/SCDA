tt = 1965:1998;

ff = figure(101);
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.4 0.3]);
set(gcf,'defaulttextinterpreter','latex');


subplot(1,2,1)
plot(tt, y_ss)
set(gca,'TickLabelInterpreter','latex');
xlim([1965,1998])
xlabel('Year','FontSize',11)
title('Lapwings: observed census data','interpreter','latex','FontSize',11)

subplot(1,2,2)
plot(tt, f_ss)
hold on
plot(tt, mean(f_ss) + 0*f_ss,'k')
hold off
set(gca,'TickLabelInterpreter','latex');
xlim([1965,1998])
xlabel('Year','FontSize',11)
title('Frost days (normalised)','interpreter','latex','FontSize',11)
 

name = ['Figures/Lapwings_data.eps'];
set(gcf,'PaperPositionMode','auto');
print(ff,name,'-depsc','-r0')