% PF_temp_dependence_independentT.m
%
% Sets up and solves inversion to get T dependence and drift from Pinon
% Flat calibration data.
%

load('../calibrations/PinonFlat/PFdata.mat','flipInfoAll','dataDec100')

% +X1 and -X
a_X1=flipInfoAll.gCal([34:3:58,61:5:end]);
t_X1=flipInfoAll.t([34:3:58,61:5:end])-flipInfoAll.t(61);
T_X1=flipInfoAll.T([34:3:58,61:5:end]);
a_negX=-flipInfoAll.gCal(65:5:end);
t_negX=flipInfoAll.t(65:5:end)-flipInfoAll.t(61);
T_negX=flipInfoAll.T(65:5:end);

G_X=[[t_X1;zeros(size(t_negX))],[zeros(size(t_X1));t_negX],[T_X1;T_negX],...
    [ones(size(T_X1));zeros(size(T_negX))],[zeros(size(T_X1));ones(size(T_negX))]];

m_X=(G_X'*G_X)\G_X'*[a_X1;a_negX];
a_X_star=G_X*m_X;

figure
hold on
plot(t_X1+flipInfoAll.t(61),a_X1-mean(a_X1),'or','markersize',20)
plot(t_X1+flipInfoAll.t(61),a_X_star(1:length(t_X1))-mean(a_X1),'*r','markersize',20)
plot(t_X1+flipInfoAll.t(61),a_X1-a_X_star(1:length(t_X1)),'xk','markersize',20)
text(flipInfoAll.t(end-10),-2.25e-4,['drift = ' num2str(round(m_X(1)*365*10^5,1)) ' ug/y'],'fontsize',15)
text(flipInfoAll.t(end-10),-3.25e-4,['Tdep = ' num2str(round(m_X(3)*10^5,1)) ' ug/C'],'fontsize',15)
set(gca,'fontsize',15)
datetick('x',6,'keeplimits')
lim_x=xlim;
xtickangle(45)
ylabel('Accel (m/s^2)')
title('+X Calibrations')
ylim([-7 6]*10^-4)
yyaxis right
plot(dataDec100.t,dataDec100.T,'linewidth',1)
ylabel(['Temperature (' char(176) 'C)'])
legend('calibration - g','model - g','residual','Temperature','location','southeast')
box on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../calibrations/PinonFlat/T_dependence/X1','-dtiff','-r300')
saveas(gcf,'../calibrations/PinonFlat/T_dependence/X1.fig')

figure
hold on
plot(t_negX+flipInfoAll.t(61),a_negX-mean(a_negX),'or','markersize',20)
plot(t_negX+flipInfoAll.t(61),a_X_star(length(t_X1)+1:end)-mean(a_negX),'*r','markersize',20)
plot(t_negX+flipInfoAll.t(61),a_negX-a_X_star(length(t_X1)+1:end),'xk','markersize',20)
text(flipInfoAll.t(end-10),-2.25e-4,['drift = ' num2str(round(m_X(2)*365*10^5,1)) ' ug/y'],'fontsize',15)
text(flipInfoAll.t(end-10),-3.25e-4,['Tdep = ' num2str(round(m_X(3)*10^5,1)) ' ug/C'],'fontsize',15)
set(gca,'fontsize',15)
xlim([lim_x(1) lim_x(2)])
datetick('x',6,'keeplimits')
xtickangle(45)
ylabel('Accel (m/s^2)')
title('-X Calibrations')
ylim([-7 6]*10^-4)
yyaxis right
plot(dataDec100.t,dataDec100.T,'linewidth',1)
ylabel(['Temperature (' char(176) 'C)'])
legend('calibration + g','model + g','residual','Temperature','location','southeast')
box on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../calibrations/PinonFlat/T_dependence/negX','-dtiff','-r300')
saveas(gcf,'../calibrations/PinonFlat/T_dependence/negX.fig')
save('../calibrations/PinonFlat/T_dependence/X1','a_X1','t_X1','T_X1','G_X','m_X','a_X_star','a_negX','t_negX','T_negX')

% +Y and -Y
a_Y=flipInfoAll.gCal([35:3:59,62:5:end]);
t_Y=flipInfoAll.t([35:3:59,62:5:end])-flipInfoAll.t(61);
T_Y=flipInfoAll.T([35:3:59,62:5:end]);
a_negY=-flipInfoAll.gCal(63:5:end);
t_negY=flipInfoAll.t(63:5:end)-flipInfoAll.t(61);
T_negY=flipInfoAll.T(63:5:end);

G_Y=[[t_Y;zeros(size(t_negY))],[zeros(size(t_Y));t_negY],[T_Y;T_negY],...
    [ones(size(T_Y));zeros(size(T_negY))],[zeros(size(T_Y));ones(size(T_negY))]];

m_Y=(G_Y'*G_Y)\G_Y'*[a_Y;a_negY];
a_Y_star=G_Y*m_Y;

figure
hold on
plot(t_Y+flipInfoAll.t(61),a_Y-mean(a_Y),'or','markersize',20)
plot(t_Y+flipInfoAll.t(61),a_Y_star(1:length(t_Y))-mean(a_Y),'*r','markersize',20)
plot(t_Y+flipInfoAll.t(61),a_Y-a_Y_star(1:length(t_Y)),'xk','markersize',20)
text(flipInfoAll.t(end-10),-2.25e-4,['drift = ' num2str(round(m_Y(1)*365*10^5,1)) ' ug/y'],'fontsize',15)
text(flipInfoAll.t(end-10),-3.25e-4,['Tdep = ' num2str(round(m_Y(3)*10^5,1)) ' ug/C'],'fontsize',15)
set(gca,'fontsize',15)
datetick('x',6,'keeplimits')
lim_x=xlim;
xtickangle(45)
ylabel('Accel (m/s^2)')
title('+Y Calibrations')
ylim([-7 6]*10^-4)
yyaxis right
plot(dataDec100.t,dataDec100.T,'linewidth',1)
ylabel(['Temperature (' char(176) 'C)'])
legend('calibration - g','model - g','residual','Temperature','location','southeast')
box on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../calibrations/PinonFlat/T_dependence/Y','-dtiff','-r300')
saveas(gcf,'../calibrations/PinonFlat/T_dependence/Y.fig')

figure
hold on
plot(t_negY+flipInfoAll.t(61),a_negY-mean(a_negY),'or','markersize',20)
plot(t_negY+flipInfoAll.t(61),a_Y_star(length(t_Y)+1:end)-mean(a_negY),'*r','markersize',20)
plot(t_negY+flipInfoAll.t(61),a_negY-a_Y_star(length(t_Y)+1:end),'xk','markersize',20)
text(flipInfoAll.t(end-10),-2.25e-4,['drift = ' num2str(round(m_Y(2)*365*10^5,1)) ' ug/y'],'fontsize',15)
text(flipInfoAll.t(end-10),-3.25e-4,['Tdep = ' num2str(round(m_Y(3)*10^5,1)) ' ug/C'],'fontsize',15)
set(gca,'fontsize',15)
xlim([lim_x(1) lim_x(2)])
datetick('x',6,'keeplimits')
xtickangle(45)
ylabel('Accel (m/s^2)')
title('-Y Calibrations')
ylim([-7 6]*10^-4)
yyaxis right
plot(dataDec100.t,dataDec100.T,'linewidth',1)
ylabel(['Temperature (' char(176) 'C)'])
legend('calibration + g','model + g','residual','Temperature','location','southeast')
box on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../calibrations/PinonFlat/T_dependence/negY','-dtiff','-r300')
saveas(gcf,'../calibrations/PinonFlat/T_dependence/negY.fig')
save('../calibrations/PinonFlat/T_dependence/Y','a_Y','t_Y','T_Y','G_Y','m_Y','a_Y_star','a_negY','t_negY','T_negY')