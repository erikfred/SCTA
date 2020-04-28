% PF_temp_dependence.m
%
% Sets up and solves inversion to get T dependence and drift from Pinon
% Flat calibration data.
%

load('../calibrations/PinonFlat/PFdata.mat','flipInfoAll')

% +X1
a_X1=flipInfoAll.gCal(61:5:end);
t_X1=flipInfoAll.t(61:5:end)-flipInfoAll.t(61);
T_X1=flipInfoAll.T(61:5:end);

G_X1=[t_X1,T_X1,ones(size(T_X1))];

m_X1=(G_X1'*G_X1)\G_X1'*a_X1;
a_X1_star=G_X1*m_X1;

figure
hold on
plot(flipInfoAll.t(61:5:end),a_X1-mean(a_X1),'ob','markersize',10)
plot(flipInfoAll.t(61:5:end),a_X1_star-mean(a_X1),'*b','markersize',10)
plot(flipInfoAll.t(61:5:end),a_X1-a_X1_star,'xk','markersize',10)
text(flipInfoAll.t(end-5),-0.5e-4,['drift = ' num2str(round(m_X1(1)*365*10^5,1)) ' ug/y'],'fontsize',15)
text(flipInfoAll.t(end-5),-0.75e-4,['Tdep = ' num2str(round(m_X1(2)*10^5,1)) ' ug/C'],'fontsize',15)
text(flipInfoAll.t(end-5),-1e-4,['gmod = ' num2str(round(m_X1(3),3)) ' m/s^2'],'fontsize',15)
set(gca,'fontsize',15)
datetick('keeplimits')
xtickangle(45)
ylabel('Accel (m/s^2)')
title('First +X Calibrations')
ylim([-1.5 1.5]*10^-4)
legend('calibration - g','model - g','residual')
box on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../calibrations/PinonFlat/T_dependence/X1','-dtiff','-r300')
saveas(gcf,'../calibrations/PinonFlat/T_dependence/X1.fig')
save('../calibrations/PinonFlat/T_dependence/X1','a_X1','t_X1','T_X1','G_X1','m_X1','a_X1_star')

% +X2
a_X2=flipInfoAll.gCal(64:5:end);
t_X2=flipInfoAll.t(64:5:end)-flipInfoAll.t(64);
T_X2=flipInfoAll.T(64:5:end);

G_X2=[t_X2,T_X2,ones(size(T_X2))];

m_X2=(G_X2'*G_X2)\G_X2'*a_X2;
a_X2_star=G_X2*m_X2;

figure
hold on
plot(flipInfoAll.t(61:5:end),a_X2-mean(a_X2),'ob','markersize',10)
plot(flipInfoAll.t(61:5:end),a_X2_star-mean(a_X2),'*b','markersize',10)
plot(flipInfoAll.t(61:5:end),a_X2-a_X2_star,'xk','markersize',10)
text(flipInfoAll.t(end-5),-0.5e-4,['drift = ' num2str(round(m_X2(1)*365*10^5,1)) ' ug/y'],'fontsize',15)
text(flipInfoAll.t(end-5),-0.75e-4,['Tdep = ' num2str(round(m_X2(2)*10^5,1)) ' ug/C'],'fontsize',15)
text(flipInfoAll.t(end-5),-1e-4,['gmod = ' num2str(round(m_X2(3),3)) ' m/s^2'],'fontsize',15)
set(gca,'fontsize',15)
datetick('keeplimits')
xtickangle(45)
ylabel('Accel (m/s^2)')
title('Second +X Calibrations')
ylim([-1.5 1.5]*10^-4)
legend('calibration - g','model - g','residual')
box on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../calibrations/PinonFlat/T_dependence/X2','-dtiff','-r300')
saveas(gcf,'../calibrations/PinonFlat/T_dependence/X2.fig')
save('../calibrations/PinonFlat/T_dependence/X2','a_X2','t_X2','T_X2','G_X2','m_X2','a_X2_star')

% -X
a_negX=flipInfoAll.gCal(65:5:end);
t_negX=flipInfoAll.t(65:5:end)-flipInfoAll.t(65);
T_negX=flipInfoAll.T(65:5:end);

G_negX=[t_negX,T_negX,ones(size(T_negX))];

m_negX=(G_negX'*G_negX)\G_negX'*a_negX;
a_negX_star=G_negX*m_negX;

figure
hold on
plot(flipInfoAll.t(61:5:end),a_negX-mean(a_negX),'ob','markersize',10)
plot(flipInfoAll.t(61:5:end),a_negX_star-mean(a_negX),'*b','markersize',10)
plot(flipInfoAll.t(61:5:end),a_negX-a_negX_star,'xk','markersize',10)
text(flipInfoAll.t(61),-0.5e-4,['drift = ' num2str(round(m_negX(1)*365*10^5,1)) ' ug/y'],'fontsize',15)
text(flipInfoAll.t(61),-0.75e-4,['Tdep = ' num2str(round(m_negX(2)*10^5,1)) ' ug/C'],'fontsize',15)
text(flipInfoAll.t(61),-1e-4,['gmod = ' num2str(round(m_negX(3),3)) ' m/s^2'],'fontsize',15)
set(gca,'fontsize',15)
datetick('keeplimits')
xtickangle(45)
ylabel('Accel (m/s^2)')
title('-X Calibrations')
ylim([-1.5 1.5]*10^-4)
legend('calibration - g','model - g','residual')
box on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../calibrations/PinonFlat/T_dependence/negX','-dtiff','-r300')
saveas(gcf,'../calibrations/PinonFlat/T_dependence/negX.fig')
save('../calibrations/PinonFlat/T_dependence/negX','a_negX','t_negX','T_negX','G_negX','m_negX','a_negX_star')

% +Y
a_Y=flipInfoAll.gCal(62:5:end);
t_Y=flipInfoAll.t(62:5:end)-flipInfoAll.t(62);
T_Y=flipInfoAll.T(62:5:end);

G_Y=[t_Y,T_Y,ones(size(T_Y))];

m_Y=(G_Y'*G_Y)\G_Y'*a_Y;
a_Y_star=G_Y*m_Y;

figure
hold on
plot(flipInfoAll.t(61:5:end),a_Y-mean(a_Y),'ob','markersize',10)
plot(flipInfoAll.t(61:5:end),a_Y_star-mean(a_Y),'*b','markersize',10)
plot(flipInfoAll.t(61:5:end),a_Y-a_Y_star,'xk','markersize',10)
text(flipInfoAll.t(end-5),-0.5e-4,['drift = ' num2str(round(m_Y(1)*365*10^5,1)) ' ug/y'],'fontsize',15)
text(flipInfoAll.t(end-5),-0.75e-4,['Tdep = ' num2str(round(m_Y(2)*10^5,1)) ' ug/C'],'fontsize',15)
text(flipInfoAll.t(end-5),-1e-4,['gmod = ' num2str(round(m_Y(3),3)) ' m/s^2'],'fontsize',15)
set(gca,'fontsize',15)
datetick('keeplimits')
xtickangle(45)
ylabel('Accel (m/s^2)')
title('+Y Calibrations')
ylim([-1.5 1.5]*10^-4)
legend('calibration - g','model - g','residual')
box on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../calibrations/PinonFlat/T_dependence/Y','-dtiff','-r300')
saveas(gcf,'../calibrations/PinonFlat/T_dependence/Y.fig')
save('../calibrations/PinonFlat/T_dependence/Y','a_Y','t_Y','T_Y','G_Y','m_Y','a_Y_star')

% -Y
a_negY=flipInfoAll.gCal(63:5:end);
t_negY=flipInfoAll.t(63:5:end)-flipInfoAll.t(63);
T_negY=flipInfoAll.T(63:5:end);

G_negY=[t_negY,T_negY,ones(size(T_negY))];

m_negY=(G_negY'*G_negY)\G_negY'*a_negY;
a_negY_star=G_negY*m_negY;

figure
hold on
plot(flipInfoAll.t(61:5:end),a_negY-mean(a_negY),'ob','markersize',10)
plot(flipInfoAll.t(61:5:end),a_negY_star-mean(a_negY),'*b','markersize',10)
plot(flipInfoAll.t(61:5:end),a_negY-a_negY_star,'xk','markersize',10)
text(flipInfoAll.t(61),-0.5e-4,['drift = ' num2str(round(m_negY(1)*365*10^5,1)) ' ug/y'],'fontsize',15)
text(flipInfoAll.t(61),-0.75e-4,['Tdep = ' num2str(round(m_negY(2)*10^5,1)) ' ug/C'],'fontsize',15)
text(flipInfoAll.t(61),-1e-4,['gmod = ' num2str(round(m_negY(3),3)) ' m/s^2'],'fontsize',15)
set(gca,'fontsize',15)
datetick('keeplimits')
xtickangle(45)
ylabel('Accel (m/s^2)')
title('-Y Calibrations')
ylim([-1.5 1.5]*10^-4)
legend('calibration - g','model - g','residual')
box on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../calibrations/PinonFlat/T_dependence/negY','-dtiff','-r300')
saveas(gcf,'../calibrations/PinonFlat/T_dependence/negY.fig')
save('../calibrations/PinonFlat/T_dependence/negY','a_negY','t_negY','T_negY','G_negY','m_negY','a_negY_star')