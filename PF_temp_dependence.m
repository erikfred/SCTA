% PF_temp_dependence.m
%
% Sets up and solves inversion to get T and dT/dt dependence and drift from
% Pinon Flat calibration data. Constrains dependencies to be the same for 
% +X1, +X2, and -X, and same for +Y and -Y. Includes option for including 
% or excluding dT/dt dependence.
%
% FUTURE ADDITIONS:
% - Option for weighting data by interval between calibrations, to mitigate
% influence of calibration interval
% - Drift may change after instrument disturbance (e.g. power loss), so may
% build in option to make linear component independent between these
% events
%

clear; close all

dT_dep=true; %false will exclude dT/dt from inversion

load('../calibrations/PinonFlat/PFdata.mat','flipInfoAll','dataDec100')

% +X1, +X2, and -X
a_X1=flipInfoAll.gCal([1:3:58,61:5:end]);
t_X1=flipInfoAll.t([1:3:58,61:5:end])-flipInfoAll.t(1);
T_X1=flipInfoAll.T([1:3:58,61:5:end]);
a_X2a=flipInfoAll.gCal(3:3:60); % need to split up +X2 to before and after
t_X2a=flipInfoAll.t(3:3:60)-flipInfoAll.t(1);
T_X2a=flipInfoAll.T(3:3:60);
a_X2b=flipInfoAll.gCal(64:5:end);
t_X2b=flipInfoAll.t(64:5:end)-flipInfoAll.t(1);
T_X2b=flipInfoAll.T(64:5:end);
a_negX=-flipInfoAll.gCal(65:5:end);
t_negX=flipInfoAll.t(65:5:end)-flipInfoAll.t(1);
T_negX=flipInfoAll.T(65:5:end);

if dT_dep
    % interpolate to daily values
    td_X1=[t_X1(1):t_X1(end)+0.1]';
    ad_X1=interp1(t_X1,a_X1,td_X1,'pchip');
    Td_X1=interp1(t_X1,T_X1,td_X1,'pchip');
    dTd_X1=diff(Td_X1); dTd_X1=[dTd_X1(1);dTd_X1];
    td_X2a=[t_X2a(1):t_X2a(end)+0.1]';
    ad_X2a=interp1(t_X2a,a_X2a,td_X2a,'pchip');
    Td_X2a=interp1(t_X2a,T_X2a,td_X2a,'pchip');
    dTd_X2a=diff(Td_X2a); dTd_X2a=[dTd_X2a(1);dTd_X2a];
    td_X2b=[t_X2b(1):t_X2b(end)+0.1]';
    ad_X2b=interp1(t_X2b,a_X2b,td_X2b,'pchip');
    Td_X2b=interp1(t_X2b,T_X2b,td_X2b,'pchip');
    dTd_X2b=diff(Td_X2b); dTd_X2b=[dTd_X2b(1);dTd_X2b];
    td_negX=[t_negX(1):t_negX(end)+0.1]';
    ad_negX=interp1(t_negX,a_negX,td_negX,'pchip');
    Td_negX=interp1(t_negX,T_negX,td_negX,'pchip');
    dTd_negX=diff(Td_negX); dTd_negX=[dTd_negX(1);dTd_negX];
    
    % inversion of interpolated data
    Gd_X=[[td_X1;zeros(size(td_X2a));zeros(size(td_X2b));zeros(size(td_negX))],...
        [zeros(size(td_X1));td_X2a;zeros(size(td_X2b));zeros(size(td_negX))],...
        [zeros(size(td_X1));zeros(size(td_X2a));td_X2b;zeros(size(td_negX))],...
        [zeros(size(td_X1));zeros(size(td_X2a));zeros(size(td_X2b));td_negX],...
        [Td_X1;Td_X2a;Td_X2b;Td_negX],...
        [dTd_X1;dTd_X2a;dTd_X2b;dTd_negX],...
        [ones(size(Td_X1));zeros(size(td_X2a));zeros(size(td_X2b));zeros(size(Td_negX))],...
        [zeros(size(Td_X1));ones(size(td_X2a));zeros(size(td_X2b));zeros(size(Td_negX))],...
        [zeros(size(Td_X1));zeros(size(td_X2a));ones(size(td_X2b));zeros(size(Td_negX))],...
        [zeros(size(Td_X1));zeros(size(td_X2a));zeros(size(td_X2b));ones(size(Td_negX))]];
    
    m_X=(Gd_X'*Gd_X)\Gd_X'*[ad_X1;ad_X2a;ad_X2b;ad_negX];
    ad_X_star=Gd_X*m_X;
    
    % resample to original times
    [~,ix1,~]=intersect(round(td_X1),round(t_X1));
    a_X_star=ad_X_star(ix1);
    dT_X1=dTd_X1(ix1);
    [~,ix2a,~]=intersect(ceil(td_X2a),ceil(t_X2a));
    a_X_star=[a_X_star;ad_X_star(length(td_X1)+ix2a)];
    dT_X2a=dTd_X2a(ix2a);
    [~,ix2b,~]=intersect(ceil(td_X2b),ceil(t_X2b));
    a_X_star=[a_X_star;ad_X_star(length(td_X1)+length(td_X2a)+ix2b)];
    dT_X2b=dTd_X2b(ix2b);
    [~,inegx,~]=intersect(ceil(td_negX),ceil(t_negX));
    a_X_star=[a_X_star;ad_X_star(length(td_X1)+length(td_X2a)+length(td_X2b)+inegx)];
    dT_negX=dTd_negX(inegx);
else
    % inversion of original data only
    G_X=[[t_X1;zeros(size(t_X2a));zeros(size(t_X2b));zeros(size(t_negX))],...
        [zeros(size(t_X1));t_X2a;zeros(size(t_X2b));zeros(size(t_negX))],...
        [zeros(size(t_X1));zeros(size(t_X2a));t_X2b;zeros(size(t_negX))],...
        [zeros(size(t_X1));zeros(size(t_X2a));zeros(size(t_X2b));t_negX],...
        [T_X1;T_X2a;T_X2b;T_negX],...
        [ones(size(T_X1));zeros(size(t_X2a));zeros(size(t_X2b));zeros(size(T_negX))],...
        [zeros(size(T_X1));ones(size(t_X2a));zeros(size(t_X2b));zeros(size(T_negX))],...
        [zeros(size(T_X1));zeros(size(t_X2a));ones(size(t_X2b));zeros(size(T_negX))],...
        [zeros(size(T_X1));zeros(size(t_X2a));zeros(size(t_X2b));ones(size(T_negX))]];
    
    m_X=(G_X'*G_X)\G_X'*[a_X1;a_X2a;a_X2b;a_negX];
    a_X_star=G_X*m_X;
end

% % PLOTS TO VERIFY EVERYTHING LOOKS OK
% figure
% plot(td_X1,ad_X1,'o')
% hold on
% plot(td_X1,ad_X_star(1:length(td_X1)),'x')
% yyaxis right
% % plot(td_X1,Td_X1) %temperature
% plot(td_X1,ad_X1-ad_X_star(1:length(td_X1)),'x') %residual
% 
% figure
% plot(td_X2a,ad_X2a,'o')
% hold on
% plot(td_X2a,ad_X_star(length(td_X1)+1:length(td_X1)+length(td_X2a)),'x')
% yyaxis right
% % plot(td_X2a,Td_X2a)
% plot(td_X2a,ad_X2a-ad_X_star(length(td_X1)+1:length(td_X1)+length(td_X2a)),'x')
% 
% figure
% plot(td_X2b,ad_X2b,'o')
% hold on
% plot(td_X2b,ad_X_star(length(td_X1)+length(td_X2a)+1:length(td_X1)+length(td_X2a)+length(td_X2b)),'x')
% yyaxis right
% % plot(td_X2b,Td_X2b)
% plot(td_X2b,ad_X2b-ad_X_star(length(td_X1)+length(td_X2a)+1:length(td_X1)+length(td_X2a)+length(td_X2b)),'x')
% 
% figure
% plot(td_negX,ad_negX,'o')
% hold on
% plot(td_negX,ad_X_star(length(td_X1)+length(td_X2a)+length(td_X2b)+1:length(td_X1)+length(td_X2a)+length(td_X2b)+length(td_negX)),'x')
% yyaxis right
% % plot(td_negX,Td_negX)
% plot(td_negX,ad_negX-ad_X_star(length(td_X1)+length(td_X2a)+length(td_X2b)+1:length(td_X1)+length(td_X2a)+length(td_X2b)+length(td_negX)),'x')

% X PLOTTING
% +X1
figure
hold on
plot(t_X1+flipInfoAll.t(61),a_X1-mean(a_X1),'or','markersize',20)
plot(t_X1+flipInfoAll.t(61),a_X_star(1:length(t_X1))-mean(a_X1),'*r','markersize',20)
plot(t_X1+flipInfoAll.t(61),a_X1-a_X_star(1:length(t_X1)),'xk','markersize',20)
text(flipInfoAll.t(end-10),-2.25e-4,['drift = ' num2str(round(m_X(1)*365*10^5,1)) ' ug/y'],'fontsize',15)
text(flipInfoAll.t(end-10),-3.25e-4,['Tdep = ' num2str(round(m_X(5)*10^5,1)) ' ug/C'],'fontsize',15)
if dT_dep
    text(flipInfoAll.t(end-10),-4.25e-4,['dTdep = ' num2str(round(m_X(6)*10^5,1)) ' ug/(C/d)'],'fontsize',15)
end
set(gca,'fontsize',15)
datetick('x',6,'keeplimits')
lim_x=xlim;
xtickangle(45)
ylabel('Accel (m/s^2)')
title('+X1 Calibrations')
ylim([-7 6]*10^-4)
yyaxis right
plot(dataDec100.t,dataDec100.T,'linewidth',1)
ylabel(['Temperature (' char(176) 'C)'])
legend('calibration - g','model - g','residual','Temperature','location','southeast')
box on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
if dT_dep
    print('../calibrations/PinonFlat/T_dependence/X1_dT','-dtiff','-r300')
    saveas(gcf,'../calibrations/PinonFlat/T_dependence/X1_dT.fig')
else
    print('../calibrations/PinonFlat/T_dependence/X1','-dtiff','-r300')
    saveas(gcf,'../calibrations/PinonFlat/T_dependence/X1.fig')
end

% +X2
figure
hold on
plot([t_X2a;t_X2b]+flipInfoAll.t(61),[a_X2a;a_X2b]-mean([a_X2a;a_X2b]),'or','markersize',20)
plot([t_X2a;t_X2b]+flipInfoAll.t(61),a_X_star(length(t_X1)+1:length([t_X1;t_X2a;t_X2b]))-mean([a_X2a;a_X2b]),'*r','markersize',20)
plot([t_X2a;t_X2b]+flipInfoAll.t(61),[a_X2a;a_X2b]-a_X_star(length(t_X1)+1:length([t_X1;t_X2a;t_X2b])),'xk','markersize',20)
text(flipInfoAll.t(end-10),-2.25e-4,['drift = ' num2str(round(m_X(2)*365*10^5,1)) ', ' num2str(round(m_X(3)*365*10^5,1)) ' ug/y'],'fontsize',15)
text(flipInfoAll.t(end-10),-3.25e-4,['Tdep = ' num2str(round(m_X(5)*10^5,1)) ' ug/C'],'fontsize',15)
if dT_dep
    text(flipInfoAll.t(end-10),-4.25e-4,['dTdep = ' num2str(round(m_X(6)*10^5,1)) ' ug/(C/d)'],'fontsize',15)
end
set(gca,'fontsize',15)
datetick('x',6,'keeplimits')
lim_x=xlim;
xtickangle(45)
ylabel('Accel (m/s^2)')
title('+X2 Calibrations')
ylim([-7 6]*10^-4)
yyaxis right
plot(dataDec100.t,dataDec100.T,'linewidth',1)
ylabel(['Temperature (' char(176) 'C)'])
legend('calibration - g','model - g','residual','Temperature','location','southeast')
box on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
if dT_dep
    print('../calibrations/PinonFlat/T_dependence/X2_dT','-dtiff','-r300')
    saveas(gcf,'../calibrations/PinonFlat/T_dependence/X2_dT.fig')
else
    print('../calibrations/PinonFlat/T_dependence/X2','-dtiff','-r300')
    saveas(gcf,'../calibrations/PinonFlat/T_dependence/X2.fig')
end

% -X
figure
hold on
plot(t_negX+flipInfoAll.t(61),a_negX-mean(a_negX),'or','markersize',20)
plot(t_negX+flipInfoAll.t(61),a_X_star(length([t_X1;t_X2a;t_X2b])+1:end)-mean(a_negX),'*r','markersize',20)
plot(t_negX+flipInfoAll.t(61),a_negX-a_X_star(length([t_X1;t_X2a;t_X2b])+1:end),'xk','markersize',20)
text(flipInfoAll.t(end-10),-2.25e-4,['drift = ' num2str(round(m_X(4)*365*10^5,1)) ' ug/y'],'fontsize',15)
text(flipInfoAll.t(end-10),-3.25e-4,['Tdep = ' num2str(round(m_X(5)*10^5,1)) ' ug/C'],'fontsize',15)
if dT_dep
    text(flipInfoAll.t(end-10),-4.25e-4,['dTdep = ' num2str(round(m_X(6)*10^5,1)) ' ug/(C/d)'],'fontsize',15)
end
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
if dT_dep
    print('../calibrations/PinonFlat/T_dependence/negX_dT','-dtiff','-r300')
    saveas(gcf,'../calibrations/PinonFlat/T_dependence/negX_dT.fig')
    save('../calibrations/PinonFlat/T_dependence/X_dT','m_X','a_X_star')
else
    print('../calibrations/PinonFlat/T_dependence/negX','-dtiff','-r300')
    saveas(gcf,'../calibrations/PinonFlat/T_dependence/negX.fig')
    save('../calibrations/PinonFlat/T_dependence/X_T','m_X','a_X_star')
end

% +Y and -Y
a_Y=flipInfoAll.gCal([2:3:59,62:5:end]);
t_Y=flipInfoAll.t([2:3:59,62:5:end])-flipInfoAll.t(1);
T_Y=flipInfoAll.T([2:3:59,62:5:end]);
a_negY=-flipInfoAll.gCal(63:5:end);
t_negY=flipInfoAll.t(63:5:end)-flipInfoAll.t(1);
T_negY=flipInfoAll.T(63:5:end);

if dT_dep
    % interpolate to daily values
    td_Y=[t_Y(1):t_Y(end)+0.1]';
    ad_Y=interp1(t_Y,a_Y,td_Y,'pchip');
    Td_Y=interp1(t_Y,T_Y,td_Y,'pchip');
    dTd_Y=diff(Td_Y); dTd_Y=[dTd_Y(1);dTd_Y];
    td_negY=[t_negY(1):t_negY(end)+0.1]';
    ad_negY=interp1(t_negY,a_negY,td_negY,'pchip');
    Td_negY=interp1(t_negY,T_negY,td_negY,'pchip');
    dTd_negY=diff(Td_negY); dTd_negY=[dTd_negY(1);dTd_negY];
    
    % inversion of interpolated data
    Gd_Y=[[td_Y;zeros(size(td_negY))],...
        [zeros(size(td_Y));td_negY],...
        [Td_Y;Td_negY],...
        [dTd_Y;dTd_negY],...
        [ones(size(Td_Y));zeros(size(Td_negY))],...
        [zeros(size(Td_Y));ones(size(Td_negY))]];
    
    m_Y=(Gd_Y'*Gd_Y)\Gd_Y'*[ad_Y;ad_negY];
    ad_Y_star=Gd_Y*m_Y;
    
    % resample to original times
    [~,iy,~]=intersect(round(td_Y),round(t_Y));
    a_Y_star=ad_Y_star(iy);
    dT_Y=dTd_Y(iy);
    [~,inegy,~]=intersect(ceil(td_negY),ceil(t_negY));
    a_Y_star=[a_Y_star;ad_Y_star(length(td_Y)+inegy)];
    dT_negY=dTd_negY(inegy);
else
    % inversion of original data only
    G_Y=[[t_Y;zeros(size(t_negY))],...
        [zeros(size(t_Y));t_negY],...
        [T_Y;T_negY],...
        [ones(size(T_Y));zeros(size(T_negY))],...
        [zeros(size(T_Y));ones(size(T_negY))]];
    
    m_Y=(G_Y'*G_Y)\G_Y'*[a_Y;a_negY];
    a_Y_star=G_Y*m_Y;
end

% % PLOTS TO VERIFY EVERYTHING LOOKS OK
% figure
% plot(td_Y,ad_Y,'o')
% hold on
% plot(td_Y,ad_Y_star(1:length(td_Y)),'x')
% yyaxis right
% plot(td_Y,Td_Y) %temperature
% % plot(td_Y,ad_Y-ad_X_star(1:length(td_Y)),'x') %residual
% 
% figure
% plot(td_negY,ad_negY,'o')
% hold on
% plot(td_negY,ad_Y_star(length(td_Y)+1:length(td_Y)+length(td_negY)),'x')
% yyaxis right
% plot(td_negY,Td_negY)
% % plot(td_negY,ad_negY-ad_Y_star(length(td_Y)+1:length(td_Y)+length(td_negY)),'x')

% PLOTTING
figure
hold on
plot(t_Y+flipInfoAll.t(61),a_Y-mean(a_Y),'or','markersize',20)
plot(t_Y+flipInfoAll.t(61),a_Y_star(1:length(t_Y))-mean(a_Y),'*r','markersize',20)
plot(t_Y+flipInfoAll.t(61),a_Y-a_Y_star(1:length(t_Y)),'xk','markersize',20)
text(flipInfoAll.t(end-10),-2.25e-4,['drift = ' num2str(round(m_Y(1)*365*10^5,1)) ' ug/y'],'fontsize',15)
text(flipInfoAll.t(end-10),-3.25e-4,['Tdep = ' num2str(round(m_Y(3)*10^5,1)) ' ug/C'],'fontsize',15)
if dT_dep
    text(flipInfoAll.t(end-10),-4.25e-4,['dTdep = ' num2str(round(m_Y(6)*10^5,1)) ' ug/(C/d)'],'fontsize',15)
end
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
if dT_dep
    print('../calibrations/PinonFlat/T_dependence/Y_dT','-dtiff','-r300')
    saveas(gcf,'../calibrations/PinonFlat/T_dependence/Y_dT.fig')
else
    print('../calibrations/PinonFlat/T_dependence/Y','-dtiff','-r300')
    saveas(gcf,'../calibrations/PinonFlat/T_dependence/Y.fig')
end

figure
hold on
plot(t_negY+flipInfoAll.t(61),a_negY-mean(a_negY),'or','markersize',20)
plot(t_negY+flipInfoAll.t(61),a_Y_star(length(t_Y)+1:end)-mean(a_negY),'*r','markersize',20)
plot(t_negY+flipInfoAll.t(61),a_negY-a_Y_star(length(t_Y)+1:end),'xk','markersize',20)
text(flipInfoAll.t(end-10),-2.25e-4,['drift = ' num2str(round(m_Y(2)*365*10^5,1)) ' ug/y'],'fontsize',15)
text(flipInfoAll.t(end-10),-3.25e-4,['Tdep = ' num2str(round(m_Y(3)*10^5,1)) ' ug/C'],'fontsize',15)
if dT_dep
    text(flipInfoAll.t(end-10),-4.25e-4,['dTdep = ' num2str(round(m_Y(6)*10^5,1)) ' ug/(C/d)'],'fontsize',15)
end
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
if dT_dep
    print('../calibrations/PinonFlat/T_dependence/negY_dT','-dtiff','-r300')
    saveas(gcf,'../calibrations/PinonFlat/T_dependence/negY_dT.fig')
    save('../calibrations/PinonFlat/T_dependence/Y_dT','m_Y','a_Y_star')
else
    print('../calibrations/PinonFlat/T_dependence/negY','-dtiff','-r300')
    saveas(gcf,'../calibrations/PinonFlat/T_dependence/negY.fig')
    save('../calibrations/PinonFlat/T_dependence/Y_T','m_Y','a_Y_star')
end