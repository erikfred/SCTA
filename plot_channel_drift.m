% plot_channel_drift.m
%
% Automated process for removing significant outliers (3*MAD method) and
% then quantifying linear drift rates and RMSE of the calibrations to those
% rates. Generates and saves plots throughout.
%

clear; close all

% AXIAL, INITIAL LOCATION

% AXIAL, NEW LOCATION

load('../calibrations/Axial/axialdata_newloc','flipInfoAll')
% identify and separate each of the calibrations
i_x=find(flipInfoAll.orientation==1);
i_xb=find(flipInfoAll.orientation==1 & flipInfoAll.t>datenum(2019,08,13));
i_y=find(flipInfoAll.orientation==2);
i_yb=find(flipInfoAll.orientation==2 & flipInfoAll.t>datenum(2019,08,13));
i_negx=find(flipInfoAll.orientation==-1);
i_negy=find(flipInfoAll.orientation==-2);

figure(55); clf
%----- +X1
% start w/ all cals
subplot(311); hold on
t_x1=flipInfoAll.t(i_x(1:2:end));
a_x1=flipInfoAll.gCal(i_x(1:2:end));
% iterative linear fits to remove outliers
ix1=true;
cor_x1=a_x1;
ni=length(cor_x1);
while any(ix1)
    P=polyfit(t_x1-t_x1(1),cor_x1,1);
    lin_x1=polyval(P,t_x1-t_x1(1));
    temp_x1=cor_x1-lin_x1;
    [~,ix1]=rmoutliers(temp_x1);
    cor_x1=cor_x1(~ix1);
    t_x1=t_x1(~ix1);
end
nf=length(cor_x1);
plot(t_x1,cor_x1,'ko')
plot(t_x1,lin_x1,'k-')
drift_x1=round(P(1)*365*10^5,2);
text(t_x1(10),cor_x1(35),['drift = ' num2str(drift_x1) ' \mug/yr'])
text(t_x1(10),cor_x1(35)-0.00001,['\sigma = ' num2str(round(rms(cor_x1-lin_x1)*10^5,2)) ' \mug'])
datetick('x',3)
xtickangle(45)
ylabel('+X1 (m/s^2)')
set(gca,'fontsize',12)
ylim([9.81185 9.81205])
xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/20,mean(yl)-[0.00005 -0.00005],...
    '-k','linewidth',2); text(xl(1)+diff(xl)/15,mean(yl),'10 \mug','fontsize',12)
%----- +X2
% start w/ all cals
yyaxis right; hold on
set(gca,'YColor','r')
t_x2=flipInfoAll.t(i_x(2:2:end));
a_x2=flipInfoAll.gCal(i_x(2:2:end));
% iterative linear fits to remove outliers
ix2=true;
cor_x2=a_x2;
ni=length(cor_x2);
while any(ix2)
    P=polyfit(t_x2-t_x2(1),cor_x2,1);
    lin_x2=polyval(P,t_x2-t_x2(1));
    hold off
    temp_x2=cor_x2-lin_x2;
    [~,ix2]=rmoutliers(temp_x2);
    cor_x2=cor_x2(~ix2);
    t_x2=t_x2(~ix2);
end
nf=length(cor_x2);
hold on
plot(t_x2,cor_x2,'ro')
plot(t_x2,lin_x2,'r-')
drift_x2=round(P(1)*365*10^5,2);
text(t_x2(35),cor_x2(15),['drift = ' num2str(drift_x2) ' \mug/yr'],'color','r')
text(t_x2(35),cor_x2(15)-0.00001,['\sigma = ' num2str(round(rms(cor_x2-lin_x2)*10^5,2)) ' \mug'],'color','r')
datetick('x',3)
xtickangle(45)
ylabel('+X2 (m/s^2)')
set(gca,'fontsize',12)
title({'Axial SCTA X Calibrations',[datestr(t_x1(1),'mmm dd, yyyy') ' - ' datestr(t_x1(end),'mmm dd, yyyy')]})
ylim([9.81215 9.81235])
box on; grid on
%----- -X
% start w/ all cals
subplot(312); hold on
t_xn=flipInfoAll.t(i_negx(1:end));
a_xn=-1*flipInfoAll.gCal(i_negx(1:end));
% iterative linear fits to remove outliers
ixn=true;
cor_xn=a_xn;
ni=length(cor_xn);
while any(ixn)
    P=polyfit(t_xn-t_xn(1),cor_xn,1);
    lin_xn=polyval(P,t_xn-t_xn(1));
    hold off
    temp_xn=cor_xn-lin_xn;
    [~,ixn]=rmoutliers(temp_xn);
    cor_xn=cor_xn(~ixn);
    t_xn=t_xn(~ixn);
end
nf=length(cor_xn);
hold on
plot(t_xn,cor_xn,'ks')
plot(t_xn,lin_xn,'k-')
drift_xn=round(P(1)*365*10^5,2);
text(t_xn(10),cor_xn(35),['drift = ' num2str(drift_xn) ' \mug/yr'])
text(t_xn(10),cor_xn(35)-0.00001,['\sigma = ' num2str(round(rms(cor_xn-lin_xn)*10^5,2)) ' \mug'])
datetick('x',3)
xtickangle(45)
ylabel('-X (m/s^2)')
set(gca,'fontsize',12)
ylim([-9.809 -9.8088])
xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/20,mean(yl)-[0.00005 -0.00005],...
    '-k','linewidth',2); text(xl(1)+diff(xl)/15,mean(yl),'10 \mug','fontsize',12)
box on; grid on
%----- X1 span
% start w/ all cals
subplot(313); hold on
t_xs1=flipInfoAll.t(i_negx(1:end));
a_xs1=flipInfoAll.gCal(i_xb(1:2:end))+flipInfoAll.gCal(i_negx(1:end));
% iterative linear fits to remove outliers
ixs1=true;
cor_xs1=a_xs1;
ni=length(cor_xs1);
while any(ixs1)
    P=polyfit(t_xs1-t_xs1(1),cor_xs1,1);
    lin_xs1=polyval(P,t_xs1-t_xs1(1));
    hold off
    temp_xs1=cor_xs1-lin_xs1;
    [~,ixs1]=rmoutliers(temp_xs1);
    cor_xs1=cor_xs1(~ixs1);
    t_xs1=t_xs1(~ixs1);
end
nf=length(cor_xs1);
hold on
plot(t_xs1,cor_xs1-cor_xs1(1),'k^')
plot(t_xs1,lin_xs1-cor_xs1(1),'k-')
drift_xs1=round(P(1)*365*10^5,2);
text(t_xs1(10),0.8e-5,['drift = ' num2str(drift_xs1) ' \mug/yr'])
text(t_xs1(10),0.7e-5,['\sigma = ' num2str(round(rms(cor_xs1-lin_xs1)*10^5,2)) ' \mug'])
datetick('x',3)
xtickangle(45)
ylabel('\Delta X1 span (m/s^2)')
set(gca,'fontsize',12)
ylim([-0.00001 0.00001])
xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/20,mean(yl)-[0.000005 -0.000005],...
    '-k','linewidth',2); text(xl(1)+diff(xl)/15,mean(yl),'1 \mug','fontsize',12)
box on; grid on
%----- X2 span
% start w/ all cals
yyaxis right; hold on
set(gca,'YColor','r')
t_xs2=flipInfoAll.t(i_negx(1:end));
a_xs2=flipInfoAll.gCal(i_xb(2:2:end))+flipInfoAll.gCal(i_negx);
% iterative linear fits to remove outliers
ixs2=true;
cor_xs2=a_xs2;
ni=length(cor_xs2);
while any(ixs2)
    P=polyfit(t_xs2-t_xs2(1),cor_xs2,1);
    lin_xs2=polyval(P,t_xs2-t_xs2(1));
    hold off
    temp_xs2=cor_xs2-lin_xs2;
    [~,ixs2]=rmoutliers(temp_xs2);
    cor_xs2=cor_xs2(~ixs2);
    t_xs2=t_xs2(~ixs2);
end
nf=length(cor_xs2);
hold on
plot(t_xs2,cor_xs2-cor_xs2(1),'r^')
plot(t_xs2,lin_xs2-cor_xs2(1),'r-')
drift_xs2=round(P(1)*365*10^5,2);
text(t_xs2(25),0.8e-5,['drift = ' num2str(drift_xs2) ' \mug/yr'],'color','r')
text(t_xs2(25),0.7e-5,['\sigma = ' num2str(round(rms(cor_xs2-lin_xs2)*10^5,2)) ' \mug'],'color','r')
datetick('x',3)
xtickangle(45)
ylabel('\Delta X2 span (m/s^2)')
set(gca,'fontsize',12)
ylim([-0.00001 0.00001])
xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/20,mean(yl)-[0.000005 -0.000005],...
    '-k','linewidth',2); text(xl(1)+diff(xl)/15,mean(yl),'1 \mug','fontsize',12)
box on; grid on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print('../calibrations/Axial/drift_rates/x_channels','-dtiff','-r300')

figure(56); clf
%----- +Y
% start w/ all cals
subplot(311); hold on
t_y1=flipInfoAll.t(i_y(1:end));
a_y1=flipInfoAll.gCal(i_y(1:end));
% iterative linear fits to remove outliers
iy1=true;
cor_y1=a_y1;
ni=length(cor_y1);
while any(iy1)
    P=polyfit(t_y1-t_y1(1),cor_y1,1);
    lin_y1=polyval(P,t_y1-t_y1(1));
    temp_y1=cor_y1-lin_y1;
    [~,iy1]=rmoutliers(temp_y1);
    cor_y1=cor_y1(~iy1);
    t_y1=t_y1(~iy1);
end
nf=length(cor_y1);
plot(t_y1,cor_y1,'ko')
plot(t_y1,lin_y1,'k-')
drift_y1=round(P(1)*365*10^5,2);
text(t_y1(10),cor_y1(35),['drift = ' num2str(drift_y1) ' \mug/yr'])
text(t_y1(10),cor_y1(35)-0.00002,['\sigma = ' num2str(round(rms(cor_y1-lin_y1)*10^5,2)) ' \mug'])
datetick('x',3)
xtickangle(45)
ylabel('+Y (m/s^2)')
set(gca,'fontsize',12)
ylim([9.8111 9.8115])
xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/20,mean(yl)-[0.00005 -0.00005],...
    '-k','linewidth',2); text(xl(1)+diff(xl)/15,mean(yl),'10 \mug','fontsize',12)
box on; grid on
%----- -Y
% start w/ all cals
subplot(312); hold on
t_yn=flipInfoAll.t(i_negy(1:end));
a_yn=-1*flipInfoAll.gCal(i_negy(1:end));
% iterative linear fits to remove outliers
iyn=true;
cor_yn=a_yn;
ni=length(cor_yn);
while any(iyn)
    P=polyfit(t_yn-t_yn(1),cor_yn,1);
    lin_yn=polyval(P,t_yn-t_yn(1));
    hold off
    temp_yn=cor_yn-lin_yn;
    [~,iyn]=rmoutliers(temp_yn);
    cor_yn=cor_yn(~iyn);
    t_yn=t_yn(~iyn);
end
nf=length(cor_yn);
hold on
plot(t_yn,cor_yn,'ks')
plot(t_yn,lin_yn,'k-')
drift_yn=round(P(1)*365*10^5,2);
text(t_yn(10),cor_yn(35),['drift = ' num2str(drift_yn) ' \mug/yr'])
text(t_yn(10),cor_yn(35)-0.00002,['\sigma = ' num2str(round(rms(cor_yn-lin_yn)*10^5,2)) ' \mug'])
datetick('x',3)
xtickangle(45)
ylabel('-Y (m/s^2)')
set(gca,'fontsize',12)
ylim([-9.8111 -9.8107])
xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/20,mean(yl)-[0.00005 -0.00005],...
    '-k','linewidth',2); text(xl(1)+diff(xl)/15,mean(yl),'10 \mug','fontsize',12)
box on; grid on
%----- Y span
% start w/ all cals
subplot(313); hold on
t_ys=flipInfoAll.t(i_negy(1:end));
a_ys=flipInfoAll.gCal(i_y(1:end))+flipInfoAll.gCal(i_negy(1:end));
% iterative linear fits to remove outliers
iys=true;
cor_ys=a_ys;
ni=length(cor_ys);
while any(iys)
    P=polyfit(t_ys-t_ys(1),cor_ys,1);
    lin_ys=polyval(P,t_ys-t_ys(1));
    hold off
    temp_ys=cor_ys-lin_ys;
    [~,iys]=rmoutliers(temp_ys);
    cor_ys=cor_ys(~iys);
    t_ys=t_ys(~iys);
end
nf=length(cor_ys);
hold on
plot(t_ys,cor_ys-cor_ys(1),'k^')
plot(t_ys,lin_ys-cor_ys(1),'k-')
drift_ys=round(P(1)*365*10^5,2);
text(t_ys(10),3e-5,['drift = ' num2str(drift_ys) ' \mug/yr'])
text(t_ys(10),2.8e-5,['\sigma = ' num2str(round(rms(cor_ys-lin_ys)*10^5,2)) ' \mug'])
datetick('x',3)
xtickangle(45)
ylabel('\Delta Y span (m/s^2)')
set(gca,'fontsize',12)
ylim([0 0.00004])
xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/20,mean(yl)-[0.000005 -0.000005],...
    '-k','linewidth',2); text(xl(1)+diff(xl)/15,mean(yl),'1 \mug','fontsize',12)
box on; grid on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print('../calibrations/Axial/drift_rates/y_channels','-dtiff','-r300')