% drift_spandrift.m
%
% Calculates the drift and span drift of each channel
%

close all; clear

Axial=true;
PF=false;
apply_T_cor=false;
exclude_dec=true;

%% Axial
if Axial
    load ../calibrations/Axial/axialdata flipInfoAll
    
    if exclude_dec
        sbset_x1=[76:5:106,121:5:length(flipInfoAll.t)];
        sbset_x2=[77:5:107,122:5:length(flipInfoAll.t)];
        sbset_y1=[74:5:104,119:5:length(flipInfoAll.t)];
        sbset_y2=[75:5:105,120:5:length(flipInfoAll.t)];
    else
        sbset_x1=76:5:length(flipInfoAll.t);
        sbset_x2=77:5:length(flipInfoAll.t);
        sbset_y1=74:5:length(flipInfoAll.t);
        sbset_y2=75:5:length(flipInfoAll.t);
    end
    
    % calculate drift
    px1=polyfit(flipInfoAll.t(sbset_x1)-flipInfoAll.t(73),...
        flipInfoAll.gCalTCor(sbset_x1)-mean(flipInfoAll.gCalTCor(sbset_x1)),1);
    x1drift=px1(1)*365;
    x1_m=px1(1)*(flipInfoAll.t(sbset_x1)-flipInfoAll.t(73))+px1(2);
    px2=polyfit(flipInfoAll.t(sbset_x2)-flipInfoAll.t(74),...
        flipInfoAll.gCalTCor(sbset_x2)-mean(flipInfoAll.gCalTCor(sbset_x2)),1);
    x2drift=px2(1)*365;
    x2_m=px2(1)*(flipInfoAll.t(sbset_x2)-flipInfoAll.t(74))+px2(2);
    pxs=polyfit(flipInfoAll.t(sbset_x2)-flipInfoAll.t(74),(flipInfoAll.gCalTCor(sbset_x1)+flipInfoAll.gCalTCor(sbset_x2))/2 ...
        -mean((flipInfoAll.gCalTCor(sbset_x1)+flipInfoAll.gCalTCor(sbset_x2))/2),1);
    xspandrift=pxs(1)*365;
    xs_m=pxs(1)*(flipInfoAll.t(sbset_x2)-flipInfoAll.t(74))+pxs(2);
    
    py1=polyfit(flipInfoAll.t(sbset_y1)-flipInfoAll.t(71),...
        flipInfoAll.gCalTCor(sbset_y1)-mean(flipInfoAll.gCalTCor(sbset_y1)),1);
    y1drift=py1(1)*365;
    y1_m=py1(1)*(flipInfoAll.t(sbset_y1)-flipInfoAll.t(71))+py1(2);
    py2=polyfit(flipInfoAll.t(sbset_y2)-flipInfoAll.t(72),...
        flipInfoAll.gCalTCor(sbset_y2)-mean(flipInfoAll.gCalTCor(sbset_y2)),1);
    y2drift=py2(1)*365;
    y2_m=py2(1)*(flipInfoAll.t(sbset_y2)-flipInfoAll.t(72))+py2(2);
    pys=polyfit(flipInfoAll.t(sbset_y2)-flipInfoAll.t(72),(flipInfoAll.gCalTCor(sbset_y1)+flipInfoAll.gCalTCor(sbset_y2))/2 ...
        -mean((flipInfoAll.gCalTCor(sbset_y1)+flipInfoAll.gCalTCor(sbset_y2))/2),1);
    yspandrift=pys(1)*365;
    ys_m=pys(1)*(flipInfoAll.t(sbset_y2)-flipInfoAll.t(72))+pys(2);
    
    % plotting
    figure(51)
    clf
    hold on
    plot(flipInfoAll.t(sbset_x1),flipInfoAll.gCalTCor(sbset_x1)-mean(flipInfoAll.gCalTCor(sbset_x1)),'xk','markersize',18);
    plot(flipInfoAll.t(sbset_x2),flipInfoAll.gCalTCor(sbset_x2)-mean(flipInfoAll.gCalTCor(sbset_x2)),'+k','markersize',18);
    plot(flipInfoAll.t(sbset_x2),(flipInfoAll.gCalTCor(sbset_x1)+flipInfoAll.gCalTCor(sbset_x2))/2 ...
        -mean((flipInfoAll.gCalTCor(sbset_x1)+flipInfoAll.gCalTCor(sbset_x2))/2),'^k','markerfacecolor','k','markersize',18);
    plot(flipInfoAll.t(sbset_x1),x1_m,'--k')
    plot(flipInfoAll.t(sbset_x2),x2_m,'--k')
    plot(flipInfoAll.t(sbset_x2),xs_m,'k','linewidth',1)
    text(flipInfoAll.t(end)-40,x1_m(end)+10^-6,[num2str(round(x1drift*10^5,2)) ' \mug/yr'],'fontsize',14)
    text(flipInfoAll.t(end)-40,x2_m(end)-10^-6,[num2str(round(x2drift*10^5,2)) ' \mug/yr'],'fontsize',14)
    text(flipInfoAll.t(end)-40,xs_m(end)+10^-6,[num2str(round(xspandrift*10^5,2)) ' \mug/yr'],'fontsize',14)
%     xl = xlim; yl = ylim;
%     plot([0 0]+xl(1)+diff(xl)/10,mean(yl)+[0 0.00001],'-k')
%     text(xl(1)+diff(xl)/9,mean(yl)+0.000005,'10^{-6} g')
    legend('+X calibration','-X calibration','X span','location','northwest')
    datetick
    title({'Axial SCTA X Calibrations',[datestr(flipInfoAll.t(73),'mmm dd, yyyy') ...
        ' - ' datestr(flipInfoAll.t(end),'mmm dd, yyyy')]})
    ylabel('Acceleration (m/s^2)')
    set(gca,'fontsize',18)
    xtickangle(45)
    box on
    
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print -djpeg ../calibrations/Axial/span/xspandrift_noDecJan.jpeg
    print -dtiff ../calibrations/Axial/span/xspandrift_noDecJan.tiff -r300
%     !open ../calibrations/Axial/span/xspandrift3.jpeg
    
    figure(52)
    clf
    hold on
    plot(flipInfoAll.t(sbset_y1),flipInfoAll.gCalTCor(sbset_y1)-mean(flipInfoAll.gCalTCor(sbset_y1)),'xb','markersize',18);
    plot(flipInfoAll.t(sbset_y2),flipInfoAll.gCalTCor(sbset_y2)-mean(flipInfoAll.gCalTCor(sbset_y2)),'+b','markersize',18);
    plot(flipInfoAll.t(sbset_y2),(flipInfoAll.gCalTCor(sbset_y1)+flipInfoAll.gCalTCor(sbset_y2))/2 ...
        -mean((flipInfoAll.gCalTCor(sbset_y1)+flipInfoAll.gCalTCor(sbset_y2))/2),'^b','markerfacecolor','b','markersize',18);
    plot(flipInfoAll.t(sbset_y1),y1_m,'--b')
    plot(flipInfoAll.t(sbset_y2),y2_m,'--b')
    plot(flipInfoAll.t(sbset_y2),ys_m,'b','linewidth',1)
    text(flipInfoAll.t(end)-40,y1_m(end)+10^-6,[num2str(round(y1drift*10^5,2)) ' \mug/yr'],'fontsize',14)
    text(flipInfoAll.t(end)-40,y2_m(end)-10^-6,[num2str(round(y2drift*10^5,2)) ' \mug/yr'],'fontsize',14)
    text(flipInfoAll.t(end)-40,ys_m(end)+10^-6,[num2str(round(yspandrift*10^5,2)) ' \mug/yr'],'fontsize',14)
%     xl = xlim; yl = ylim;
%     plot([0 0]+xl(1)+diff(xl)/10,mean(yl)+[0 0.00001],'-k')
%     text(xl(1)+diff(xl)/9,mean(yl)+0.000005,'10^{-6} g')
    legend('+Y calibration','-Y calibration','Y span','location','northwest')
    datetick
    title({'Axial SCTA Y Calibrations',[datestr(flipInfoAll.t(73),'mmm dd, yyyy') ...
        ' - ' datestr(flipInfoAll.t(end),'mmm dd, yyyy')]})
    ylabel('Acceleration (m/s^2)')
    set(gca,'fontsize',18)
    xtickangle(45)
    box on
    
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print -djpeg ../calibrations/Axial/span/yspandrift_noDecJan.jpeg
    print -dtiff ../calibrations/Axial/span/yspandrift_noDecJan.tiff -r300
%     !open ../calibrations/Axial/span/yspandrift3.jpeg
end

%% Pinon Flat
if PF
    load ../calibrations/PinonFlat/PFdata
    
    t_X2=flipInfoAll.t(64:5:end)-flipInfoAll.t(64);
    a_X2=flipInfoAll.gCal(64:5:end);
    t_negX=flipInfoAll.t(65:5:end)-flipInfoAll.t(65);
    a_negX=flipInfoAll.gCal(65:5:end);
    t_Y=flipInfoAll.t(62:5:end)-flipInfoAll.t(62);
    a_Y=flipInfoAll.gCal(62:5:end);
    t_negY=flipInfoAll.t(63:5:end)-flipInfoAll.t(63);
    a_negY=flipInfoAll.gCal(63:5:end);
    if apply_T_cor
        load('../calibrations/PinonFlat/T_dependence/X1.mat','m_X')
        load('../calibrations/PinonFlat/T_dependence/Y.mat','m_Y')
        a_X2=a_X2-m_X(3)*flipInfoAll.T(64:5:end);
        a_negX=a_negX+m_X(3)*flipInfoAll.T(65:5:end);
        a_Y=a_Y-m_Y(3)*flipInfoAll.T(62:5:end);
        a_negY=a_negY+m_Y(3)*flipInfoAll.T(63:5:end);
    end
    
    % calculate drift
    px1=polyfit(t_X2,a_X2-mean(a_X2),1);
    x1drift=px1(1)*365;
    x1_m=px1(1)*(t_X2)+px1(2);
    px2=polyfit(t_negX,a_negX-mean(a_negX),1);
    x2drift=px2(1)*365;
    x2_m=px2(1)*(t_negX)+px2(2);
    pxs=polyfit(t_negX,(a_X2+a_negX)/2-mean((a_X2+a_negX)/2),1);
    xspandrift=pxs(1)*365;
    xs_m=pxs(1)*(t_negX)+pxs(2);
    
    py1=polyfit(t_Y,a_Y-mean(a_Y),1);
    y1drift=py1(1)*365;
    y1_m=py1(1)*(t_Y)+py1(2);
    py2=polyfit(t_negY,a_negY-mean(a_negY),1);
    y2drift=py2(1)*365;
    y2_m=py2(1)*(t_negY)+py2(2);
    pys=polyfit(t_negY,(a_Y+a_negY)/2-mean((a_Y+a_negY)/2),1);
    yspandrift=pys(1)*365;
    ys_m=pys(1)*(t_negY)+pys(2);
    
    % plotting
    figure(51)
    clf
    hold on
    plot(flipInfoAll.t(64:5:end),a_X2-mean(a_X2),'ob','markersize',18);
    plot(flipInfoAll.t(65:5:end),a_negX-mean(a_negX),'sb','markersize',18);
    plot(flipInfoAll.t(65:5:end),(a_X2+a_negX)/2 ...
        -mean((a_X2+a_negX)/2),'^k','markerfacecolor','k','markersize',18);
    plot(flipInfoAll.t(64:5:end),x1_m,'b')
    plot(flipInfoAll.t(65:5:end),x2_m,'b')
    plot(flipInfoAll.t(65:5:end),xs_m,'k','linewidth',1)
    text(flipInfoAll.t(end)-2,x1_m(end)+10^-6,[num2str(round(x1drift*10^5,2)) ' \mug/yr'],'fontsize',15)
    text(flipInfoAll.t(end)-2,x2_m(end)-10^-6,[num2str(round(x2drift*10^5,2)) ' \mug/yr'],'fontsize',15)
    text(flipInfoAll.t(end)-2,xs_m(end)+10^-6,[num2str(round(xspandrift*10^5,2)) ' \mug/yr'],'fontsize',15)
%     xl = xlim; yl = ylim;
%     plot([0 0]+xl(1)+diff(xl)/10,mean(yl)+[0 0.00001],'-k')
%     text(xl(1)+diff(xl)/9,mean(yl)+0.000005,'10^{-6} g')
    legend('+X','-X','X span','location','best')
    datetick
    title({'PF SCTA X Calibrations',[datestr(flipInfoAll.t(64),'mmm dd, yyyy') ...
        ' - ' datestr(flipInfoAll.t(end),'mmm dd, yyyy')]})
    xlabel('Date')
    ylabel('Calibration, m/s^2')
    set(gca,'fontsize',15)
    
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print -djpeg ../calibrations/PinonFlat/span/xspandrift.jpeg
    print -dtiff ../calibrations/PinonFlat/span/xspandrift.tiff -r300
    
    figure(52)
    clf
    hold on
    plot(flipInfoAll.t(62:5:end),a_Y-mean(a_Y),'or','markersize',18);
    plot(flipInfoAll.t(63:5:end),a_negY-mean(a_negY),'sr','markersize',18);
    plot(flipInfoAll.t(63:5:end),(a_Y+a_negY)/2 ...
        -mean((a_Y+a_negY)/2),'^k','markerfacecolor','k','markersize',18);
    plot(flipInfoAll.t(62:5:end),y1_m,'r')
    plot(flipInfoAll.t(63:5:end),y2_m,'r')
    plot(flipInfoAll.t(63:5:end),ys_m,'k','linewidth',1)
    text(flipInfoAll.t(end)-2,y1_m(end)+10^-6,[num2str(round(y1drift*10^5,2)) ' \mug/yr'],'fontsize',15)
    text(flipInfoAll.t(end)-2,y2_m(end)-10^-6,[num2str(round(y2drift*10^5,2)) ' \mug/yr'],'fontsize',15)
    text(flipInfoAll.t(end)-2,ys_m(end)+10^-6,[num2str(round(yspandrift*10^5,2)) ' \mug/yr'],'fontsize',15)
%     xl = xlim; yl = ylim;
%     plot([0 0]+xl(1)+diff(xl)/10,mean(yl)+[0 0.00001],'-k')
%     text(xl(1)+diff(xl)/9,mean(yl)+0.000005,'10^{-6} g')
    legend('+Y','-Y','Y span','location','best')
    datetick
    title({'PF SCTA Y Calibrations',[datestr(flipInfoAll.t(64),'mmm dd, yyyy') ...
        ' - ' datestr(flipInfoAll.t(end),'mmm dd, yyyy')]})
    xlabel('Date')
    ylabel('Calibration, m/s^2')
    set(gca,'fontsize',15)
    
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print -djpeg ../calibrations/PinonFlat/span/yspandrift.jpeg
    print -dtiff ../calibrations/PinonFlat/span/yspandrift.tiff -r300
end