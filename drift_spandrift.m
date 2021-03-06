% drift_spandrift.m
%
% Calculates the drift and span drift of each channel
%

close all; clear

Axial=true;
PF=false;
apply_T_cor=false;

%% Axial
if Axial
    load('../calibrations/Axial/detailed_flipInfo');
    load('../calibrations/Axial/axialdata','flipInfoAll');
    
    % identify bad calibrations
    bad_x1=[datenum(2019,12,01) datenum(2020,01,08) datenum(2020,01,15) datenum(2020,01,22) datenum(2020,01,29) ...
        datenum(2020,09,09)];
    bad_y=[datenum(2019,12,01) datenum(2020,01,08) datenum(2020,06,03)];
    bad_negy=[datenum(2019,12,01) datenum(2020,01,08) datenum(2020,01,15) datenum(2020,01,22) datenum(2020,02,19) ...
        datenum(2020,03,18) datenum(2020,06,03) datenum(2020,06,17) datenum(2020,09,09)];
    bad_x2=[datenum(2019,12,01) datenum(2020,01,08) datenum(2020,01,15) datenum(2020,04,29) ...
        datenum(2020,09,02)];
    bad_negx=[datenum(2019,11,01) datenum(2019,12,01) datenum(2020,09,02)];
    bads={bad_x1, bad_y, bad_negy, bad_x2, bad_negx};
    ubads=cat(2,bad_x1,bad_y,bad_negy,bad_x2,bad_negx); ubads=unique(ubads);
    save('../calibrations/Axial/badcaldates','bad_x1','bad_x2','bad_negx','bad_y','bad_negy','bads','ubads')
    
    % find linear best fits (and drifts) to only good calibrations
    x1_t=[floor(flipInfoAll.t(1:3:73));floor(flipInfoSome.t(1:5:end))];
    x1=[flipInfoAll.gCalTCor(1:3:73);flipInfoSome.gCalTCor(1:5:end)];
    [~,ix1]=setdiff(x1_t,bad_x1);
    % all cals
    px1=polyfit(x1_t(ix1)-x1_t(ix1(1)),x1(ix1),1);
    x1drift=px1(1)*365;
    x1_m=px1(1)*(x1_t-x1_t(1))+px1(2);
    % only cals under +/- flipping scheme
    px1=polyfit(x1_t(ix1)-x1_t(ix1(1)),x1(ix1),1);
    x1drift=px1(1)*365;
    x1_m=px1(1)*(x1_t-x1_t(1))+px1(2);
    
    y_t=[floor(flipInfoAll.t(2:3:74));floor(flipInfoSome.t(2:5:end))];
    y=[flipInfoAll.gCalTCor(2:3:74);flipInfoSome.gCalTCor(2:5:end)];
    [~,iy]=setdiff(y_t,bad_y);
    py=polyfit(y_t(iy)-y_t(iy(1)),y(iy),1);
    ydrift=py(1)*365;
    y_m=py(1)*(y_t-y_t(1))+py(2);
    
    [~,inegy]=setdiff(floor(flipInfoSome.t),bad_negy); inegy=inegy+2;
    pnegy=polyfit(flipInfoSome.t(inegy)-flipInfoSome.t(inegy(1)),flipInfoSome.gCalTCor(inegy),1);
    negydrift=pnegy(1)*365;
    negy_m=pnegy(1)*(flipInfoSome.t(3:5:end)-flipInfoSome.t(3))+pnegy(2);
    
    [~,ix2]=setdiff(floor(flipInfoSome.t),bad_x2); ix2=ix2+3;
    px2=polyfit(flipInfoSome.t(ix2)-flipInfoSome.t(ix2(1)),flipInfoSome.gCalTCor(ix2),1);
    x2drift=px2(1)*365;
    x2_m=px2(1)*(flipInfoSome.t(4:5:end)-flipInfoSome.t(4))+px2(2);
    
    [~,inegx]=setdiff(floor(flipInfoSome.t),bad_negx); inegx=inegx+4;
    pnegx=polyfit(flipInfoSome.t(inegx)-flipInfoSome.t(inegx(1)),flipInfoSome.gCalTCor(inegx),1);
    negxdrift=pnegx(1)*365;
    negx_m=pnegx(1)*(flipInfoSome.t(5:5:end)-flipInfoSome.t(5))+pnegx(2);
    
    [~,ixs1]=setdiff(floor(flipInfoSome.t),[bad_x1,bad_negx]);
    pxs1=polyfit(flipInfoSome.t(ixs1)-flipInfoSome.t(ixs1(1)),...
        (flipInfoSome.gCalTCor(ixs1)+flipInfoSome.gCalTCor(ixs1+4))/2,1);
    xspandrift1=pxs1(1)*365;
    xs1_m=pxs1(1)*(flipInfoSome.t(1:5:end)-flipInfoSome.t(1))+pxs1(2);
    
    [~,ixs2]=setdiff(floor(flipInfoSome.t),[bad_x2,bad_negx]); ixs2=ixs2+3;
    pxs2=polyfit(flipInfoSome.t(ixs2)-flipInfoSome.t(ixs2(1)),...
        (flipInfoSome.gCalTCor(ixs2)+flipInfoSome.gCalTCor(ixs2+1))/2,1);
    xspandrift2=pxs2(1)*365;
    xs2_m=pxs2(1)*(flipInfoSome.t(4:5:end)-flipInfoSome.t(4))+pxs2(2);
    
    [~,iys]=setdiff(floor(flipInfoSome.t),[bad_y,bad_negy]); iys=iys+1;
    pys=polyfit(flipInfoSome.t(iys)-flipInfoSome.t(iys(1)),...
        (flipInfoSome.gCalTCor(iys)+flipInfoSome.gCalTCor(iys+1))/2,1);
    yspandrift=pys(1)*365;
    ys_m=pys(1)*(flipInfoSome.t(2:5:end)-flipInfoSome.t(2))+pys(2);
    
    % plotting
    figure(50)
    clf
    hold on
    plot(x1_t,x1-mean(x1_m),'+r','markersize',18,'linewidth',2);
    h1=plot(x1_t(ix1),x1(ix1)-mean(x1_m),'+k','markersize',18,'linewidth',2);
    plot(flipInfoSome.t(5:5:end),flipInfoSome.gCalTCor(5:5:end)-mean(negx_m),'xr','markersize',18,'linewidth',2);
    h2=plot(flipInfoSome.t(inegx),flipInfoSome.gCalTCor(inegx)-mean(negx_m),'xk','markersize',18,'linewidth',2);
    h3=plot(flipInfoSome.t(ixs1),(flipInfoSome.gCalTCor(ixs1)+flipInfoSome.gCalTCor(ixs1+4))/2 ...
        -mean(xs1_m),'^k','markerfacecolor','k','markersize',18);
    plot(x1_t,x1_m-mean(x1_m),'--k')
    plot(flipInfoSome.t(5:5:end),negx_m-mean(negx_m),'--k')
    plot(flipInfoSome.t(1:5:end),xs1_m-mean(xs1_m),'k','linewidth',1)
    text(flipInfoSome.t(end)-40,x1_m(end)-mean(x1_m)+10^-5,[num2str(round(x1drift*10^5,2)) ' \mug/yr'],'fontsize',14)
    text(flipInfoSome.t(end)-40,negx_m(end)-mean(negx_m)-10^-5,[num2str(round(negxdrift*10^5,2)) ' \mug/yr'],'fontsize',14)
    text(flipInfoSome.t(end)-40,xs1_m(end)-mean(xs1_m)+10^-5,[num2str(round(xspandrift1*10^5,2)) ' \mug/yr'],'fontsize',14)
    legend([h1,h2,h3],'+X calibration','-X calibration','X span','location','northwest')
    datetick('x',6)
    title({'Axial SCTA X Calibrations',[datestr(x1_t(1),'mmm dd, yyyy') ...
        ' - ' datestr(flipInfoSome.t(end),'mmm dd, yyyy')]})
    ylabel('Acceleration (m/s^2)')
    set(gca,'fontsize',18)
    xtickangle(45)
    box on
    
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print -djpeg ../calibrations/Axial/span/x1spandrift.jpeg
    print -dtiff ../calibrations/Axial/span/x1spandrift.tiff -r300
%     !open ../calibrations/Axial/span/x1spandrift.jpeg

    figure(51)
    clf
    hold on
    plot(flipInfoSome.t(4:5:end),flipInfoSome.gCalTCor(4:5:end)-mean(x2_m),'+r','markersize',18,'linewidth',2);
    h1=plot(flipInfoSome.t(ix2),flipInfoSome.gCalTCor(ix2)-mean(x2_m),'+k','markersize',18,'linewidth',2);
    plot(flipInfoSome.t(5:5:end),flipInfoSome.gCalTCor(5:5:end)-mean(negx_m),'xr','markersize',18,'linewidth',2);
    h2=plot(flipInfoSome.t(inegx),flipInfoSome.gCalTCor(inegx)-mean(negx_m),'xk','markersize',18,'linewidth',2);
    h3=plot(flipInfoSome.t(ixs2),(flipInfoSome.gCalTCor(ixs2)+flipInfoSome.gCalTCor(ixs2+1))/2 ...
        -mean(xs2_m),'^k','markerfacecolor','k','markersize',18);
    plot(flipInfoSome.t(4:5:end),x2_m-mean(x2_m),'--k')
    plot(flipInfoSome.t(5:5:end),negx_m-mean(negx_m),'--k')
    plot(flipInfoSome.t(4:5:end),xs2_m-mean(xs2_m),'k','linewidth',1)
    text(flipInfoSome.t(end)-40,x2_m(end)-mean(x2_m)+10^-5,[num2str(round(x2drift*10^5,2)) ' \mug/yr'],'fontsize',14)
    text(flipInfoSome.t(end)-40,negx_m(end)-mean(negx_m)-10^-5,[num2str(round(negxdrift*10^5,2)) ' \mug/yr'],'fontsize',14)
    text(flipInfoSome.t(end)-40,xs2_m(end)-mean(xs2_m)+10^-5,[num2str(round(xspandrift2*10^5,2)) ' \mug/yr'],'fontsize',14)
    legend([h1,h2,h3],'+X calibration','-X calibration','X span','location','northwest')
    datetick('x',6)
    title({'Axial SCTA X Calibrations',[datestr(flipInfoSome.t(1),'mmm dd, yyyy') ...
        ' - ' datestr(flipInfoSome.t(end),'mmm dd, yyyy')]})
    ylabel('Acceleration (m/s^2)')
    set(gca,'fontsize',18)
    xtickangle(45)
    box on
    
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print -djpeg ../calibrations/Axial/span/x2spandrift.jpeg
    print -dtiff ../calibrations/Axial/span/x2spandrift.tiff -r300
%     !open ../calibrations/Axial/span/x2spandrift.jpeg
    
    figure(52)
    clf
    hold on
    plot(y_t,y-mean(y_m),'+r','markersize',18,'linewidth',2);
    h1=plot(y_t(iy),y(iy)-mean(y_m),'+k','markersize',18,'linewidth',2);
    plot(flipInfoSome.t(3:5:end),flipInfoSome.gCalTCor(3:5:end)-mean(negy_m),'xr','markersize',18,'linewidth',2);
    h2=plot(flipInfoSome.t(inegy),flipInfoSome.gCalTCor(inegy)-mean(negy_m),'xk','markersize',18,'linewidth',2);
    h3=plot(flipInfoSome.t(iys),(flipInfoSome.gCalTCor(iys)+flipInfoSome.gCalTCor(iys+1))/2 ...
        -mean(ys_m),'^k','markerfacecolor','k','markersize',18);
    plot(y_t,y_m-mean(y_m),'--k')
    plot(flipInfoSome.t(3:5:end),negy_m-mean(negy_m),'--k')
    plot(flipInfoSome.t(2:5:end),ys_m-mean(ys_m),'k','linewidth',1)
    text(flipInfoSome.t(end)-40,y_m(end)-mean(y_m)+10^-5,[num2str(round(ydrift*10^5,2)) ' \mug/yr'],'fontsize',14)
    text(flipInfoSome.t(end)-40,negy_m(end)-mean(negy_m)-10^-5,[num2str(round(negydrift*10^5,2)) ' \mug/yr'],'fontsize',14)
    text(flipInfoSome.t(end)-40,ys_m(end)-mean(ys_m)+10^-5,[num2str(round(yspandrift*10^5,2)) ' \mug/yr'],'fontsize',14)
    legend([h1,h2,h3],'+Y calibration','-Y calibration','Y span','location','northwest')
    datetick('x',6)
    title({'Axial SCTA Y Calibrations',[datestr(y_t(1),'mmm dd, yyyy') ...
        ' - ' datestr(flipInfoSome.t(end),'mmm dd, yyyy')]})
    ylabel('Acceleration (m/s^2)')
    set(gca,'fontsize',18)
    xtickangle(45)
    box on
    
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print -djpeg ../calibrations/Axial/span/yspandrift.jpeg
    print -dtiff ../calibrations/Axial/span/yspandrift.tiff -r300
%     !open ../calibrations/Axial/span/yspandrift.jpeg
end

%% Pinon Flat
if PF
    load ../calibrations/PinonFlat/PFdata
    
    t_X2=flipInfoSome.t(64:5:end)-flipInfoSome.t(64);
    a_X2=flipInfoSome.gCal(64:5:end);
    t_negX=flipInfoSome.t(65:5:end)-flipInfoSome.t(65);
    a_negX=flipInfoSome.gCal(65:5:end);
    t_Y=flipInfoSome.t(62:5:end)-flipInfoSome.t(62);
    a_Y=flipInfoSome.gCal(62:5:end);
    t_negY=flipInfoSome.t(63:5:end)-flipInfoSome.t(63);
    a_negY=flipInfoSome.gCal(63:5:end);
    if apply_T_cor
        load('../calibrations/PinonFlat/T_dependence/X1.mat','m_X')
        load('../calibrations/PinonFlat/T_dependence/Y.mat','m_Y')
        a_X2=a_X2-m_X(3)*flipInfoSome.T(64:5:end);
        a_negX=a_negX+m_X(3)*flipInfoSome.T(65:5:end);
        a_Y=a_Y-m_Y(3)*flipInfoSome.T(62:5:end);
        a_negY=a_negY+m_Y(3)*flipInfoSome.T(63:5:end);
    end
    
    % calculate drift
    px1=polyfit(t_X2,a_X2-mean(a_X2),1);
    x1drift=px1(1)*365;
    x1_m=px1(1)*(t_X2)+px1(2);
    px2=polyfit(t_negX,a_negX-mean(a_negX),1);
    x2drift=px2(1)*365;
    x2_m=px2(1)*(t_negX)+px2(2);
    pxs2=polyfit(t_negX,(a_X2+a_negX)/2-mean((a_X2+a_negX)/2),1);
    xspandrift2=pxs2(1)*365;
    xs2_m=pxs2(1)*(t_negX)+pxs2(2);
    
    py=polyfit(t_Y,a_Y-mean(a_Y),1);
    ydrift=py(1)*365;
    y1_m=py(1)*(t_Y)+py(2);
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
    plot(flipInfoSome.t(64:5:end),a_X2-mean(a_X2),'ob','markersize',18);
    plot(flipInfoSome.t(65:5:end),a_negX-mean(a_negX),'sb','markersize',18);
    plot(flipInfoSome.t(65:5:end),(a_X2+a_negX)/2 ...
        -mean((a_X2+a_negX)/2),'^k','markerfacecolor','k','markersize',18);
    plot(flipInfoSome.t(64:5:end),x1_m,'b')
    plot(flipInfoSome.t(65:5:end),x2_m,'b')
    plot(flipInfoSome.t(65:5:end),xs2_m,'k','linewidth',1)
    text(flipInfoSome.t(end)-2,x1_m(end)+10^-6,[num2str(round(x1drift*10^5,2)) ' \mug/yr'],'fontsize',15)
    text(flipInfoSome.t(end)-2,x2_m(end)-10^-6,[num2str(round(x2drift*10^5,2)) ' \mug/yr'],'fontsize',15)
    text(flipInfoSome.t(end)-2,xs2_m(end)+10^-6,[num2str(round(xspandrift2*10^5,2)) ' \mug/yr'],'fontsize',15)
%     xl = xlim; yl = ylim;
%     plot([0 0]+xl(1)+diff(xl)/10,mean(yl)+[0 0.00001],'-k')
%     text(xl(1)+diff(xl)/9,mean(yl)+0.000005,'10^{-6} g')
    legend('+X','-X','X span','location','best')
    datetick
    title({'PF SCTA X Calibrations',[datestr(flipInfoSome.t(64),'mmm dd, yyyy') ...
        ' - ' datestr(flipInfoSome.t(end),'mmm dd, yyyy')]})
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
    plot(flipInfoSome.t(62:5:end),a_Y-mean(a_Y),'or','markersize',18);
    plot(flipInfoSome.t(63:5:end),a_negY-mean(a_negY),'sr','markersize',18);
    plot(flipInfoSome.t(63:5:end),(a_Y+a_negY)/2 ...
        -mean((a_Y+a_negY)/2),'^k','markerfacecolor','k','markersize',18);
    plot(flipInfoSome.t(62:5:end),y1_m,'r')
    plot(flipInfoSome.t(63:5:end),y2_m,'r')
    plot(flipInfoSome.t(63:5:end),ys_m,'k','linewidth',1)
    text(flipInfoSome.t(end)-2,y1_m(end)+10^-6,[num2str(round(ydrift*10^5,2)) ' \mug/yr'],'fontsize',15)
    text(flipInfoSome.t(end)-2,y2_m(end)-10^-6,[num2str(round(y2drift*10^5,2)) ' \mug/yr'],'fontsize',15)
    text(flipInfoSome.t(end)-2,ys_m(end)+10^-6,[num2str(round(yspandrift*10^5,2)) ' \mug/yr'],'fontsize',15)
%     xl = xlim; yl = ylim;
%     plot([0 0]+xl(1)+diff(xl)/10,mean(yl)+[0 0.00001],'-k')
%     text(xl(1)+diff(xl)/9,mean(yl)+0.000005,'10^{-6} g')
    legend('+Y','-Y','Y span','location','best')
    datetick
    title({'PF SCTA Y Calibrations',[datestr(flipInfoSome.t(64),'mmm dd, yyyy') ...
        ' - ' datestr(flipInfoSome.t(end),'mmm dd, yyyy')]})
    xlabel('Date')
    ylabel('Calibration, m/s^2')
    set(gca,'fontsize',15)
    
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print -djpeg ../calibrations/PinonFlat/span/yspandrift.jpeg
    print -dtiff ../calibrations/PinonFlat/span/yspandrift.tiff -r300
end