load ../calibrations/PinonFlat/PFdata.mat

% Plot PF calibrations w/out T correction
figure
clf
plot(flipInfoAll.t(1:3:60),flipInfoAll.gCal(1:3:60),'ok','markersize',18);
hold on
plot(flipInfoAll.t(2:3:60),flipInfoAll.gCal(2:3:60),'or','markersize',18);
plot(flipInfoAll.t(3:3:60),flipInfoAll.gCal(3:3:60),'ob','markersize',18);
plot(flipInfoAll.t(65:5:end),flipInfoAll.gCal(65:5:end),'sk','markersize',18);
plot(flipInfoAll.t(63:5:end),flipInfoAll.gCal(63:5:end),'sr','markersize',18);
plot(flipInfoAll.t(65:5:end),(flipInfoAll.gCal(64:5:end)+flipInfoAll.gCal(65:5:end))/2,'^k','markersize',18);
plot(flipInfoAll.t(63:5:end),(flipInfoAll.gCal(62:5:end)+flipInfoAll.gCal(63:5:end))/2,'^r','markersize',18);
plot(flipInfoAll.t(61:5:end),flipInfoAll.gCal(61:5:end),'ok','markersize',18);
plot(flipInfoAll.t(62:5:end),flipInfoAll.gCal(62:5:end),'or','markersize',18);
plot(flipInfoAll.t(64:5:end),flipInfoAll.gCal(64:5:end),'ob','markersize',18);
xl = xlim; yl = ylim;
plot([0 0]+xl(1)+diff(xl)/10,mean(yl)+[0 0.0001],'-k')
text(xl(1)+diff(xl)/9,mean(yl)+0.00005,'10^{-5} g')
legend('1st X','Y','2nd X','-X','-Y','X span','Y span','location','best')
datetick
xtickangle(45)
title({'Pinon Flat SCTA Calibrations',[datestr(flipInfoAll.t(1),'mmm dd, yyyy')...
    ' - ' datestr(flipInfoAll.t(end),'mmm dd, yyyy')]})
ylabel('Calibration, m/s^2')
set(gca,'FontSize',15)
yyaxis right
plot(dataDec100.t,dataDec100.T,'linewidth',1)
ylabel(['Temperature (' char(176) 'C)'])

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print -dtiff ../calibrations/PinonFlat/PF_cals_and_T.tiff -r300

% % now w/ T correction
% figure
% clf
% plot(flipInfoAll.t(1:3:60),flipInfoAll.gCal(1:3:60),'ok','markersize',18);
% hold on
% plot(flipInfoAll.t(2:3:60),flipInfoAll.gCal(2:3:60),'or','markersize',18);
% plot(flipInfoAll.t(3:3:60),flipInfoAll.gCal(3:3:60),'ob','markersize',18);
% plot(flipInfoAll.t(65:5:end),flipInfoAll.gCal(65:5:end),'sk','markersize',18);
% plot(flipInfoAll.t(63:5:end),flipInfoAll.gCal(63:5:end),'sr','markersize',18);
% plot(flipInfoAll.t(65:5:end),(flipInfoAll.gCal(64:5:end)+flipInfoAll.gCal(65:5:end))/2,'^k','markersize',18);
% plot(flipInfoAll.t(63:5:end),(flipInfoAll.gCal(62:5:end)+flipInfoAll.gCal(63:5:end))/2,'^r','markersize',18);
% plot(flipInfoAll.t(61:5:end),flipInfoAll.gCal(61:5:end),'ok','markersize',18);
% plot(flipInfoAll.t(62:5:end),flipInfoAll.gCal(62:5:end),'or','markersize',18);
% plot(flipInfoAll.t(64:5:end),flipInfoAll.gCal(64:5:end),'ob','markersize',18);
% xl = xlim; yl = ylim;
% plot([0 0]+xl(1)+diff(xl)/10,mean(yl)+[0 0.0001],'-k')
% text(xl(1)+diff(xl)/9,mean(yl)+0.00005,'10^{-5} g')
% legend('1st X','Y','2nd X','-X','-Y','X span','Y span','location','best')
% datetick
% xtickangle(45)
% title({'Pinon Flat SCTA Calibrations',[datestr(flipInfoAll.t(1),'mmm dd, yyyy')...
%     ' - ' datestr(flipInfoAll.t(end),'mmm dd, yyyy')]})
% ylabel('Calibration, m/s^2')
% set(gca,'FontSize',15)
% yyaxis right
% plot(dataDec100.t,dataDec100.T,'linewidth',1)
% ylabel(['Temperature (' char(176) 'C)'])
% 
% fh=gcf;
% fh.PaperUnits='inches';
% fh.PaperPosition=[0 0 11 8.5];
% print -dtiff ../calibrations/PinonFlat/PF_cals_and_T.tiff -r300