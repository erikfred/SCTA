figure
clf
plot(flipInfoAll.t(1:3:69),flipInfoAll.gCalTCor(1:3:69),'xr','markersize',18);
hold on
plot(flipInfoAll.t(2:3:69),flipInfoAll.gCalTCor(2:3:69),'xb','markersize',18);
plot(flipInfoAll.t(3:3:69),flipInfoAll.gCalTCor(3:3:69),'xk','markersize',18);
plot(flipInfoAll.t(74:5:end-5),flipInfoAll.gCalTCor(74:5:end-5),'+k','markersize',18);
plot(flipInfoAll.t(72:5:end-5),flipInfoAll.gCalTCor(72:5:end-5),'+b','markersize',18);
plot(flipInfoAll.t(70:5:end-5),flipInfoAll.gCalTCor(70:5:end-5),'xr','markersize',18);
plot(flipInfoAll.t(71:5:end-5),flipInfoAll.gCalTCor(71:5:end-5),'xb','markersize',18);
plot(flipInfoAll.t(73:5:end-5),flipInfoAll.gCalTCor(73:5:end-5),'xk','markersize',18);
xl = xlim; yl = ylim;
plot([0 0]+xl(1)+diff(xl)/10,mean(yl)+[0 0.0001],'-k')
text(xl(1)+diff(xl)/9,mean(yl)+0.00005,'10 \mug','fontsize',14)
legend('1st X','Y','2nd X','-X','-Y','location','best')
datetick
title({'Axial SCTA Calibrations August 2018 - November 2019'})
ylabel('Calibration, m/s^2')
set(gca,'fontsize',18)
xtickangle(45)
box on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../calibrations/Axial/process_Axial_AGU','-dtiff')