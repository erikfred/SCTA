% daily_compplots.m
%
% loads AXCC1 and SCTA data day by day and plots them on the same axis for
% comparison.
%

%%%%%%%%%%%%CONFIG%%%%%%%%%%%%%%
cha1='HNZ';
cha2='BNZ';

fac1=314252;
fac2=10^7;
%%%%%%%%END CONFIG%%%%%%%%%%%%%%

for i=1:30
    t0=datenum(2019,4,1)+i-1;
    t0_s=datestr(t0,31); t0_s=t0_s(1:10);
    
    AXCC1=rdmseed(['../tiltcompare/AXCC1/AXCC1_' cha1 '_' t0_s '.miniseed']);
    t1=cat(1,AXCC1.t);
    z1=cat(1,AXCC1.d);
    
    AXCC2=rdmseed(['../tiltcompare/AXCC2/AXCC2_' cha2 '_' t0_s '.miniseed']);
    t2=cat(1,AXCC2.t);
    z2=cat(1,AXCC2.d);
    
    figure(77)
    clf
    plot(t1,(z1-mean(z1))/fac1,'linewidth',1)
    hold on
    plot(t2,(z2-mean(z2))/fac2,'linewidth',1)
    set(gca,'fontsize',18)
    ylabel('Acceleration (m/s^2)')
    datetick('x')
    title(t0_s)
    
    orient landscape
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print(['../tiltcompare/dailycomps/' t0_s],'-dtiff')
end