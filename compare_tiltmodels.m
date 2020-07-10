% compare_tiltmodels.m
%
% Compare tidal tilts at Pinon Flat from Mark Zumberge against tidal tilts
% predicted by spotl to ensure that spotl produces expected result.
%

clear; close all

% load predicted tilts
fid=fopen('../tidal_comp/PinonFlat/ew.tilt.tides.ascii','r');
PF.e=fscanf(fid,'%f'); % [nrad]
fclose(fid);
fid=fopen('../tidal_comp/PinonFlat/ns.tilt.tides.ascii','r');
PF.n=fscanf(fid,'%f'); % [nrad]
fclose(fid);
PF.t=linspace(datenum(2018,10,09),datenum(2020,03,05,23,55,00),length(PF.e))';

% calculate spotl tilts
path(path,genpath('/Users/erikfred/Documents/spotl/EKF_spotl/'))
[spotl.n,spotl.e,spotl.t,~]=tide_tilt(33.610199,-116.455498,...
    -datenum(2018,10,09),-datenum(2020,03,05),3600,1290,'osu','continent','datenum');

% PF SCTA data
load('../calibrations/PinonFlat/PFstitch_min_temp.mat','stitch_min')
[~,i1]=min(abs(stitch_min.t-datenum(2019,07,08)));
[~,i2]=min(abs(stitch_min.t-datenum(2019,08,04)));
SCTA.t=stitch_min.t(i1:i2);
% note orientation such that y:west and x:north
SCTA.n=-1*fillmissing(stitch_min.LAX(i1:i2),'linear')*10^3; % [nrad]
SCTA.e=fillmissing(stitch_min.LAY(i1:i2),'linear')*10^3; % [nrad]
[b,a]=butter(4,2*(1/60/60/48)/(1/60),'high');
SCTA.e=filtfilt(b,a,SCTA.e);
SCTA.n=filtfilt(b,a,SCTA.n);

% plot and compare
figure
subplot(211)
hold on
plot(PF.t,PF.e,'linewidth',1)
plot(spotl.t,spotl.e,'linewidth',1)
plot(SCTA.t,SCTA.e)
xlim([datenum(2019,07,08) datenum(2019,08,04)])
datetick('x','keeplimits')
ylabel('E tilt (nrad)')
legend('Zumberge','SPOTL','SCTA')
set(gca,'fontsize',12)
box on
subplot(212)
hold on
plot(PF.t,PF.n,'linewidth',1)
plot(spotl.t,spotl.n,'linewidth',1)
plot(SCTA.t,SCTA.n)
xlim([datenum(2019,07,08) datenum(2019,08,04)])
datetick('x','keeplimits')
ylabel('N tilt (nrad)')
legend('Zumberge','SPOTL','SCTA')
set(gca,'fontsize',12)
box on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../tidal_comp/PinonFlat/tiltcompare','-dtiff','-r300')