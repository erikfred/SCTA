% demonstrate_data.m
%
% For each of Axial and PF, generates example time series plots of 1) a full
% calibration sequence and 2) a non-calibration day.
%

%% Pinon Flat

load('../../calibrations/PinonFlat/PFdata.mat','flipInfoAll')

% calibration examples
dayPF=floor(flipInfoAll.t(91));
dataPF=get_sctaDay('/Volumes/GoogleDrive/My Drive/Oceanography/SCTA-Share/OOI-PF/SCTA-PF/ParsedData',dayPF);

% full calibration sequence
figure(50); clf
plot(dataPF.t,dataPF.a(:,1),'linewidth',1)
hold on
plot(dataPF.t,dataPF.a(:,2),'linewidth',1)
xlim([dayPF+datenum(0000,00,00,21,00,00) dayPF+datenum(0000,00,00,21,12,00)])
ylim([-10 10])
datetick('x','keeplimits')
xtickangle(45)
legend('X','Y')
ylabel('Acceleration (m/s^2)')
title(['Example Calibration Sequence at Pinon Flat on ' datestr(dayPF)])
set(gca,'fontsize',14)
box on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../../paper_figures/cal_sequence_PF','-dtiff','-r300')

% one-orientation interval
figure(51); clf
plot(dataPF.t,dataPF.as,'linewidth',1)
xlim([dayPF+datenum(0000,00,00,21,00,50) dayPF+datenum(0000,00,00,21,02,35)])
datetick('x','keeplimits')
xtickangle(45)
lim_y=ylim;
hold on
hp=patch([dayPF+datenum(0000,00,00,21,01,50) dayPF+datenum(0000,00,00,21,02,20) ...
    dayPF+datenum(0000,00,00,21,02,20) dayPF+datenum(0000,00,00,21,01,50)],...
    [lim_y(1) lim_y(1) lim_y(2) lim_y(2)],[0.8 0.8 0.8]);
hp.EdgeColor='none';
hp.FaceVertexAlphaData=0.2;
hp.FaceAlpha='flat';
ylim(lim_y)
ylabel('Acceleration (m/s^2)')
title(['Example Calibration +X calibration at Pinon Flat on ' datestr(dayPF)])
set(gca,'fontsize',14)
box on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../../paper_figures/x_cal_PF','-dtiff','-r300')

% non-calibration example
dayPF=dayPF-1;
dataPF=get_sctaDay('/Volumes/GoogleDrive/My Drive/Oceanography/SCTA-Share/OOI-PF/SCTA-PF/ParsedData',dayPF);

dataPF.x=movmean(dataPF.a(:,1),40*60);
dataPF.y=movmean(dataPF.a(:,2),40*60);

% time series of day of data
figure(52); clf
plot(dataPF.t,dataPF.x,'linewidth',1)
hold on
plot(dataPF.t,dataPF.y+0.03349,'linewidth',1)
xlim([dayPF dayPF+1])
datetick('x','keeplimits')
xtickangle(45)
set(gca,'ytick',-0.0047:0.000002:-0.004686)
set(gca,'yticklabel',0:0.2:1.4)
legend('X','Y')
ylabel('Acceleration (\mug)')
title(['Example Day at Pinon Flat on ' datestr(dayPF)])
set(gca,'fontsize',14)
box on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../../paper_figures/example_day_time_PF','-dtiff','-r300')

% spectra of data day
[dataPF.pxx,dataPF.fx]=pwelch(dataPF.a(:,1)-mean(dataPF.a(:,1)),40*60*60,40*60*60*0.9,2^18,40);
[dataPF.pyy,dataPF.fy]=pwelch(dataPF.a(:,2)-mean(dataPF.a(:,2)),40*60*60,40*60*60*0.9,2^18,40);

figure(53); clf
loglog(dataPF.fx,dataPF.pxx,'linewidth',1)
hold on
loglog(dataPF.fy,dataPF.pyy,'linewidth',1)
legend('X','Y')
ylabel('Amplitude')
xlabel('Frequency (Hz)')
title(['Example Day at Pinon Flat on ' datestr(dayPF)])
set(gca,'fontsize',14)
box on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../../paper_figures/example_day_spectral_PF','-dtiff','-r300')

%% Axial

load('../../calibrations/Axial/axialdata.mat','flipInfoAll')

% calibration examples
dayA=datenum(2020,10,14);
dataA=get_sctaDay('/Volumes/GoogleDrive/My Drive/Oceanography/SCTA-Share/OOI-SCTA/ParsedData',dayA);

% full calibration sequence
figure(60); clf
plot(dataA.t,dataA.a(:,1),'linewidth',1)
hold on
plot(dataA.t,dataA.a(:,2),'linewidth',1)
xlim([dayA+datenum(0000,00,00,20,59,00) dayA+datenum(0000,00,00,21,11,00)])
ylim([-10 10])
datetick('x','keeplimits')
xtickangle(45)
legend('X','Y')
ylabel('Acceleration (m/s^2)')
title(['Example Calibration Sequence at Axial on ' datestr(dayA)])
set(gca,'fontsize',14)
box on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../../paper_figures/cal_sequence_Axial','-dtiff','-r300')

% one-orientation interval
figure(61); clf
plot(dataA.t,dataA.as,'linewidth',1)
xlim([dayA+datenum(0000,00,00,20,59,55) dayA+datenum(0000,00,00,21,01,40)])
datetick('x','keeplimits')
xtickangle(45)
lim_y=ylim;
hold on
hp=patch([dayA+datenum(0000,00,00,21,00,55) dayA+datenum(0000,00,00,21,01,25) ...
    dayA+datenum(0000,00,00,21,01,25) dayA+datenum(0000,00,00,21,00,55)],...
    [lim_y(1) lim_y(1) lim_y(2) lim_y(2)],[0.8 0.8 0.8]);
hp.EdgeColor='none';
hp.FaceVertexAlphaData=0.2;
hp.FaceAlpha='flat';
ylim(lim_y)
ylabel('Acceleration (m/s^2)')
title(['Example Calibration +X calibration at Axial on ' datestr(dayA)])
set(gca,'fontsize',14)
box on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../../paper_figures/x_cal_Axial','-dtiff','-r300')

% non-calibration example
dayA=dayA-1;
dataA=get_sctaDay('/Volumes/GoogleDrive/My Drive/Oceanography/SCTA-Share/OOI-SCTA/ParsedData',dayA);

dataA.x=movmean(dataA.a(:,1),40*60);
dataA.y=movmean(dataA.a(:,2),40*60);

% time series of day of data
figure(62); clf
plot(dataA.t,dataA.x,'linewidth',1)
hold on
plot(dataA.t,dataA.y+0.17375,'linewidth',1)
xlim([dayA dayA+1])
datetick('x','keeplimits')
xtickangle(45)
set(gca,'ytick',0.0307:0.00001:0.03079)
set(gca,'yticklabel',0:1:9)
legend('X','Y')
ylabel('Acceleration (\mug)')
title(['Example Day at Axial on ' datestr(dayA)])
set(gca,'fontsize',14)
box on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../../paper_figures/example_day_time_Axial','-dtiff','-r300')

% spectra of data day
[dataA.pxx,dataA.fx]=pwelch(dataA.a(:,1)-mean(dataA.a(:,1)),40*60*60,40*60*60*0.9,2^18,40);
[dataA.pyy,dataA.fy]=pwelch(dataA.a(:,2)-mean(dataA.a(:,2)),40*60*60,40*60*60*0.9,2^18,40);

figure(63); clf
loglog(dataA.fx,dataA.pxx,'linewidth',1)
hold on
loglog(dataA.fy,dataA.pyy,'linewidth',1)
legend('X','Y')
ylabel('Amplitude')
xlabel('Frequency (Hz)')
title(['Example Day at Axial on ' datestr(dayA)])
set(gca,'fontsize',14)
box on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../../paper_figures/example_day_spectral_Axial','-dtiff','-r300')