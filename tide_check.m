% tide_check.m
%
% Compares tidal signals between Axial SCTA and BOPT sensor at AXCC1. Also
% generates plots for the other BOPT locations. Will eventually include
% comparison to a tidal model. Should allow us to resolve
% orientation/polarity ambiguity between instruments.
%

clear; close all

%% CONFIG

t0=datenum(2018,11,01);
tf=datenum(2019,12,01);

t0_str=datestr(t0,30);
tf_str=datestr(tf,30);
svstr=[t0_str(1:8) '-' tf_str(1:8)];

%load data
load('../calibrations/Axial/axialstitch_min_temp.mat','stitch_min')
load('../tiltcompare/SCTA_Lily_comp/AXCC1data_min.mat','lily_min')

load('../compass_directions.mat')
CCMP=CCMP([4 3 2 1 5]);

bopt={'AXCC1','AXEC2','AXID1'};

% t1=t0;
% AXEC2.time=[];AXEC2.LAX=[];AXEC2.LAY=[];AXEC2.BDO=[];
% AXID1.time=[];AXID1.LAX=[];AXID1.LAY=[];AXID1.BDO=[];
% ASHES.time=[];ASHES.LAX=[];ASHES.LAY=[];ASHES.BDO=[];
%% PULL DATA FOR OTHER AXIAL INSTRUMENTS
% while t1<=tf
%     t1_s=datestr(t1,31); t1_s=t1_s(1:10);
%     
%     % grab all BOPT data
%     for i=1:length(sta2)
%         for j=1:length(cha2)
%             if ~exist(['../tiltcompare/' sta2{i} '/' sta2{i} '_' cha2{j} '_' datestr(t1,29) '.miniseed'],'file')
%                 IRIS_data_pull(sta2{i},cha2{j},'11',t1,t1+1);
%             end
%             temp=rdmseed(['../tiltcompare/' sta2{i} '/' sta2{i} '_' cha2{j} '_' datestr(t1,29) '.miniseed']);
%             ttemp=cat(1,temp.t);
%             dtemp=double(cat(1,temp.d));
%             % decimate to 1 sample/min
%             dtemp1=decimate(dtemp,6,'fir');
%             data.a(:,j)=decimate(dtemp1,10,'fir');
%             data.t=linspace(ttemp(1),ttemp(1)+1439/1440,length(data.a))';
%         end
%         % add pressure
%         if ~exist(['../tiltcompare/' sta2{i} '/' sta2{i} '_BDO_' datestr(t1,29) '.miniseed'],'file')
%             IRIS_data_pull(sta2{i},'BDO','11',t1,t1+1);
%         end
%         temp=rdmseed(['../tiltcompare/' sta2{i} '/' sta2{i} '_BDO_' datestr(t1,29) '.miniseed']);
%         dtemp=double(cat(1,temp.d));
%         % decimate to 1 sample/min
%         dtemp1=decimate(dtemp,10,'fir');
%         dtemp2=decimate(dtemp1,12,'fir');
%         data.p=decimate(dtemp2,10,'fir');
%         % append to appropriate structure
%         eval([sta2{i} '.time=[' sta2{i} '.time; data.t];']);
%         eval([sta2{i} '.LAX=[' sta2{i} '.LAX; data.a(:,1)];']);
%         eval([sta2{i} '.LAY=[' sta2{i} '.LAY; data.a(:,2)];']);
%         eval([sta2{i} '.BDO=[' sta2{i} '.BDO; data.p];']);
%         clear data
%     end
% end

%% PROCESSING
% trim data
SCTA.t=stitch_min.t(stitch_min.t>t0 & stitch_min.t<tf);
SCTA.LAX=stitch_min.LAX(stitch_min.t>t0 & stitch_min.t<tf);
SCTA.LAY=stitch_min.LAY(stitch_min.t>t0 & stitch_min.t<tf);

AXCC1.t=lily_min.t(lily_min.t>t0 & lily_min.t<tf);
AXCC1.LAX=lily_min.LAX(lily_min.t>t0 & lily_min.t<tf);
AXCC1.LAY=lily_min.LAY(lily_min.t>t0 & lily_min.t<tf);
AXCC1.BDO=lily_min.BDO(lily_min.t>t0 & lily_min.t<tf);

% BOPT compass correction
for l=1%:length(bopt)
    crd_rot=[cosd(CCMP(l).plus_y) sind(CCMP(l).plus_y); -sind(CCMP(l).plus_y) cosd(CCMP(l).plus_y)];
    eval(['temp=crd_rot*[' bopt{l} '.LAX'';' bopt{l} '.LAY''];']);
    eval([bopt{l} '.LAX=temp(1,:)'';']);
    eval([bopt{l} '.LAY=temp(2,:)'';']);
end

% SCTA orientation correction
crd_rot=[cosd(CCMP(5).plus_y) sind(CCMP(5).plus_y); -sind(CCMP(5).plus_y) cosd(CCMP(5).plus_y)];
temp=crd_rot*[SCTA.LAX';SCTA.LAY'];
SCTA.LAX=temp(1,:)';
SCTA.LAY=temp(2,:)';

% downsample to hourly, ensure even sampling
[~,SCTA.LAX,~,~]=downsample_uneven5(SCTA.t,SCTA.LAX,1/24);
[SCTA.t,SCTA.LAY,~,~]=downsample_uneven5(SCTA.t,SCTA.LAY,1/24);

[~,AXCC1.BDO,~,~]=downsample_uneven5(AXCC1.t,AXCC1.BDO,1/24);
[~,AXCC1.LAX,~,~]=downsample_uneven5(AXCC1.t,AXCC1.LAX,1/24);
[AXCC1.t,AXCC1.LAY,~,~]=downsample_uneven5(AXCC1.t,AXCC1.LAY,1/24);

% interpolate NaNs
SCTA.LAX=fillgaps(SCTA.LAX,24*7);
SCTA.LAY=fillgaps(SCTA.LAY,24*7);

AXCC1.LAX=fillgaps(AXCC1.LAX,24*7);
AXCC1.LAY=fillgaps(AXCC1.LAY,24*7);

% bandpass 1-2 sample/day
[b,a]=butter(4,[2*(1/60/60/12.5)/(1/3600) 2*(1/60/60/12.3711)/(1/3600)]);
SCTA.LAX=filtfilt(b,a,SCTA.LAX);
SCTA.LAY=filtfilt(b,a,SCTA.LAY);
AXCC1.LAX=filtfilt(b,a,AXCC1.LAX);
AXCC1.LAY=filtfilt(b,a,AXCC1.LAY);
% AXEC2.LAX=filtfilt(b,a,AXEC2.LAX);
% AXEC2.LAY=filtfilt(b,a,AXEC2.LAY);
% AXID1.LAX=filtfilt(b,a,AXID1.LAX);
% AXID1.LAY=filtfilt(b,a,AXID1.LAY);

%% TIDAL TILT MODEL
path(path,genpath('/Users/erikfred/Documents/spotl/EKF_spotl/'))

[AXCCsptl.LAY,AXCCsptl.LAX,AXCCsptl.t,~]=tide_tilt(45.955002,-130.008743,...
    -t0,-tf,3600,-1529.9,'osu','ocean','datenum');
AXCCsptl.LAY=AXCCsptl.LAY/10^3; % convert from nrad to urad
AXCCsptl.LAX=AXCCsptl.LAX/10^3;

[AXECsptl.LAY,AXECsptl.LAX,AXECsptl.t,~]=tide_tilt(45.939671,-129.973801,...
    -t0,-tf,3600,-1519.0,'osu','ocean','datenum');
AXECsptl.LAY=AXECsptl.LAY/10^3; % convert from nrad to urad
AXECsptl.LAX=AXECsptl.LAX/10^3;

[AXIDsptl.LAY,AXIDsptl.LAX,AXIDsptl.t,~]=tide_tilt(45.925732,-129.977997,...
    -t0,-tf,3600,-1527.5,'osu','ocean','datenum');
AXIDsptl.LAY=AXIDsptl.LAY/10^3; % convert from nrad to urad
AXIDsptl.LAX=AXIDsptl.LAX/10^3;

% bandpass 1-2 sample/day
[b,a]=butter(4,[2*(1/60/60/12.5)/(1/3600) 2*(1/60/60/12.3711)/(1/3600)]);
AXCCsptl.LAX=filtfilt(b,a,AXCCsptl.LAX);
AXCCsptl.LAY=filtfilt(b,a,AXCCsptl.LAY);
AXECsptl.LAX=filtfilt(b,a,AXECsptl.LAX);
AXECsptl.LAY=filtfilt(b,a,AXECsptl.LAY);
AXIDsptl.LAX=filtfilt(b,a,AXIDsptl.LAX);
AXIDsptl.LAY=filtfilt(b,a,AXIDsptl.LAY);

%% TIME SERIES PLOTTING

figure
% Central Caldera
subplot(211)
plot(AXCC1.t,detrend(AXCC1.LAX-mean(AXCC1.LAX)),'linewidth',1)
hold on
plot(SCTA.t,detrend(SCTA.LAX-mean(SCTA.LAX)),'linewidth',1)
plot(AXCCsptl.t,AXCCsptl.LAX,'k','linewidth',1)
xlim([0 30]+round(mean(t0,tf)))
datetick('x','keeplimits')
set(gca,'xticklabel',[])
ylabel('x-tilt (\murad)')
title('Central Caldera')
legend('BOPT','SCTA','SPOTL','location','northeast')
subplot(212)
plot(AXCC1.t,detrend(AXCC1.LAY-mean(AXCC1.LAY)),'linewidth',1)
hold on
plot(SCTA.t,detrend(SCTA.LAY-mean(SCTA.LAY)),'linewidth',1)
plot(AXCCsptl.t,AXCCsptl.LAY,'k','linewidth',1)
xlim([0 30]+round(mean(t0,tf)))
datetick('x','keeplimits')
ylabel('y-tilt (\murad)')
legend('BOPT','SCTA','SPOTL','location','northeast')

% % Eastern Caldera
% subplot(422)
% plot(AXEC2.t,detrend(AXEC2.LAX-mean(AXEC2.LAX)),'linewidth',1)
% hold on
% plot(AXECsptl.t,AXECsptl.LAX*4,'k','linewidth',1)
% datetick('x','keeplimits')
% set(gca,'xticklabel',[])
% ylabel('x-tilt (\murad)')
% title('Eastern Caldera')
% legend('BOPT','SPOTL','location','northeast')
% subplot(424)
% plot(AXEC2.t,detrend(AXEC2.LAY-mean(AXEC2.LAY)),'linewidth',1)
% hold on
% plot(AXECsptl.t,AXECsptl.LAY*4,'k','linewidth',1)
% datetick('x','keeplimits')
% ylabel('y-tilt (\murad)')
% legend('BOPT','SPOTL','location','northeast')
% 
% % ASHES
% subplot(425)
% datetick('x','keeplimits')
% set(gca,'xticklabel',[])
% ylabel('x-tilt (\murad)')
% title('Ashes Vent Field')
% subplot(427)
% datetick('x','keeplimits')
% ylabel('y-tilt (\murad)')
% 
% % International District
% subplot(426)
% plot(AXID1.t,detrend(AXID1.LAX-mean(AXID1.LAX)),'linewidth',1)
% hold on
% plot(AXIDsptl.t,AXIDsptl.LAX*4,'k','linewidth',1)
% datetick('x','keeplimits')
% set(gca,'xticklabel',[])
% ylabel('x-tilt (\murad)')
% title('International District')
% legend('BOPT','SPOTL','location','northeast')
% subplot(428)
% plot(AXID1.t,detrend(AXID1.LAY-mean(AXID1.LAY)),'linewidth',1)
% hold on
% plot(AXIDsptl.t,AXIDsptl.LAY*4,'k','linewidth',1)
% datetick('x','keeplimits')
% ylabel('y-tilt (\murad)')
% legend('BOPT','SPOTL','location','northeast')
    
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print(['../tidal_comp/Axial/tilt_' svstr],'-dtiff','-r300')

% compare AXCC1 pressure to tidal amplitude
[h,th]=oceantide(45.955002,-130.008743,-t0,-tf,3600,'osu','datenum');

figure
hold on
yyaxis left
plot(AXCC1.t,(AXCC1.BDO-nanmean(AXCC1.BDO))/10065)
ylabel('bopt pressure (m)')
yyaxis right
plot(th,h)
xlim([0 30]+round(mean(t0,tf)))
datetick('x','keeplimits')
ylabel('spotl tidal amplitude (m)')
title('central caldera')
set(gca,'fontsize',14)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print(['../tidal_comp/Axial/AXCC1_pressureVspotl_' svstr],'-dtiff','-r300')

%% Spectral Plots
% power spectra
[SCTA.pxx,SCTA.fx]=pwelch(SCTA.LAX-mean(SCTA.LAX),24*30,12*10,2^16,1/3600);
[SCTA.pyy,SCTA.fy]=pwelch(SCTA.LAY-mean(SCTA.LAY),24*30,12*10,2^16,1/3600);

[AXCC1.pxx,AXCC1.fx]=pwelch(AXCC1.LAX-mean(AXCC1.LAX),24*30,12*10,2^16,1/3600);
[AXCC1.pyy,AXCC1.fy]=pwelch(AXCC1.LAY-mean(AXCC1.LAY),24*30,12*10,2^16,1/3600);

[AXCCsptl.pxx,AXCCsptl.fx]=pwelch(AXCCsptl.LAX-mean(AXCCsptl.LAX),24*30,12*10,2^16,1/3600);
[AXCCsptl.pyy,AXCCsptl.fy]=pwelch(AXCCsptl.LAY-mean(AXCCsptl.LAY),24*30,12*10,2^16,1/3600);

figure
subplot(211)
loglog(AXCC1.fx*60*60*24,AXCC1.pxx,'linewidth',1)
hold on
loglog(SCTA.fx*60*60*24,SCTA.pxx,'linewidth',1)
loglog(AXCCsptl.fx*60*60*24,AXCCsptl.pxx,'k','linewidth',1)
ylabel('PSD (x)')
xlabel('Frequency (Hz)')
legend('BOPT','SCTA','SPOTL','location','northeast')
set(gca,'fontsize',14)
subplot(212)
loglog(AXCC1.fy*60*60*24,AXCC1.pyy,'linewidth',1)
hold on
loglog(SCTA.fy*60*60*24,SCTA.pyy,'linewidth',1)
loglog(AXCCsptl.fy*60*60*24,AXCCsptl.pyy,'k','linewidth',1)
ylabel('PSD (y)')
xlabel('Frequency (Hz)')
legend('BOPT','SCTA','SPOTL','location','northeast')
set(gca,'fontsize',14)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print(['../tidal_comp/Axial/PSD_' svstr],'-dtiff','-r300')

% phase amplitude and angle
[SCTA.Fv,SCTA.amp_ftsx,SCTA.phs_ftsx]=amplitudespectrum(SCTA.t,SCTA.LAX-mean(SCTA.LAX));
[SCTA.Fv,SCTA.amp_ftsy,SCTA.phs_ftsy]=amplitudespectrum(SCTA.t,SCTA.LAY-mean(SCTA.LAY));

[AXCC1.Fv,AXCC1.amp_ftsx,AXCC1.phs_ftsx]=amplitudespectrum(AXCC1.t,AXCC1.LAX-mean(AXCC1.LAX));
[AXCC1.Fv,AXCC1.amp_ftsy,AXCC1.phs_ftsy]=amplitudespectrum(AXCC1.t,AXCC1.LAY-mean(AXCC1.LAY));

[AXCCsptl.Fv,AXCCsptl.amp_ftsx,AXCCsptl.phs_ftsx]=amplitudespectrum(AXCCsptl.t,AXCCsptl.LAX-mean(AXCCsptl.LAX));
[AXCCsptl.Fv,AXCCsptl.amp_ftsy,AXCCsptl.phs_ftsy]=amplitudespectrum(AXCCsptl.t,AXCCsptl.LAY-mean(AXCCsptl.LAY));

figure
subplot(211)
loglog(AXCC1.Fv,AXCC1.amp_ftsx,'linewidth',1)
hold on
loglog(SCTA.Fv,SCTA.amp_ftsx,'linewidth',1)
loglog(AXCCsptl.Fv,AXCCsptl.amp_ftsx,'k','linewidth',1)
ylabel('X FFT')
xlabel('Frequency (1/day)')
legend('BOPT','SCTA','SPOTL','location','northeast')
set(gca,'fontsize',14)
subplot(212)
loglog(AXCC1.Fv,AXCC1.amp_ftsy,'linewidth',1)
hold on
loglog(SCTA.Fv,SCTA.amp_ftsy,'linewidth',1)
loglog(AXCCsptl.Fv,AXCCsptl.amp_ftsy,'k','linewidth',1)
ylabel('Y FFT')
xlabel('Frequency (1/day)')
legend('BOPT','SCTA','SPOTL','location','northeast')
set(gca,'fontsize',14)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print(['../tidal_comp/Axial/fft_' svstr],'-dtiff','-r300')

figure
subplot(211)
semilogx(AXCC1.Fv,AXCC1.phs_ftsx,'linewidth',1)
hold on
semilogx(SCTA.Fv,SCTA.phs_ftsx,'linewidth',1)
semilogx(AXCCsptl.Fv,AXCCsptl.phs_ftsx,'k','linewidth',1)
ylabel('X Phase Angle (rad)')
xlabel('Frequency (1/day)')
legend('BOPT','SCTA','SPOTL','location','northeast')
set(gca,'fontsize',14)
subplot(212)
semilogx(AXCC1.Fv,AXCC1.phs_ftsy,'linewidth',1)
hold on
semilogx(SCTA.Fv,SCTA.phs_ftsy,'linewidth',1)
semilogx(AXCCsptl.Fv,AXCCsptl.phs_ftsy,'k','linewidth',1)
ylabel('Y Phase Angle (rad)')
xlabel('Frequency (1/day)')
legend('BOPT','SCTA','SPOTL','location','northeast')
set(gca,'fontsize',14)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print(['../tidal_comp/Axial/phase_' svstr],'-dtiff','-r300')

%% Cross Correlations
% resample model to 1 sample/minute
xtemp1=interp1(AXCCsptl.t,AXCCsptl.LAX,AXCC1.t);
ytemp1=interp1(AXCCsptl.t,AXCCsptl.LAY,AXCC1.t);
[cross1x,lag1x]=xcorr(AXCC1.LAX,xtemp1,6,'normalized');
[cross1y,lag1y]=xcorr(AXCC1.LAY,ytemp1,6,'normalized');

xtemp2=interp1(AXCCsptl.t,AXCCsptl.LAX,SCTA.t);
ytemp2=interp1(AXCCsptl.t,AXCCsptl.LAY,SCTA.t);
[cross2x,lag2x]=xcorr(-SCTA.LAX,xtemp2,6,'normalized');
[cross2y,lag2y]=xcorr(-SCTA.LAY,ytemp2,6,'normalized');

figure
subplot(211)
plot(lag1x,cross1x,'linewidth',1)
hold on
plot(lag2x,cross2x,'linewidth',1)
text(0,-0.7,[num2str(round(cross1x(7),2)) ' @ 0 lag'],'color','b','fontsize',12)
text(0,0.7,[num2str(round(cross2x(7),2)) ' @ 0 lag'],'color','r','fontsize',12)
ylim([-1 1])
ylabel('X')
set(gca,'xticklabel',[])
title('Cross Correlation Amplitudes')
legend('BOPT','SCTA','location','northeast')
set(gca,'fontsize',14)
subplot(212)
plot(lag1y,cross1y,'linewidth',1)
hold on
plot(lag2y,cross2y,'linewidth',1)
text(0,-0.5,[num2str(round(cross1y(7),2)) ' @ 0 lag'],'color','b','fontsize',12)
text(0,0.5,[num2str(round(cross2y(7),2)) ' @ 0 lag'],'color','r','fontsize',12)
ylim([-1 1])
ylabel('Y')
xlabel('Lags (hours)')
legend('BOPT','SCTA','location','northeast')
set(gca,'fontsize',14)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print(['../tidal_comp/Axial/xcorr_' svstr],'-dtiff','-r300')