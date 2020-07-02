% tide_check.m
%
% Compares tidal signals between Axial SCTA and BOPT sensor at AXCC1. Also
% generates plots for the other BOPT locations. Will eventually include
% comparison to a tidal model. Should allow us to resolve
% orientation/polarity ambiguity between instruments.
%

clear; close all

%% CONFIG
load('../compass_directions.mat')
CCMP=CCMP([4 3 2 1]);

t0=datenum(2019,01,15);
tf=datenum(2019,01,21);

t0_str=datestr(t0,30);
tf_str=datestr(tf,30);
svstr=[t0_str(1:8) '-' tf_str(7:8)];

%load data
sta1={'AXCC2'};
cha1={'MNE','MNN'};
sta2={'AXCC1','AXEC2','AXID1'};
cha2={'LAX','LAY'};

t1=t0;
AXCC1.time=[];AXCC1.LAX=[];AXCC1.LAY=[];AXCC1.BDO=[];
AXEC2.time=[];AXEC2.LAX=[];AXEC2.LAY=[];AXEC2.BDO=[];
AXID1.time=[];AXID1.LAX=[];AXID1.LAY=[];AXID1.BDO=[];
ASHES.time=[];ASHES.LAX=[];ASHES.LAY=[];ASHES.BDO=[];
AXCC2.time=[];AXCC2.MNE=[];AXCC2.MNN=[];
%% PULL DATA
while t1<=tf
    t1_s=datestr(t1,31); t1_s=t1_s(1:10);
    
    % grab all BOPT data
    for i=1:length(sta2)
        for j=1:length(cha2)
            if ~exist(['../tiltcompare/' sta2{i} '/' sta2{i} '_' cha2{j} '_' datestr(t1,29) '.miniseed'],'file')
                IRIS_data_pull(sta2{i},cha2{j},'11',t1,t1+1);
            end
            temp=rdmseed(['../tiltcompare/' sta2{i} '/' sta2{i} '_' cha2{j} '_' datestr(t1,29) '.miniseed']);
            ttemp=cat(1,temp.t);
            dtemp=double(cat(1,temp.d));
            % decimate to 1 sample/min
            dtemp1=decimate(dtemp,6,'fir');
            data.a(:,j)=decimate(dtemp1,10,'fir');
            data.t=linspace(ttemp(1),ttemp(1)+1439/1440,length(data.a))';
        end
        % add pressure
        if ~exist(['../tiltcompare/' sta2{i} '/' sta2{i} '_BDO_' datestr(t1,29) '.miniseed'],'file')
            IRIS_data_pull(sta2{i},'BDO','11',t1,t1+1);
        end
        temp=rdmseed(['../tiltcompare/' sta2{i} '/' sta2{i} '_BDO_' datestr(t1,29) '.miniseed']);
        dtemp=double(cat(1,temp.d));
        % decimate to 1 sample/min
        dtemp1=decimate(dtemp,10,'fir');
        dtemp2=decimate(dtemp1,12,'fir');
        data.p=decimate(dtemp2,10,'fir');
        % append to appropriate structure
        eval([sta2{i} '.time=[' sta2{i} '.time; data.t];']);
        eval([sta2{i} '.LAX=[' sta2{i} '.LAX; data.a(:,1)];']);
        eval([sta2{i} '.LAY=[' sta2{i} '.LAY; data.a(:,2)];']);
        eval([sta2{i} '.BDO=[' sta2{i} '.BDO; data.p];']);
        clear data
    end
    
    % grab AXCC2 (SCTA) data
    for k=1:length(cha1)
        if ~exist(['../tiltcompare/AXCC2/AXCC2_' cha1{k} '_' datestr(t1,29) '.miniseed'],'file')
            IRIS_data_pull('AXCC2',cha1{k},'--',t1,t1+1);
        end
        temp=rdmseed(['../tiltcompare/AXCC2/AXCC2_' cha1{k} '_' datestr(t1,29) '.miniseed']);
        ttemp=cat(1,temp.t);
        dtemp=double(cat(1,temp.d))/10^7;
        % decimate to 1 sample/min
        dtemp1=decimate(dtemp,8,'fir');
        dtemp2=decimate(dtemp1,6,'fir');
        data.a(:,k)=decimate(dtemp2,10,'fir');
        data.t=linspace(ttemp(1),ttemp(1)+1439/1440,length(data.a))';
    end
    % append to AXCC2
    AXCC2.time=[AXCC2.time; data.t];
    AXCC2.MNE=[AXCC2.MNE; data.a(:,1)];
    AXCC2.MNN=[AXCC2.MNN; data.a(:,2)];
    clear data
    
    t1=t1+1;
end

% check to ensure time series have same length
keyboard

%% PROCESSING
% BOPT compass correction
for l=1:length(sta2)
    crd_rot=[cosd(CCMP(l).plus_y) sind(CCMP(l).plus_y); -sind(CCMP(l).plus_y) cosd(CCMP(l).plus_y)];
    eval(['temp=crd_rot*[' sta2{l} '.LAX'';' sta2{l} '.LAY''];']);
    eval([sta2{l} '.LAX=temp(1,:)'';']);
    eval([sta2{l} '.LAY=temp(2,:)'';']);
end

% convert AXCC2 acceleration to tilt
AXCC2.LAX=asin(AXCC2.MNE/9.81)*10^6;
AXCC2.LAY=asin(AXCC2.MNN/9.81)*10^6;

% 2 sample/day low pass filter
[b,a]=butter(4,2*(1/60/60/12)/(1/60),'low');
AXCC2.LAX=filtfilt(b,a,AXCC2.LAX);
AXCC2.LAY=filtfilt(b,a,AXCC2.LAY);
AXCC1.LAX=filtfilt(b,a,AXCC1.LAX);
AXCC1.LAY=filtfilt(b,a,AXCC1.LAY);
AXEC2.LAX=filtfilt(b,a,AXEC2.LAX);
AXEC2.LAY=filtfilt(b,a,AXEC2.LAY);
AXID1.LAX=filtfilt(b,a,AXID1.LAX);
AXID1.LAY=filtfilt(b,a,AXID1.LAY);

% 1 sample/day high pass filter
[b,a]=butter(4,2*(1/60/60/24)/(1/60),'high');
AXCC2.LAX=filtfilt(b,a,AXCC2.LAX);
AXCC2.LAY=filtfilt(b,a,AXCC2.LAY);
AXCC1.LAX=filtfilt(b,a,AXCC1.LAX);
AXCC1.LAY=filtfilt(b,a,AXCC1.LAY);
AXEC2.LAX=filtfilt(b,a,AXEC2.LAX);
AXEC2.LAY=filtfilt(b,a,AXEC2.LAY);
AXID1.LAX=filtfilt(b,a,AXID1.LAX);
AXID1.LAY=filtfilt(b,a,AXID1.LAY);

%% TIDAL TILT MODEL
path(path,genpath('/Users/erikfred/Documents/spotl/EKF_spotl/'))

[AXCCsptl.LAY,AXCCsptl.LAX,AXCCsptl.time,~]=tide_tilt(45.955002,-130.008743,...
    -t0,-t1,3600,-1529.9,'osu','ocean','datenum');
AXCCsptl.LAY=AXCCsptl.LAY/10^3; % convert from nrad to urad
AXCCsptl.LAX=AXCCsptl.LAX/10^3;

[AXECsptl.LAY,AXECsptl.LAX,AXECsptl.time,~]=tide_tilt(45.939671,-129.973801,...
    -t0,-t1,3600,-1519.0,'osu','ocean','datenum');
AXECsptl.LAY=AXECsptl.LAY/10^3; % convert from nrad to urad
AXECsptl.LAX=AXECsptl.LAX/10^3;

[AXIDsptl.LAY,AXIDsptl.LAX,AXIDsptl.time,~]=tide_tilt(45.925732,-129.977997,...
    -t0,-t1,3600,-1527.5,'osu','ocean','datenum');
AXIDsptl.LAY=AXIDsptl.LAY/10^3; % convert from nrad to urad
AXIDsptl.LAX=AXIDsptl.LAX/10^3;

% 2 sample/day low pass filter
[b,a]=butter(4,2*(1/3600/12)/(1/3600),'low');
AXCCsptl.LAX=filtfilt(b,a,AXCCsptl.LAX);
AXCCsptl.LAY=filtfilt(b,a,AXCCsptl.LAY);
AXECsptl.LAX=filtfilt(b,a,AXECsptl.LAX);
AXECsptl.LAY=filtfilt(b,a,AXECsptl.LAY);
AXIDsptl.LAX=filtfilt(b,a,AXIDsptl.LAX);
AXIDsptl.LAY=filtfilt(b,a,AXIDsptl.LAY);

% 1 sample/day high pass filter
[b,a]=butter(4,2*(1/3600/24)/(1/3600),'high');
AXCCsptl.LAX=filtfilt(b,a,AXCCsptl.LAX);
AXCCsptl.LAY=filtfilt(b,a,AXCCsptl.LAY);
AXECsptl.LAX=filtfilt(b,a,AXECsptl.LAX);
AXECsptl.LAY=filtfilt(b,a,AXECsptl.LAY);
AXIDsptl.LAX=filtfilt(b,a,AXIDsptl.LAX);
AXIDsptl.LAY=filtfilt(b,a,AXIDsptl.LAY);

%% TIME SERIES PLOTTING

figure
% Central Caldera
subplot(421)
plot(AXCC1.time,detrend(AXCC1.LAX-mean(AXCC1.LAX)),'linewidth',1)
hold on
plot(AXCC2.time,-detrend(AXCC2.LAX-mean(AXCC2.LAX)),'linewidth',1)
plot(AXCCsptl.time,AXCCsptl.LAX*4,'k','linewidth',1)
datetick('x','keeplimits')
set(gca,'xticklabel',[])
ylabel('x-tilt (\murad)')
title('Central Caldera')
legend('BOPT','SCTA','SPOTL','location','northeast')
subplot(423)
plot(AXCC2.time,detrend(AXCC2.LAY-mean(AXCC2.LAY)),'linewidth',1)
hold on
plot(AXCC1.time,-detrend(AXCC1.LAY-mean(AXCC1.LAY)),'linewidth',1)
plot(AXCCsptl.time,AXCCsptl.LAY*4,'k','linewidth',1)
datetick('x','keeplimits')
ylabel('y-tilt (\murad)')
legend('BOPT','SCTA','SPOTL','location','northeast')

% Eastern Caldera
subplot(422)
plot(AXEC2.time,detrend(AXEC2.LAX-mean(AXEC2.LAX)),'linewidth',1)
hold on
plot(AXECsptl.time,AXECsptl.LAX*4,'k','linewidth',1)
datetick('x','keeplimits')
set(gca,'xticklabel',[])
ylabel('x-tilt (\murad)')
title('Eastern Caldera')
legend('BOPT','SPOTL','location','northeast')
subplot(424)
plot(AXEC2.time,detrend(AXEC2.LAY-mean(AXEC2.LAY)),'linewidth',1)
hold on
plot(AXECsptl.time,AXECsptl.LAY*4,'k','linewidth',1)
datetick('x','keeplimits')
ylabel('y-tilt (\murad)')
legend('BOPT','SPOTL','location','northeast')

% ASHES
subplot(425)
datetick('x','keeplimits')
set(gca,'xticklabel',[])
ylabel('x-tilt (\murad)')
title('Ashes Vent Field')
subplot(427)
datetick('x','keeplimits')
ylabel('y-tilt (\murad)')

% International District
subplot(426)
plot(AXID1.time,detrend(AXID1.LAX-mean(AXID1.LAX)),'linewidth',1)
hold on
plot(AXIDsptl.time,AXIDsptl.LAX*4,'k','linewidth',1)
datetick('x','keeplimits')
set(gca,'xticklabel',[])
ylabel('x-tilt (\murad)')
title('International District')
legend('BOPT','SPOTL','location','northeast')
subplot(428)
plot(AXID1.time,detrend(AXID1.LAY-mean(AXID1.LAY)),'linewidth',1)
hold on
plot(AXIDsptl.time,AXIDsptl.LAY*4,'k','linewidth',1)
datetick('x','keeplimits')
ylabel('y-tilt (\murad)')
legend('BOPT','SPOTL','location','northeast')
    
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print(['../tidal_comp/Axial/tilt_' svstr],'-dtiff','-r300')

% compare AXCC1 pressure to tidal amplitude
[h,th]=oceantide(45.955002,-130.008743,-t0,-t1,3600,'osu','datenum');

figure
hold on
yyaxis left
plot(AXCC1.time,(AXCC1.BDO-mean(AXCC1.BDO))/1.45038E-4/10065)
ylabel('bopt pressure (m)')
yyaxis right
plot(th,h)
datetick('x','keeplimits')
ylabel('spotl tidal amplitude (m)')
title('central caldera')
set(gca,'fontsize',14)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print(['../tidal_comp/Axial/AXCC1_pressureVspotl_' svstr],'-dtiff','-r300')

%% COMPARISON METRICS
% power spectra
[AXCC2.fx,AXCC2.ax]=amplitudespectrum(AXCC2.time,AXCC2.LAX-mean(AXCC2.LAX));
[AXCC2.fy,AXCC2.ay]=amplitudespectrum(AXCC2.time,AXCC2.LAY-mean(AXCC2.LAX));

[AXCC1.fx,AXCC1.ax]=amplitudespectrum(AXCC1.time,AXCC1.LAX-mean(AXCC1.LAX));
[AXCC1.fy,AXCC1.ay]=amplitudespectrum(AXCC1.time,AXCC1.LAY-mean(AXCC1.LAX));

[AXCCsptl.fx,AXCCsptl.ax]=amplitudespectrum(AXCCsptl.time,AXCCsptl.LAX-mean(AXCCsptl.LAX));
[AXCCsptl.fy,AXCCsptl.ay]=amplitudespectrum(AXCCsptl.time,AXCCsptl.LAY-mean(AXCCsptl.LAX));

figure
loglog(AXCC2.fx,AXCC2.ax,'linewidth',1)
hold on
loglog(AXCC1.fx,AXCC1.ax,'linewidth',1)
loglog(AXCCsptl.fx,AXCCsptl.ax,'linewidth',1)
ylabel('Amplitude')
xlabel('Frequency (1/day)')
legend('BOPT','SCTA','SPOTL','location','northeast')
set(gca,'fontsize',14)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print(['../tidal_comp/Axial/spectra_' svstr],'-dtiff','-r300')

% cross correlations
% resample model to 1 sample/minute
xtemp1=interp1(AXCCsptl.time,AXCCsptl.LAX,AXCC1.time);
ytemp1=interp1(AXCCsptl.time,AXCCsptl.LAY,AXCC1.time);
[cross1x,lag1x]=xcorr(AXCC1.LAX,xtemp1,360,'normalized');
[cross1y,lag1y]=xcorr(AXCC1.LAY,ytemp1,360,'normalized');

xtemp2=interp1(AXCCsptl.time,AXCCsptl.LAX,AXCC2.time);
ytemp2=interp1(AXCCsptl.time,AXCCsptl.LAY,AXCC2.time);
[cross2x,lag2x]=xcorr(-AXCC2.LAX,xtemp2,360,'normalized');
[cross2y,lag2y]=xcorr(-AXCC2.LAY,ytemp2,360,'normalized');

figure
subplot(211)
plot(lag1x/60,cross1x,'linewidth',1)
hold on
plot(lag2x/60,cross2x,'linewidth',1)
text(3,0.1,num2str(round(max(abs(cross1x)),2)),'color','b','fontsize',12)
text(3,-0.1,num2str(round(max(abs(cross2x)),2)),'color','r','fontsize',12)
ylabel('X')
set(gca,'xticklabel',[])
title('Cross Correlation Amplitudes')
legend('BOPT','SCTA','location','northeast')
set(gca,'fontsize',14)
subplot(212)
plot(lag1y/60,cross1y,'linewidth',1)
hold on
plot(lag2y/60,cross2y,'linewidth',1)
text(3,0.1,num2str(round(max(abs(cross1y)),2)),'color','b','fontsize',12)
text(3,-0.1,num2str(round(max(abs(cross2y)),2)),'color','r','fontsize',12)
ylabel('Y')
xlabel('Lags (hours)')
legend('BOPT','SCTA','location','northeast')
set(gca,'fontsize',14)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print(['../tidal_comp/Axial/xcorr_' svstr],'-dtiff','-r300')