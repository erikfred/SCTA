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

t0=datenum(2019,07,10);
tf=datenum(2019,07,16);

%load data
sta1={'AXCC2'};
cha1={'MNE','MNN'};
sta2={'AXCC1','AXEC2','AXID1'};
cha2={'LAX','LAY'};

t1=t0;
AXCC1.time=[];AXCC1.LAX=[];AXCC1.LAY=[];
AXEC2.time=[];AXEC2.LAX=[];AXEC2.LAY=[];
AXID1.time=[];AXID1.LAX=[];AXID1.LAY=[];
ASHES.time=[];ASHES.LAX=[];ASHES.LAY=[];
AXCC2.time=[];AXCC2.MNE=[];AXCC2.MNN=[];
%% PULL DATA
while t1<=tf
    t1_s=datestr(t1,31); t1_s=t1_s(1:10);
    
    % grab all BOPT data
    for i=1:length(sta2)
        for j=1:length(cha2)
            IRIS_data_pull(sta2{i},cha2{j},'11',t1,t1+1)
            if ~exist(['../tiltcompare/' sta2{i} '/' sta2{i} '_' cha2{j} '_' datestr(t1,29) '.miniseed'],'file')
                if ~exist(['../tiltcompare/' sta2{i} '/' sta2{i} '_' cha2{j} '_' datestr(t1,29) '.miniseed'],'file')
                    IRIS_data_pull(sta2{i},cha2{j},'11',t1,t1+1);
                end
            end
            temp=rdmseed(['../tiltcompare/' sta2{i} '/' sta2{i} '_' cha2{j} '_' datestr(t1,29) '.miniseed']);
            ttemp=cat(1,temp.t);
            dtemp=double(cat(1,temp.d));
            % decimate to 1 sample/min
            dtemp1=decimate(dtemp,6,'fir');
            data.a(:,j)=decimate(dtemp1,10,'fir');
            data.t=linspace(ttemp(1),ttemp(end),length(data.a))';
        end
        % append to appropriate structure
        eval([sta2{i} '.time=[' sta2{i} '.time; data.t];']);
        eval([sta2{i} '.LAX=[' sta2{i} '.LAX; data.a(:,1)];']);
        eval([sta2{i} '.LAY=[' sta2{i} '.LAY; data.a(:,2)];']);
        clear data
    end
    
    % grab AXCC2 (SCTA) data
    for k=1:length(cha1)
        IRIS_data_pull('AXCC2',cha1{k},'--',t1,t1+1)
        if ~exist(['../tiltcompare/AXCC2/AXCC2_' cha1{k} '_' datestr(t1,29) '.miniseed'],'file')
            if ~exist(['../tiltcompare/AXCC2/AXCC2_' cha1{k} '_' datestr(t1,29) '.miniseed'],'file')
                IRIS_data_pull('AXCC2',cha1{k},'--',t1,t1+1);
            end
        end
        temp=rdmseed(['../tiltcompare/AXCC2/AXCC2_' cha1{k} '_' datestr(t1,29) '.miniseed']);
        ttemp=cat(1,temp.t);
        dtemp=double(cat(1,temp.d))/10^7;
        % decimate to 1 sample/min
        dtemp1=decimate(dtemp,8,'fir');
        dtemp2=decimate(dtemp1,6,'fir');
        data.a(:,k)=decimate(dtemp2,10,'fir');
        data.t=linspace(ttemp(1),ttemp(end),length(data.a))';
    end
    % append to AXCC2
    AXCC2.time=[AXCC2.time; data.t];
    AXCC2.MNE=[AXCC2.MNE; data.a(:,1)];
    AXCC2.MNN=[AXCC2.MNN; data.a(:,2)];
    clear data
    
    t1=t1+1;
end

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

% 6 sample/day low pass filter
[b,a]=butter(4,2*(1/60/60/12)/(1/60),'low');
AXCC2.LAX=filtfilt(b,a,AXCC2.LAX);
AXCC2.LAY=filtfilt(b,a,AXCC2.LAY);
AXCC1.LAX=filtfilt(b,a,AXCC1.LAX);
AXCC1.LAY=filtfilt(b,a,AXCC1.LAY);
AXEC2.LAX=filtfilt(b,a,AXEC2.LAX);
AXEC2.LAY=filtfilt(b,a,AXEC2.LAY);
AXID1.LAX=filtfilt(b,a,AXID1.LAX);
AXID1.LAY=filtfilt(b,a,AXID1.LAY);

%% PLOTTING

figure
% Central Caldera
subplot(421)
plot(AXCC1.time,detrend(AXCC1.LAX-mean(AXCC1.LAX)),'linewidth',1)
hold on
plot(AXCC2.time,-detrend(AXCC2.LAX-mean(AXCC2.LAX)),'linewidth',1)
datetick('x','keeplimits')
set(gca,'xticklabel',[])
ylabel('x-tilt (\murad)')
title('Central Caldera')
legend('BOPT','SCTA','location','northeast')
subplot(423)
plot(AXCC2.time,detrend(AXCC2.LAY-mean(AXCC2.LAY)),'linewidth',1)
hold on
plot(AXCC1.time,-detrend(AXCC1.LAY-mean(AXCC1.LAY)),'linewidth',1)
datetick('x','keeplimits')
ylabel('y-tilt (\murad)')
legend('BOPT','SCTA','location','northeast')
% Eastern Caldera
subplot(422)
plot(AXEC2.time,detrend(AXEC2.LAX-mean(AXEC2.LAX)),'linewidth',1)
datetick('x','keeplimits')
set(gca,'xticklabel',[])
ylabel('x-tilt (\murad)')
title('Eastern Caldera')
subplot(424)
plot(AXEC2.time,detrend(AXEC2.LAY-mean(AXEC2.LAY)),'linewidth',1)
datetick('x','keeplimits')
ylabel('y-tilt (\murad)')
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
datetick('x','keeplimits')
set(gca,'xticklabel',[])
ylabel('x-tilt (\murad)')
title('International District')
subplot(428)
plot(AXID1.time,detrend(AXID1.LAY-mean(AXID1.LAY)),'linewidth',1)
datetick('x','keeplimits')
ylabel('y-tilt (\murad)')
    
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print(['../test'],'-dtiff','-r300')