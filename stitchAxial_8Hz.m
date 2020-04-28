% stitchAxial_8Hz.m
%
% Stitches together 8Hz data stream and performs calibrations. Options to
% compare against other tilt instruments at Axial.
%
% Future fixes:
%   - address plateaus that occur when NaNs at beginning/end of calibration
%   segments
%       - similarly, having each segment start where the last left off is
%       overly simple, since I am cutting out a few hours of data during
%       and after each calibration. Should fit a linear trend and calculate
%       where next segment expected to start.
%   - options to remove offsets from glitches/fishbumps/etc.
%

clear; close all

molist=['Dec';'Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec2'];

%startdate=datenum(2018,12,01);
startdate=datenum(2019,09,30);

load('../calibrations/Axial/axialdata.mat','flipInfoAll')
load('../compass_directions.mat')

for j=1:size(molist,1)
    load(['../monthly_plots/' molist(j,:) '/AXCC2.mat'])
    if ~exist('data','var')
        data=AXCC2;
    else
        data=merge_oneElementStructure(data, AXCC2, 'iflip');
    end
end

tflip=flipInfoAll.t([2:3:68, 71:5:end]); %times of y-flips (middle of calibration)
tflip(tflip<startdate)=[]; %cuts out early, short data segments

stitch.t=[];
stitch.x1=[];
stitch.x2=[];
stitch.y=[];
for i=1:length(tflip)
    %isolate segment between flips
    [~,it]=min(abs(data.time-tflip(i)));
    if i==length(tflip)
        ft=length(data.time);
    else
        [~,ft]=min(abs(data.time-tflip(i+1)));
    end
    
    ttemp=data.time(it+60:ft-10);
    etemp=data.MNE(it+60:ft-10);
    ntemp=data.MNN(it+60:ft-10);
    
    if i~=length(tflip)
        %determine drift rate between calibrations
        iy1=find(flipInfoAll.t==tflip(i));
        iy2=find(flipInfoAll.t==tflip(i+1));
        xcal1=flipInfoAll.gCalTCor(iy2-1)-flipInfoAll.gCalTCor(iy1-1);
        xcal2=flipInfoAll.gCalTCor(iy2+1)-flipInfoAll.gCalTCor(iy1+1);
        ycal=flipInfoAll.gCalTCor(iy2)-flipInfoAll.gCalTCor(iy1);
    end
    
    %apply drift correction
    x1lin=linspace(0,xcal1,length(etemp))';
    x2lin=linspace(0,xcal2,length(etemp))';
    ylin=linspace(0,ycal,length(ntemp))';
    
    etemp1=etemp-x1lin;
    etemp2=etemp-x2lin;
    ntemp=ntemp-ylin;

    %remove offset after calibration
    if ~isempty(stitch.x1)
        ix1a=find(~isnan(etemp1),1);
        ix1b=find(~isnan(stitch.x1),1,'last');
        ix2a=find(~isnan(etemp2),1);
        ix2b=find(~isnan(stitch.x2),1,'last');
        iya=find(~isnan(ntemp),1);
        iyb=find(~isnan(stitch.y),1,'last');
        
        etemp1=etemp1-(etemp1(ix1a)-stitch.x1(ix1b));
        etemp2=etemp2-(etemp2(ix2a)-stitch.x2(ix2b));
        ntemp=ntemp-(ntemp(iya)-stitch.y(iyb));
    end
    
    stitch.t=[stitch.t;ttemp];
    stitch.x1=[stitch.x1;etemp1];
    stitch.x2=[stitch.x2;etemp2];
    stitch.y=[stitch.y;ntemp];
end

%convert accel to tilt in microrad
stitch.LAX1=asin(stitch.x1/9.81)*10^6;stitch.LAX1=stitch.LAX1-stitch.LAX1(1);
stitch.LAX2=asin(stitch.x2/9.81)*10^6;stitch.LAX2=stitch.LAX2-stitch.LAX2(1);
stitch.LAY=asin(stitch.y/9.81)*10^6;stitch.LAY=stitch.LAY-stitch.LAY(1);

figure(67)
clf
hold on
%can't combine 2nd X calibrations from different orders
% plot(stitch.t,stitch.LAX2,'linewidth',2,'color',[0.8 0.8 0.8])
% plot(stitch.t,-stitch.LAX1,'linewidth',1,'color','k')
plot(stitch.t,-stitch.LAX1,'linewidth',1,'color','k')
legend('SCTA','location','northeast')
datetick
set(gca,'fontsize',18)
ylabel('East Tilt (\murad)')

figure(68)
clf
hold on
plot(stitch.t,-stitch.LAY,'r','linewidth',1)
legend('SCTA','location','northeast')
datetick
set(gca,'fontsize',18)
ylabel('North Tilt (\murad)')

% %% compare to AXCC1 coarse "iris" tilt
% 
% %AXCC1 compass correction
% crd_rot=[cosd(CCMP(4).plus_y) sind(CCMP(4).plus_y); -sind(CCMP(4).plus_y) cosd(CCMP(4).plus_y)];
% temp=crd_rot*[coarse.x';coarse.y']; coarse.x_rot=temp(1,:)'; coarse.y_rot=temp(2,:)';
% 
% %start from 0
% coarse.x_rot=coarse.x_rot-coarse.x_rot(1);
% coarse.y_rot=coarse.y_rot-coarse.y_rot(1);
% 
% figure(67)
% plot(coarse.t,coarse.x_rot,'b','linewidth',1)
% legend('1st cal.','2nd cal.','Iris')
% 
% figure(68)
% plot(coarse.t,coarse.y_rot,'b','linewidth',1)
% legend('SCTA','Iris')

%% compare to AXCC1 Lily tiltmeter

t0=startdate;
tf=floor(data.time(end));

%load data
sta='AXCC1';
cha1={'LAX','LAY'};

t1=t0;
AXCC1.time=[];AXCC1.LAX=[];AXCC1.LAY=[];
while t1<=tf
    t1_s=datestr(t1,31); t1_s=t1_s(1:10);
    for i=1:2
        file_string=['AXCC1_' cha1{i} '_' t1_s '.miniseed'];
        %attempt download if file not found
        if ~exist(['../tiltcompare/AXCC1/' file_string],'file')
            IRIS_data_pull('AXCC1',cha1{i},'11',t1,t1+1);
        end
        %some dates have no data (power failure, etc.)
        if exist(['../tiltcompare/AXCC1/' file_string],'file')
            temp=rdmseed(['../tiltcompare/AXCC1/' file_string]);
            dtemp=double(cat(1,temp.d));
            if length(dtemp)<3600 %short timeseries will throw errors
                t1=t1+1;
                continue
            end
            dtempd=decimate(dtemp,6,'fir');
            dtempd=decimate(dtempd,10,'fir');
            eval(['AXCC1.' cha1{i} '=[AXCC1.' cha1{i} ';dtempd];']);
            if i==1
                ttemp=cat(1,temp.t);
                ttempd=decimate(ttemp,6,'fir');
                ttempd=decimate(ttempd,10,'fir');
                AXCC1.time=[AXCC1.time;ttempd];
            end
        end
    end
    t1=t1+1;
end

%AXCC1 compass correction
crd_rot=[cosd(CCMP(4).plus_y) sind(CCMP(4).plus_y); -sind(CCMP(4).plus_y) cosd(CCMP(4).plus_y)];
temp=crd_rot*[AXCC1.LAX';AXCC1.LAY']; AXCC1.LAX_rot=temp(1,:)'; AXCC1.LAY_rot=temp(2,:)';

%start from 0
AXCC1.LAX_rot=AXCC1.LAX_rot-AXCC1.LAX_rot(1);
AXCC1.LAY_rot=AXCC1.LAY_rot-AXCC1.LAY_rot(1);

figure(67)
plot(AXCC1.time,-AXCC1.LAX_rot,'b','linewidth',1)
legend('SCTA','LILY','location','northeast')

figure(68)
plot(AXCC1.time,AXCC1.LAY_rot,'b','linewidth',1)
legend('SCTA','LILY','location','northeast')


%% placeholder code  to be incorporated later
% intended to mark where calibrations are made to check that stitching is
% properly lining up segments

keyboard % to ensure this doesn't run until I get it in the proper place

[~,ical,~]=unique(round(flipInfoAll.t)); % finds calibration days
t_flip=flipInfoAll.t(ical); % finds time of calibration

for j=1:length(t_flip)
[~,im(j)]=min(abs(stitch.t-t_flip(j)));
end

plot(t_flip,stitch.LAX1(im),'ok','markersize',10)

%% save figures, variables

figure(67)
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../tiltcompare/SCTA_Lily_comp/east_Jul-Nov','-dtiff')

figure(68)
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../tiltcompare/SCTA_Lily_comp/north_Jul-Nov','-dtiff')

save('../tiltcompare/SCTA_Lily_comp/data','data','AXCC1','stitch')