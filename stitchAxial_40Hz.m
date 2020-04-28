% stitchAxial_40Hz.m
%
% Stitches together 40Hz data and performs calibrations
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

load ../calibrations/Axial/axialdata.mat

startdate=datenum(2018,08,01);

tflip=flipInfoAll.t(2:3:end); %times of y-flips (middle of calibration)
tflip(tflip<startdate)=[]; %cuts out early, short data segments

stitch.t=[];
stitch.x1=[];
stitch.x2=[];
stitch.y=[];
for i=1:length(tflip)
    %isolate segment between flips
    [~,it]=min(abs(dataDec100.t-tflip(i)));
    if i==length(tflip)
        ft=length(dataDec100.t);
    else
        [~,ft]=min(abs(dataDec100.t-tflip(i+1)));
    end
    
    ttemp=dataDec100.t(it+60:ft-10);
    etemp=dataDec100.a(it+60:ft-10,1);
    ntemp=dataDec100.a(it+60:ft-10,2);
    
    if i~=length(tflip)
        %determine drift rate between calibrations
        iy1=find(flipInfoAll.t==tflip(i));
        iy2=find(flipInfoAll.t==tflip(i+1));
        xcal1=flipInfoAll.gCalTCor(iy1-1)-flipInfoAll.gCalTCor(iy2-1);
        xcal2=flipInfoAll.gCalTCor(iy1+1)-flipInfoAll.gCalTCor(iy2+1);
        ycal=flipInfoAll.gCalTCor(iy1)-flipInfoAll.gCalTCor(iy2);
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
hold on
plot(stitch.t,stitch.LAX1,'linewidth',1)
plot(stitch.t,stitch.LAX2,'linewidth',1)
legend('1st cal.','2nd cal.')
datetick
set(gca,'fontsize',15)

figure(68)
hold on
plot(stitch.t,stitch.LAY,'linewidth',1)
datetick
set(gca,'fontsize',15)