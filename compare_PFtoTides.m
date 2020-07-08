% Script to compare the Pinon Flat Tilt data to tidal model
%
% For Paroscientific Tilt Sensors 
%   Channel X was pointing north
%   Channel Y was pointing west
%   Their convention is to report +ve values when the sensor is tilted up (Z accelerations are positive)
%
% The above appears to be the orientation of the SCTA based on matching data
%
% For SPOTL and the Tidal model provided by Mark Zumberge the convention is
%   East tilt positive downwards
%   North tile postive downwards

% Load the data we collected - created by process_PinonFlat_ver1.m
load pfData

% Load tilts from Mark Zumberge
load ../zumberge/ew.tilt.tides.ascii
ewZum = ew_tilt_tides;
load ../zumberge/ns.tilt.tides.ascii
nsZum = ns_tilt_tides;
t0 = datenum('10/09/2018 00:00:00');
t1 = datenum('03/05/2020 23:55:00');
tzum = (t0:5/24/60:t1)';
i = tzum>datenum('2018-10-17') & tzum<datenum('2019-02-15');
tzum= tzum(i);
ewZum = ewZum(i);
nsZum = nsZum(i);

% Acceleration of gravity (Unknown for PF)
g = 9.81;

% First tiltmeter
t1 = tilt1Dec100.t;
i = t1>datenum('2018-10-17') & t1<datenum('2019-02-15');
t1 = t1(i);
x1 = tilt1Dec100.a(i,1);
y1 = tilt1Dec100.a(i,2);

% Second tiltmeter
t2 = tilt2Dec100.t;
i = t2>datenum('2018-10-17') & t2<datenum('2019-02-15');
t2 = t2(i);
x2 = tilt2Dec100.a(i,1);
y2 = tilt2Dec100.a(i,2);

% Acceletometer
ta = sctaDec100.t;
i = ta>datenum('2018-10-17') & ta<datenum('2019-02-15');
ta = ta(i);
xa = sctaDec100.a(i,1);
ya = sctaDec100.a(i,2);

% A simple detrending of tiltmeter 1 based on daily average
tbar = floor(t1(1))+0.5:ceil(t1(end));
xbar = nan(size(tbar));
ybar = nan(size(tbar));
j = 0;
for t = tbar
  j = j+1;
  i = t1>t-0.5 & t1<t+0.5;
  xbar(j) = nanmean(x1(i));
  ybar(j) = nanmean(y1(i));
end
x1 = x1 - interp1(tbar,xbar,t1,'spline','extrap');
y1 = y1 - interp1(tbar,ybar,t1,'spline','extrap');

% A simple detrending of tiltmeter 2 based on daily average
tbar = floor(t2(1))+0.5:ceil(t2(end));
xbar = nan(size(tbar));
ybar = nan(size(tbar));
j = 0;
for t = tbar
  j = j+1;
  i = t2>t-0.5 & t2<t+0.5;
  xbar(j) = nanmean(x2(i));
  ybar(j) = nanmean(y2(i));
end
x2 = x2 - interp1(tbar,xbar,t2,'spline','extrap');
y2 = y2 - interp1(tbar,ybar,t2,'spline','extrap');

% A simple detrending of accelerometer based on daily median (to ignore
% calibrations)
tbar = floor(ta(1))+0.5:ceil(ta(end));
xbar = nan(size(tbar));
ybar = nan(size(tbar));
j = 0;
for t = tbar
  j = j+1;
  i = ta>t-0.5 & ta<t+0.5;
  xbar(j) = nanmedian(xa(i));
  ybar(j) = nanmedian(ya(i));
end
xa = xa - interp1(tbar,xbar,ta,'spline','extrap');
ya = ya - interp1(tbar,ybar,ta,'spline','extrap');

thour1 = floor(t1(1)*24)/24;
thour2 = ceil(t1(end)*24)/24;

% SPOTL tilts
[nTiltPF,eTiltPF,timeTilt] = tide_tilt(33.609,-116.455,-thour1,-thour2,3600,1280,[],[],'datenum');

% East-West Plot
figure(1)
clf
plot(ta,ya/9.81*1e9,'g','linewidth',0.1)
hold on
plot(t1,y1/9.81*1e9,'y',t2,y2/9.81*1e9,'c')
plot(timeTilt,eTiltPF,'r',tzum,ewZum,'k','linewidth',1)
ylim([-150 150])
xlim([737404 737430])
legend('SCTA','Parosci. Tilt 1','Parosci. Tilt 2','SPOTL','Theoretical','location','southeast')
datetick('keeplimits')
ylabel('Tilt, nanoradians')
title('Pinon Flat East-West Tilt (+ve down in East Direction')
print -djpeg compare_PFtoTides_EW.jpg

xlim([737415 737417])
datetick('x',6,'keeplimits')
print -djpeg compare_PFtoTides_EWzoom.jpg

% North-South Plot
figure(2)
clf
plot(ta,-xa/9.81*1e9,'g','linewidth',0.1)
hold on
plot(t1,-x1/9.81*1e9,'y',t2,-x2/9.81*1e9,'c')
plot(timeTilt,nTiltPF,'r',tzum,nsZum,'k','linewidth',1)
ylim([-150 150])
xlim([737404 737430])
legend('SCTA','Parosci. Tilt 1','Parosci. Tilt 2','SPOTL','Theoretical','location','southeast')
datetick('keeplimits')
ylabel('Tilt, nanoradians')
title('Pinon Flat North-South Tilt (+ve down in North Direction')
print -djpeg compare_PFtoTides_NS.jpg

xlim([737415 737417])
datetick('x',6,'keeplimits')
print -djpeg compare_PFtoTides_NSzoom.jpg

