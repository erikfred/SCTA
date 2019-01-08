% Script to process Pinon Flat Flips
% This is the first try and is very simple just aimed at making sure data 
% acquisition during flips is okay and that tilt data plots well

getData = false;

p.dadT=[NaN NaN NaN];
p.TRef = 30;
% 
% Days with flips
startDate = datenum('2018-10-17');
waitTiltDays = 4;

if getData 
  flipInfoAll = [];
  tilt1Dec1 = [];
  tilt2Dec1 = [];
  tilt1Dec100 = [];
  tilt2Dec100 = [];
  sctaDec1 = [];
  sctaDec100 = [];

  % Load Data
  for dn = startDate:now
    testDate = datestr(dn);
    scta = get_sctaDay('/Users/wilcock/Mydrive/APL/SCTA-Share/OOI-PF/SCTA-PF/ParsedData',testDate);
    tilt1 = get_tiltDay('/Users/wilcock/Mydrive/APL/SCTA-Share/OOI-PF/SCTA-PF/ParsedData','SCTA-Tilt_20_19_',testDate);
    tilt2 = get_tiltDay('/Users/wilcock/Mydrive/APL/SCTA-Share/OOI-PF/SCTA-PF/ParsedData','SCTA-Tilt_22_25_',testDate);

    if isempty(scta.t)

      fprintf(['No data on ' testDate '\n\n'])

    else

      % Decimate the data
      [scta1DayDec] = decimate_SCTA(scta,1);
      [sctaDec1] = decimate_SCTA(scta,1,sctaDec1);
      [sctaDec100] = decimate_SCTA(scta,100,sctaDec100);
      if dn-startDate>=waitTiltDays
        [tilt1Dec100] = decimate_SCTA(tilt1,100,tilt1Dec100);
        [tilt2Dec100] = decimate_SCTA(tilt2,100,tilt2Dec100);
        if isempty(tilt1Dec1) 
          tilt1Dec1 = tilt1;
          tilt2Dec1 = tilt2;
        else
          tilt1Dec1 = merge_oneElementStructure(tilt1Dec1, tilt1);
          tilt2Dec1 = merge_oneElementStructure(tilt2Dec1, tilt2);
        end      
      end

      % Find Flips
      p.cosThreshVert = 0.99;         % 0.99 = 8 degrees from vertical - threshold for a flip
      p.minTime4Flip = 100;           % Minimum duration in seconds for a flip to be counted
      p.complexRange80 = 1e-2;        % If 90% value - 10% value is greater than this then not a simple flip into one orientation
      p.cosThreshNorm = 0.9996;       % Threshold for normal (one channel vertical) orientation (2°)
      p.tBufferNorm = 120;            % Make non-normal any sample within this of a non-normal orientation
      p.nMadNorm = 6;                 % If total acceleration in normal position is this many MADs from the median, make it non-normal
      [flipInfo,lNormOrt] = find_flip(scta1DayDec.t,scta1DayDec.a,scta1DayDec.as,p);

      if isempty(flipInfo(1).t)

        fprintf(['No flips found on ' testDate '\n\n'])

      elseif length(flipInfo(1).t)~=3

        warning(['Unextpected number of flips found on ' testDate])
        keyboard

      else

        % Process Flips
        p.daMax = 1e-4;                 % During a calibration successive samples will not change by more than this
        p.tCalLim = [60 90];            % Time limits for calibration in seconds since start of stable output
        flipInfo2 = analyze_flips(scta1DayDec,flipInfo,p,1);
  %         disp(['Paused for ' testDates{i}])

        if isempty(flipInfoAll)
          flipInfoAll = flipInfo2(1);
          for i=2:length(flipInfo2)
            flipInfoAll = merge_oneElementStructure(flipInfoAll, flipInfo2(i), 'normalOrientation');
          end
        else
          for i=1:length(flipInfo2)
            flipInfoAll = merge_oneElementStructure(flipInfoAll, flipInfo2(i), 'normalOrientation');
          end
        end
      end
    end
  end

  if length(tilt1Dec100.t) > length(tilt2Dec100.t)
    n = length(tilt2Dec100.t);
    tilt1Dec100.t = tilt1Dec100.t(1:n);
    tilt1Dec100.a = tilt1Dec100.a(1:n,:);
    tilt1Dec100.T = tilt1Dec100.T(1:n);
    tilt1Dec100.n = tilt1Dec100.n(1:n);
  elseif length(tilt1Dec100.t) < length(tilt2Dec100.t)
    n = length(tilt1Dec100.t);
    tilt2Dec100.t = tilt2Dec100.t(1:n);
    tilt2Dec100.a = tilt2Dec100.a(1:n,:);
    tilt2Dec100.T = tilt2Dec100.T(1:n);
    tilt2Dec100.n = tilt2Dec100.n(1:n);
  end

  save pfdata tilt1Dec1 tilt2Dec1 tilt1Dec100 tilt2Dec100 sctaDec1 sctaDec100 p flipInfoAll

else
  load pfdata
end
  
  % Consistency
disp('Uncorrected (1st - 2nd) - mean(1st - 2nd)')
(flipInfoAll.gCal(1:3:7) - flipInfoAll.gCal(3:3:9)) - mean(flipInfoAll.gCal(1:3:7) - flipInfoAll.gCal(3:3:9))
disp('T Corrected (1st - 2nd) - mean(1st - 2nd)')
(flipInfoAll.gCalTCor(1:3:7) - flipInfoAll.gCalTCor(3:3:9)) - mean(flipInfoAll.gCalTCor(1:3:7) - flipInfoAll.gCalTCor(3:3:9))

% Plot calibrations
figure
clf
plot(flipInfoAll.t(1:3:end),flipInfoAll.gCal(1:3:end),'ok',flipInfoAll.t(1:3:end),flipInfoAll.gCalTCor(1:3:end),'xk','markersize',18);
hold on
plot(flipInfoAll.t(2:3:end),flipInfoAll.gCal(2:3:end),'or',flipInfoAll.t(2:3:end),flipInfoAll.gCalTCor(2:3:end),'xr','markersize',18);
plot(flipInfoAll.t(3:3:end),flipInfoAll.gCal(3:3:end),'ob',flipInfoAll.t(3:3:end),flipInfoAll.gCalTCor(3:3:end),'xb','markersize',18);
legend('1st X','1st X (T Corrected)','Y','Y (T Corrected)','2nd X','2nd X (T Corrected)','location','east')
datetick
title({'Pinon Flat SCTA Calibrations',[datestr(startDate,'mmm dd, yyyy') ' - ' datestr(floor(now),'mmm dd, yyyy')]})
xlabel('Date')
ylabel('Calibration, m/s^2')
xl = xlim; yl = ylim;
plot([0 0]+xl(1)+diff(xl)/10,mean(yl)+[0 0.0001],'-k')
text(xl(1)+diff(xl)/9,mean(yl)+0.00005,'10^{-5} g')
print -djpeg SCTA_PF_ver1.jpeg
print -dtiff SCTA_PF_ver1.tiff -r300
!open SCTA_PF_ver1.jpeg

%% Plot tiltdata X
figure
clf
px1 = polyfit(tilt1Dec100.t - tilt1Dec100.t(1),tilt1Dec100.a(:,1),3);
subplot(311)
yyaxis left
plot(tilt1Dec100.t,tilt1Dec100.a(:,1),tilt1Dec100.t,polyval(px1,tilt1Dec100.t - tilt1Dec100.t(1)))
ylabel('Acceleration, m/s^{2}')
yyaxis right
plot(tilt1Dec100.t,tilt1Dec100.T)
ylabel('Temperature, °C')
title('X - Axis, Tiltmeter 20-19');
datetick
xlim([startDate+waitTiltDays now])

px2 = polyfit(tilt2Dec100.t - tilt2Dec100.t(1),tilt2Dec100.a(:,1),3);
subplot(312)
yyaxis left
plot(tilt2Dec100.t,tilt2Dec100.a(:,1),tilt2Dec100.t,polyval(px2,tilt2Dec100.t - tilt2Dec100.t(1)))
ylabel('Acceleration, m/s^{2}')
yyaxis right
plot(tilt1Dec100.t,tilt2Dec100.T)
ylabel('Temperature, °C')
title('X - Axis, Tiltmeter 22-25');
datetick
xlim([startDate+waitTiltDays now])

subplot(313)
yyaxis left
plot(tilt1Dec100.t,tilt1Dec100.a(:,1)-polyval(px1,tilt1Dec100.t - tilt1Dec100.t(1)),'b')
hold on
plot(tilt2Dec100.t,tilt2Dec100.a(:,1)-polyval(px2,tilt2Dec100.t - tilt2Dec100.t(1)),'-g')
ylabel('Acceleration, m/s^{2}')
yyaxis right
plot(tilt1Dec100.t,tilt2Dec100.T)
legend('20-19','22-25','Temperature')
datetick
title('X Detrended')
xlim([startDate+waitTiltDays now])
print -djpeg tiltX_PF_ver1.jpeg
!open tiltX_PF_ver1.jpeg

%% Plot tiltdata Y
figure
clf
py1 = polyfit(tilt1Dec100.t - tilt1Dec100.t(1),tilt1Dec100.a(:,2),3);
subplot(311)
yyaxis left
plot(tilt1Dec100.t,tilt1Dec100.a(:,2),tilt1Dec100.t,polyval(py1,tilt1Dec100.t - tilt1Dec100.t(1)))
ylabel('Acceleration, m/s^{2}')
yyaxis right
plot(tilt1Dec100.t,tilt1Dec100.T)
ylabel('Temperature, °C')
title('Y - Axis, Tiltmeter 20-19');
datetick
xlim([startDate+waitTiltDays now])

py2 = polyfit(tilt2Dec100.t - tilt2Dec100.t(1),tilt2Dec100.a(:,2),3);
subplot(312)
yyaxis left
plot(tilt2Dec100.t,tilt2Dec100.a(:,2),tilt2Dec100.t,polyval(py2,tilt2Dec100.t - tilt2Dec100.t(1)))
ylabel('Acceleration, m/s^{2}')
yyaxis right
plot(tilt1Dec100.t,tilt2Dec100.T)
ylabel('Temperature, °C')
title('Y - Axis, Tiltmeter 22-25');
datetick
xlim([startDate+waitTiltDays now])

subplot(313)
yyaxis left
plot(tilt1Dec100.t,tilt1Dec100.a(:,2)-polyval(py1,tilt1Dec100.t - tilt1Dec100.t(1)),'b')
hold on
plot(tilt2Dec100.t,tilt2Dec100.a(:,2)-polyval(py2,tilt2Dec100.t - tilt2Dec100.t(1)),'-g')
ylabel('Acceleration, m/s^{2}')
yyaxis right
plot(tilt1Dec100.t,tilt2Dec100.T)
legend('20-19','22-25','Temperature','location','south')
datetick
title('Y Detrended')
xlim([startDate+waitTiltDays now])
print -djpeg tiltY_PF_ver1.jpeg
!open tiltY_PF_ver1.jpeg

%% Plot Detrended Tilt
figure
clf

subplot(211)
yyaxis left
plot(tilt1Dec100.t,tilt1Dec100.a(:,1)-polyval(px1,tilt1Dec100.t - tilt1Dec100.t(1)),'b')
hold on
plot(tilt2Dec100.t,tilt2Dec100.a(:,1)-polyval(px2,tilt2Dec100.t - tilt2Dec100.t(1)),'-g')
ylabel('Acceleration, m/s^{2}')
yyaxis right
plot(tilt1Dec100.t,tilt2Dec100.T)
legend('20-19','22-25','Temperature')
datetick
title('X Detrended')
xlim([startDate+waitTiltDays now])

subplot(212)
yyaxis left
plot(tilt1Dec100.t,tilt1Dec100.a(:,2)-polyval(py1,tilt1Dec100.t - tilt1Dec100.t(1)),'b')
hold on
plot(tilt2Dec100.t,tilt2Dec100.a(:,2)-polyval(py2,tilt2Dec100.t - tilt2Dec100.t(1)),'-g')
ylabel('Acceleration, m/s^{2}')
yyaxis right
plot(tilt1Dec100.t,tilt2Dec100.T)
legend('20-19','22-25','Temperature')
datetick
title('Y Detrended')
xlim([startDate+waitTiltDays now])
print -djpeg tiltDetrended_ver1.jpeg
!open tiltDetrended_ver1.jpeg

%% Plot tiltdata on Accelerometer
figure

tCorrect0 = flipInfoAll.t(end-3);
tCorrect1 = flipInfoAll.t(end);
Xcorrect1 = flipInfoAll.gCal(end-2) - flipInfoAll.gCal(end-5);
Ycorrect1 = flipInfoAll.gCal(end-1) - flipInfoAll.gCal(end-4);

Date1 = datenum('Nov 12, 2018 21:15');
Date1 = datenum('Nov 13, 2018');
Date2 = datenum('Dec 6, 2018 20:50');
i1 = (sctaDec1.t > Date1 & sctaDec1.t < Date2);
i100 = (sctaDec100.t > Date1 & sctaDec100.t < Date2);

Xc = (sctaDec100.t(i100)-tCorrect0)/(tCorrect1-tCorrect0)*Xcorrect1;
Yc = (sctaDec100.t(i100)-tCorrect0)/(tCorrect1-tCorrect0)*Ycorrect1;

clf
subplot(211)
yyaxis left
h1 = plot(sctaDec100.t(i100),sctaDec100.a(i100,1),'b');
hold on
xp = xlim;
xp = xp(1)*2/3 + xp(2)/3;
plot([xp xp],[0.003681 0.003686],'-k','linewidth',2)
text(xp+0.5,0.0036835,'0.5 \murad')
% plot(sctaDec100.t(i100),sctaDec100.a(i100,1)+Xc,'g','linewidth',2)
ylabel('Acceleration, m/s^{2}')
ylim([0.00368 0.00370]);
yyaxis right
h2 = plot(sctaDec1.t(i1),sctaDec1.T(i1),'linewidth',1)
ylabel('Temperature, °C')
% legend('SCTA 0.01 Hz','SCTA 0.01Hz - no drift','Temperature','location','east')
legend([h1 h2],'SCTA 0.01 Hz','Temperature','location','east')
datetick
title('X SCTA')
xlim([Date1 now])

subplot(212)
yyaxis left
h1 = plot(sctaDec100.t(i100),sctaDec100.a(i100,2),'b');
hold on
xp = xlim;
xp = xp(1)/2 + xp(2)/2;
plot([xp xp],[-0.039185 -0.039175],'-k','linewidth',2);
text(xp+0.5,-0.03918,'1 \murad')
% plot(sctaDec100.t(i100),sctaDec100.a(i100,2)+Yc,'g','linewidth',2)
ylabel('Acceleration, m/s^{2}')
ylim([-0.03922 -0.03917]);
yyaxis right
h2 = plot(sctaDec1.t(i1),sctaDec1.T(i1),'linewidth',1)
ylabel('Temperature, °C')
% legend('SCTA 0.01 Hz','SCTA 0.01Hz - no drift','Temperature','location','east')
legend([h1 h2],'SCTA 0.01 Hz','Temperature','location','east')
datetick
title('Y SCTA')
xlim([Date1 now])
print -djpeg sctaTilt_ver1.jpeg
print -dtiff sctaTilt_ver1.tiff -r300
!open sctaTilt_ver1.jpeg


%% Plot 2nd tiltmeter for AGU
%% Plot tiltdata X
figure
clf

i = tilt2Dec100.t<datenum('Nov 16, 2018');
subplot(211)
px2 = polyfit(tilt2Dec100.t(i) - tilt2Dec100.t(1),tilt2Dec100.a(i,1),3);
yyaxis left
plot(tilt2Dec100.t(i),tilt2Dec100.a(i,1),'linewidth',1)
hold on
plot(tilt2Dec100.t(i),polyval(px2,tilt2Dec100.t(i) - tilt2Dec100.t(1)))
ylabel('Acceleration, m/s^{2}')
yyaxis right
plot(tilt1Dec100.t(i),tilt2Dec100.T(i),'linewidth',1)
ylabel('Temperature, °C')
title('X - Axis, Tiltmeter 22-25');
datetick
xlim([startDate+waitTiltDays datenum('11/16/18')])

py2 = polyfit(tilt2Dec100.t(i) - tilt2Dec100.t(1),tilt2Dec100.a(i,2),3);
subplot(212)
yyaxis left
plot(tilt2Dec100.t(i),tilt2Dec100.a(i,2),'linewidth',1)
hold on
plot(tilt2Dec100.t(i),polyval(py2,tilt2Dec100.t(i) - tilt2Dec100.t(1)))
ylabel('Acceleration, m/s^{2}')
yyaxis right
plot(tilt1Dec100.t(i),tilt2Dec100.T(i),'linewidth',1)
ylabel('Temperature, °C')
title('Y - Axis, Tiltmeter 22-25');
datetick
xlim([startDate+waitTiltDays datenum('11/16/18')])

print -djpeg tiltAGU2018_PF_ver1.jpeg
print -dtiff tiltAGU2018_PF_ver1.tiff -r300

