% Script to process Axial Seamount Flips
% This is the first try and is very simple just aimed at making sure data
% acquisition during flips is okay

p.dadT=[6.2023e-5 2.9562e-5 NaN];
p.TRef = 5.6;

% Days with flips
testDates = {'8/8/18' '8/9/18' '8/10/18' '8/11/18' '8/12/18' '8/13/18' '8/14/18' '8/15/18' '8/16/18' '8/17/18' '8/18/18' ...
              '8/19/18' '8/20/18' '8/27/18' '9/3/18' '9/10/18' '9/17/18' '9/24/18' '10/1/18' '10/8/18' '10/15/18' '10/22/18'};

flipInfoAll = [];

% Load Data
for i = 1 :length(testDates)
  if datenum(testDates{i})+1<now
    data = get_sctaDay('/users/wilcock/MyDrive/APL/SCTA-Share/OOI-SCTA/ParsedData',testDates{i});

    if isempty(data.t)
      
      fprintf(['No data on ' testDates{i} '\n\n'])
      
    else

      % Decimate the data
      [dataDec] = decimate_SCTA(data,1);

      % Find Flips
      p.cosThreshVert = 0.99;         % 0.99 = 8 degrees from vertical - threshold for a flip
      p.minTime4Flip = 100;           % Minimum duration in seconds for a flip to be counted
      p.complexRange80 = 1e-2;        % If 90% value - 10% value is greater than this then not a simple flip into one orientation
      p.cosThreshNorm = 0.9996;       % Threshold for normal (one channel vertical) orientation (2°)
      p.tBufferNorm = 120;            % Make non-normal any sample within this of a non-normal orientation
      p.nMadNorm = 6;                 % If total acceleration in normal position is this many MADs from the median, make it non-normal
      [flipInfo,lNormOrt] = find_flip(dataDec.t,dataDec.a,dataDec.as,p);
      
      if isempty(flipInfo.t)
  
        fprintf(['No flips found on ' testDates{i} '\n\n'])
      
      elseif length(flipInfo.t)~=3
        
        warning(['Odd number of flips found on ' testDates{i}])
        keyboard

      else
        
        % Process Flips
        p.daMax = 1e-4;                 % During a calibration successive samples will not change by more than this
        p.tCalLim = [60 90];            % Time limits for calibration in seconds since start of stable output
        flipInfo2 = analyze_flips(dataDec,flipInfo,p,1);
%         disp(['Paused for ' testDates{i}])
        
        if isempty(flipInfoAll)
          flipInfoAll = flipInfo2;
        else
          flipInfoAll = merge_oneElementStructure(flipInfoAll, flipInfo2, 'normalOrientation');
        end
      end
    end
  end
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
title({'Axial SCTA Calibrations',[datestr(testDates{1}) ' - ' datestr(floor(now-1))]})
xlabel('Date')
ylabel('Calibration, m/s^2')
xl = xlim; yl = ylim;
plot([0 0]+xl(1)+diff(xl)/10,mean(yl)+[0 0.0001],'-k')
text(xl(1)+diff(xl)/9,mean(yl)+0.00005,'1 part in 10^{5}')
print -djpeg process_Axial_ver1.jpeg
!open process_Axial_ver1.jpeg

