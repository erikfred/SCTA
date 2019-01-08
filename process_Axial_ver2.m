% Script to process Axial Seamount Flips
% 
% This version processes all days and creates a continuous record of flips

%% Parameters
% dataLoaded: 0 - need to load from raw files; 1 need to get from matlab file; 2 - already in memory
dataLoaded = 0;

% Start and end date
startDate = datenum('8/6/18');
% startDate = datenum('10/8/18');
endDate = floor(now-1);

% Temperature sensitivity parameters
p.dadT=[6.2023e-5 2.9562e-5 NaN];
p.TRef = 5.6;

% Find flip parameters
p.cosThreshVert = 0.99;         % 0.99 = 8 degrees from vertical - threshold for a flip
p.minTime4Flip = 100;           % Minimum duration in seconds for a flip to be counted
p.complexRange80 = 1e-2;        % If 90% value - 10% value is greater than this then not a simple flip into one orientation
p.cosThreshNorm = 0.9996;       % Threshold for normal (one channel vertical) orientation (2°)
p.tBufferNorm = 120;            % Make non-normal any sample within this of a non-normal orientation
p.nMadNorm = 6;                 % If total acceleration in normal position is this many MADs from the median, make it non-normal

% Process flips parameters
p.daMax = 1e-4;                 % During a calibration successive samples will not change by more than this
p.tCalLim = [60 90];            % Time limits for calibration in seconds since start of stable output

%% Load Data

% Empty matrices
flipInfoAll = [];
dataDec1 = [];
dataDec100 = [];

if dataLoaded == 0
  for day = startDate:endDate
    data = get_sctaDay('/users/wilcock/MyDrive/APL/SCTA-Share/OOI-SCTA/ParsedData',day);

    if isempty(data.t)

      fprintf(['No data on ' datestr(day) '\n\n'])

    else

      % Decimate the data
      [data1DayDec] = decimate_SCTA(data,1);
      [dataDec1] = decimate_SCTA(data,1,dataDec1);
      [dataDec100] = decimate_SCTA(data,100,dataDec100);

      % Find flips
      [flipInfo,lNormOrt] = find_flip(data1DayDec.t,data1DayDec.a,data1DayDec.as,p);

      if isempty(flipInfo.t)

        fprintf(['No flips found on ' datestr(day) '\n\n'])

      elseif length(flipInfo.t)~=3

        warning(['Peculiar number of flips found on ' datestr(day)])
        keyboard

      else

        % Process Flips
        flipInfo2 = analyze_flips(data1DayDec,flipInfo,p,1);

        if isempty(flipInfoAll)
          flipInfoAll = flipInfo2;
        else
          flipInfoAll = merge_oneElementStructure(flipInfoAll, flipInfo2, 'normalOrientation');
        end
      end
    end
  end

  dataDec1 = NANgap_scta(dataDec1);
  dataDec100 = NANgap_scta(dataDec100);
  save axialdata dataDec1 dataDec100 flipInfoAll
  
elseif dataLoaded==1
  load axialdata 

end

% Consistency
disp('Uncorrected (1st - 2nd) - mean(1st - 2nd)')
(flipInfoAll.gCal(1:3:7) - flipInfoAll.gCal(3:3:9)) - mean(flipInfoAll.gCal(1:3:7) - flipInfoAll.gCal(3:3:9));
disp('T Corrected (1st - 2nd) - mean(1st - 2nd)')
(flipInfoAll.gCalTCor(1:3:7) - flipInfoAll.gCalTCor(3:3:9)) - mean(flipInfoAll.gCalTCor(1:3:7) - flipInfoAll.gCalTCor(3:3:9));

% Plot calibrations
figure
clf
plot(flipInfoAll.t(1:3:end),flipInfoAll.gCal(1:3:end),'ok',flipInfoAll.t(1:3:end),flipInfoAll.gCalTCor(1:3:end),'xk','markersize',18);
hold on
plot(flipInfoAll.t(2:3:end),flipInfoAll.gCal(2:3:end),'or',flipInfoAll.t(2:3:end),flipInfoAll.gCalTCor(2:3:end),'xr','markersize',18);
plot(flipInfoAll.t(3:3:end),flipInfoAll.gCal(3:3:end),'ob',flipInfoAll.t(3:3:end),flipInfoAll.gCalTCor(3:3:end),'xb','markersize',18);
legend('1st X','1st X (T Corrected)','Y','Y (T Corrected)','2nd X','2nd X (T Corrected)','location','east')
datetick
title({'Axial SCTA Calibrations',[datestr(startDate,'mmm dd, yyyy') ' - ' datestr(endDate,'mmm dd, yyyy')]})
xlabel('Date')
ylabel('Calibration, m/s^2')
xl = xlim; yl = ylim;
plot([0 0]+xl(1)+diff(xl)/10,mean(yl)+[0 0.0001],'-k')
text(xl(1)+diff(xl)/9,mean(yl)+0.00005,'10^{-5} g')
print -djpeg process_Axial_ver1.jpeg
print -dtiff process_Axial_ver1.tiff -r300
!open process_Axial_ver2.jpeg

