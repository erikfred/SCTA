% compare_calibrations_40v8.m
%
% For input date range, calculates and compares calibrations under
% different sample rates and methods. Goal is to figure out why -X drift is
% not consistently linear.
%

dayi=datenum(2019,11,01);
dayf=dayi; %same as dayi if inspecting single calibration

clear flipInfo40 flipInfo8

% flip & calibration conditions
p.dadT=[6.2023e-5 2.9562e-5 NaN];
p.TRef = 5.7;

p.cosThreshVert = 0.99;         % 0.99 = 8 degrees from vertical - threshold for a flip
p.minTime4Flip = 100;           % Minimum duration in seconds for a flip to be counted
p.complexRange80 = 1e-2;        % If 90% value - 10% value is greater than this then not a simple flip into one orientation
p.cosThreshNorm = 0.9996;       % Threshold for normal (one channel vertical) orientation (2°)
p.tBufferNorm = 120;            % Make non-normal any sample within this of a non-normal orientation
p.nMadNorm = 6;                 % If total acceleration in normal position is this many MADs from the median, make it non-normal

p.daMax = 1e-4;                 % During a calibration successive samples will not change by more than this
p.tCalLim = [60 90];            % Time limits for calibration in seconds since start of stable output

% load original calibration structure
if ~exist('flipInfoAll','var')
    load ../calibrations/Axial/axialdata
end
flipInfoOrig=flipInfoAll;

for i=dayi:dayf
    
    % download and load necessary data at 40 and 8 Hz
    data40=get_sctaDay('/Volumes/GoogleDrive/My Drive/Oceanography/SCTA-Share/OOI-SCTA/ParsedData',dayi);
    
    data8=[];
    cha={'MNE','MNN','MNZ','MKA'};
    chastr={'a(:,1)','a(:,2)','a(:,3)','T'};
    for m=1:length(cha)
        IRIS_data_pull('AXCC2',cha{m},'--',dayi,dayi+1)
        temp=rdmseed(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayi,29) '.miniseed']);
        data8.t=cat(1,temp.t);
        eval(['data8.' chastr{m} '=cat(1,temp.d)/10^7;']);
    end
    data8.as=sqrt(data8.a(:,1).^2+data8.a(:,2).^2+data8.a(:,3).^2);
    
    % decimate, find flips, calibrate
    dataDec40=decimate_SCTA(data40,1);
    [flipInfo,~]=find_flip(dataDec40.t,dataDec40.a,dataDec40.as,p);
    
    if isempty(flipInfo.t)
        fprintf(['No flips found on ' datestr(dayi) '\n\n'])
    elseif length(flipInfo.t)~=3 && length(flipInfo.t)~=5
        warning(['Peculiar number of flips found on ' datestr(dayi)])
        keyboard
    else
        flipInfo40=analyze_flips(dataDec40,flipInfo,p,1);
    end
    
    dataDec8=decimate_SCTA(data8,1);
    [flipInfo,~]=find_flip(dataDec8.t,dataDec8.a,dataDec8.as,p);
    
    if isempty(flipInfo.t)
        fprintf(['No flips found on ' datestr(dayi) '\n\n'])
    elseif length(flipInfo.t)~=3 && length(flipInfo.t)~=5
        warning(['Peculiar number of flips found on ' datestr(dayi)])
        keyboard
    else
        flipInfo8=analyze_flips(dataDec8,flipInfo,p,1);
    end
    
    % compare calibrations
    disp(' ')
    disp('40Hz Calibration')
    disp(flipInfo40.gCal(end))
    disp(' ')
    disp('40Hz Calibration (TCor)')
    disp(flipInfo40.gCalTCor(end))
    disp(' ')
    disp('diff')
    disp(flipInfo40.gCalTCor(end)-flipInfo40.gCal(end))
    disp(' ')
    disp('8Hz Calibration')
    disp(flipInfo8.gCal(end))
    disp(' ')
    disp('8Hz Calibration (TCor)')
    disp(flipInfo8.gCalTCor(end))
    disp(' ')
    disp('diff')
    disp(flipInfo8.gCalTCor(end)-flipInfo8.gCal(end))
    disp(' ')
    disp('8Hz - 40Hz')
    disp(flipInfo8.gCal(end)-flipInfo40.gCal(end))
    disp(' ')
    disp('8Hz - 40Hz (TCor)')
    disp(flipInfo8.gCalTCor(end)-flipInfo40.gCalTCor(end))
    
end