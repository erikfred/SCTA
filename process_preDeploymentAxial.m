% Script to process 24 hour preDeployment flip data test for Axial SCTA

% Load Data
data = get_sctaDay('/Volumes/GoogleDrive/My Drive/Oceanography/SCTA-Share/PredeploymentTesting/SCTA-Accel_2018-05-24.nc');

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

% Process Flips
p.daMax = 1e-4;                 % During a calibration successive samples will not change by more than this
p.tCalLim = [60 90];            % Time limits for calibration in seconds since start of stable output
flipInfo2 = analyze_flips(dataDec,flipInfo,p,0);

% Plot calibrations
figure
clf
plot(flipInfo2.gCal,'o');
title('PreDeployment Axial SCTA Test')
xlabel('Flip')
ylabel('Calibration, m/s^2')