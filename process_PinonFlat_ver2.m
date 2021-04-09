% Script to process Pinon Flat Flips
%
% This version processes all days and creates a continuous record of flips

%% Parameters
% dataLoaded: 0 - need to load from raw files; 1 need to append to existing matlab file; 2 - already in memory
dataLoaded = 0;

% Start and end date
startDate = datenum('10/18/18'); %daily flips begin 10/18/18
endDate = floor(now-1);

% Temperature sensitivity parameters
p.dadT=[2.9562e-5 6.2023e-5 NaN]; % [dxdT dydT dzdT]
p.TRef=30;

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
        data = get_sctaDay('/Volumes/GoogleDrive/My Drive/Oceanography/SCTA-Share/OOI-PF/SCTA-PF/ParsedData',day);
        
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
                
            elseif length(flipInfo.t)~=3 && length(flipInfo.t)~=5
                
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
    save ../calibrations/PinonFlat/PFdata dataDec1 dataDec100 flipInfoAll
    
elseif dataLoaded==1
    load ../calibrations/PinonFlat/PFdata
    startDate2=floor(dataDec1.t(end));
    
    for dayn = startDate2:endDate
        data = get_sctaDay('/Volumes/GoogleDrive/My Drive/Oceanography/SCTA-Share/OOI-PF/SCTA-PF/ParsedData',dayn);
        
        if dayn==datenum(2020,02,06)
            keyboard
        end
        
        if isempty(data.t)
            
            fprintf(['No data on ' datestr(dayn) '\n\n'])
            
        else
            
            % Decimate the data
            [data1DayDec] = decimate_SCTA(data,1);
            [dataDec1] = decimate_SCTA(data,1,dataDec1);
            [dataDec100] = decimate_SCTA(data,100,dataDec100);
            
            % Find flips
            [flipInfo,lNormOrt] = find_flip(data1DayDec.t,data1DayDec.a,data1DayDec.as,p);
            
            if isempty(flipInfo.t)
                
                fprintf(['No flips found on ' datestr(dayn) '\n\n'])
                
            elseif length(flipInfo.t)~=3 && length(flipInfo.t)~=5
                
                warning(['Peculiar number of flips found on ' datestr(dayn)])
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
    save ../calibrations/PinonFlat/PFdata dataDec1 dataDec100 flipInfoAll -v7.3
    
end

% Plot calibrations
figure
clf
plot(flipInfoAll.t(1:3:60),flipInfoAll.gCal(1:3:60),'ok',flipInfoAll.t(1:3:60),flipInfoAll.gCalTCor(1:3:60),'xk','markersize',18);
hold on
plot(flipInfoAll.t(2:3:60),flipInfoAll.gCal(2:3:60),'or',flipInfoAll.t(2:3:60),flipInfoAll.gCalTCor(2:3:60),'xr','markersize',18);
plot(flipInfoAll.t(3:3:60),flipInfoAll.gCal(3:3:60),'ob',flipInfoAll.t(3:3:60),flipInfoAll.gCalTCor(3:3:60),'xb','markersize',18);
plot(flipInfoAll.t(65:5:end),flipInfoAll.gCal(65:5:end),'sk',flipInfoAll.t(65:5:end),flipInfoAll.gCalTCor(65:5:end),'+k','markersize',18);
plot(flipInfoAll.t(63:5:end),flipInfoAll.gCal(63:5:end),'sr',flipInfoAll.t(63:5:end),flipInfoAll.gCalTCor(63:5:end),'+r','markersize',18);
plot(flipInfoAll.t(65:5:end),(flipInfoAll.gCal(64:5:end)+flipInfoAll.gCal(65:5:end))/2,'^k',...
    flipInfoAll.t(65:5:end),(flipInfoAll.gCalTCor(64:5:end)+flipInfoAll.gCalTCor(65:5:end))/2,'+k','markersize',18);
plot(flipInfoAll.t(63:5:end),(flipInfoAll.gCal(62:5:end)+flipInfoAll.gCal(63:5:end))/2,'^r',...
    flipInfoAll.t(63:5:end),(flipInfoAll.gCalTCor(62:5:end)+flipInfoAll.gCalTCor(63:5:end))/2,'+r','markersize',18);
plot(flipInfoAll.t(61:5:end),flipInfoAll.gCal(61:5:end),'ok',flipInfoAll.t(61:5:end),flipInfoAll.gCalTCor(61:5:end),'xk','markersize',18);
plot(flipInfoAll.t(62:5:end),flipInfoAll.gCal(62:5:end),'or',flipInfoAll.t(62:5:end),flipInfoAll.gCalTCor(62:5:end),'xr','markersize',18);
plot(flipInfoAll.t(64:5:end),flipInfoAll.gCal(64:5:end),'ob',flipInfoAll.t(64:5:end),flipInfoAll.gCalTCor(64:5:end),'xb','markersize',18);
xl = xlim; yl = ylim;
plot([0 0]+xl(1)+diff(xl)/10,mean(yl)+[0 0.0001],'-k')
text(xl(1)+diff(xl)/9,mean(yl)+0.00005,'10^{-5} g')
legend('1st X','1st X (T Corrected)','Y','Y (T Corrected)','2nd X','2nd X (T Corrected)','-X','-X (T Corrected)',...
    '-Y','-Y (T Corrected)','X span','X span (T corrected)','Y span','Y span (T corrected)','location','best')
datetick
title({'Pinon Flat SCTA Calibrations',[datestr(startDate,'mmm dd, yyyy') ' - ' datestr(endDate,'mmm dd, yyyy')]})
xlabel('Date')
ylabel('Calibration, m/s^2')
set(gca,'FontSize',16)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print -djpeg ../calibrations/PinonFlat/process_PinonFlat.jpeg
print -dtiff ../calibrations/PinonFlat/process_PinonFlat.tiff -r300
!open ../calibrations/PinonFlat/process_PinonFlat.jpeg

