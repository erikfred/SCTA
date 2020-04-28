% Script to process Axial Seamount Flips
% 
% This version processes all days and creates a continuous record of flips

%% Parameters
% dataLoaded: 0 - need to load from raw files; 1 need to append to existing matlab file; 2 - already in memory
dataLoaded = 1;

% Start and end date
startDate = datenum('8/1/18');
% startDate = datenum('10/8/18');
endDate = floor(now-1);

% Temperature sensitivity parameters
p.dadT=[6.2023e-5 2.9562e-5 NaN];
p.TRef = 5.7;

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
    for dayn = startDate:endDate
        data = get_sctaDay('/Volumes/GoogleDrive/My Drive/Oceanography/SCTA-Share/OOI-SCTA/ParsedData',dayn);
        
        if isempty(data.t) && dayn<datenum(2019,08,13) %temporary fix, should apply to entire series
            
            fprintf(['No data on ' datestr(dayn) '\n\n'])
            
        else
            if isempty(data.t) && dayn>=datenum(2019,08,13)
                data=[];
                cha={'MNE','MNN','MNZ','MKA'};
                chastr={'a(:,1)','a(:,2)','a(:,3)','T'};
                for m=1:length(cha)
                    IRIS_data_pull('AXCC2',cha{m},'--',dayn,dayn+1)
                    temp=rdmseed(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayn,29) '.miniseed']);
                    data.t=cat(1,temp.t);
                    eval(['data.' chastr{m} '=cat(1,temp.d)/10^7;']);
                end
                data.as=sqrt(data.a(:,1).^2+data.a(:,2).^2+data.a(:,3).^2);
            end
            
            % Decimate the data
            [data1DayDec] = decimate_SCTA(data,1);
            [dataDec1] = decimate_SCTA(data,1,dataDec1);
            [dataDec100] = decimate_SCTA(data,100,dataDec100);
            
            % Find flips using undecimated data
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
            
            if isempty(flipInfo.t) && dayn>datenum(2018,12,15) && str2double(datestr(dayn,'dd'))==1
                data=[];
                cha={'MNE','MNN','MNZ','MKA'};
                chastr={'a(:,1)','a(:,2)','a(:,3)','T'};
                for m=1:length(cha)
                    IRIS_data_pull('AXCC2',cha{m},'--',dayn,dayn+1)
                    temp=rdmseed(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayn,29) '.miniseed']);
                    data.t=cat(1,temp.t);
                    eval(['data.' chastr{m} '=cat(1,temp.d)/10^7;']);
                end
                data.as=sqrt(data.a(:,1).^2+data.a(:,2).^2+data.a(:,3).^2);
                
                [data1DayDec] = decimate_SCTA(data,1);
                [flipInfo,lNormOrt] = find_flip(data1DayDec.t,data1DayDec.a,data1DayDec.as,p);
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
    save ../calibrations/Axial/axialdata dataDec1 dataDec100 flipInfoAll -v7.3
    
elseif dataLoaded==1
    load ../calibrations/Axial/axialdata
    startDate2=floor(dataDec1.t(end))+1;
    
    for dayn = startDate2:endDate
        data = get_sctaDay('/Volumes/GoogleDrive/My Drive/Oceanography/SCTA-Share/OOI-SCTA/ParsedData',dayn);
        
        if isempty(data.t) && dayn<datenum(2019,08,13) %temporary fix, should apply to entire series
            
            fprintf(['No data on ' datestr(dayn) '\n\n'])
            
        else
            
            fprintf(['No 40Hz data on ' datestr(dayn) '\n\n'])
            
            if isempty(data.t) && dayn>=datenum(2019,08,13)
                data=[];
                cha={'MNE','MNN','MNZ','MKA'};
                chastr={'a(:,1)','a(:,2)','a(:,3)','T'};
                for m=1:length(cha)
                    IRIS_data_pull('AXCC2',cha{m},'--',dayn,dayn+1)
                    temp=rdmseed(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayn,29) '.miniseed']);
                    data.t=cat(1,temp.t);
                    eval(['data.' chastr{m} '=cat(1,temp.d)/10^7;']);
                end
                data.as=sqrt(data.a(:,1).^2+data.a(:,2).^2+data.a(:,3).^2);
            end
            
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
            
            if (isempty(flipInfo.t) && dayn<datenum(2019,8,13) && str2double(datestr(dayn,'dd'))==1) || ...
                    (isempty(flipInfo.t) && dayn>=datenum(2019,8,13) && dayn<datenum(2019,11,30) && strcmp(datestr(dayn,'ddd'),'Tue')) || ...
                    (isempty(flipInfo.t) && dayn>=datenum(2019,11,30) && dayn<datenum(2020,01,08) && str2double(datestr(dayn,'dd'))==1) ||...
                    (isempty(flipInfo.t) && dayn>=datenum(2020,01,08) && strcmp(datestr(dayn,'ddd'),'Wed'))
                data=[];
                cha={'MNE','MNN','MNZ','MKA'};
                chastr={'a(:,1)','a(:,2)','a(:,3)','T'};
                for m=1:length(cha)
                    IRIS_data_pull('AXCC2',cha{m},'--',dayn,dayn+1)
                    temp=rdmseed(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayn,29) '.miniseed']);
                    data.t=cat(1,temp.t);
                    eval(['data.' chastr{m} '=cat(1,temp.d)/10^7;']);
                end
                data.as=sqrt(data.a(:,1).^2+data.a(:,2).^2+data.a(:,3).^2);
                
                [data1DayDec] = decimate_SCTA(data,1);
                [flipInfo,lNormOrt] = find_flip(data1DayDec.t,data1DayDec.a,data1DayDec.as,p);
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
    save ../calibrations/Axial/axialdata dataDec1 dataDec100 flipInfoAll -v7.3
    
end

% % Consistency
% disp('Uncorrected (1st - 2nd) - mean(1st - 2nd)')
% (flipInfoAll.gCal(1:3:7) - flipInfoAll.gCal(3:3:9)) - mean(flipInfoAll.gCal(1:3:7) - flipInfoAll.gCal(3:3:9));
% disp('T Corrected (1st - 2nd) - mean(1st - 2nd)')
% (flipInfoAll.gCalTCor(1:3:7) - flipInfoAll.gCalTCor(3:3:9)) - mean(flipInfoAll.gCalTCor(1:3:7) - flipInfoAll.gCalTCor(3:3:9));

% Plot calibrations
figure
clf
plot(flipInfoAll.t(1:3:72),flipInfoAll.gCal(1:3:72),'ok',flipInfoAll.t(1:3:72),flipInfoAll.gCalTCor(1:3:72),'xk','markersize',18);
hold on
plot(flipInfoAll.t(2:3:72),flipInfoAll.gCal(2:3:72),'or',flipInfoAll.t(2:3:72),flipInfoAll.gCalTCor(2:3:72),'xr','markersize',18);
plot(flipInfoAll.t(3:3:72),flipInfoAll.gCal(3:3:72),'ob',flipInfoAll.t(3:3:72),flipInfoAll.gCalTCor(3:3:72),'xb','markersize',18);
plot(flipInfoAll.t(77:5:end),flipInfoAll.gCal(77:5:end),'sk',flipInfoAll.t(77:5:end),flipInfoAll.gCalTCor(77:5:end),'+k','markersize',18);
plot(flipInfoAll.t(75:5:end),flipInfoAll.gCal(75:5:end),'sr',flipInfoAll.t(75:5:end),flipInfoAll.gCalTCor(75:5:end),'+r','markersize',18);
plot(flipInfoAll.t(77:5:end),(flipInfoAll.gCal(76:5:end)+flipInfoAll.gCal(77:5:end))/2,'^k',...
    flipInfoAll.t(77:5:end),(flipInfoAll.gCalTCor(76:5:end)+flipInfoAll.gCalTCor(77:5:end))/2,'+k','markersize',18);
plot(flipInfoAll.t(75:5:end),(flipInfoAll.gCal(74:5:end)+flipInfoAll.gCal(75:5:end))/2,'^r',...
    flipInfoAll.t(75:5:end),(flipInfoAll.gCalTCor(74:5:end)+flipInfoAll.gCalTCor(75:5:end))/2,'+r','markersize',18);
plot(flipInfoAll.t(73:5:end),flipInfoAll.gCal(73:5:end),'ok',flipInfoAll.t(73:5:end),flipInfoAll.gCalTCor(73:5:end),'xk','markersize',18);
plot(flipInfoAll.t(74:5:end),flipInfoAll.gCal(74:5:end),'or',flipInfoAll.t(74:5:end),flipInfoAll.gCalTCor(74:5:end),'xr','markersize',18);
plot(flipInfoAll.t(76:5:end),flipInfoAll.gCal(76:5:end),'ob',flipInfoAll.t(76:5:end),flipInfoAll.gCalTCor(76:5:end),'xb','markersize',18);
xl = xlim; yl = ylim;
plot([0 0]+xl(1)+diff(xl)/10,mean(yl)+[0 0.0001],'-k')
text(xl(1)+diff(xl)/9,mean(yl)+0.00005,'10^{-5} g')
legend('1st X','1st X (T Corrected)','Y','Y (T Corrected)','2nd X','2nd X (T Corrected)','-X','-X (T Corrected)',...
    '-Y','-Y (T Corrected)','X span','X span (T corrected)','Y span','Y span (T corrected)','location','northwest')
datetick
title({'Axial SCTA Calibrations',[datestr(startDate,'mmm dd, yyyy') ' - ' datestr(endDate,'mmm dd, yyyy')]})
xtickangle(45)
ylabel('Calibration, m/s^2')
set(gca,'fontsize',15)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print -djpeg ../calibrations/Axial/process_Axial_ver2b.jpeg
print -dtiff ../calibrations/Axial/process_Axial_ver2b.tiff -r300
!open ../calibrations/Axial/process_Axial_ver2b.jpeg

