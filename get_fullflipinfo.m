% get_fullflipinfo.m
%
% Goes through all flips in flipInfoAll and extracts more details about the
% calibration (values for all channels, ???)
%

% dataLoaded: 0 - need to load from raw files; 1 - append to existing matlab file
dataLoaded = 0;

% interval: 0 - before relocation; 1 - after relocation
interval = 1;
if interval==0
    suffix='';
else
    suffix='_newloc';
end

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

load(['../calibrations/Axial/axialdata' suffix],'flipInfoAll')
[~,id,~]=unique(floor(flipInfoAll.t));
daylist=floor(flipInfoAll.t(id));
daylist(daylist<datenum(2019,08,13))=[]; % exclude 3-position calibrations

flipInfoSome=[];
dataDec1=[];
dataDec100=[];

if dataLoaded==1
    load(['../calibrations/Axial/detailed_flipInfo' suffix])
    
    id2=find(floor(flipInfoAll.t)>floor(flipInfoSome.t(end)),1);
    daylist=unique(floor(flipInfoAll.t(id2:5:end))); % all the new calibrations
end

for dayn = daylist'
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
            flipInfo2 = analyze_flips_ver2(data1DayDec,flipInfo,p,1);
            
            if isempty(flipInfoSome)
                flipInfoSome = flipInfo2;
            else
                flipInfoSome = merge_oneElementStructure(flipInfoSome, flipInfo2, 'normalOrientation');
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
            flipInfo2 = analyze_flips_ver2(data1DayDec,flipInfo,p,1);
            if isempty(flipInfoSome)
                flipInfoSome = flipInfo2;
            else
                flipInfoSome = merge_oneElementStructure(flipInfoSome, flipInfo2, 'normalOrientation');
            end
        end
    end
end
save(['../calibrations/Axial/detailed_flipInfo' suffix],'flipInfoSome')

% % write to matrix
% calmat=floor(flipInfoSome.t(1:5:end));
% calmat(:,2)=flipInfoSome.gCal(1:5:end);
% calmat(:,3)=flipInfoSome.xCal(1:5:end);
% calmat(:,4)=flipInfoSome.yCal(1:5:end);
% calmat(:,5)=flipInfoSome.zCal(1:5:end);
% calmat(:,6)=flipInfoSome.gCal(2:5:end);
% calmat(:,7)=flipInfoSome.xCal(2:5:end);
% calmat(:,8)=flipInfoSome.yCal(2:5:end);
% calmat(:,9)=flipInfoSome.zCal(2:5:end);
% calmat(:,10)=flipInfoSome.gCal(3:5:end);
% calmat(:,11)=flipInfoSome.xCal(3:5:end);
% calmat(:,12)=flipInfoSome.yCal(3:5:end);
% calmat(:,13)=flipInfoSome.zCal(3:5:end);
% calmat(:,14)=flipInfoSome.gCal(4:5:end);
% calmat(:,15)=flipInfoSome.xCal(4:5:end);
% calmat(:,16)=flipInfoSome.yCal(4:5:end);
% calmat(:,17)=flipInfoSome.zCal(4:5:end);
% calmat(:,18)=flipInfoSome.gCal(5:5:end);
% calmat(:,19)=flipInfoSome.xCal(5:5:end);
% calmat(:,20)=flipInfoSome.yCal(5:5:end);
% calmat(:,21)=flipInfoSome.zCal(5:5:end);
% 
% writematrix(calmat,'../calibrations/Axial/detailed_flipInfo.xls')