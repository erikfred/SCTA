% Script to process Axial Seamount Flips
% 
% This version processes all days and creates a continuous record of flips

%% Parameters
% dataLoaded: 0 - need to load from raw files; 1 - append to existing matlab file
dataLoaded = 1;

% Start and end date
% startDate = datenum('8/1/18');
startDate = datenum('09/11/20');
endDate = datenum('08/27/21'); % SCTA recovered 8/27/2021
if startDate==datenum('09/11/20')
    suffix='_newloc';
else
    suffix='';
end

% Temperature sensitivity parameters
p.dadT=[2.9562e-5 6.2023e-5 NaN]; % [dxdT dydT dzdT]
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
        data = get_sctaDay('/Users/erikfred/Google Drive/My Drive/Oceanography/SCTA-Share/OOI-SCTA/ParsedData',dayn);
        
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
            
            if dayn==datenum(2020,09,11) % remove first half day of data to account for instrument move
                cond=data.t>datenum(2020,09,11,12,0,0);
                data.a=data.a(cond,:);
                data.as=data.as(cond);
                data.T=data.T(cond);
                data.t=data.t(cond);
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
            
            if isempty(flipInfo.t) && dayn>datenum(2018,12,15) && str2double(datestr(dayn,'dd'))==1 || ...
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
                
                if dayn==datenum(2020,09,11) % remove first half day of data to account for instrument move
                    cond=data.t>datenum(2020,09,11,12,0,0);
                    data.a=data.a(cond);
                    data.as=data.as(cond);
                    data.T=data.T(cond);
                    data.t=data.t(cond);
                end
                
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
    save(['../calibrations/Axial/axialdata' suffix],'dataDec1','dataDec100','flipInfoAll','-v7.3')
    
elseif dataLoaded==1
    load(['../calibrations/Axial/axialdata' suffix])
    startDate2=floor(dataDec1.t(end))+1;
    
    for dayn = startDate2:endDate
        data = get_sctaDay('/Volumes/GoogleDrive/My Drive/Oceanography/SCTA-Share/OOI-SCTA/ParsedData',dayn);
        
        if isempty(data.t) && dayn<datenum(2019,08,13) %temporary fix, should apply to entire series
            
            fprintf(['No data on ' datestr(dayn) '\n\n'])
            
        else
            
            fprintf(['No 40Hz data on ' datestr(dayn) '\n\n'])
            
            if dayn>datenum(2020,05,09) && dayn<=datenum(2020,06,01)
                % datastream cut during this interval due to cable issues
                continue
            elseif dayn>datenum(2021,01,13) && dayn<=datenum(2021,01,17)
                % datastream cut during this interval due to cable issues
                continue
            elseif isempty(data.t) && dayn>=datenum(2019,08,13)
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
    save(['../calibrations/Axial/axialdata' suffix],'dataDec1','dataDec100','flipInfoAll','-v7.3')
    
end

% identify and separate each of the calibrations
i_x=find(flipInfoAll.orientation==1);
i_xb=find(flipInfoAll.orientation==1 & flipInfoAll.t>datenum(2019,08,13));
i_y=find(flipInfoAll.orientation==2);
i_yb=find(flipInfoAll.orientation==2 & flipInfoAll.t>datenum(2019,08,13));
i_negx=find(flipInfoAll.orientation==-1);
i_negy=find(flipInfoAll.orientation==-2);

% Plot calibrations
figure
clf
plot(flipInfoAll.t(i_x(1:2:end)),flipInfoAll.gCal(i_x(1:2:end)),'ok',flipInfoAll.t(i_x(1:2:end)),flipInfoAll.gCalTCor(i_x(1:2:end)),'xk','markersize',18);
hold on
plot(flipInfoAll.t(i_x(2:2:end)),flipInfoAll.gCal(i_x(2:2:end)),'ob',flipInfoAll.t(i_x(2:2:end)),flipInfoAll.gCalTCor(i_x(2:2:end)),'xb','markersize',18);
plot(flipInfoAll.t(i_negx),flipInfoAll.gCal(i_negx),'sk',flipInfoAll.t(i_negx),flipInfoAll.gCalTCor(i_negx),'+k','markersize',18);
plot(flipInfoAll.t(i_y),flipInfoAll.gCal(i_y),'or',flipInfoAll.t(i_y),flipInfoAll.gCalTCor(i_y),'xr','markersize',18);
plot(flipInfoAll.t(i_negy),flipInfoAll.gCal(i_negy),'sr',flipInfoAll.t(i_negy),flipInfoAll.gCalTCor(i_negy),'+r','markersize',18);
plot(flipInfoAll.t(i_xb(2:2:end)),(flipInfoAll.gCal(i_xb(2:2:end))+flipInfoAll.gCal(i_negx))/2,'^k',...
    flipInfoAll.t(i_xb(2:2:end)),(flipInfoAll.gCalTCor(i_xb(2:2:end))+flipInfoAll.gCalTCor(i_negx))/2,'+k','markersize',18);
plot(flipInfoAll.t(i_yb),(flipInfoAll.gCal(i_yb)+flipInfoAll.gCal(i_negy))/2,'^r',...
    flipInfoAll.t(i_yb),(flipInfoAll.gCalTCor(i_yb)+flipInfoAll.gCalTCor(i_negy))/2,'+r','markersize',18);
xl = xlim; yl = ylim;
plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.0005 0.0004],'-k','linewidth',2)
text(xl(1)+diff(xl)/9,mean(yl)-0.00045,'10 \mug','fontsize',12)
legend('1st X','1st X (T Corrected)','2nd X','2nd X (T Corrected)','-X','-X (T Corrected)','Y','Y (T Corrected)',...
    '-Y','-Y (T Corrected)','X span','X span (T corrected)','Y span','Y span (T corrected)','location','northwest')
datetick('x',3)
title({'Axial SCTA Calibrations',[datestr(startDate,'mmm dd, yyyy') ' - ' datestr(endDate,'mmm dd, yyyy')]})
xtickangle(45)
ylabel('Calibration, m/s^2')
set(gca,'fontsize',15)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print(['../calibrations/Axial/process_Axial' suffix],'-dtiff','-r300')

% Plot calibrations on separate axes
if strcmp(suffix,'_newloc')
    i_x1=i_x([1:2:end-6 end-3:2:end]); i_x1b=i_x(116);
    i_x2=i_x(2:2:end);
    i_xneg=i_negx([1:end-13 end-11:end]); i_xnegb=i_negx(48);
    
    figure
    subplot(311); hold on
    plot(flipInfoAll.t(i_x1),flipInfoAll.gCal(i_x1),'ko','markersize',12);
    datetick('x',3)
    xtickangle(45)
    ylabel('+X1 (m/s^2)')
    set(gca,'fontsize',12)
    ylim([9.81185 9.81205])
%     xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.00005 -0.00005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'10 \mug','fontsize',12)
    yyaxis right
    plot(flipInfoAll.t(i_x2),flipInfoAll.gCal(i_x2),'ro','markersize',12);
%     title({'Axial SCTA X Calibrations',[datestr(startDate,'mmm dd, yyyy') ' - ' datestr(endDate,'mmm dd, yyyy')]})
    datetick('x',3)
    xtickangle(45)
    ylabel('+X2 (m/s^2)')
    set(gca,'fontsize',12)
    set(gca,'YColor','r')
    ylim([9.81215 9.81215+0.0002])
    xline(flipInfoAll.t(i_x1b),'k:','linewidth',2)
    box on; grid on
    subplot(312); hold on
    plot(flipInfoAll.t(i_xneg),-1*flipInfoAll.gCal(i_xneg),'ks','markersize',12);
    datetick('x',3)
    xtickangle(45)
    ylabel('-X (m/s^2)')
    set(gca,'fontsize',12)
    ylim([-9.809 -9.8088])
    xline(flipInfoAll.t(i_xnegb),'k:','linewidth',2)
%     xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.00005 -0.00005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'10 \mug','fontsize',12)
    box on; grid on
    subplot(313); hold on
    xspan1=(flipInfoAll.gCal(intersect(i_x1,i_xneg-4))+flipInfoAll.gCal(intersect(i_x1+4,i_xneg)));
    plot(flipInfoAll.t(intersect(i_x1,i_xneg-4)),xspan1-xspan1(1),'k^','markersize',12);
    datetick('x',3)
    xtickangle(45)
    ylabel('\Delta X1 span (m/s^2)')
    set(gca,'fontsize',12)
    ylim([0 0.00004])
    yyaxis right
    xspan2=(flipInfoAll.gCal(intersect(i_x2,i_xneg-1))+flipInfoAll.gCal(intersect(i_x2+1,i_xneg)));
    plot(flipInfoAll.t(intersect(i_x2,i_xneg-1)),xspan2-xspan2(1),'r^','markersize',12);
    datetick('x',3)
    xtickangle(45)
    ylabel('\Delta X2 span (m/s^2)')
    set(gca,'fontsize',12)
    set(gca,'YColor','r')
    ylim([0 0.00004])
    xline(flipInfoAll.t(i_x1b),'k:','linewidth',2)
    xline(flipInfoAll.t(i_xnegb),'k:','linewidth',2)
%     xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)+[0 0.00001],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl)+0.00000,'1 \mug','fontsize',12)
    box on; grid on
    
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print(['../calibrations/Axial/process_Axial' suffix '_manuscript_x'],'-dtiff','-r300')
    print(['../calibrations/Axial/process_Axial' suffix '_manuscript_x'],'-depsc','-painters')
    
    i_y1=i_y;
    i_yneg=i_negy([1:36 38:end]); i_ynegb=i_negy(37);
    
    figure
    subplot(311); hold on
    plot(flipInfoAll.t(i_y1),flipInfoAll.gCal(i_y1),'ko','markersize',12);
%     title({'Axial SCTA Y Calibrations',[datestr(startDate,'mmm dd, yyyy') ' - ' datestr(endDate,'mmm dd, yyyy')]})
    datetick('x',3)
    xtickangle(45)
    ylabel('+Y (m/s^2)')
    set(gca,'fontsize',12)
    ylim([9.8111 9.8115])
%     xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.00005 -0.00005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'10 \mug','fontsize',12)
    box on; grid on
    subplot(312); hold on
    plot(flipInfoAll.t(i_yneg),-1*flipInfoAll.gCal(i_yneg),'ks','markersize',12);
    datetick('x',3)
    xtickangle(45)
    ylabel('-Y (m/s^2)')
    set(gca,'fontsize',12)
    ylim([-9.8111 -9.8107])
    xline(flipInfoAll.t(i_ynegb),'k:','linewidth',2)
%     xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.00005 -0.00005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'10 \mug','fontsize',12)
    box on; grid on
    subplot(313); hold on
    yspan=(flipInfoAll.gCal(intersect(i_y1,i_yneg-1))+flipInfoAll.gCal(intersect(i_y1+1,i_yneg)));
    plot(flipInfoAll.t(intersect(i_y1,i_yneg-1)),yspan-yspan(1),'k^','markersize',12);
    datetick('x',3)
    xtickangle(45)
    ylabel('\Delta Y span (m/s^2)')
    set(gca,'fontsize',12)
    ylim([0 0.00004])
    xline(flipInfoAll.t(i_ynegb),'k:','linewidth',2)
%     xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0 -0.00001],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'1 \mug','fontsize',12)
    box on; grid on
    
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print(['../calibrations/Axial/process_Axial' suffix '_manuscript_y'],'-dtiff','-r300')
    print(['../calibrations/Axial/process_Axial' suffix '_manuscript_y'],'-depsc','-painters')
elseif strcmp(suffix,'')
    figure
    subplot(311); hold on
    plot(flipInfoAll.t(i_x([1:2:end-6 end-3:2:end])),flipInfoAll.gCal(i_x([1:2:end-6 end-3:2:end])),'ko','markersize',12);
    datetick('x',3)
    xtickangle(45)
    ylabel('+X1 (m/s^2)')
    set(gca,'fontsize',12)
    ylim([9.8111 9.8121])
    xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.00005 -0.00005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'10 \mug','fontsize',12)
    yyaxis right
    plot(flipInfoAll.t(i_x(2:2:end)),flipInfoAll.gCal(i_x(2:2:end)),'ro','markersize',12);
    title({'Axial SCTA X Calibrations',[datestr(startDate,'mmm dd, yyyy') ' - ' datestr(endDate,'mmm dd, yyyy')]})
    datetick('x',3)
    xtickangle(45)
    ylabel('+X2 (m/s^2)')
    set(gca,'fontsize',12)
    set(gca,'YColor','r')
    ylim([9.8113 9.8123])
    box on; grid on
    subplot(312); hold on
    plot(flipInfoAll.t(i_negx([1:end-13 end-11:end])),-1*flipInfoAll.gCal(i_negx([1:end-13 end-11:end])),'ks','markersize',12);
    xlim(xl)
    datetick('x',3,'keeplimits')
    xtickangle(45)
    ylabel('-X (m/s^2)')
    set(gca,'fontsize',12)
    ylim([-9.8099 -9.8089])
    xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.00005 -0.00005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'10 \mug','fontsize',12)
    box on; grid on
    subplot(313); hold on
    xspan1=(flipInfoAll.gCal(i_xb([1:2:end-26 end-23:2:end-6 end-3:2:end]))+flipInfoAll.gCal(i_negx([1:end-13 end-11:end-3 end-1:end])));
    plot(flipInfoAll.t(i_xb([1:2:end-26 end-23:2:end-6 end-3:2:end])),xspan1-xspan1(1),'k^','markersize',12);
    xlim(xl)
    datetick('x',3,'keeplimits')
    xtickangle(45)
    ylabel('\Delta X1 span (m/s^2)')
    set(gca,'fontsize',12)
    ylim([0 0.0001])
    yyaxis right
    xspan2=(flipInfoAll.gCal(i_xb(2:2:end))+flipInfoAll.gCal(i_negx));
    plot(flipInfoAll.t(i_xb(2:2:end)),xspan2-xspan2(1),'r^','markersize',12);
    xlim(xl)
    datetick('x',3,'keeplimits')
    xtickangle(45)
    ylabel('\Delta X2 span (m/s^2)')
    set(gca,'fontsize',12)
    set(gca,'YColor','r')
    ylim([0 0.0001])
    xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)+[0 0.00001],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl)+0.00000,'1 \mug','fontsize',12)
    box on; grid on
    
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print(['../calibrations/Axial/process_Axial' suffix '_manuscript_x'],'-dtiff','-r300')
    print(['../calibrations/Axial/process_Axial' suffix '_manuscript_x'],'-depsc','-painters')
    
    figure
    subplot(311); hold on
    plot(flipInfoAll.t(i_y),flipInfoAll.gCal(i_y),'ko','markersize',12);
    title({'Axial SCTA Y Calibrations',[datestr(startDate,'mmm dd, yyyy') ' - ' datestr(endDate,'mmm dd, yyyy')]})
    datetick('x',3)
    xtickangle(45)
    ylabel('+Y (m/s^2)')
    set(gca,'fontsize',12)
    ylim([9.8098 9.8113])
    xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.00005 -0.00005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'10 \mug','fontsize',12)
    box on; grid on
    subplot(312); hold on
    plot(flipInfoAll.t(i_negy([1:36 38:end])),-1*flipInfoAll.gCal(i_negy([1:36 38:end])),'ks','markersize',12);
    xlim(xl)
    datetick('x',3,'keeplimits')
    xtickangle(45)
    ylabel('-Y (m/s^2)')
    set(gca,'fontsize',12)
    ylim([-9.8125 -9.8110])
    xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.00005 -0.00005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'10 \mug','fontsize',12)
    box on; grid on
    subplot(313); hold on
    yspan=(flipInfoAll.gCal(i_yb([1:36 38:end]))+flipInfoAll.gCal(i_negy([1:36 38:end])));
    plot(flipInfoAll.t(i_yb([1:36 38:end])),yspan-yspan(1),'k^','markersize',12);
    xlim(xl)
    datetick('x',3,'keeplimits')
    xtickangle(45)
    ylabel('\Delta Y span (m/s^2)')
    set(gca,'fontsize',12)
    ylim([0 0.00015])
    xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.000005 -0.000005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'1 \mug','fontsize',12)
    box on; grid on
    
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print(['../calibrations/Axial/process_Axial' suffix '_manuscipt_y'],'-dtiff','-r300')
    print(['../calibrations/Axial/process_Axial' suffix '_manuscipt_y'],'-depsc','-painters')
end