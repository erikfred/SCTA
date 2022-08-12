% Script to process Pinon Flat Flips
%
% This version processes all days and creates a continuous record of flips

%% Parameters
% dataLoaded: 0 - need to load from raw files; 1 need to append to existing matlab file; 2 - already in memory
dataLoaded = 2;

% Start and end date
startDate = datenum('10/18/18'); % daily flips begin 10/18/18
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

if dataLoaded == 0
    % Empty matrices
    flipInfoAll = [];
    dataDec1 = [];
    dataDec100 = [];

    for day = startDate:endDate
        data = get_sctaDay('/Users/erikfred/Google Drive/My Drive/Oceanography/SCTA-Share/OOI-PF/SCTA-PF/ParsedData',day);
        
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
        data = get_sctaDay('/Users/erikfred/Google Drive/My Drive/Oceanography/SCTA-Share/OOI-PF/SCTA-PF/ParsedData',dayn);
        
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

% identify and separate each of the calibrations
inta=find(flipInfoAll.t<=datenum(2019,08,09));
intb=find(flipInfoAll.t>datenum(2019,08,09));
i_x=find(flipInfoAll.orientation==1);
i_x1=i_x(1:2:end); i_x1a=intersect(i_x1,inta); i_x1b=intersect(i_x1,intb);
i_x2=i_x(2:2:end); i_x2a=intersect(i_x2,inta); i_x2b=intersect(i_x2,intb);
i_y=find(flipInfoAll.orientation==2);
i_ya=intersect(i_y,inta); i_yb=intersect(i_y,intb);
i_negx=find(flipInfoAll.orientation==-1);
i_negy=find(flipInfoAll.orientation==-2);

% Plot calibrations
figure; clf; hold on
plot(flipInfoAll.t(i_x1),flipInfoAll.gCal(i_x1),'ok',flipInfoAll.t(i_x1),flipInfoAll.gCalTCor(i_x1),'xk','markersize',18);
plot(flipInfoAll.t(i_x2),flipInfoAll.gCal(i_x2),'ob',flipInfoAll.t(i_x2),flipInfoAll.gCalTCor(i_x2),'xb','markersize',18);
plot(flipInfoAll.t(i_y),flipInfoAll.gCal(i_y),'or',flipInfoAll.t(i_y),flipInfoAll.gCalTCor(i_y),'xr','markersize',18);
plot(flipInfoAll.t(i_negx),flipInfoAll.gCal(i_negx),'sk',flipInfoAll.t(i_negx),flipInfoAll.gCalTCor(i_negx),'+k','markersize',18);
plot(flipInfoAll.t(i_negy),flipInfoAll.gCal(i_negy),'sr',flipInfoAll.t(i_negy),flipInfoAll.gCalTCor(i_negy),'+r','markersize',18);
plot(flipInfoAll.t(i_x1b),(flipInfoAll.gCal(i_x1b)+flipInfoAll.gCal(i_negx))/2,'^k',...
    flipInfoAll.t(i_x1b),(flipInfoAll.gCalTCor(i_x1b)+flipInfoAll.gCalTCor(i_negx))/2,'+k','markersize',18);
plot(flipInfoAll.t(i_x2b),(flipInfoAll.gCal(i_x2b)+flipInfoAll.gCal(i_negx))/2,'^b',...
    flipInfoAll.t(i_x2b),(flipInfoAll.gCalTCor(i_x2b)+flipInfoAll.gCalTCor(i_negx))/2,'+b','markersize',18);
plot(flipInfoAll.t(i_yb),(flipInfoAll.gCal(i_yb)+flipInfoAll.gCal(i_negy))/2,'^r',...
    flipInfoAll.t(i_yb),(flipInfoAll.gCalTCor(i_yb)+flipInfoAll.gCalTCor(i_negy))/2,'+r','markersize',18);
xl = xlim; yl = ylim;
plot([0 0]+xl(1)+diff(xl)/10,mean(yl)+[0 0.0001],'-k')
text(xl(1)+diff(xl)/9,mean(yl)+0.00005,'10^{-5} g')
legend('1st X','1st X (T Corrected)','2nd X','2nd X (T Corrected)','Y','Y (T Corrected)','-X','-X (T Corrected)',...
    '-Y','-Y (T Corrected)','X1 span','X1 span (T corrected)','X2 span','X2 span (T corrected)','Y span','Y span (T corrected)','location','best')
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

% Plot calibrations on separate axes
figure
subplot(311); hold on
plot(flipInfoAll.t(i_x1),flipInfoAll.gCal(i_x1),'ko','markersize',12);
plot(flipInfoAll.t(i_x2),flipInfoAll.gCal(i_x2),'ro','markersize',12);
datetick('x',3)
xtickangle(45)
ylabel('+X (m/s^2)')
title({'Pinon Flat SCTA X Calibrations',[datestr(startDate,'mmm dd, yyyy') ' - ' datestr(endDate,'mmm dd, yyyy')]})
set(gca,'fontsize',12)
ylim([9.7927 9.7942])
xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.0005 -0.0005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'100 \mug','fontsize',12)
yyaxis right
plot(dataDec100.t(dataDec100.t>=flipInfoAll.t(i_x1(1)) & dataDec100.t<=flipInfoAll.t(i_x1(end))),...
    dataDec100.T(dataDec100.t>=flipInfoAll.t(i_x1(1)) & dataDec100.t<=flipInfoAll.t(i_x1(end))),'c','linewidth',1)
xlim(xl)
datetick('x',3,'keeplimits')
xtickangle(45)
ylabel(['T (' char(176) 'C)'])
legend('X1','X2','location','northwest')
set(gca,'fontsize',12)
set(gca,'YColor','c')
ylim([23.5 31])
box on; grid on
subplot(312); hold on
plot(flipInfoAll.t(i_negx),-1*flipInfoAll.gCal(i_negx),'ks','markersize',12);
xlim(xl)
datetick('x',3,'keeplimits')
xtickangle(45)
ylabel('-X (m/s^2)')
set(gca,'fontsize',12)
ylim([-9.7933 -9.7918])
xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.0005 -0.0005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'100 \mug','fontsize',12)
yyaxis right
plot(dataDec100.t(dataDec100.t>=flipInfoAll.t(i_negx(1)) & dataDec100.t<=flipInfoAll.t(i_negx(end))),...
    dataDec100.T(dataDec100.t>=flipInfoAll.t(i_negx(1)) & dataDec100.t<=flipInfoAll.t(i_negx(end))),'c','linewidth',1)
xlim(xl)
datetick('x',3,'keeplimits')
xtickangle(45)
ylabel(['T (' char(176) 'C)'])
set(gca,'fontsize',12)
set(gca,'YColor','c')
ylim([23.5 31])
box on; grid on
subplot(313); hold on
xspan1=(flipInfoAll.gCal(i_x1b)+flipInfoAll.gCal(i_negx));
plot(flipInfoAll.t(i_x1b),xspan1-xspan1(1),'k^','markersize',12);
xspan2=(flipInfoAll.gCal(i_x2b)+flipInfoAll.gCal(i_negx));
plot(flipInfoAll.t(i_x2b),xspan2-xspan2(1),'r^','markersize',12);
xlim(xl)
datetick('x',3,'keeplimits')
xtickangle(45)
ylabel('\Delta X1 span (m/s^2)')
set(gca,'fontsize',12)
ylim([-0.00005 0.0001])
xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.00005 -0.00005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'10 \mug','fontsize',12)
yyaxis right
plot(dataDec100.t(dataDec100.t>=flipInfoAll.t(i_negx(1)) & dataDec100.t<=flipInfoAll.t(i_negx(end))),...
    dataDec100.T(dataDec100.t>=flipInfoAll.t(i_negx(1)) & dataDec100.t<=flipInfoAll.t(i_negx(end))),'c','linewidth',1)
xlim(xl)
datetick('x',3,'keeplimits')
xtickangle(45)
ylabel(['T (' char(176) 'C)'])
legend('X1','X2','location','northwest')
set(gca,'fontsize',12)
set(gca,'YColor','c')
ylim([27.5 30])
box on; grid on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print('../calibrations/PinonFlat/process_PinonFlat_x_alt2','-dtiff','-r300')

figure
subplot(311); hold on
plot(flipInfoAll.t(i_y),flipInfoAll.gCal(i_y),'ko','markersize',12);
title({'Pinon Flat SCTA Y Calibrations',[datestr(startDate,'mmm dd, yyyy') ' - ' datestr(endDate,'mmm dd, yyyy')]})
datetick('x',3)
xtickangle(45)
ylabel('+Y (m/s^2)')
set(gca,'fontsize',12)
ylim([9.7938 9.7953])
xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.0005 -0.0005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'100 \mug','fontsize',12)
yyaxis right
plot(dataDec100.t(dataDec100.t>=flipInfoAll.t(i_y(1)) & dataDec100.t<=flipInfoAll.t(i_y(end))),...
    dataDec100.T(dataDec100.t>=flipInfoAll.t(i_y(1)) & dataDec100.t<=flipInfoAll.t(i_y(end))),'c','linewidth',1)
xlim(xl)
datetick('x',3,'keeplimits')
xtickangle(45)
ylabel(['T (' char(176) 'C)'])
set(gca,'fontsize',12)
set(gca,'YColor','c')
ylim([23.5 31])
box on; grid on
subplot(312); hold on
plot(flipInfoAll.t(i_negy),-1*flipInfoAll.gCal(i_negy),'ks','markersize',12);
xlim(xl)
datetick('x',3,'keeplimits')
xtickangle(45)
ylabel('-Y (m/s^2)')
set(gca,'fontsize',12)
ylim([-9.7920 -9.7905])
xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.0005 -0.0005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'100 \mug','fontsize',12)
yyaxis right
plot(dataDec100.t(dataDec100.t>=flipInfoAll.t(i_negy(1)) & dataDec100.t<=flipInfoAll.t(i_negy(end))),...
    dataDec100.T(dataDec100.t>=flipInfoAll.t(i_negy(1)) & dataDec100.t<=flipInfoAll.t(i_negy(end))),'c','linewidth',1)
xlim(xl)
datetick('x',3,'keeplimits')
xtickangle(45)
ylabel(['T (' char(176) 'C)'])
set(gca,'fontsize',12)
set(gca,'YColor','c')
ylim([23.5 31])
box on; grid on
subplot(313); hold on
yspan=(flipInfoAll.gCal(i_yb)+flipInfoAll.gCal(i_negy));
plot(flipInfoAll.t(i_yb),yspan-yspan(1),'k^','markersize',12);
xlim(xl)
datetick('x',3,'keeplimits')
xtickangle(45)
ylabel('\Delta Y span (m/s^2)')
set(gca,'fontsize',12)
ylim([-0.00005 0.0001])
xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.00005 -0.00005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'10 \mug','fontsize',12)
yyaxis right
plot(dataDec100.t(dataDec100.t>=flipInfoAll.t(i_negy(1)) & dataDec100.t<=flipInfoAll.t(i_negy(end))),...
    dataDec100.T(dataDec100.t>=flipInfoAll.t(i_negy(1)) & dataDec100.t<=flipInfoAll.t(i_negy(end))),'c','linewidth',1)
xlim(xl)
datetick('x',3,'keeplimits')
xtickangle(45)
ylabel(['T (' char(176) 'C)'])
set(gca,'fontsize',12)
set(gca,'YColor','c')
ylim([27.5 30])
box on; grid on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print('../calibrations/PinonFlat/process_PinonFlat_y_alt2','-dtiff','-r300')