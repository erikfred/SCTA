% Script to process Axial Seamount Flips
% 
% This version processes all days and creates a continuous record of flips

%% Parameters
% dataLoaded: 0 - need to load from raw files; 1 - append to existing matlab file
dataLoaded = 0;

% Start and end date
% startDate = datenum('8/06/18');
% endDate = datenum('09/10/20'); % SCTA moved 9/11/2020
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
%     % Location 1 manually determined calibrations
%     cal_list=sort([737280;737281;737282;737290;737313;737320;737327;737341;737348;...
%         737355;737362;737369;737376;737383;737390;737397;737426;737457;...
%         737485;737516;737546;737577;737607;737638;737650;737657;737664;...
%         737671;737678;737713;737730;737760;737798;737805;737812;737819;...
%         737826;737833;737840;737847;737854;737861;737868;737875;737882;...
%         737889;737896;737903;737910;737917;737945;737952;737959;737966;...
%         737973;737980;737987;737994;738001;738008;738015;738022;738029;...
%         738036;738043])';
%     for dayn = cal_list
%     % Location 2 manually determined calibrations
%     cal_list=sort([738050;738057;738064;738071;738078;738085;738092;738099;738106;...
%         738113;738120;738134;738141;738148;738155;738176;738183;738190;...
%         738197;738204;738211;738218;738225;738232;738239;738246;738253;...
%         738260;738267;738274;738281;738288;738295;738302;738309;738316;...
%         738323;738366;738367;738368;738369;738371;738372;738373;738374;...
%         738375;738376;738377;738378;738379;738380;738381;738382;738383;...
%         738385;738386;738387;738391;738392;738393;738127;738162;738370;...
%         738384;738388;738389;738390;738394])';
    for dayn = startDate:endDate
        data = get_sctaDay('/Users/erikfred/Google Drive/My Drive/Oceanography/SCTA-Share/OOI-SCTA/ParsedData',dayn);
        
        % intervals for which we don't expect data
        if dayn>datenum(2020,05,09) && dayn<=datenum(2020,06,01)
            % datastream cut during this interval due to cable issues
            fprintf(['Data outage on ' datestr(dayn) '\n\n'])
            continue
        elseif dayn>datenum(2021,01,13) && dayn<=datenum(2021,01,17)
            % datastream cut during this interval due to cable issues
            fprintf(['Data outage on ' datestr(dayn) '\n\n'])
            continue
        end
        
        % parsing issue with some early APL data
        if dayn==737282
            icut=2991877:length(data.t);
            data.t(icut)=[]; data.a(icut,:)=[]; data.as(icut)=[]; data.T(icut)=[];
            dt=data.t(2)-data.t(1);
            data.t=[data.t(1:end-1);(data.t(end):dt:dayn+datenum(0,0,0,22,0,0))'];
            dn=length(data.t)-length(data.as);
            atemp=repmat(data.a(end,:),dn+1,1);
            data.a(end:end+dn,:)=atemp;
            data.as(end:end+dn)=data.as(end);
            data.T(end:end+dn)=data.T(end);
        end
        
        % load 8 Hz data, as needed
        if isempty(data.t)
            fprintf(['No 40Hz data on ' datestr(dayn) '\n\n'])
            data=[];
            cha={'MNE','MNN','MNZ','MXG','MKA'};
            chastr={'a(:,1)','a(:,2)','a(:,3)','as','T'};
            for m=1:length(cha)
                IRIS_data_pull('AXCC2',cha{m},'--',dayn,dayn+1)
                temp=rdmseed(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayn,29) '.miniseed']);
                data.t=cat(1,temp.t);
                eval(['data.' chastr{m} '=cat(1,temp.d)/10^7;']);
            end
        else
            t_base=data.t-floor(data.t(1)); % dateless time
            t_cal=t_base(t_base>datenum(0,0,0,20,0,0) & t_base<datenum(0,0,0,22,0,0)); % hour before and after calibration
            count=round(length(t_cal)/40); % ensures calibration interval contains expected # samples
            if count<7200
                fprintf(['No 40Hz data on ' datestr(dayn) '\n\n'])
                data=[];
                cha={'MNE','MNN','MNZ','MXG','MKA'};
                chastr={'a(:,1)','a(:,2)','a(:,3)','as','T'};
                for m=1:length(cha)
                    IRIS_data_pull('AXCC2',cha{m},'--',dayn,dayn+1)
                    temp=rdmseed(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayn,29) '.miniseed']);
                    data.t=cat(1,temp.t);
                    eval(['data.' chastr{m} '=cat(1,temp.d)/10^7;']);
                end
            end
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
    end
    
    dataDec1 = NANgap_scta(dataDec1);
    dataDec100 = NANgap_scta(dataDec100);
    save(['../calibrations/Axial/axialdata' suffix],'dataDec1','dataDec100','flipInfoAll','-v7.3')
    
elseif dataLoaded==1
    load(['../calibrations/Axial/axialdata' suffix])
    startDate2=floor(dataDec1.t(end))+1;
    
    for dayn = startDate2:endDate
        data = get_sctaDay('/Users/erikfred/Google Drive/My Drive/Oceanography/SCTA-Share/OOI-SCTA/ParsedData',dayn);
        
        % intervals for which we don't expect data
        if dayn>datenum(2020,05,09) && dayn<=datenum(2020,06,01)
            % datastream cut during this interval due to cable issues
            fprintf(['Data outage on ' datestr(dayn) '\n\n'])
            continue
        elseif dayn>datenum(2021,01,13) && dayn<=datenum(2021,01,17)
            % datastream cut during this interval due to cable issues
            fprintf(['Data outage on ' datestr(dayn) '\n\n'])
            continue
        end
        
        % load 8 Hz data, as needed
        if isempty(data.t)
            fprintf(['No 40Hz data on ' datestr(dayn) '\n\n'])
            data=[];
            cha={'MNE','MNN','MNZ','MXG','MKA'};
            chastr={'a(:,1)','a(:,2)','a(:,3)','as','T'};
            for m=1:length(cha)
                IRIS_data_pull('AXCC2',cha{m},'--',dayn,dayn+1)
                temp=rdmseed(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayn,29) '.miniseed']);
                data.t=cat(1,temp.t);
                eval(['data.' chastr{m} '=cat(1,temp.d)/10^7;']);
            end
        else
            t_base=data.t-floor(data.t(1)); % dateless time
            t_cal=t_base(t_base>datenum(0,0,0,21,0,0) & t_base<datenum(0,0,0,23,0,0)); % hour before and after calibration
            count=round(length(t_cal)/40); % ensures calibration interval contains expected # samples
            if count<7200
                fprintf(['No 40Hz data on ' datestr(dayn) '\n\n'])
                data=[];
                cha={'MNE','MNN','MNZ','MXG','MKA'};
                chastr={'a(:,1)','a(:,2)','a(:,3)','as','T'};
                for m=1:length(cha)
                    IRIS_data_pull('AXCC2',cha{m},'--',dayn,dayn+1)
                    temp=rdmseed(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayn,29) '.miniseed']);
                    data.t=cat(1,temp.t);
                    eval(['data.' chastr{m} '=cat(1,temp.d)/10^7;']);
                end
            end
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
    end
    
    dataDec1 = NANgap_scta(dataDec1);
    dataDec100 = NANgap_scta(dataDec100);
    save(['../calibrations/Axial/axialdata' suffix],'dataDec1','dataDec100','flipInfoAll','-v7.3')
    
end

% identify and separate each of the calibrations
i_x=find(flipInfoAll.orientation==1);
i_xb=find(flipInfoAll.orientation==1 & flipInfoAll.t>datenum(2019,08,13));
i_xa=setxor(i_x,i_xb);
i_y=find(flipInfoAll.orientation==2);
i_yb=find(flipInfoAll.orientation==2 & flipInfoAll.t>datenum(2019,08,13));
i_ya=setxor(i_y,i_yb);
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
    i_x1=i_x([1:2:113 117:2:end]); i_x1b=i_x(115);
    i_x2=i_x(2:2:end);
    i_xneg=i_negx([1:47 49:end]); i_xnegb=i_negx(48);
    
    p_x1=polyfit(flipInfoAll.t(i_x1)-flipInfoAll.t(i_x1(1)),flipInfoAll.gCal(i_x1),1);
    m_x1=polyval(p_x1,flipInfoAll.t(i_x1)-flipInfoAll.t(i_x1(1)));
    disp(['+X1 Slope = ' num2str(p_x1(1)*365*10^5) ' ug/yr'])
    disp(['+X1 Misfit = ' num2str(std(flipInfoAll.gCal(i_x1)-m_x1)*10^5) ' ug'])
        
    p_x2=polyfit(flipInfoAll.t(i_x2)-flipInfoAll.t(i_x2(1)),flipInfoAll.gCal(i_x2),1);
    m_x2=polyval(p_x2,flipInfoAll.t(i_x2)-flipInfoAll.t(i_x2(1)));
    disp(['+X2 Slope = ' num2str(p_x2(1)*365*10^5) ' ug/yr'])
    disp(['+X2 Misfit = ' num2str(std(flipInfoAll.gCal(i_x2)-m_x2)*10^5) ' ug'])
    
    p_xneg=polyfit(flipInfoAll.t(i_xneg)-flipInfoAll.t(i_xneg(1)),flipInfoAll.gCal(i_xneg),1);
    m_xneg=polyval(p_xneg,flipInfoAll.t(i_xneg)-flipInfoAll.t(i_xneg(1)));
    disp(['-X Slope = ' num2str(p_xneg(1)*365*10^5) ' ug/yr'])
    disp(['-X Misfit = ' num2str(std(flipInfoAll.gCal(i_xneg)-m_xneg)*10^5) ' ug'])
    
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
    ylim([9.81215 9.81235])
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
    xspan1b=flipInfoAll.gCal([i_x1b,i_xnegb-4])+flipInfoAll.gCal([i_x1b+4,i_xnegb]);
    plot(flipInfoAll.t(intersect(i_x1,i_xneg-4)),xspan1-xspan1(1),'k^','markersize',12);
    datetick('x',3)
    xtickangle(45)
    ylabel('\Delta X1 span (m/s^2)')
    set(gca,'fontsize',12)
    ylim([0 0.00004])
    yyaxis right
    xspan2=(flipInfoAll.gCal(intersect(i_x2,i_xneg-1))+flipInfoAll.gCal(intersect(i_x2+1,i_xneg)));
    xspan2b=flipInfoAll.gCal(i_xnegb-1)+flipInfoAll.gCal(i_xnegb);
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
    
    tspan1=flipInfoAll.t(intersect(i_x1,i_xneg-4));
    p_xspan1=polyfit(tspan1-tspan1(1),xspan1,1);
    m_xspan1=polyval(p_xspan1,tspan1-tspan1(1));
    disp(['X1 Span Slope = ' num2str(p_xspan1(1)*365*10^5) ' ug/yr'])
    disp(['X1 Span Misfit = ' num2str(std(xspan1-m_xspan1)*10^5) ' ug'])
    
    tspan2=flipInfoAll.t(intersect(i_x2,i_xneg-1));
    p_xspan2=polyfit(tspan2-tspan2(1),xspan2,1);
    m_xspan2=polyval(p_xspan2,tspan2-tspan2(1));
    disp(['X2 Span Slope = ' num2str(p_xspan2(1)*365*10^5) ' ug/yr'])
    disp(['X2 Span Misfit = ' num2str(std(xspan2-m_xspan2)*10^5) ' ug'])
    
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print(['../calibrations/Axial/process_Axial' suffix '_manuscript_x'],'-dtiff','-r300')
    print(['../calibrations/Axial/process_Axial' suffix '_manuscript_x'],'-depsc','-painters')
    
    i_y1=i_y;
    i_yneg=i_negy([1:36 38:end]); i_ynegb=i_negy(37);
    
    p_y1=polyfit(flipInfoAll.t(i_y1)-flipInfoAll.t(i_y1(1)),flipInfoAll.gCal(i_y1),1);
    m_y1=polyval(p_y1,flipInfoAll.t(i_y1)-flipInfoAll.t(i_y1(1)));
    disp(['+Y Slope = ' num2str(p_y1(1)*365*10^5) ' ug/yr'])
    disp(['+Y Misfit = ' num2str(std(flipInfoAll.gCal(i_y1)-m_y1)*10^5) ' ug'])
    
    p_yneg=polyfit(flipInfoAll.t(i_yneg)-flipInfoAll.t(i_yneg(1)),flipInfoAll.gCal(i_yneg),1);
    m_yneg=polyval(p_yneg,flipInfoAll.t(i_yneg)-flipInfoAll.t(i_yneg(1)));
    disp(['-Y Slope = ' num2str(p_yneg(1)*365*10^5) ' ug/yr'])
    disp(['-Y Misfit = ' num2str(std(flipInfoAll.gCal(i_yneg)-m_yneg)*10^5) ' ug'])
    
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
    yspan=flipInfoAll.gCal(intersect(i_y1,i_yneg-1))+flipInfoAll.gCal(intersect(i_y1+1,i_yneg));
    yspanb=flipInfoAll.gCal(i_ynegb-1)+flipInfoAll.gCal(i_ynegb);
    plot(flipInfoAll.t(intersect(i_y1,i_yneg-1)),yspan-yspan(1),'k^','markersize',12);
    datetick('x',3)
    xtickangle(45)
    ylabel('\Delta Y span (m/s^2)')
    set(gca,'fontsize',12)
    ylim([0 0.00004])
    xline(flipInfoAll.t(i_ynegb),'k:','linewidth',2)
%     xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0 -0.00001],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'1 \mug','fontsize',12)
    box on; grid on
    
    tspan=flipInfoAll.t(intersect(i_y1,i_yneg-1));
    p_yspan=polyfit(tspan-tspan(1),yspan,1);
    m_yspan=polyval(p_yspan,tspan-tspan(1));
    disp(['Y Span Slope = ' num2str(p_yspan(1)*365*10^5) ' ug/yr'])
    disp(['Y Span Misfit = ' num2str(std(yspan-m_yspan)*10^5) ' ug'])
    
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print(['../calibrations/Axial/process_Axial' suffix '_manuscript_y'],'-dtiff','-r300')
    print(['../calibrations/Axial/process_Axial' suffix '_manuscript_y'],'-depsc','-painters')
elseif strcmp(suffix,'')
    i_x1=i_x([1:2:61,73:2:127]); i_x1b=i_x([63,65,67,69,71,129]);
    i_x2=i_x([2:2:62,68:2:96,100:2:126,130]); i_x2b=i_x([64,66,98,128]);
    i_xneg=i_negx([1:6,9:39,41]); i_xnegb=i_negx([7,8,40]);
    
    p_x1=polyfit(flipInfoAll.t(i_x1)-flipInfoAll.t(i_x1(1)),flipInfoAll.gCal(i_x1),1);
    m_x1=polyval(p_x1,flipInfoAll.t(i_x1)-flipInfoAll.t(i_x1(1)));
    disp(['+X1 Slope = ' num2str(p_x1(1)*365*10^5) ' ug/yr'])
    disp(['+X1 Misfit = ' num2str(std(flipInfoAll.gCal(i_x1)-m_x1)*10^5) ' ug'])
    
    ttest=flipInfoAll.t(i_x1)-flipInfoAll.t(i_x1(1));
    xtest=flipInfoAll.gCal(i_x1)-m_x1;
    f_x1=fit(ttest(1:31),xtest(1:31)-xtest(31),'exp1','Lower',[-1 -1],'Upper',[0 0],'StartPoint',[-0.00007 -0.01]);
    p_x1a=polyfit(ttest(1:31),flipInfoAll.gCal(i_x1(1:31))-f_x1.a*exp(f_x1.b*(ttest(1:31))),1);
    m_x1a=polyval(p_x1a,ttest(1:31));
    tempa=flipInfoAll.gCal(i_x1(1:31))-f_x1.a*exp(f_x1.b*(ttest(1:31)))-m_x1a;
    disp(['+X1a Slope = ' num2str(p_x1a(1)*365*10^5) ' ug/yr'])
    p_x1b=polyfit(ttest(32:end),flipInfoAll.gCal(i_x1(32:end))-f_x1.a*exp(f_x1.b*(ttest(32:end))),1);
    m_x1b=polyval(p_x1b,ttest(32:end));
    tempb=flipInfoAll.gCal(i_x1(32:end))-f_x1.a*exp(f_x1.b*(ttest(32:end)))-m_x1b;
    disp(['+X1b Slope = ' num2str(p_x1b(1)*365*10^5) ' ug/yr'])
    disp(['+X1 Misfit = ' num2str(std([tempa;tempb])*10^5) ' ug'])
    
    p_x2=polyfit(flipInfoAll.t(i_x2)-flipInfoAll.t(i_x2(1)),flipInfoAll.gCal(i_x2),1);
    m_x2=polyval(p_x2,flipInfoAll.t(i_x2)-flipInfoAll.t(i_x2(1)));
    disp(['+X2 Slope = ' num2str(p_x2(1)*365*10^5) ' ug/yr'])
    disp(['+X2 Misfit = ' num2str(std(flipInfoAll.gCal(i_x2)-m_x2)*10^5) ' ug'])
    
    p_xneg=polyfit(flipInfoAll.t(i_xneg)-flipInfoAll.t(i_xneg(1)),flipInfoAll.gCal(i_xneg),1);
    m_xneg=polyval(p_xneg,flipInfoAll.t(i_xneg)-flipInfoAll.t(i_xneg(1)));
    disp(['-X Slope = ' num2str(p_xneg(1)*365*10^5) ' ug/yr'])
    disp(['-X Misfit = ' num2str(std(flipInfoAll.gCal(i_xneg)-m_xneg)*10^5) ' ug'])
    
    figure
    subplot(311); hold on
    plot(flipInfoAll.t(i_x1),flipInfoAll.gCal(i_x1),'ko','markersize',12);
    datetick('x',3)
    xtickangle(45)
    ylabel('+X1 (m/s^2)')
    set(gca,'fontsize',12)
    ylim([9.8111 9.8121])
    for l=1:length(i_x1b)
        xline(flipInfoAll.t(i_x1b(l)),'k:','linewidth',1)
    end
    xl=xlim;
%     xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.00005 -0.00005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'10 \mug','fontsize',12)
    yyaxis right
    plot(flipInfoAll.t(i_x2),flipInfoAll.gCal(i_x2),'ro','markersize',12);
%     title({'Axial SCTA X Calibrations',[datestr(startDate,'mmm dd, yyyy') ' - ' datestr(endDate,'mmm dd, yyyy')]})
    datetick('x',3)
    xtickangle(45)
    ylabel('+X2 (m/s^2)')
    set(gca,'fontsize',12)
    set(gca,'YColor','r')
    ylim([9.8113 9.8123])
    for l=1:length(i_x2b)
        xline(flipInfoAll.t(i_x2b(l)),'r--','linewidth',0.5)
    end
    box on; grid on
    subplot(312); hold on
    plot(flipInfoAll.t(i_xneg),-1*flipInfoAll.gCal(i_xneg),'ks','markersize',12);
    xlim(xl)
    datetick('x',3,'keeplimits')
    xtickangle(45)
    ylabel('-X (m/s^2)')
    set(gca,'fontsize',12)
    ylim([-9.80960 -9.80860])
    for l=1:length(i_xnegb)
        xline(flipInfoAll.t(i_xnegb(l)),'k:','linewidth',1)
    end
%     xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.00005 -0.00005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'10 \mug','fontsize',12)
    box on; grid on
    subplot(313); hold on
    xspan1=flipInfoAll.gCal(intersect(i_x1,i_xneg-4))+flipInfoAll.gCal(intersect(i_x1+4,i_xneg));
    plot(flipInfoAll.t(intersect(i_x1,i_xneg-4)),xspan1-xspan1(1),'k^','markersize',12);
    xlim(xl)
    datetick('x',3,'keeplimits')
    xtickangle(45)
    ylabel('\Delta X1 span (m/s^2)')
    set(gca,'fontsize',12)
    ylim([0 0.0001])
    yyaxis right
    xspan2=(flipInfoAll.gCal(intersect(i_x2,i_xneg-1))+flipInfoAll.gCal(intersect(i_x2+1,i_xneg)));
    plot(flipInfoAll.t(intersect(i_x2,i_xneg-1)),xspan2-xspan2(1),'r^','markersize',12);
    xlim(xl)
    datetick('x',3,'keeplimits')
    xtickangle(45)
    ylabel('\Delta X2 span (m/s^2)')
    set(gca,'fontsize',12)
    set(gca,'YColor','r')
    ylim([0 0.0001])
    i_s1b=unique([i_x1b; i_xnegb-4]); i_s2b=unique([i_x2b; i_xnegb-1]);
    xspan1b=flipInfoAll.gCal(i_s1b)+flipInfoAll.gCal(i_s1b+4);
    xspan2b=flipInfoAll.gCal(i_s2b)+flipInfoAll.gCal(i_s2b+1);
    for l=1:length(i_s1b)
        xline(flipInfoAll.t(i_s1b(l)),'k:','linewidth',1)
    end
    for l=1:length(i_s2b)
        xline(flipInfoAll.t(i_s2b(l)),'r--','linewidth',0.5)
    end
%     xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)+[0 0.00001],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl)+0.00000,'1 \mug','fontsize',12)
    box on; grid on
    
    tspan1=flipInfoAll.t(intersect(i_x1,i_xneg-4));
    p_xspan1=polyfit(tspan1-tspan1(1),xspan1,1);
    m_xspan1=polyval(p_xspan1,tspan1-tspan1(1));
    disp(['X1 Span Slope = ' num2str(p_xspan1(1)*365*10^5) ' ug/yr'])
    disp(['X1 Span Misfit = ' num2str(std(xspan1-m_xspan1)*10^5) ' ug'])
    
    tspan2=flipInfoAll.t(intersect(i_x2,i_xneg-1));
    p_xspan2=polyfit(tspan2-tspan2(1),xspan2,1);
    m_xspan2=polyval(p_xspan2,tspan2-tspan2(1));
    disp(['X2 Span Slope = ' num2str(p_xspan2(1)*365*10^5) ' ug/yr'])
    disp(['X2 Span Misfit = ' num2str(std(xspan2-m_xspan2)*10^5) ' ug'])
    
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print(['../calibrations/Axial/process_Axial' suffix '_manuscript_x'],'-dtiff','-r300')
    print(['../calibrations/Axial/process_Axial' suffix '_manuscript_x'],'-depsc','-painters')
    
    i_y1=i_y([1:31,34:50,52:end]); i_y1b=i_y([32,33,51]);
    i_yneg=i_negy([1:7,12:14,16:18,20:28,30:40]); i_ynegb=i_negy([8,9,10,11,15,19,29,41]);
    
    p_y1=polyfit(flipInfoAll.t(i_y1)-flipInfoAll.t(i_y1(1)),flipInfoAll.gCal(i_y1),1);
    m_y1=polyval(p_y1,flipInfoAll.t(i_y1)-flipInfoAll.t(i_y1(1)));
    disp(['+Y Slope = ' num2str(p_y1(1)*365*10^5) ' ug/yr'])
    disp(['+Y Misfit = ' num2str(std(flipInfoAll.gCal(i_y1)-m_y1)*10^5) ' ug'])
    
    p_yneg=polyfit(flipInfoAll.t(i_yneg)-flipInfoAll.t(i_yneg(1)),flipInfoAll.gCal(i_yneg),1);
    m_yneg=polyval(p_yneg,flipInfoAll.t(i_yneg)-flipInfoAll.t(i_yneg(1)));
    disp(['-Y Slope = ' num2str(p_yneg(1)*365*10^5) ' ug/yr'])
    disp(['-Y Misfit = ' num2str(std(flipInfoAll.gCal(i_yneg)-m_yneg)*10^5) ' ug'])
    
    figure
    subplot(311); hold on
    plot(flipInfoAll.t(i_y1),flipInfoAll.gCal(i_y1),'ko','markersize',12);
%     title({'Axial SCTA Y Calibrations',[datestr(startDate,'mmm dd, yyyy') ' - ' datestr(endDate,'mmm dd, yyyy')]})
    datetick('x',3)
    xtickangle(45)
    ylabel('+Y (m/s^2)')
    set(gca,'fontsize',12)
    ylim([9.80980 9.81130])
    for l=1:length(i_y1b)
        xline(flipInfoAll.t(i_y1b(l)),'k:','linewidth',1)
    end
    xl=xlim;
%     xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.00005 -0.00005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'10 \mug','fontsize',12)
    box on; grid on
    subplot(312); hold on
    plot(flipInfoAll.t(i_yneg),-1*flipInfoAll.gCal(i_yneg),'ks','markersize',12);
    xlim(xl)
    datetick('x',3,'keeplimits')
    xtickangle(45)
    ylabel('-Y (m/s^2)')
    set(gca,'fontsize',12)
    ylim([-9.8125 -9.8110])
    for l=1:length(i_ynegb)
        xline(flipInfoAll.t(i_ynegb(l)),'k:','linewidth',1)
    end
%     xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.00005 -0.00005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'10 \mug','fontsize',12)
    box on; grid on
    subplot(313); hold on
    yspan=flipInfoAll.gCal(intersect(i_y1,i_yneg-1))+flipInfoAll.gCal(intersect(i_y1+1,i_yneg));
    plot(flipInfoAll.t(intersect(i_y1,i_yneg-1)),yspan-yspan(1),'k^','markersize',12);
    xlim(xl)
    datetick('x',3,'keeplimits')
    xtickangle(45)
    ylabel('\Delta Y span (m/s^2)')
    set(gca,'fontsize',12)
    ylim([-0.00002 0.00008])
    i_s1b=unique([i_y1b; i_ynegb-1]);
    yspanb=flipInfoAll.gCal(i_s1b)+flipInfoAll.gCal(i_s1b+1);
    for l=1:length(i_s1b)
        xline(flipInfoAll.t(i_s1b(l)),'k:','linewidth',1)
    end
%     xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.000005 -0.000005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'1 \mug','fontsize',12)
    box on; grid on
    
    tspan=flipInfoAll.t(intersect(i_y1,i_yneg-1));
    p_yspan=polyfit(tspan-tspan(1),yspan,1);
    m_yspan=polyval(p_yspan,tspan-tspan(1));
    disp(['Y Span Slope = ' num2str(p_yspan(1)*365*10^5) ' ug/yr'])
    disp(['Y Span Misfit = ' num2str(std(yspan-m_yspan)*10^5) ' ug'])
    
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print(['../calibrations/Axial/process_Axial' suffix '_manuscript_y'],'-dtiff','-r300')
    print(['../calibrations/Axial/process_Axial' suffix '_manuscript_y'],'-depsc','-painters')
end