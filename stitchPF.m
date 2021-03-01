% stitchPF.m
%
% Stitches together inter-calibration segments from decimated 40Hz data
% stream, removing offsets and applying calibrations.
%
% Future fixes:
%   - comparison against tilt model and/or collocated sensor
%   - troubleshoot hourly version
%   - options to remove offsets from glitches/fishbumps/etc.
%

clear; close all

%%%%% CONFIG %%%%%

load('../calibrations/PinonFlat/PFdata.mat','flipInfoAll')
load('../calibrations/PinonFlat/PFdata_hr.mat')
load('../calibrations/PinonFlat/PFdata_min.mat')

load('../calibrations/PinonFlat/T_dependence/X_dT')
load('../calibrations/PinonFlat/T_dependence/Y_dT')

ix=[1:3:58,61:5:length(flipInfoAll.t)]';
xcals=flipInfoAll.gCal(ix);
iy=[2:3:59,62:5:length(flipInfoAll.t)]';
ycals=flipInfoAll.gCal(iy);

lindrift=false; % if false, calculates drift independently each interval
interpolate=0; % 1 to interpolate between calibrations, 0 to extrapolate data out

stitch_method=1; % 0: interpolation and best fit, 1: inversion

%%% END CONFIG %%%

% pick up data structure from previous run, or start from scratch
if exist('../calibrations/PinonFlat/PFstitch_hr.mat','file') && ...
        exist('../calibrations/PinonFlat/PFstitch_min.mat','file')
    load('../calibrations/PinonFlat/PFstitch_hr.mat')
    load('../calibrations/PinonFlat/PFstitch_min.mat')
    
    i0=length(stitch_hr)+1;
    
    stitch_hr.t=[stitch_hr.t; data_hr.t(i0:end)];
    stitch_hr.MNE=[stitch_hr.MNE; data_hr.MNE(i0:end)];
    stitch_hr.MNN=[stitch_hr.MNN; data_hr.MNN(i0:end)];
    stitch_hr.MNZ=[stitch_hr.MNZ; data_hr.MNZ(i0:end)];
    stitch_hr.MKA=[stitch_hr.MKA; data_hr.MKA(i0:end)];
    stitch_hr.iflip=[stitch_hr.iflip; data_hr.iflip(i0:end)];
    
    stitch_min.t=[stitch_min.t; data_min.t(i0:end)];
    stitch_min.MNE=[stitch_min.MNE; data_min.MNE(i0:end)];
    stitch_min.MNN=[stitch_min.MNN; data_min.MNN(i0:end)];
    stitch_min.MNZ=[stitch_min.MNZ; data_min.MNZ(i0:end)];
    stitch_min.MKA=[stitch_min.MKA; data_min.MKA(i0:end)];
    stitch_min.iflip=[stitch_min.iflip; data_min.iflip(i0:end)];
else
    i0=1;
    
    stitch_hr=data_hr;
    stitch_min=data_min;
    
    % remove interval when communication was lost and data recording stop
    [~,i1]=min(abs(datenum(2019,02,14)-stitch_min.t));
    [~,i2]=min(abs(datenum(2019,04,03)-stitch_min.t));
    stitch_min.t(i1:i2)=NaN;
    stitch_min.MNE(i1:i2)=NaN;
    stitch_min.MNN(i1:i2)=NaN;
    stitch_min.MNZ(i1:i2)=NaN;
    stitch_min.MKA(i1:i2)=NaN;
end

%% 1 sample/minute version
% identify indices of intervals
iflipstart_min=ischange(stitch_min.t(stitch_min.iflip)); iflipstart_min(1)=true;
flipstart_min=stitch_min.iflip(iflipstart_min);
flipstart_min=[flipstart_min(1:15);178253;flipstart_min(16:end)]; % manually add segment start after connection loss
flipdate_min=stitch_min.t(flipstart_min); %datenums of flips for reference

iend=length(stitch_min.t);

cal_log=[];
transients=[];
for i=1:length(flipstart_min)
    if i~=16
        %specify samples to exclude around flip
        ipre=4;
        ipost=539;
        
        transients(i).t=stitch_min.t(flipstart_min(i):flipstart_min(i)+ipost);
        transients(i).x=stitch_min.MNE(flipstart_min(i):flipstart_min(i)+ipost);
        transients(i).y=stitch_min.MNN(flipstart_min(i):flipstart_min(i)+ipost);
        
        %substitute points during calibration and recovery with NaNs
        inan=flipstart_min(i)-ipre:flipstart_min(i)+ipost;
        stitch_min.MNE(inan)=nan;
        stitch_min.MNN(inan)=nan;
        stitch_min.MNZ(inan)=nan;
        stitch_min.MKA(inan)=nan;
    end
    
    if i==1
        eint=stitch_min.MNE(1:flipstart_min(i)-(ipre+1));
        nint=stitch_min.MNN(1:flipstart_min(i)-(ipre+1));
        Tint=stitch_min.MKA(1:flipstart_min(i)-(ipre+1));
        tint=stitch_min.t(1:flipstart_min(i)-(ipre+1));
    elseif i==16
        %substitute limited intervening points with NaNs
        inan=171395:178252;
        stitch_min.MNE(inan)=nan;
        stitch_min.MNN(inan)=nan;
        stitch_min.MNZ(inan)=nan;
        stitch_min.MKA(inan)=nan;
        
        eint=stitch_min.MNE(flipstart_min(i-1)+(ipost+1):171394);
        nint=stitch_min.MNN(flipstart_min(i-1)+(ipost+1):171394);
        Tint=stitch_min.MKA(flipstart_min(i-1)+(ipost+1):171394);
        tint=stitch_min.t(flipstart_min(i-1)+(ipost+1):171394);
        
        transients(i).t=zeros(size(transients(i-1).t));
        transients(i).x=zeros(size(transients(i-1).t));
        transients(i).y=zeros(size(transients(i-1).t));
    elseif i==17
        eint=stitch_min.MNE(178253:flipstart_min(i)-(ipre+1));
        nint=stitch_min.MNN(178253:flipstart_min(i)-(ipre+1));
        Tint=stitch_min.MKA(178253:flipstart_min(i)-(ipre+1));
        tint=stitch_min.t(178253:flipstart_min(i)-(ipre+1));
    else
        eint=stitch_min.MNE(flipstart_min(i-1)+(ipost+1):flipstart_min(i)-(ipre+1));
        nint=stitch_min.MNN(flipstart_min(i-1)+(ipost+1):flipstart_min(i)-(ipre+1));
        Tint=stitch_min.MKA(flipstart_min(i-1)+(ipost+1):flipstart_min(i)-(ipre+1));
        tint=stitch_min.t(flipstart_min(i-1)+(ipost+1):flipstart_min(i)-(ipre+1));
    end
    
    if lindrift
        xdrift=m_X(1)/24/60;
        ydrift=m_Y(1)/24/60;
        disp(['X drift rate = ' num2str(xdrift*60*24*365*10^5) ' \mug/yr'])
        disp(['Y drift rate = ' num2str(ydrift*60*24*365*10^5) ' \mug/yr'])
    elseif interpolate==1
        
        % pull calibration value, special handling of certain segments
        if i==1 % assume no drift for first seg
            cal2.x_plus=xcals(i);
            cal2.y_plus=ycals(i);
            cal1=cal2;
            continue
        elseif i>16 % i=16 is a gap, not a calibration
            cal2.x_plus=xcals(i-1);
            cal2.y_plus=ycals(i-1);
        else
            cal2.x_plus=xcals(i);
            cal2.y_plus=ycals(i);
        end
        
        % calculate drift
        if i==16 || i==17 % i=16 is a gap, not a calibration
            xdrift=(cal2.x_plus-cal1.x_plus)/(flipdate_min(17)-flipdate_min(15))/24/60; %per minute
            ydrift=(cal2.y_plus-cal1.y_plus)/(flipdate_min(17)-flipdate_min(15))/24/60; %per minute
            disp(['X drift rate = ' num2str(xdrift*60*24*365*10^5) ' \mug/yr'])
            disp(['Y drift rate = ' num2str(ydrift*60*24*365*10^5) ' \mug/yr'])
        else
            xdrift=(cal2.x_plus-cal1.x_plus)/(flipdate_min(i)-flipdate_min(i-1))/24/60; %per minute
            ydrift=(cal2.y_plus-cal1.y_plus)/(flipdate_min(i)-flipdate_min(i-1))/24/60; %per minute
            disp(['X drift rate = ' num2str(xdrift*60*24*365*10^5) ' \mug/yr'])
            disp(['Y drift rate = ' num2str(ydrift*60*24*365*10^5) ' \mug/yr'])
        end
        
        % current cal gets stored as previous cal for next iteration
        if i~=16 % i=16 is a gap, not a calibration
            cal1=cal2;
        end
    else
        % pull calibration value, special handling of certain segments
        if i==1 % assume no drift for first seg
            cal2.x_plus=xcals(i);
            cal2.y_plus=ycals(i);
            cal1=cal2;
            cal_log(i,1)=0; cal_log(i,2)=0;
            continue
        elseif i>16 % i=16 is a gap, not a calibration
            cal2.x_plus=xcals(i-1);
            cal2.y_plus=ycals(i-1);
        else
            cal2.x_plus=xcals(i);
            cal2.y_plus=ycals(i);
        end
        
        % calculate drift
        if i==16 || i==17 % i=16 is a gap, not a calibration
            xdrift=(cal2.x_plus-cal1.x_plus)/(flipdate_min(17)-flipdate_min(15))/24/60; %per minute
            ydrift=(cal2.y_plus-cal1.y_plus)/(flipdate_min(17)-flipdate_min(15))/24/60; %per minute
            disp(['X drift rate = ' num2str(xdrift*60*24*365*10^5) ' \mug/yr'])
            disp(['Y drift rate = ' num2str(ydrift*60*24*365*10^5) ' \mug/yr'])
            cal_log(i,1)=0; cal_log(i,2)=0;
        else
            delta_ax=0; % need to sort out how to determine this
            cal1.ax=(eint(1)+delta_ax)*(tint(end)-flipdate_min(i-1))/(tint(end)-tint(1))-...
                eint(end)*(tint(1)-flipdate_min(i-1))/(tint(end)-tint(1));
            cal2.ax=eint(end)*(flipdate_min(i)-tint(1))/(tint(end)-tint(1))-...
                (eint(1)+delta_ax)*(flipdate_min(i)-tint(end))/(tint(end)-tint(1));
            cal2.x_dif_cal=cal2.ax-cal1.ax-(cal2.x_plus-cal1.x_plus);
            
            delta_ay=0; % need to sort out how to determine this
            cal1.ay=(nint(1)+delta_ay)*(tint(end)-flipdate_min(i-1))/(tint(end)-tint(1))-...
                nint(end)*(tint(1)-flipdate_min(i-1))/(tint(end)-tint(1));
            cal2.ay=nint(end)*(flipdate_min(i)-tint(1))/(tint(end)-tint(1))-...
                (nint(1)+delta_ay)*(flipdate_min(i)-tint(end))/(tint(end)-tint(1));
            cal2.y_dif_cal=cal2.ay-cal1.ay-(cal2.y_plus-cal1.y_plus);
            
            xdrift=(cal2.x_plus-cal1.x_plus)/(flipdate_min(i)-flipdate_min(i-1))/24/60; %per minute
            ydrift=(cal2.y_plus-cal1.y_plus)/(flipdate_min(i)-flipdate_min(i-1))/24/60; %per minute
            disp(['X drift rate = ' num2str(xdrift*60*24*365*10^5) ' \mug/yr'])
            disp(['Y drift rate = ' num2str(ydrift*60*24*365*10^5) ' \mug/yr'])
            
            figure(12)
            plot(tint,eint,'b','linewidth',1)
            hold on
            plot([flipdate_min(i-1) flipdate_min(i)],[cal1.ax cal2.ax],'bo:')
            datetick
            title('X segments')
            figure(13)
            plot(tint,nint,'b','linewidth',1)
            hold on
            plot([flipdate_min(i-1) flipdate_min(i)],[cal1.ay cal2.ay],'bo:')
            datetick
            title('Y segments')
            
            cal_log(i,1)=cal2.x_dif_cal; cal_log(i,2)=cal2.y_dif_cal;
        end
        
        % current cal gets stored as previous cal for next iteration
        if i~=16 % i=16 is a gap, not a calibration
            cal1=cal2;
        end
    end
    
    if stitch_method==0
        % apply drift and temperature correction to segment
        elin=linspace(0,xdrift*length(eint),length(eint))';
        nlin=linspace(0,ydrift*length(nint),length(nint))';
        etemp=Tint*m_X(5);%2.363636e-5;
        ntemp=Tint*m_Y(3);
        eint_cal=eint-elin;
        nint_cal=nint-nlin;
        
        % plots to confirm calibrations are working
        %     if i>1
        %         figure(44); clf;
        %         subplot(211); hold on
        %         plot((1:length(eint))/60/24,eint-nanmean(eint),'b','linewidth',1)
        %         plot((1:length(eint))/60/24,(eint-elin)-nanmean(eint-elin),'k','linewidth',1)
        %         legend('Uncorrected','Corrected','location','northeast')
        %         ylabel('Acceleration (m/s^2)')
        %         xlabel('Time (d)')
        %         title(['X-tilt ' datestr(flipstart_min(i),'dd-mmm')])
        % %         yyaxis right
        % %         plot((1:length(eint))/60/24,Tint)
        %         set(gca,'fontsize',14)
        %         box on
        %         subplot(212); hold on
        %         plot(flipdate_min,xcals,'o')
        %         plot(flipdate_min(i-1:i),xcals(i-1:i),'ok','markerfacecolor','k')
        %         text(flipdate_min(end-10),mean(xcals),['cal dif = ' num2str(diff(xcals(i-1:i)))],'fontsize',12)
        %         ylabel('Acceleration (m/s^2)')
        %         datetick('x',6,'keeplimits')
        %         set(gca,'fontsize',14)
        %         box on
        %         fh=gcf;
        %         fh.PaperUnits='inches';
        %         fh.PaperPosition=[0 0 11 8.5];
        %         print(['../longterm_tilt/PinonFlat/correction_tests/x' num2str(i)],'-djpeg','-r100')
        %
        %
        %         figure(44); clf;
        %         subplot(211); hold on
        %         plot((1:length(nint))/60/24,nint-nanmean(nint),'b','linewidth',1)
        %         plot((1:length(nint))/60/24,(nint-nlin)-nanmean(nint-nlin),'k','linewidth',1)
        %         legend('Uncorrected','Corrected','location','northeast')
        %         ylabel('Acceleration (m/s^2)')
        %         xlabel('Time (d)')
        %         title(['Y-tilt ' datestr(flipstart_min(i),'dd-mmm')])
        % %         yyaxis right
        % %         plot((1:length(eint))/60/24,Tint)
        %         set(gca,'fontsize',16)
        %         box on
        %         subplot(212); hold on
        %         plot(flipdate_min,ycals,'o')
        %         plot(flipdate_min(i-1:i),ycals(i-1:i),'ok','markerfacecolor','k')
        %         text(flipdate_min(end-10),mean(ycals),['cal dif = ' num2str(diff(ycals(i-1:i)))],'fontsize',12)
        %         ylabel('Acceleration (m/s^2)')
        %         datetick('x',6,'keeplimits')
        %         set(gca,'fontsize',14)
        %         box on
        %         fh=gcf;
        %         fh.PaperUnits='inches';
        %         fh.PaperPosition=[0 0 11 8.5];
        %         print(['../longterm_tilt/PinonFlat/correction_tests/y' num2str(i)],'-djpeg','-r100')
        %     end
        
        %minimize offset between segments by aligning end and start of segments
        if i==1
            eint_cor=eint_cal;
            nint_cor=nint_cal;
            
            %replace data segment in 'stitch_min'
            stitch_min.MNE(1:flipstart_min(i)-(ipre+1))=eint_cor;
            stitch_min.MNN(1:flipstart_min(i)-(ipre+1))=nint_cor;
        elseif i==17 % accounts for non-flip at i=16
            %use end segment from previous interval to calculate linear
            %trend, fit current segment to that trend
            lfit=5000; %length of data on either side of flip to be aligned
            ltime=(stitch_min.t(178254)-stitch_min.t(171394))*24*60;
            
            p_e1=polyfit([-lfit:0]',stitch_min.MNE(171394-lfit:171394),1);
            lin_e2=polyval(p_e1,[ltime+1:ltime+lfit]');
            e2_dif=lin_e2-eint_cal(1:lfit);
            Ge=ones(size(e2_dif));
            e_offset=inv(Ge'*Ge)*Ge'*e2_dif;
            
            p_n1=polyfit([-lfit:0]',stitch_min.MNN(171394-lfit:171394),1);
            lin_n2=polyval(p_n1,[ltime+1:ltime+lfit]');
            n2_dif=lin_n2-nint_cal(1:lfit);
            Gn=ones(size(n2_dif));
            n_offset=inv(Ge'*Ge)*Ge'*n2_dif;
            
            eint_cor=eint_cal+e_offset;
            nint_cor=nint_cal+n_offset;
            
            stitch_min.MNE(flipstart_min(i-1):flipstart_min(i)-(ipre+1))=eint_cor;
            stitch_min.MNN(flipstart_min(i-1):flipstart_min(i)-(ipre+1))=nint_cor;
        else
            %use end segment from previous interval to calculate linear
            %trend, fit current segment to that trend
            lfit=870; %length of data on either side of flip to be aligned
            
            p_e1=polyfit([flipstart_min(i-1)-(lfit+ipre):flipstart_min(i-1)-(ipre+1)]'-flipstart_min(i-1),stitch_min.MNE(flipstart_min(i-1)-(lfit+ipre):flipstart_min(i-1)-(ipre+1)),1);
            lin_e2=polyval(p_e1,[flipstart_min(i-1)+(ipost+1):flipstart_min(i-1)+(lfit+ipost)]'-flipstart_min(i-1));
            e2_dif=lin_e2-eint_cal(1:lfit);
            Ge=ones(size(e2_dif));
            e_offset=inv(Ge'*Ge)*Ge'*e2_dif;
            
            p_n1=polyfit([flipstart_min(i-1)-(lfit+ipre):flipstart_min(i-1)-(ipre+1)]'-flipstart_min(i-1),stitch_min.MNN(flipstart_min(i-1)-(lfit+ipre):flipstart_min(i-1)-(ipre+1)),1);
            lin_n2=polyval(p_n1,[flipstart_min(i-1)+(ipost+1):flipstart_min(i-1)+(lfit+ipost)]'-flipstart_min(i-1));
            n2_dif=lin_n2-nint_cal(1:lfit);
            Gn=ones(size(n2_dif));
            n_offset=inv(Ge'*Ge)*Ge'*n2_dif;
            
            eint_cor=eint_cal+e_offset;
            nint_cor=nint_cal+n_offset;
            
            %replace data segment in 'stitch_min'
            if i==16
                stitch_min.MNE(flipstart_min(i-1)+(ipost+1):171394)=eint_cor;
                stitch_min.MNN(flipstart_min(i-1)+(ipost+1):171394)=nint_cor;
            else
                stitch_min.MNE(flipstart_min(i-1)+(ipost+1):flipstart_min(i)-(ipre+1))=eint_cor;
                stitch_min.MNN(flipstart_min(i-1)+(ipost+1):flipstart_min(i)-(ipre+1))=nint_cor;
            end
        end
        
        segs(i).e=eint_cor;
        segs(i).n=nint_cor;
        segs(i).t=tint;
        segs(i).T=Tint;
        
    elseif stitch_method==1
        % apply drift correction to segment
        elin=linspace(0,xdrift*length(eint),length(eint))';
        nlin=linspace(0,ydrift*length(nint),length(nint))';
        eint_cal=eint-elin;
        nint_cal=nint-nlin;
%         eint_cal=eint-mean(elin);
%         nint_cal=nint-mean(nlin);
        eint_cor=eint_cal;
        nint_cor=nint_cal;
        
        segs(i).e=eint_cal;
        segs(i).n=nint_cal;
        segs(i).t=tint;
        segs(i).T=Tint;
    end
    
    %% piecewise plots to differentiate segments
    if i==1
        figure(45); clf; hold on
        figure(46); clf; hold on
    end
    figure(45); hold on
    plot(tint,eint_cor,'linewidth',1)
    if i==length(flipstart_min)
        ylabel('Acceleration (m/s^2)')
        xlabel('Time (d)')
        title('X-tilt ')
        datetick('x',6,'keeplimits')
%         yyaxis right
%         plot((1:length(eint))/60/24,Tint)
        set(gca,'fontsize',14)
        box on
        fh=gcf;
        fh.PaperUnits='inches';
        fh.PaperPosition=[0 0 11 8.5];
%         print('../longterm_tilt/PinonFlat/correction_tests/x_seg','-djpeg','-r100')
%         saveas(gcf,'../longterm_tilt/PinonFlat/correction_tests/x_seg.fig')
    end
    
    figure(46); hold on
    plot(tint,nint_cor,'linewidth',1)
    if i==length(flipstart_min)
        ylabel('Acceleration (m/s^2)')
        xlabel('Time (d)')
        title('Y-tilt')
        datetick('x',6,'keeplimits')
%         yyaxis right
%         plot((1:length(eint))/60/24,Tint)
        set(gca,'fontsize',14)
        box on
        fh=gcf;
        fh.PaperUnits='inches';
        fh.PaperPosition=[0 0 11 8.5];
%         print('../longterm_tilt/PinonFlat/correction_tests/y_seg','-djpeg','-r100')
%         saveas(gcf,'../longterm_tilt/PinonFlat/correction_tests/y_seg.fig')
    end
end

% tack on remaining data with drift from previous interval
eint=stitch_min.MNE(flipstart_min(i)+ipost+1:end);
nint=stitch_min.MNN(flipstart_min(i)+ipost+1:end);
tint=stitch_min.t(flipstart_min(i)+ipost+1:end);
Tint=stitch_min.MKA(flipstart_min(i)+ipost+1:end);
       
elin=linspace(0,xdrift*length(eint),length(eint))';
nlin=linspace(0,ydrift*length(nint),length(nint))';
eint_cal=eint-elin;
nint_cal=nint-nlin;

if stitch_method==0
    p_e1=polyfit([flipstart_min(i)-(lfit+ipre):flipstart_min(i)-(ipre+1)]'-flipstart_min(i),stitch_min.MNE(flipstart_min(i)-(lfit+ipre):flipstart_min(i)-(ipre+1)),1);
    lin_e2=polyval(p_e1,[flipstart_min(i)+(ipost+1):flipstart_min(i)+(lfit+ipost)]'-flipstart_min(i));
    e2_dif=lin_e2-eint_cal(1:lfit);
    Ge=ones(size(e2_dif));
    e_offset=inv(Ge'*Ge)*Ge'*e2_dif;
    
    p_n1=polyfit([flipstart_min(i)-(lfit+ipre):flipstart_min(i)-(ipre+1)]'-flipstart_min(i),stitch_min.MNN(flipstart_min(i)-(lfit+ipre):flipstart_min(i)-(ipre+1)),1);
    lin_n2=polyval(p_n1,[flipstart_min(i)+(ipost+1):flipstart_min(i)+(lfit+ipost)]'-flipstart_min(i));
    n2_dif=lin_n2-nint_cal(1:lfit);
    Gn=ones(size(n2_dif));
    n_offset=inv(Ge'*Ge)*Ge'*n2_dif;
    
    eint_cor=eint_cal+e_offset;
    nint_cor=nint_cal+n_offset;
    
    stitch_min.MNE(flipstart_min(i)+ipost+1:end)=eint_cor;
    stitch_min.MNN(flipstart_min(i)+ipost+1:end)=nint_cor;
    
    segs(i+1).e=eint_cor;
    segs(i+1).n=nint_cor;
    segs(i+1).t=tint;
    segs(i+1).T=Tint;
else
    segs(i+1).e=eint_cal;
    segs(i+1).n=nint_cal;
    segs(i+1).t=tint;
    segs(i+1).T=Tint;
end

%% manual stitching

if stitch_method==0
    for i=2:length(segs)
        inan=isnan(segs(i).e);
        segs(i).t(inan)=[];
        segs(i).T(inan)=[];
        segs(i).e(inan)=[];
        segs(i).n(inan)=[];
    end
    
    % Inversion matrix
    Tp1=cat(1,segs(2:16).T);
    Tp2=cat(1,segs(17:end).T);
    G=[[Tp1;Tp2],...
        [ones(length(Tp1),1);zeros(length(Tp2),1)],...
        [zeros(length(Tp1),1);ones(length(Tp2),1)]];
    
    m_east=inv(G'*G)*G'*cat(1,segs(2:end).e);
    e_star=G*m_east;
    
    m_north=inv(G'*G)*G'*cat(1,segs(2:end).n);
    n_star=G*m_north;
    
    % plotting
    figure(6); clf
    i2=0;
    for j=2:length(segs)
        subplot(211); hold on
        plot(segs(j).t,segs(j).e,'linewidth',1)
        datetick('x')
        ylabel('X (m/s^2)')
        title('Observed and predicted acceleration')
        set(gca,'fontsize',14)
        box on
        subplot(212); hold on
        i1=1+i2;
        i2=length(cat(1,segs(2:j).t));
        edif=segs(j).e-e_star(i1:i2);
        plot(segs(j).t,edif,'linewidth',1)
        datetick('x')
        ylabel('X (m/s^2)')
        title('Misfit')
        set(gca,'fontsize',14)
        box on
    end
    subplot(211)
    plot(cat(1,segs(2:end).t),e_star,'.k')
    text(datenum(2018,11,01),4.4e-3,['dX/dT = ' num2str(m_east(1)) ' (m/s^2)/C'],'fontsize',12)
    
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print('../longterm_tilt/PinonFlat/correction_tests/man_corrected_X','-dtiff','-r300')
    
    figure(7); clf
    i2=0;
    for j=2:length(segs)
        subplot(211); hold on
        plot(segs(j).t,segs(j).n,'linewidth',1)
        datetick('x')
        ylabel('Y (m/s^2)')
        title('Observed and predicted acceleration')
        set(gca,'fontsize',14)
        box on
        subplot(212); hold on
        i1=1+i2;
        i2=length(cat(1,segs(2:j).t));
        ndif=segs(j).n-n_star(i1:i2);
        plot(segs(j).t,ndif,'linewidth',1)
        datetick('x')
        ylabel('Y (m/s^2)')
        title('Misfit')
        set(gca,'fontsize',14)
        box on
    end
    subplot(211)
    plot(cat(1,segs(2:end).t),n_star,'.k')
    text(datenum(2018,11,01),-0.0386,['dY/dT = ' num2str(m_north(1)) ' (m/s^2)/C'],'fontsize',12)
    
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print('../longterm_tilt/PinonFlat/correction_tests/man_corrected_Y','-dtiff','-r300')

%% inversion stitching
else
    % set up inversion to align segmented data
    etest=cat(1,segs(2:end).e);
    ntest=cat(1,segs(2:end).n);
    ttest=cat(1,segs(2:end).t);
    Ttest=cat(1,segs(2:end).T);
    eye_ish=[];
    for k=2:length(segs)
        eye_k=zeros(length(segs(k).e),length(segs)-1);
        eye_k(:,k-1)=1;
        eye_ish=[eye_ish;eye_k];
    end
    Gtest=[Ttest,eye_ish];
    mtest_e=(Gtest'*Gtest)\Gtest'*etest;
    etest_star=Gtest*mtest_e;
    mtest_n=(Gtest'*Gtest)\Gtest'*ntest;
    ntest_star=Gtest*mtest_n;
    
    figure(4); clf
    for j=2:length(segs)
        subplot(212); hold on
        plot(segs(j).t,segs(j).e-mtest_e(j)-mtest_e(1)*segs(j).T,'linewidth',1)
        datetick('x')
        ylabel('X (m/s^2)')
        title('Offset- and temperature-corrected acceleration')
        set(gca,'fontsize',14)
        box on
        subplot(211); hold on
        h1(j)=plot(segs(j).t,segs(j).e-mtest_e(j),'linewidth',1);
        h2(j)=plot(segs(j).t,mtest_e(1)*segs(j).T,'k');
        datetick('x')
        ylabel('X (m/s^2)')
        title('Observed acceleration with inverted offsets removed')
        set(gca,'fontsize',14)
        box on
    end
    legend([h1(2) h2(2)],'Tilt segments','Inverted T dependence')
    text(datenum(2018,11,01),-7.5e-4,['dX/dT = ' num2str(mtest_e(1)) ' (m/s^2)/C'],'fontsize',12)
    
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print('../longterm_tilt/PinonFlat/correction_tests/inv_corrected_X','-dtiff','-r300')
    
    figure(5); clf
    for j=2:length(segs)
        subplot(212); hold on
        plot(segs(j).t,segs(j).n-mtest_n(j)-mtest_n(1)*segs(j).T,'linewidth',1)
        datetick('x')
        ylabel('Y (m/s^2)')
        title('Offset- and temperature-corrected acceleration')
        set(gca,'fontsize',14)
        box on
        subplot(211); hold on
        h1(j)=plot(segs(j).t,segs(j).n-mtest_n(j),'linewidth',1);
        h1(j)=plot(segs(j).t,mtest_n(1)*segs(j).T,'k');
        datetick('x')
        ylabel('Y (m/s^2)')
        title('Observed acceleration with inverted offsets removed')
        set(gca,'fontsize',14)
        box on
    end
    legend([h1(2) h2(2)],'Tilt segments','Inverted T dependence')
    text(datenum(2018,11,01),-2.7e-3,['dY/dT = ' num2str(mtest_n(1)) ' (m/s^2)/C'],'fontsize',12)
    
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print('../longterm_tilt/PinonFlat/correction_tests/inv_corrected_Y','-dtiff','-r300')
end
%% resume original code

keyboard;

%% Manually correct connection loss offset (come up with better solution)

p_e1=polyfit(stitch_min.t(i1-5000:i1-1)-stitch_min.t(i1-5000),stitch_min.MNE(i1-5000:i1-1),1);
lin_e2=polyval(p_e1,stitch_min.t(i2+1:i2+5000)-stitch_min.t(i1-5000));
e2_dif=lin_e2-stitch_min.MNE(i2+1:i2+5000);
Ge=ones(size(e2_dif));
e_offset=inv(Ge'*Ge)*Ge'*e2_dif;
stitch_min.MNE(i2+1:end)=stitch_min.MNE(i2+1:end)+e_offset;

p_n1=polyfit(stitch_min.t(i1-5000:i1-1)-stitch_min.t(i1-5000),stitch_min.MNN(i1-5000:i1-1),1);
lin_n2=polyval(p_n1,stitch_min.t(i2+1:i2+5000)-stitch_min.t(i1-5000));
n2_dif=lin_n2-stitch_min.MNN(i2+1:i2+5000);
Gn=ones(size(n2_dif));
n_offset=inv(Gn'*Gn)*Gn'*n2_dif;
stitch_min.MNN(i2+1:end)=stitch_min.MNN(i2+1:end)+n_offset;

%% convert accel to tilt in microrad
stitch_min.LAX=asin(stitch_min.MNE/9.81)*10^6;stitch_min.LAX=stitch_min.LAX-stitch_min.LAX(1);
stitch_min.LAY=asin(stitch_min.MNN/9.81)*10^6;stitch_min.LAY=stitch_min.LAY-stitch_min.LAY(1);
stitch_hr.LAX=asin(stitch_hr.MNE/9.81)*10^6;stitch_hr.LAX=stitch_hr.LAX-stitch_hr.LAX(1);
stitch_hr.LAY=asin(stitch_hr.MNN/9.81)*10^6;stitch_hr.LAY=stitch_hr.LAY-stitch_hr.LAY(1);

%save figures, variables
% save('../calibrations/Axial/axialstitch_hr.mat','stitch_hr')
save('../calibrations/PinonFlat/PFstitch_min_temp2.mat','stitch_min')