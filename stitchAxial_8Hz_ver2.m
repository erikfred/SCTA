% stitchAxial_8Hz_ver2.m
%
% Stitches together inter-calibration segments from decimated 8Hz data
% stream, removing offsets and applying calibrations. Options to compare
% against other tilt instruments at Axial.
%
% Future fixes:
%   - all calibrations after Dec 2019 have been excluded
%   - troubleshoot hourly version
%   - options to remove offsets from glitches/fishbumps/etc.
%

clear; close all

comptilts=true;

load('../calibrations/Axial/axialdata.mat','flipInfoAll')
load('../calibrations/Axial/axialdata_hr.mat')
load('../calibrations/Axial/axialdata_min.mat')
load('../compass_directions.mat')

% pick up data structure from previous run, or start from scratch
if exist('../calibrations/Axial/axialstitch_hr.mat','file') && ...
        exist('../calibrations/Axial/axialstitch_min.mat','file')
    load('../calibrations/Axial/axialstitch_hr.mat')
    load('../calibrations/Axial/axialstitch_min.mat')
    
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
end

%% 1 sample/minute version
% identify indices of intervals
iflipstart_min=ischange(stitch_min.t(stitch_min.iflip)); iflipstart_min(1)=true;
flipstart_min=stitch_min.iflip(iflipstart_min);
flipdate_min=stitch_min.t(flipstart_min); %datenums of flips for reference

iend=length(stitch_min.t);
% figure % for in-line testing
% hold on
for i=1:length(flipstart_min)
    %substitute points during calibration and recovery with NaNs
    inan=flipstart_min(i)-4:flipstart_min(i)+119;
    stitch_min.MNE(inan)=nan;
    stitch_min.MNN(inan)=nan;
    stitch_min.MNZ(inan)=nan;
    
    %store calibration values
    icals=find(ceil(stitch_min.t(flipstart_min(i)))==ceil(flipInfoAll.t));
    if stitch_min.t(flipstart_min(i))<datenum(2019,08,13)
        cal2.x_plus=flipInfoAll.gCalTCor(icals(1));
        cal2.y_plus=flipInfoAll.gCalTCor(icals(2));
    else
        cal2.x_plus=flipInfoAll.gCalTCor(icals(1));
        cal2.y_plus=flipInfoAll.gCalTCor(icals(2));
        cal2.y_minus=flipInfoAll.gCalTCor(icals(3));
        cal2.x2_plus=flipInfoAll.gCalTCor(icals(4));
        cal2.x_minus=flipInfoAll.gCalTCor(icals(5));
    end
    
    if i==1
        cal1=cal2;
        continue
    end
    
    % apply calibration correction to interval
    eint=stitch_min.MNE(flipstart_min(i-1)+120:flipstart_min(i)-5);
    nint=stitch_min.MNN(flipstart_min(i-1)+120:flipstart_min(i)-5);
    
    if i<24 %|| i>28 %bad calibrations from Dec-Feb, use previous drift rates
        xdrift=(cal2.x_plus-cal1.x_plus)/(flipdate_min(i)-flipdate_min(i-1))/24/60; %per minute
        ydrift=(cal2.y_plus-cal1.y_plus)/(flipdate_min(i)-flipdate_min(i-1))/24/60; %per minute
    end
    disp(['X drift rate = ' num2str(xdrift*60*24*365*10^5) ' \mug/yr'])
    disp(['Y drift rate = ' num2str(ydrift*60*24*365*10^5) ' \mug/yr'])
    
    elin=linspace(0,xdrift*length(eint),length(eint))';
    nlin=linspace(0,ydrift*length(nint),length(nint))';
    
    eint_cal=eint-elin;
    nint_cal=nint-nlin;
    
    %minimize offset between segments by aligning last and first day
    p_e1=polyfit([flipstart_min(i-1)-1444:flipstart_min(i-1)-5]',stitch_min.MNE(flipstart_min(i-1)-1444:flipstart_min(i-1)-5),1);
    lin_e2=polyval(p_e1,[flipstart_min(i-1)+120:flipstart_min(i-1)+1559]');
    e2_dif=lin_e2-eint_cal(1:1440);
    Ge=ones(size(e2_dif));
    e_offset=inv(Ge'*Ge)*Ge'*e2_dif;
    
    p_n1=polyfit([flipstart_min(i-1)-1444:flipstart_min(i-1)-5]',stitch_min.MNN(flipstart_min(i-1)-1444:flipstart_min(i-1)-5),1);
    lin_n2=polyval(p_n1,[flipstart_min(i-1)+120:flipstart_min(i-1)+1559]');
    n2_dif=lin_n2-nint_cal(1:1440);
    Gn=ones(size(n2_dif));
    n_offset=inv(Ge'*Ge)*Ge'*n2_dif;
    
    eint_cor=eint_cal+e_offset;
    nint_cor=nint_cal+n_offset;
    
    %replace data segment in 'stitch_min'
    stitch_min.MNE(flipstart_min(i-1)+120:flipstart_min(i)-5)=eint_cor;
    stitch_min.MNN(flipstart_min(i-1)+120:flipstart_min(i)-5)=nint_cor;
    
    % current cal gets stored as previous cal for next iteration
%     plot(stitch_min.t,stitch_min.MNE,'linewidth',1) %for in-line testing
    cal1=cal2;
end

%tack on remaining data without detrending
stitch_min.MNE(flipstart_min(i)+120:end)=stitch_min.MNE(flipstart_min(i)+120:end)-...
    (stitch_min.MNE(flipstart_min(i)+120)-stitch_min.MNE(flipstart_min(i)-5));
stitch_min.MNN(flipstart_min(i)+120:end)=stitch_min.MNN(flipstart_min(i)+120:end)-...
    (stitch_min.MNN(flipstart_min(i)+120)-stitch_min.MNN(flipstart_min(i)-5));

% %% 1 sample/hr version
% % identify indices of intervals
% iflipstart_hr=ischange(stitch_hr.t(stitch_hr.iflip)); iflipstart_hr(1)=true;
% flipstart_hr=stitch_hr.iflip(iflipstart_hr);
% flipdate_hr=stitch_min.t(flipstart_hr); %datenums of flips for reference
% 
% iend=length(stitch_hr.t);
% figure
% hold on
% for j=1:length(flipstart_hr)
%     %substitute points during calibration and recovery with NaNs
%     inan=flipstart_hr(j)-1:flipstart_hr(j)+4;
%     stitch_hr.MNE(inan)=nan;
%     stitch_hr.MNN(inan)=nan;
%     stitch_hr.MNZ(inan)=nan;
%     
%     %store calibration values
%     icals=find(ceil(stitch_hr.t(flipstart_hr(j)))==ceil(flipInfoAll.t));
%     if stitch_hr.t(flipstart_hr(1))<datenum(2019,08,13)
%         cal2.x_plus=flipInfoAll.gCalTCor(icals(1));
%         cal2.y_plus=flipInfoAll.gCalTCor(icals(2));
%     else
%         cal2.x_plus=flipInfoAll.gCalTCor(icals(1));
%         cal2.y_plus=flipInfoAll.gCalTCor(icals(2));
%         cal2.y_minus=flipInfoAll.gCalTCor(icals(3));
%         cal2.x2_plus=flipInfoAll.gCalTCor(icals(4));
%         cal2.x_minus=flipInfoAll.gCalTCor(icals(5));
%     end
%     
%     if j==1
%         cal1=cal2;
%         continue
%     end
%     
%     % apply calibration correction to interval
%     eint=stitch_hr.MNE(flipstart_hr(j-1)+5:flipstart_hr(j)-2);
%     nint=stitch_hr.MNN(flipstart_hr(j-1)+5:flipstart_hr(j)-2);
%     
%     if i<24 || i>28 %bad calibrations from Dec-Feb, use previous drift rates
%         xdrift=(cal2.x_plus-cal1.x_plus)/(flipdate_min(i)-flipdate_min(i-1))/24; %per hour
%         ydrift=(cal2.y_plus-cal1.y_plus)/(flipdate_min(i)-flipdate_min(i-1))/24; %per hour
%     end
%     disp(['X drift rate = ' num2str(xdrift*24*365*10^5) ' \mug/yr'])
%     disp(['Y drift rate = ' num2str(ydrift*24*365*10^5) ' \mug/yr'])
%     
%     elin=linspace(0,xdrift*length(eint),length(eint))';
%     nlin=linspace(0,ydrift*length(nint),length(nint))';
%     
%     eint_cal=eint-elin;
%     nint_cal=nint-nlin;
%     
%     %minimize offset between segments by aligning last and first day
%     p_e1=polyfit([flipstart_hr(i-1)-25:flipstart_hr(i-1)-2]',stitch_hr.MNE(flipstart_hr(i-1)-25:flipstart_hr(i-1)-2),1);
%     lin_e2=polyval(p_e1,[flipstart_hr(i-1)+5:flipstart_hr(i-1)+28]');
%     e2_dif=lin_e2-eint_cal(1:24);
%     Ge=ones(size(e2_dif));
%     e_offset=inv(Ge'*Ge)*Ge'*e2_dif;
%     
%     p_n1=polyfit([flipstart_hr(i-1)-25:flipstart_hr(i-1)-2]',stitch_hr.MNN(flipstart_hr(i-1)-25:flipstart_hr(i-1)-2),1);
%     lin_n2=polyval(p_n1,[flipstart_hr(i-1)+5:flipstart_hr(i-1)+28]');
%     n2_dif=lin_n2-nint_cal(1:24);
%     Gn=ones(size(n2_dif));
%     n_offset=inv(Ge'*Ge)*Ge'*n2_dif;
%     
%     eint_cor=eint_cal+e_offset;
%     nint_cor=nint_cal+n_offset;
%     
%     %replace data segment in 'stitch_hr'
%     stitch_hr.MNE(flipstart_hr(j-1)+5:flipstart_hr(j)-2)=eint_cor;
%     stitch_hr.MNN(flipstart_hr(j-1)+5:flipstart_hr(j)-2)=nint_cor;
%     
%     % current cal gets stored as previous cal for next iteration
%     plot(stitch_hr.t,stitch_hr.MNE,'linewidth',1)
%     cal1=cal2;
% end

%convert accel to tilt in microrad
stitch_min.LAX=asin(stitch_min.MNE/9.81)*10^6;stitch_min.LAX=stitch_min.LAX-stitch_min.LAX(1);
stitch_min.LAY=asin(stitch_min.MNN/9.81)*10^6;stitch_min.LAY=stitch_min.LAY-stitch_min.LAY(1);
stitch_hr.LAX=asin(stitch_hr.MNE/9.81)*10^6;stitch_hr.LAX=stitch_hr.LAX-stitch_hr.LAX(1);
stitch_hr.LAY=asin(stitch_hr.MNN/9.81)*10^6;stitch_hr.LAY=stitch_hr.LAY-stitch_hr.LAY(1);

%save figures, variables
% save('../calibrations/Axial/axialstitch_hr.mat','stitch_hr')
save('../calibrations/Axial/axialstitch_min_temp.mat','stitch_min')

%% plot and compare to AXCC1 Lily tiltmeter
if comptilts
    t0=floor(stitch_min.t(1));
    tf=floor(stitch_min.t(end));
    
    %load data
    sta='AXCC1';
    cha1={'LAX','LAY'};
    if exist('../calibrations/Axial/AXCC1_min.mat','file')
        load('../calibrations/Axial/AXCC1_min.mat')
        t1=ceil(AXCC1.time(end));
    else
        t1=t0;
        AXCC1.time=[];AXCC1.LAX=[];AXCC1.LAY=[];
    end
    
    while t1<=tf
        t1_s=datestr(t1,31); t1_s=t1_s(1:10);
        for i=1:2
            file_string=['AXCC1_' cha1{i} '_' t1_s '.miniseed'];
            %attempt download if file not found
            if ~exist(['../tiltcompare/AXCC1/' file_string],'file')
                IRIS_data_pull('AXCC1',cha1{i},'11',t1,t1+1);
            end
            %some dates have no data (power failure, etc.)
            if exist(['../tiltcompare/AXCC1/' file_string],'file')
                temp=rdmseed(['../tiltcompare/AXCC1/' file_string]);
                dtemp=double(cat(1,temp.d));
                if length(dtemp)<3600 %short timeseries will throw errors
                    t1=t1+1;
                    continue
                end
                dtempd=decimate(dtemp,6,'fir');
                dtempd=decimate(dtempd,10,'fir');
                eval(['AXCC1.' cha1{i} '=[AXCC1.' cha1{i} ';dtempd];']);
                if i==1
                    ttemp=cat(1,temp.t);
                    ttempd=decimate(ttemp,6,'fir');
                    ttempd=decimate(ttempd,10,'fir');
                    AXCC1.time=[AXCC1.time;ttempd];
                end
            end
        end
        t1=t1+1;
    end
    
    %AXCC1 compass correction
    crd_rot=[cosd(CCMP(4).plus_y) sind(CCMP(4).plus_y); -sind(CCMP(4).plus_y) cosd(CCMP(4).plus_y)];
    temp=crd_rot*[AXCC1.LAX';AXCC1.LAY']; AXCC1.LAX_rot=temp(1,:)'; AXCC1.LAY_rot=temp(2,:)';
    
    %start from 0
    AXCC1.LAX_rot=AXCC1.LAX_rot-AXCC1.LAX_rot(1);
    AXCC1.LAY_rot=AXCC1.LAY_rot-AXCC1.LAY_rot(1);
    
    figure(67)
    hold on
    plot(stitch_min.t,stitch_min.LAX,'k','linewidth',1)
    plot(AXCC1.time,AXCC1.LAX_rot,'b','linewidth',1)
    legend('SCTA','LILY','location','northwest')
    
    figure(68)
    hold on
    plot(stitch_min.t,stitch_min.LAY,'k','linewidth',1)
    plot(AXCC1.time,-AXCC1.LAY_rot,'b','linewidth',1)
    legend('SCTA','LILY','location','northwest')
end

save('../calibrations/Axial/AXCC1_min','AXCC1')
%% placeholder code  to be incorporated later
% intended to mark where calibrations are made to check that stitching is
% properly lining up segments

% keyboard % to ensure this doesn't run until I get it in the proper place
% 
% [~,ical,~]=unique(round(flipInfoAll.t)); % finds calibration days
% t_flip=flipInfoAll.t(ical); % finds time of calibration
% 
% for j=1:length(t_flip)
% [~,im(j)]=min(abs(stitch_min.t-t_flip(j)));
% end
% 
% plot(t_flip,stitch_min.LAX1(im),'ok','markersize',10)