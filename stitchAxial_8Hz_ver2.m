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

era=2; % 1 - pre-move era, 2 - post-move era
method=2; % 1 - cal-to-cal drift, 2 - constant drift, 3 - zero drift

if era==1
    flipfile='../calibrations/Axial/axialdata.mat';
    datafile_hr='../calibrations/Axial/axialdata_hr.mat';
    datafile_min='../calibrations/Axial/axialdata_min.mat';
    savefile_hr='../calibrations/Axial/axialstitch_hr.mat';
    savefile_min='../calibrations/Axial/axialstitch_min.mat';
elseif era==2
    flipfile='../calibrations/Axial/axialdata_newloc.mat';
    datafile_hr='../calibrations/Axial/axialdata_newloc_hr.mat';
    datafile_min='../calibrations/Axial/axialdata_newloc_min.mat';
    savefile_hr='../calibrations/Axial/axialstitch_newloc_hr.mat';
    savefile_min='../calibrations/Axial/axialstitch_newloc_min.mat';
else
    error('specify era!')
end

load(flipfile,'flipInfoAll')
load(datafile_hr)
load(datafile_min)
load('../compass_directions.mat')
load('../calibrations/Axial/badcaldates','bad_x1','bad_x2','bad_negx','bad_y','bad_negy','bads')
% will add these permanently later when I understand them fully
cat(2,bad_x1,datenum(2021,08,23))
cat(2,bad_negy,datenum(2021,06,16))

% pick up data structure from previous run, or start from scratch
if exist(savefile_hr,'file') && exist(savefile_min,'file')
    load(savefile_hr)
    load(savefile_min)
    
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
flipstart_min(44)=[]; % temp fix - why isn't this flip recorded in flipInfoAll?
flipstart_min(57)=[]; % temp fix
flipstart_min(60)=[]; % temp fix
flipstart_min(60)=[]; % temp fix
flipstart_min(60)=[]; % temp fix
flipstart_min(63:end)=[]; % temp fix
flipdate_min=stitch_min.t(flipstart_min); %datenums of flips for reference

iend=length(stitch_min.t);
for i=1:length(flipstart_min)

    %substitute points during calibration and recovery with NaNs
    inan=flipstart_min(i)-4:flipstart_min(i)+119;
    stitch_min.MNE(inan)=nan;
    stitch_min.MNN(inan)=nan;
    stitch_min.MNZ(inan)=nan;
    
    %store calibration values
    icals=find(ceil(stitch_min.t(flipstart_min(i)))==ceil(flipInfoAll.t));
    if stitch_min.t(flipstart_min(i))<datenum(2019,08,13)
        cal2(i).x_plus=flipInfoAll.gCalTCor(icals(1));
        cal2(i).y_plus=flipInfoAll.gCalTCor(icals(2));
    elseif i==12 || i==17
        cal2(i)=cal2(i-1);
    else
        cal2(i).x_plus=flipInfoAll.gCalTCor(icals(1));
        cal2(i).y_plus=flipInfoAll.gCalTCor(icals(2));
        cal2(i).y_minus=flipInfoAll.gCalTCor(icals(3));
        cal2(i).x2_plus=flipInfoAll.gCalTCor(icals(4));
        cal2(i).x_minus=flipInfoAll.gCalTCor(icals(5));
    end
    
    if i==1
        continue
    end
    
    % apply calibration correction to interval
    eint=stitch_min.MNE(flipstart_min(i-1)+120:flipstart_min(i)-5);
    nint=stitch_min.MNN(flipstart_min(i-1)+120:flipstart_min(i)-5);
    
    if ~any(floor(flipdate_min(i))==bad_x1) && ~any(floor(flipdate_min(i-1))==bad_x1) %skip drift calculation on bad calibrations
        xdrift=(cal2(i).x_plus-cal2(i-1).x_plus)/(flipdate_min(i)-flipdate_min(i-1))/24/60; %per minute
    end
    if ~any(floor(flipdate_min(i))==bad_y) && ~any(floor(flipdate_min(i-1))==bad_y) %skip drift calculation on bad calibrations
        ydrift=(cal2(i).y_plus-cal2(i-1).y_plus)/(flipdate_min(i)-flipdate_min(i-1))/24/60; %per minute
    end
    
    % overwrite individual drift calculations, if called for in config
    if method==2
        % x
        ixcals=[1:5:5*57, 5*58+1:5:length(flipInfoAll.t)]; % exclude bad cals (manual for now)
        linx=polyfit(flipInfoAll.t(ixcals),flipInfoAll.gCal(ixcals),1); % fit good cals w/linear
        xdrift=linx(1)/24/60; % convert drift to per minute
        % y
        iycals=2:5:length(flipInfoAll.t);
        liny=polyfit(flipInfoAll.t(iycals),flipInfoAll.gCal(iycals),1);
        ydrift=liny(1)/24/60;
        
        if i<16
            ydrift=9.207e-7/24/60;
        elseif i>=16
            ydrift=8.176e-7/24/60;
        end
    elseif method==3
        xdrift=0;
        ydrift=0;
    end
    
    disp(['X drift rate = ' num2str(xdrift*60*24*365*10^5) ' \mug/yr'])
    disp(['Y drift rate = ' num2str(ydrift*60*24*365*10^5) ' \mug/yr'])
    
    elin=linspace(0,xdrift*length(eint),length(eint))';
    nlin=linspace(0,ydrift*length(nint),length(nint))';
    
    eint_cal=eint-elin;
    nint_cal=nint-nlin;
    
    figure(3); hold on %for in-line testing
    plot(stitch_min.t(flipstart_min(i-1)+120:flipstart_min(i)-5),nint_cal-nanmean(nint_cal),'linewidth',1)
%     keyboard
    
    % remove anomalies (manual for now)
    if i==2
        iseg1=1:7990; seg1e=eint_cal(iseg1); seg1n=nint_cal(iseg1);
        p1e=polyfit(iseg1,seg1e,1); p1n=polyfit(iseg1,seg1n,1);
        iseg2=iseg1(end)+1:length(eint_cal); seg2e_lin=polyval(p1e,iseg2); seg2n_lin=polyval(p1n,iseg2);
        eint_cal(iseg2)=seg2e_lin; nint_cal(iseg2)=seg2n_lin;
    elseif i==8
        iseg=4639:4700; tseg=stitch_min.t(flipstart_min(i-1)+120:flipstart_min(i)-5);
        eint_cal(iseg)=interp1([tseg(1:iseg(1)-1); tseg(iseg(end)+1:end)],[eint_cal(1:iseg(1)-1); eint_cal(iseg(end)+1:end)],tseg(iseg));
%         iseg=4639:4700; eint_cal(iseg)=[];
%         tseg=stitch_min.t(flipstart_min(i-1)+120:flipstart_min(i)-5); tseg(iseg)=[];
%         tseg2=[tseg(1:iseg(1)-1); [tseg(iseg(1)-1)+1/24/60:1/24/60:tseg(iseg(end)+1)-1/24/60]'; tseg(iseg(end)+1:end)];
%         eint_cal=[eint_cal(1:iseg(1)-1); interp1(tseg,eint_cal,tseg2); eint_cal(iseg(end)+1:end)];
    elseif i==18
        t_all=stitch_min.t(flipstart_min(i-1)+120:flipstart_min(i)-5);
        iseg1=1:9341; tseg1=t_all(iseg1); seg1e=eint_cal(iseg1); seg1n=nint_cal(iseg1);
        iseg2=9850:length(eint_cal); tseg2=t_all(iseg2); seg2e=eint_cal(iseg2); seg2n=nint_cal(iseg2);
        p1e=polyfit(tseg1,seg1e,1); seg2e_lin=polyval(p1e,tseg2); p1n=polyfit(tseg1,seg1n,1); seg2n_lin=polyval(p1n,tseg2);
        e2dif=seg2e_lin-seg2e; n2dif=seg2n_lin-seg2n;
        Ge=ones(size(e2dif));
        seg2e=seg2e+inv(Ge'*Ge)*Ge'*e2dif; seg2n=seg2n+inv(Ge'*Ge)*Ge'*n2dif;
        eint_cal=[seg1e; polyval(p1e,t_all(iseg1(end)+1:iseg2(1)-1)); seg2e]; nint_cal=[seg1n; polyval(p1n,t_all(iseg1(end)+1:iseg2(1)-1)); seg2n];
    end

    %minimize offset between segments by aligning last and first day
    dur=1440; % no of points in each segment to use for alignment [min]
    int1a=dur+4; int1b=5;
    int2a=120; int2b=dur+119;
    p_e1=polyfit([flipstart_min(i-1)-int1a:flipstart_min(i-1)-int1b]',stitch_min.MNE(flipstart_min(i-1)-int1a:flipstart_min(i-1)-int1b),1);
    lin_e2=polyval(p_e1,[flipstart_min(i-1)+int2a:flipstart_min(i-1)+int2b]');
    if length(eint)>=dur
        e2_dif=lin_e2-eint_cal(1:dur);
    else
        e2_dif=lin_e2(1:length(eint))-eint_cal;
    end
    Ge=ones(size(e2_dif));
    e_offset=inv(Ge'*Ge)*Ge'*e2_dif;
    
    p_n1=polyfit([flipstart_min(i-1)-1444:flipstart_min(i-1)-5]',stitch_min.MNN(flipstart_min(i-1)-1444:flipstart_min(i-1)-5),1);
    lin_n2=polyval(p_n1,[flipstart_min(i-1)+120:flipstart_min(i-1)+1559]');
    if length(nint)>=1440
        n2_dif=lin_n2-nint_cal(1:1440);
    else
        n2_dif=lin_n2(1:length(nint))-nint_cal;
    end
    Gn=ones(size(n2_dif));
    n_offset=inv(Gn'*Gn)*Gn'*n2_dif;
    
    eint_cor=eint_cal+e_offset;
    nint_cor=nint_cal+n_offset;
    % another temporary fix to correct offset
    if i==20
        eint_cor=eint_cor+0.00001;
    elseif i==21
        nint_cor=nint_cor+2.5e-6;
    elseif i==26
        eint_cor=eint_cor+6e-6;
    elseif i==29
        eint_cor=eint_cor-7e-6;
    elseif i==33
        eint_cor=eint_cor-5e-6;
    elseif i==40
        eint_cor(33415:end)=eint_cor(33415:end)-1.25e-5;
        eint_cor(59897:end)=eint_cor(59897);
    end
    
    figure(2); clf %for in-line testing
    subplot(211); hold on
    plot(stitch_min.MNE(flipstart_min(i-1)+120:flipstart_min(i)-5)-nanmean(stitch_min.MNE(flipstart_min(i-1)+120:flipstart_min(i)-5)),'r','linewidth',1)
    plot(eint_cor-nanmean(eint_cor),'b','linewidth',1)
    subplot(212); hold on
    plot(stitch_min.MNN(flipstart_min(i-1)+120:flipstart_min(i)-5)-nanmean(stitch_min.MNN(flipstart_min(i-1)+120:flipstart_min(i)-5)),'r','linewidth',1)
    plot(nint_cor-nanmean(nint_cor),'b','linewidth',1)
%     keyboard
    
    %replace data segment in 'stitch_min'
    stitch_min.MNE(flipstart_min(i-1)+120:flipstart_min(i)-5)=eint_cor;
    stitch_min.MNN(flipstart_min(i-1)+120:flipstart_min(i)-5)=nint_cor;
    
    figure(1); hold on %for in-line testing
    plot(stitch_min.t(flipstart_min(i-1)+120:flipstart_min(i)-5),stitch_min.MNN(flipstart_min(i-1)+120:flipstart_min(i)-5),'linewidth',1)
%     plot(stitch_min.t([flipstart_min(i-1)+int2a:flipstart_min(i-1)+int2b]),lin_e2)
%     if i>15
%         keyboard
%     end
end

% tack on remaining data without detrending
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

% rotate channels into east and north
if era==1
    crd_rot=[cosd(CCMP(5).plus_y) sind(CCMP(5).plus_y); -sind(CCMP(5).plus_y) cosd(CCMP(5).plus_y)];
elseif era==2
    crd_rot=[cosd(CCMP(6).plus_y) sind(CCMP(6).plus_y); -sind(CCMP(6).plus_y) cosd(CCMP(6).plus_y)];
else
    error('specify era!')
end
temp=crd_rot*[stitch_min.MNE';stitch_min.MNN']; stitch_min.MNE_rot=temp(1,:)'; stitch_min.MNN_rot=temp(2,:)';
    temp=crd_rot*[stitch_hr.MNE';stitch_hr.MNN']; stitch_hr.MNE_rot=temp(1,:)'; stitch_hr.MNN_rot=temp(2,:)';

%convert accel to tilt in microrad, start from zero
stitch_min.LAX=asin(stitch_min.MNE_rot/9.81)*10^6;stitch_min.LAX=stitch_min.LAX-stitch_min.LAX(1);
stitch_min.LAY=asin(stitch_min.MNN_rot/9.81)*10^6;stitch_min.LAY=stitch_min.LAY-stitch_min.LAY(1);
stitch_hr.LAX=asin(stitch_hr.MNE_rot/9.81)*10^6;stitch_hr.LAX=stitch_hr.LAX-stitch_hr.LAX(1);
stitch_hr.LAY=asin(stitch_hr.MNN_rot/9.81)*10^6;stitch_hr.LAY=stitch_hr.LAY-stitch_hr.LAY(1);

% more quick fixes for poster
stitch_min.LAY(176944:end)=stitch_min.LAY(176944:end)+0.5;

%save figures, variables
% save('../calibrations/Axial/axialstitch_hr.mat','stitch_hr')
% save('../calibrations/Axial/axialstitch_min_temp.mat','stitch_min')

%% plot and compare to AXCC1 Lily tiltmeter
if comptilts
    %load data
    if exist('../tiltcompare/SCTA_Lily_comp/AXCC1data_min.mat','file')
        load('../tiltcompare/SCTA_Lily_comp/AXCC1data_min.mat')
    else
        warning('Could not find LILY data')
    end
    
    if era==1
        cond=lily_min.t<datenum(2020,09,11);
        lily_min.t=lily_min.t(cond);
        lily_min.LAX=lily_min.LAX(cond);
        lily_min.LAY=lily_min.LAY(cond);
        lily_min.BDO=lily_min.BDO(cond);
    elseif era==2
        cond=lily_min.t>datenum(2020,09,11,12,0,0);
        lily_min.t=lily_min.t(cond);
        lily_min.LAX=lily_min.LAX(cond);
        lily_min.LAY=lily_min.LAY(cond);
        lily_min.BDO=lily_min.BDO(cond);
    else
        error('specify era!')
    end
    
    %AXCC1 compass correction
    crd_rot=[cosd(CCMP(4).plus_y) sind(CCMP(4).plus_y); -sind(CCMP(4).plus_y) cosd(CCMP(4).plus_y)];
    temp=crd_rot*[lily_min.LAX';lily_min.LAY']; lily_min.LAX_rot=temp(1,:)'; lily_min.LAY_rot=temp(2,:)';
    
    %start from 0
    lily_min.LAX_rot=lily_min.LAX_rot-lily_min.LAX_rot(1);
    lily_min.LAY_rot=lily_min.LAY_rot-lily_min.LAY_rot(1);
    
    % more quick fixes for poster
    lily_min.LAX_rot(164128:end)=lily_min.LAX_rot(164128:end)+4;
    lily_min.LAY_rot(164128:end)=lily_min.LAY_rot(164128:end)-4;
    
    %interpolate onto same time step
    lily_min.LAX_rot_int=interp1(lily_min.t,lily_min.LAX_rot,stitch_min.t);
    lily_min.LAY_rot_int=interp1(lily_min.t,lily_min.LAY_rot,stitch_min.t);
    
    
    figure(67)
    clf; hold on
    plot(stitch_min.t,-stitch_min.LAX,'r','linewidth',1)
    plot(lily_min.t,lily_min.LAX_rot,'b','linewidth',1)
    plot(stitch_min.t,-stitch_min.LAX-lily_min.LAX_rot_int,'k','linewidth',1)
    plot(flipdate_min,-stitch_min.LAX(flipstart_min-5),'+k','markersize',5,'linewidth',2)
    legend('SCTA','LILY','difference','location','northwest')
    title(['East Tilt ' datestr(floor(stitch_min.t(1))) ' to ' datestr(floor(stitch_min.t(end)))])
    ylabel('Tilt (\murad)')
    datetick('x',12)
    xtickangle(45)
    set(gca,'fontsize',14)
    box on; grid on; grid minor
    
    figure(68)
    clf; hold on
    plot(stitch_min.t,-stitch_min.LAY,'r','linewidth',1)
    plot(lily_min.t,lily_min.LAY_rot,'b','linewidth',1)
    plot(stitch_min.t,-stitch_min.LAY-lily_min.LAY_rot_int,'k','linewidth',1)
    plot(flipdate_min,-stitch_min.LAY(flipstart_min-5),'+k','markersize',5,'linewidth',2)
    legend('SCTA','LILY','difference','location','northeast')
    title(['North Tilt ' datestr(floor(stitch_min.t(1))) ' to ' datestr(floor(stitch_min.t(end)))])
    ylabel('Tilt (\murad)')
    datetick('x',12)
    xtickangle(45)
    set(gca,'fontsize',14)
    box on; grid on; grid minor
end

% save('../tiltcompare/SCTA_Lily_comp/AXCC1data_min_temp','lily_min')
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

%% fixing gaps/steps in Lily data

lily_min.LAX(437922:end)=lily_min.LAX(437922:end)-0.4;
lily_min.LAX(380931:382251)=NaN; lily_min.LAX(382251:end)=lily_min.LAX(382251:end)-2.5;
lily_min.LAX(253002:end)=lily_min.LAX(253002:end)+1.25;
lily_min.LAX(236467:240372)=NaN; lily_min.LAX(240372:end)=lily_min.LAX(240372:end)-0.25;
lily_min.LAX(164126:165088)=NaN; lily_min.LAX(165088:end)=lily_min.LAX(165088:end)+2.5;
lily_min.LAX(16464:17505)=NaN; lily_min.LAX(17505:end)=lily_min.LAX(17505:end)+1.5;

lily_min.LAY(413104:413411)=NaN; lily_min.LAY(413411:end)=lily_min.LAY(413411:end)-1;
lily_min.LAY(346602:end)=lily_min.LAY(346602:end)-1;
lily_min.LAY(253002:end)=lily_min.LAY(253002:end)-1;
lily_min.LAY(236467:240372)=NaN; lily_min.LAY(240372:end)=lily_min.LAY(240372:end)-0.5;
lily_min.LAY(164126:165088)=NaN; lily_min.LAY(165088:end)=lily_min.LAY(165088:end)-1;
lily_min.LAY(16464:17505)=NaN; lily_min.LAY(17505:end)=lily_min.LAY(17505:end)-2.5;

% fill in NaNs for easier checking
lilyx=inpaint_nans(lily_min.LAX);
lilyy=inpaint_nans(lily_min.LAY);