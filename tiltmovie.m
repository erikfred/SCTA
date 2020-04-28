%tiltmovie.m
%
% Generates an animation of how tilt at Axial Seamount changes from late
% April through late July
%

clear; close all

load('../tiltcompare/inflation_reversal/jun2019/mar01-jul24.mat')
load('../calderacoordinates.mat')

%truncate data to avoid offsets at start
CC1.time=AXCC1.time(AXCC1.time>datenum(2019,04,23));
CC1.LAX_rot=AXCC1.LAX_rot(AXCC1.time>datenum(2019,04,23));
CC1.LAY_rot=AXCC1.LAY_rot(AXCC1.time>datenum(2019,04,23));
CC1.BDO=interp1(AXCC1.p_time,AXCC1.BDO,CC1.time,'pchip');

EC2.time=AXEC2.time(AXEC2.time>datenum(2019,04,23));
EC2.LAX_rot=AXEC2.LAX_rot(AXEC2.time>datenum(2019,04,23));
EC2.LAY_rot=AXEC2.LAY_rot(AXEC2.time>datenum(2019,04,23));
EC2.BDO=interp1(AXEC2.p_time,AXEC2.BDO,EC2.time,'pchip');

ASH.time=ASHES.time(ASHES.time>datenum(2019,04,23));
ASH.LAX_rot=ASHES.LAX_rot(ASHES.time>datenum(2019,04,23));
ASH.LAY_rot=ASHES.LAY_rot(ASHES.time>datenum(2019,04,23));
ASH.BDO=interp1(ASHES.p_time,ASHES.BDO,ASH.time,'pchip');

ID1.time=AXID1.time(AXID1.time>datenum(2019,04,23));
ID1.LAX_rot=AXID1.LAX_rot(AXID1.time>datenum(2019,04,23));
ID1.LAY_rot=AXID1.LAY_rot(AXID1.time>datenum(2019,04,23));
ID1.BDO=interp1(AXID1.p_time,AXID1.BDO,ID1.time,'pchip');

%interpolate onto same time array
EC2.LAX_rot=interp1(EC2.time,EC2.LAX_rot,CC1.time,'pchip');
EC2.LAY_rot=interp1(EC2.time,EC2.LAY_rot,CC1.time,'pchip');
EC2.BDO=interp1(EC2.time,EC2.BDO,CC1.time,'pchip');
EC2.time=CC1.time;

ASH.LAX_rot=interp1(ASH.time,ASH.LAX_rot,CC1.time,'pchip');
ASH.LAY_rot=interp1(ASH.time,ASH.LAY_rot,CC1.time,'pchip');
ASH.BDO=interp1(ASH.time,ASH.BDO,CC1.time,'pchip');
ASH.time=CC1.time;

ID1.LAX_rot=interp1(ID1.time,ID1.LAX_rot,CC1.time,'pchip');
ID1.LAY_rot=interp1(ID1.time,ID1.LAY_rot,CC1.time,'pchip');
ID1.BDO=interp1(ID1.time,ID1.BDO,CC1.time,'pchip');
ID1.time=CC1.time;

%smooth with 7-day lowpass filter
[b,a]=butter(3,2/(7*24*60),'low');

CC1.LAX_rot=filtfilt(b,a,CC1.LAX_rot);
CC1.LAY_rot=filtfilt(b,a,CC1.LAY_rot);
CC1.BDO=filtfilt(b,a,CC1.BDO);

EC2.LAX_rot=filtfilt(b,a,EC2.LAX_rot);
EC2.LAY_rot=filtfilt(b,a,EC2.LAY_rot);
EC2.BDO=filtfilt(b,a,EC2.BDO);

ASH.LAX_rot=filtfilt(b,a,ASH.LAX_rot);
ASH.LAY_rot=filtfilt(b,a,ASH.LAY_rot);
ASH.BDO=filtfilt(b,a,ASH.BDO);

ID1.LAX_rot=filtfilt(b,a,ID1.LAX_rot);
ID1.LAY_rot=filtfilt(b,a,ID1.LAY_rot);
ID1.BDO=filtfilt(b,a,ID1.BDO);

%prep for plots
caldera_loc=llh2local(caldera_rim',[AXCC1.lon,AXCC1.lat])';
crds=llh2local([AXCC1.lon,AXCC1.lat]',[AXCC1.lon,AXCC1.lat]);
AXCC1.x=crds(1); AXCC1.y=crds(2);
crds=llh2local([AXEC2.lon,AXEC2.lat]',[AXCC1.lon,AXCC1.lat]);
AXEC2.x=crds(1); AXEC2.y=crds(2);
crds=llh2local([AXID1.lon,AXID1.lat]',[AXCC1.lon,AXCC1.lat]);
AXID1.x=crds(1); AXID1.y=crds(2);
crds=llh2local([ASHES.lon,ASHES.lat]',[AXCC1.lon,AXCC1.lat]);
ASHES.x=crds(1); ASHES.y=crds(2);


bin=1:9*24*60:length(CC1.time);
for i=1:length(bin)-1
    int1=bin(i);
    int2=bin(i+1);
    
    %calculate tilt gradients
    temp=polyfit([int1:int2]',CC1.LAX_rot(int1:int2),1); temp2=polyval(temp,[int1:int2]');
    CC1.x=temp2(end)-temp2(1);
    temp=polyfit([int1:int2]',CC1.LAY_rot(int1:int2),1); temp2=polyval(temp,[int1:int2]');
    CC1.y=temp2(end)-temp2(1);
    temp=polyfit([int1:int2]',EC2.LAX_rot(int1:int2),1); temp2=polyval(temp,[int1:int2]');
    EC2.x=temp2(end)-temp2(1);
    temp=polyfit([int1:int2]',EC2.LAY_rot(int1:int2),1); temp2=polyval(temp,[int1:int2]');
    EC2.y=temp2(end)-temp2(1);
    temp=polyfit([int1:int2]',ASH.LAX_rot(int1:int2),1); temp2=polyval(temp,[int1:int2]');
    ASH.x=temp2(end)-temp2(1);
    temp=polyfit([int1:int2]',ASH.LAY_rot(int1:int2),1); temp2=polyval(temp,[int1:int2]');
    ASH.y=temp2(end)-temp2(1);
    temp=polyfit([int1:int2]',ID1.LAX_rot(int1:int2),1); temp2=polyval(temp,[int1:int2]');
    ID1.x=temp2(end)-temp2(1);
    temp=polyfit([int1:int2]',ID1.LAY_rot(int1:int2),1); temp2=polyval(temp,[int1:int2]');
    ID1.y=temp2(end)-temp2(1);
    
    %calculate elevation changes (includes factor of -1)
    temp=polyfit([int1:int2]',CC1.BDO(int1:int2),1); temp2=polyval(temp,[int1:int2]');
    CC1.z=temp2(1)-temp2(end);
    temp=polyfit([int1:int2]',EC2.BDO(int1:int2),1); temp2=polyval(temp,[int1:int2]');
    EC2.z=temp2(1)-temp2(end);
    temp=polyfit([int1:int2]',ASH.BDO(int1:int2),1); temp2=polyval(temp,[int1:int2]');
    ASH.z=temp2(1)-temp2(end);
    temp=polyfit([int1:int2]',ID1.BDO(int1:int2),1); temp2=polyval(temp,[int1:int2]');
    ID1.z=temp2(1)-temp2(end);
    
    %plot tilt direction changes
    figure(24)
    clf
    plot(caldera_loc(:,1),caldera_loc(:,2),'k')
    hold on
    plot(AXCC1.x,AXCC1.y,'k^','markersize',15,'linewidth',2)
    plot(AXID1.x,AXID1.y,'k^','markersize',15,'linewidth',2)
    plot(AXEC2.x,AXEC2.y,'k^','markersize',15,'linewidth',2)
    plot(ASHES.x,ASHES.y,'k^','markersize',15,'linewidth',2)
    quiver(AXCC1.x,AXCC1.y,CC1.x,CC1.y,'b','linewidth',2)
%     quiver(AXCC1.x,AXCC1.y,0,CC1.z,'k','linewidth',2)
    quiver(AXEC2.x,AXEC2.y,EC2.x,EC2.y,'b','linewidth',2)
%     quiver(AXEC2.x,AXEC2.y,0,EC2.z,'k','linewidth',2)
    quiver(AXID1.x,AXID1.y,ID1.x,ID1.y,'b','linewidth',2)
%     quiver(AXID1.x,AXID1.y,0,ID1.z,'k','linewidth',2)
    h1=quiver(ASHES.x,ASHES.y,ASH.x,ASH.y,'b','linewidth',2);
%     h2=quiver(ASHES.x,ASHES.y,0,ASH.z,'k','linewidth',2);
    axis([-6 6 -6 6])
    set(gca,'fontsize',18)
%     legend([h1 h2],'tilt','elev')
    title([datestr(round(CC1.time(int1)),'mmmdd') ' - ' datestr(round(CC1.time(int2)),'mmmdd')])
    
    M(i)=getframe(gcf);
    
    orient tall
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 5.5 5.5];
    print(['../tiltcompare/inflation_reversal/movie/mapview' num2str(i) 'c'],'-dtiff')
end

save('../tiltcompare/inflation_reversal/movie/movie','CC1','EC2','ASH','ID1','M','bin')