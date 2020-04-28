% axial_tilt_reversals.m
%
% Examines May-Jun inflation reversal at Axial, comparing tilt and pressure
% signals before and during the reversal for evidence of asymmetry
%

close all; clear

%%%%%CONFIG%%%%%%
eventnum=3; %1=12/2017 2=6/2018 3=5/2019

filtpres=true; %plots filtered pressures

strts=[datenum(2017,11,15),datenum(2018,5,1),datenum(2019,3,1),datenum(2019,6,01)];
enders=[datenum(2018,01,15),datenum(2018,7,1),datenum(2019,6,1),datenum(2019,7,24)];
fldrs={'dec2017','jun2018','may2019','jun2019'};
ints=[5000 10000 10000 0;40000 30000 28000 0;46000 37500 32000 0;62000 47500 42500 0];

t0=strts(eventnum);
tf=enders(eventnum);
fldr=fldrs{eventnum};
int1a=ints(1,eventnum);
int1b=ints(2,eventnum);
int2a=ints(3,eventnum);
int2b=ints(4,eventnum);
%%%%%END CONFIG%%%%%

load('../compass_directions.mat')
load('../calderacoordinates.mat')

pfac=6894.745/10065*100; %6.89Pa/count / 10065Pa/m * 100cm/m

%load data
sta={'AXCC1','AXEC2','AXID1'};
cha={'LAX','LAY','BDO'};
dec=[6 6 12;5 5 10;2 2 10];

%ASHES first
if ~exist(['../tiltcompare/inflation_reversal/' fldr '/ASHES.mat'],'file')
    [fl,dr]=uigetfile('../tiltcompare/ASHES/*.nc','Select tilt file');
    [fl1,dr1]=uigetfile('../tiltcompare/ASHES/pressure/*.nc','Select pressure files',...
        'multiselect','on');
    fl1=sort(fl1); %puts in chronological order
    %tilt
    temp=ncread([dr fl],'time');
    temp1=decimate(temp,dec(1,1),'fir');
    temp2=decimate(temp1,dec(2,1),'fir');
    temp3=decimate(temp2,dec(3,1),'fir');
    ASHES.time=temp3/60/60/24+datenum(1900,1,1);
    temp=ncread([dr fl],'lily_x_tilt');
    temp1=decimate(temp,dec(1,1),'fir');
    temp2=decimate(temp1,dec(2,1),'fir');
    temp3=decimate(temp2,dec(3,1),'fir');
    ASHES.LAX=temp3;
    temp=ncread([dr fl],'lily_y_tilt');
    temp1=decimate(temp,dec(1,1),'fir');
    temp2=decimate(temp1,dec(2,1),'fir');
    temp3=decimate(temp2,dec(3,1),'fir');
    ASHES.LAY=temp3;
    
    ASHES.p_time=[];
    ASHES.BDO=[];
    %pressure
    for l=1:length(fl1)
        ttemp=ncread([dr1 fl1{l}],'time');
        ttemp1=decimate(ttemp,dec(1,3),'fir');
        ttemp2=decimate(ttemp1,dec(2,3),'fir');
        ttemp3=decimate(ttemp2,dec(3,3),'fir');
        ASHES.p_time=[ASHES.p_time; ttemp3/60/60/24+datenum(1900,1,1)];
        ptemp=ncread([dr1 fl1{l}],'bottom_pressure');
        ptemp1=decimate(ptemp,dec(1,3),'fir');
        ptemp2=decimate(ptemp1,dec(2,3),'fir');
        ptemp3=decimate(ptemp2,dec(3,3),'fir');
        ASHES.BDO=[ASHES.BDO; ptemp3*pfac];
    end
    ASHES.lon=-130.01368;ASHES.lat=45.93363;ASHES.depth=1542;
else
    load(['../tiltcompare/inflation_reversal/' fldr '/ASHES.mat'])
end
%now others
if ~exist(['../tiltcompare/inflation_reversal/' fldr '/AXEC2.mat'],'file')
    for i=1:length(sta)
        eval([sta{i} '.time=[];']);
        eval([sta{i} '.p_time=[];']);
        for j=1:length(cha)
            t1=t0;
            file_str=[sta{i} '_' cha{j} '_'];
            eval([sta{i} '.' cha{j} '=[];']);
            while t1<tf
                t1_s=datestr(t1,31); t1_s=t1_s(1:10);
                file_string=[file_str t1_s '.miniseed'];
                %attempt download if file not found
                if ~exist(['../tiltcompare/' sta{i} '/' file_string],'file')
                    IRIS_data_pull(sta{i},cha{j},'11',t1,t1+1);
                end
                %some dates have no data (power failure, etc.)
                if exist(['../tiltcompare/' sta{i} '/' file_string],'file')
                    temp=rdmseed(['../tiltcompare/' sta{i} '/' file_string]);
                    dtemp=double(cat(1,temp.d));
                    if length(dtemp)<3600 %short timeseries will throw errors
                        t1=t1+1;
                        continue
                    end
                    dtempd=decimate(dtemp,dec(1,j),'fir');
                    dtempd=decimate(dtempd,dec(2,j),'fir');
                    dtempd=decimate(dtempd,dec(3,j),'fir');
                    eval([sta{i} '.' cha{j} '=[' sta{i} '.' cha{j} ';dtempd];']);
                    if j==1
                        ttemp=cat(1,temp.t);
                        ttempd=decimate(ttemp,dec(1,j),'fir');
                        ttempd=decimate(ttempd,dec(2,j),'fir');
                        ttempd=decimate(ttempd,dec(3,j),'fir');
                        eval([sta{i} '.time=[' sta{i} '.time;ttempd];']);
                    elseif j==3
                        ttemp=cat(1,temp.t);
                        ttempd=decimate(ttemp,dec(1,j),'fir');
                        ttempd=decimate(ttempd,dec(2,j),'fir');
                        ttempd=decimate(ttempd,dec(3,j),'fir');
                        eval([sta{i} '.p_time=[' sta{i} '.p_time;ttempd];']);
                    end
                end
                t1=t1+1;
            end
        end
        eval([sta{i} '.BDO=' sta{i} '.BDO*pfac;']);
    end
    AXCC1.lat=45.954681;AXCC1.lon=-130.008896;AXCC1.depth=-1528;
    AXEC2.lat=45.939671;AXEC2.lon=-129.973801;AXEC2.depth=-1519;
    AXID1.lat=45.925732;AXID1.lon=-129.977997;AXID1.depth=-1527.5;
else
    load(['../tiltcompare/inflation_reversal/' fldr '/AXEC2.mat'])
    load(['../tiltcompare/inflation_reversal/' fldr '/AXCC1.mat'])
    load(['../tiltcompare/inflation_reversal/' fldr '/AXID1.mat'])
end

%save everything before changing
save(['../tiltcompare/inflation_reversal/' fldr '/AXEC2'],'AXEC2')
save(['../tiltcompare/inflation_reversal/' fldr '/AXCC1'],'AXCC1')
save(['../tiltcompare/inflation_reversal/' fldr '/AXID1'],'AXID1')
save(['../tiltcompare/inflation_reversal/' fldr '/ASHES'],'ASHES')

%subsect for plots
if eventnum==3
    AXCC1.time=AXCC1.time(60855:end);
    AXCC1.LAX=AXCC1.LAX(60855:end);
    AXCC1.LAY=AXCC1.LAY(60855:end);
    AXCC1.p_time=AXCC1.p_time(60855:end);
    AXCC1.BDO=AXCC1.BDO(60855:end);
    
    AXEC2.time=AXEC2.time(60855:end);
    AXEC2.LAX=AXEC2.LAX(60855:end);
    AXEC2.LAY=AXEC2.LAY(60855:end);
    AXEC2.p_time=AXEC2.p_time(60855:end);
    AXEC2.BDO=AXEC2.BDO(60855:end);
    
    AXID1.time=AXID1.time(60855:end);
    AXID1.LAX=AXID1.LAX(60855:end);
    AXID1.LAY=AXID1.LAY(60855:end);
    AXID1.p_time=AXID1.p_time(60855:end);
    AXID1.BDO=AXID1.BDO(60855:end);
    
    ASHES.time=ASHES.time(65000:end);
    ASHES.LAX=ASHES.LAX(65000:end);
    ASHES.LAY=ASHES.LAY(65000:end);
    ASHES.p_time=ASHES.p_time(65000:end);
    ASHES.BDO=ASHES.BDO(65000:end);
end

if eventnum==2
    AXCC1.time=AXCC1.time(29000:end);
    AXCC1.LAX=AXCC1.LAX(29000:end);
    AXCC1.LAY=AXCC1.LAY(29000:end);
    AXCC1.BDO=AXCC1.BDO(29000:end);
    AXCC1.p_time=AXCC1.p_time(29000:end);
    
    AXEC2.time=AXEC2.time(29000:end);
    AXEC2.LAX=AXEC2.LAX(29000:end);
    AXEC2.LAY=AXEC2.LAY(29000:end);
    AXEC2.BDO=AXEC2.BDO(29000:end);
    AXEC2.p_time=AXEC2.p_time(29000:end);
    
    AXID1.time=AXID1.time(29000:end);
    AXID1.LAX=AXID1.LAX(29000:end);
    AXID1.LAY=AXID1.LAY(29000:end);
    AXID1.BDO=AXID1.BDO(29000:end);
    AXID1.p_time=AXID1.p_time(29000:end);
    
    ASHES.time=ASHES.time(29000:end);
    ASHES.LAX=ASHES.LAX(29000:end);
    ASHES.LAY=ASHES.LAY(29000:end);
    ASHES.BDO=ASHES.BDO(29000:end);
    ASHES.p_time=ASHES.p_time(29000:end);
end

%tilt magnitudes
temp=polyfit([int1a:int1b]',AXCC1.LAX(int1a:int1b),1);temp2=polyval(temp,int1a:int1b);
AXCC1.x_int1=(temp2(end)-temp2(1))/(AXCC1.time(int1b)-AXCC1.time(int1a));
temp=polyfit([int1a:int1b]',AXCC1.LAY(int1a:int1b),1);temp2=polyval(temp,int1a:int1b);
AXCC1.y_int1=(temp2(end)-temp2(1))/(AXCC1.time(int1b)-AXCC1.time(int1a));
temp=polyfit([int2a:int2b]',AXCC1.LAX(int2a:int2b),1);temp2=polyval(temp,int2a:int2b);
AXCC1.x_int2=(temp2(end)-temp2(1))/(AXCC1.time(int2b)-AXCC1.time(int2a));
temp=polyfit([int2a:int2b]',AXCC1.LAY(int2a:int2b),1);temp2=polyval(temp,int2a:int2b);
AXCC1.y_int2=(temp2(end)-temp2(1))/(AXCC1.time(int2b)-AXCC1.time(int2a));

temp=polyfit([int1a:int1b]',AXEC2.LAX(int1a:int1b),1);temp2=polyval(temp,int1a:int1b);
AXEC2.x_int1=(temp2(end)-temp2(1))/(AXEC2.time(int1b)-AXEC2.time(int1a));
temp=polyfit([int1a:int1b]',AXEC2.LAY(int1a:int1b),1);temp2=polyval(temp,int1a:int1b);
AXEC2.y_int1=(temp2(end)-temp2(1))/(AXEC2.time(int1b)-AXEC2.time(int1a));
temp=polyfit([int2a:int2b]',AXEC2.LAX(int2a:int2b),1);temp2=polyval(temp,int2a:int2b);
AXEC2.x_int2=(temp2(end)-temp2(1))/(AXEC2.time(int2b)-AXEC2.time(int2a));
temp=polyfit([int2a:int2b]',AXEC2.LAY(int2a:int2b),1);temp2=polyval(temp,int2a:int2b);
AXEC2.y_int2=(temp2(end)-temp2(1))/(AXEC2.time(int2b)-AXEC2.time(int2a));

temp=polyfit([int1a:int1b]',AXID1.LAX(int1a:int1b),1);temp2=polyval(temp,int1a:int1b);
AXID1.x_int1=(temp2(end)-temp2(1))/(AXID1.time(int1b)-AXID1.time(int1a));
temp=polyfit([int1a:int1b]',AXID1.LAY(int1a:int1b),1);temp2=polyval(temp,int1a:int1b);
AXID1.y_int1=(temp2(end)-temp2(1))/(AXID1.time(int1b)-AXID1.time(int1a));
temp=polyfit([int2a:int2b]',AXID1.LAX(int2a:int2b),1);temp2=polyval(temp,int2a:int2b);
AXID1.x_int2=(temp2(end)-temp2(1))/(AXID1.time(int2b)-AXID1.time(int2a));
temp=polyfit([int2a:int2b]',AXID1.LAY(int2a:int2b),1);temp2=polyval(temp,int2a:int2b);
AXID1.y_int2=(temp2(end)-temp2(1))/(AXID1.time(int2b)-AXID1.time(int2a));

temp=polyfit([int1a:int1b]',ASHES.LAX(int1a:int1b),1);temp2=polyval(temp,int1a:int1b);
ASHES.x_int1=(temp2(end)-temp2(1))/(ASHES.time(int1b)-ASHES.time(int1a));
temp=polyfit([int1a:int1b]',ASHES.LAY(int1a:int1b),1);temp2=polyval(temp,int1a:int1b);
ASHES.y_int1=(temp2(end)-temp2(1))/(ASHES.time(int1b)-ASHES.time(int1a));
temp=polyfit([int2a:int2b]',ASHES.LAX(int2a:int2b),1);temp2=polyval(temp,int2a:int2b);
ASHES.x_int2=(temp2(end)-temp2(1))/(ASHES.time(int2b)-ASHES.time(int2a));
temp=polyfit([int2a:int2b]',ASHES.LAY(int2a:int2b),1);temp2=polyval(temp,int2a:int2b);
ASHES.y_int2=(temp2(end)-temp2(1))/(ASHES.time(int2b)-ASHES.time(int2a));

%difference inflation relative to AXEC2, add sign reversal (decreased pressure = uplift)
AXCC1.p_dif=AXEC2.BDO-interp1(AXCC1.p_time,AXCC1.BDO,AXEC2.p_time);
AXEC2.p_dif=AXEC2.BDO*0;
AXID1.p_dif=AXEC2.BDO-interp1(AXID1.p_time,AXID1.BDO,AXEC2.p_time);
ASHES.p_dif=AXEC2.BDO-interp1(ASHES.p_time,ASHES.BDO,AXEC2.p_time);

%cut out high-amplitude signals/blips
AXCC1.p_dif(abs(AXCC1.p_dif-nanmedian(AXCC1.p_dif))>4)=nan;
AXEC2.p_dif(abs(AXEC2.p_dif-median(AXEC2.p_dif))>4)=nan;
AXID1.p_dif(abs(AXID1.p_dif-nanmedian(AXID1.p_dif))>4)=nan;
ASHES.p_dif(abs(ASHES.p_dif-nanmedian(ASHES.p_dif))>4)=nan;

%12-hr median filter to smooth remaining blips
AXCC1.p_dif=medfilt1(AXCC1.p_dif,60*12,'truncate');
AXEC2.p_dif=medfilt1(AXEC2.p_dif,60*12,'truncate');
AXID1.p_dif=medfilt1(AXID1.p_dif,60*12,'truncate');
ASHES.p_dif=medfilt1(ASHES.p_dif,60*12,'truncate');

%inflation rates
inan=~isnan(AXCC1.p_dif(int1a:int1b)); int1=int1a:int1b; int1=int1(inan);
temp=polyfit([int1]',AXCC1.p_dif(int1),1);temp2=polyval(temp,int1);
AXCC1.p_int1=(temp2(end)-temp2(1))/(AXCC1.time(int1b)-AXCC1.time(int1a));
inan=~isnan(AXCC1.p_dif(int2a:int2b)); int2=int2a:int2b; int2=int2(inan);
temp=polyfit([int2]',AXCC1.p_dif(int2),1);temp2=polyval(temp,int2);
AXCC1.p_int2=(temp2(end)-temp2(1))/(AXCC1.time(int2b)-AXCC1.time(int2a));

inan=~isnan(AXEC2.p_dif(int1a:int1b)); int1=int1a:int1b; int1=int1(inan);
temp=polyfit([int1]',AXEC2.p_dif(int1),1);temp2=polyval(temp,int1);
AXEC2.p_int1=(temp2(end)-temp2(1))/(AXEC2.time(int1b)-AXEC2.time(int1a));
inan=~isnan(AXEC2.p_dif(int2a:int2b)); int2=int2a:int2b; int2=int2(inan);
temp=polyfit([int2]',AXEC2.p_dif(int2),1);temp2=polyval(temp,int2);
AXEC2.p_int2=(temp2(end)-temp2(1))/(AXEC2.time(int2b)-AXEC2.time(int2a));

inan=~isnan(AXID1.p_dif(int1a:int1b)); int1=int1a:int1b; int1=int1(inan);
temp=polyfit([int1]',AXID1.p_dif(int1),1);temp2=polyval(temp,int1);
AXID1.p_int1=(temp2(end)-temp2(1))/(AXID1.time(int1b)-AXID1.time(int1a));
inan=~isnan(AXID1.p_dif(int2a:int2b)); int2=int2a:int2b; int2=int2(inan);
temp=polyfit([int2]',AXID1.p_dif(int2),1);temp2=polyval(temp,int2);
AXID1.p_int2=(temp2(end)-temp2(1))/(AXID1.time(int2b)-AXID1.time(int2a));

inan=~isnan(ASHES.p_dif(int1a:int1b)); int1=int1a:int1b; int1=int1(inan);
temp=polyfit([int1]',ASHES.p_dif(int1),1);temp2=polyval(temp,int1);
ASHES.p_int1=(temp2(end)-temp2(1))/(ASHES.time(int1b)-ASHES.time(int1a));
inan=~isnan(ASHES.p_dif(int2a:int2b)); int2=int2a:int2b; int2=int2(inan);
temp=polyfit([int2]',ASHES.p_dif(int2),1);temp2=polyval(temp,int2);
ASHES.p_int2=(temp2(end)-temp2(1))/(ASHES.time(int2b)-ASHES.time(int2a));

%tidally filtered, non-differenced pressures
stastr=datestr(AXCC1.p_time(1),0); stastr=[stastr(1:17) ':00'];
endstr=datestr(AXCC1.p_time(end),0); endstr=[endstr(1:17) ':00'];
ttemp=[datenum(stastr):(24*60)^-1:datenum(endstr)]';
ttemp_m=round(ttemp*24*60);
p_time_m=round(AXCC1.p_time*24*60);
[~,ia,ib]=intersect(p_time_m,ttemp_m);
nanvec=nan(size(ttemp_m));
nanvec(ib)=AXCC1.BDO(ia);
AXCC1.t_even=ttemp;
AXCC1.p_even=nanvec;
AXCC1.p_even_filt=Z_godin(nanvec,1/60);
AXCC1.p_even_filt_filled=inpaint_nans(AXCC1.p_even_filt,3);

stastr=datestr(AXEC2.p_time(1),0); stastr=[stastr(1:17) ':00'];
endstr=datestr(AXEC2.p_time(end),0); endstr=[endstr(1:17) ':00'];
ttemp=[datenum(stastr):(24*60)^-1:datenum(endstr)]';
ttemp_m=round(ttemp*24*60);
p_time_m=round(AXEC2.p_time*24*60);
[~,ia,ib]=intersect(p_time_m,ttemp_m);
nanvec=nan(size(ttemp_m));
nanvec(ib)=AXEC2.BDO(ia);
AXEC2.t_even=ttemp;
AXEC2.p_even=nanvec;
AXEC2.p_even_filt=Z_godin(nanvec,1/60);
AXEC2.p_even_filt_filled=inpaint_nans(AXEC2.p_even_filt,3);

stastr=datestr(ASHES.p_time(1),0); stastr=[stastr(1:17) ':00'];
endstr=datestr(ASHES.p_time(end),0); endstr=[endstr(1:17) ':00'];
ttemp=[datenum(stastr):(24*60)^-1:datenum(endstr)]';
ttemp_m=round(ttemp*24*60);
p_time_m=round(ASHES.p_time*24*60);
[~,ia,ib]=intersect(p_time_m,ttemp_m);
nanvec=nan(size(ttemp_m));
nanvec(ib)=ASHES.BDO(ia);
ASHES.t_even=ttemp;
ASHES.p_even=nanvec;
ASHES.p_even_filt=Z_godin(nanvec,1/60);
ASHES.p_even_filt_filled=inpaint_nans(ASHES.p_even_filt,3);

stastr=datestr(AXID1.p_time(1),0); stastr=[stastr(1:17) ':00'];
endstr=datestr(AXID1.p_time(end),0); endstr=[endstr(1:17) ':00'];
ttemp=[datenum(stastr):(24*60)^-1:datenum(endstr)]';
ttemp_m=round(ttemp*24*60);
p_time_m=round(AXID1.p_time*24*60);
[~,ia,ib]=intersect(p_time_m,ttemp_m);
nanvec=nan(size(ttemp_m));
nanvec(ib)=AXID1.BDO(ia);
AXID1.t_even=ttemp;
AXID1.p_even=nanvec;
AXID1.p_even_filt=Z_godin(nanvec,1/60);
AXID1.p_even_filt_filled=inpaint_nans(AXID1.p_even_filt,3);

%compass corrections
crd_rot=[cosd(CCMP(4).plus_y) sind(CCMP(4).plus_y); -sind(CCMP(4).plus_y) cosd(CCMP(4).plus_y)];
temp3=crd_rot*[AXCC1.x_int1;AXCC1.y_int1]; AXCC1.x_rot1=temp3(1); AXCC1.y_rot1=temp3(2);
temp4=crd_rot*[AXCC1.x_int2;AXCC1.y_int2]; AXCC1.x_rot2=temp4(1); AXCC1.y_rot2=temp4(2);
temp5=crd_rot*[AXCC1.LAX';AXCC1.LAY']; AXCC1.LAX_rot=temp5(1,:); AXCC1.LAY_rot=temp5(2,:);

crd_rot=[cosd(CCMP(3).plus_y) sind(CCMP(3).plus_y); -sind(CCMP(3).plus_y) cosd(CCMP(3).plus_y)];
temp3=crd_rot*[AXEC2.x_int1;AXEC2.y_int1]; AXEC2.x_rot1=temp3(1); AXEC2.y_rot1=temp3(2);
temp4=crd_rot*[AXEC2.x_int2;AXEC2.y_int2]; AXEC2.x_rot2=temp4(1); AXEC2.y_rot2=temp4(2);
temp5=crd_rot*[AXEC2.LAX';AXEC2.LAY']; AXEC2.LAX_rot=temp5(1,:); AXEC2.LAY_rot=temp5(2,:);

crd_rot=[cosd(CCMP(2).plus_y) sind(CCMP(2).plus_y); -sind(CCMP(2).plus_y) cosd(CCMP(2).plus_y)];
temp3=crd_rot*[AXID1.x_int1;AXID1.y_int1]; AXID1.x_rot1=temp3(1); AXID1.y_rot1=temp3(2);
temp4=crd_rot*[AXID1.x_int2;AXID1.y_int2]; AXID1.x_rot2=temp4(1); AXID1.y_rot2=temp4(2);
temp5=crd_rot*[AXID1.LAX';AXID1.LAY']; AXID1.LAX_rot=temp5(1,:); AXID1.LAY_rot=temp5(2,:);

crd_rot=[cosd(CCMP(1).plus_y) sind(CCMP(1).plus_y); -sind(CCMP(1).plus_y) cosd(CCMP(1).plus_y)];
temp3=crd_rot*[ASHES.x_int1;ASHES.y_int1]; ASHES.x_rot1=temp3(1); ASHES.y_rot1=temp3(2);
temp4=crd_rot*[ASHES.x_int2;ASHES.y_int2]; ASHES.x_rot2=temp4(1); ASHES.y_rot2=temp4(2);
temp5=crd_rot*[ASHES.LAX';ASHES.LAY']; ASHES.LAX_rot=temp5(1,:); ASHES.LAY_rot=temp5(2,:);


%% plotting

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
    
%plot tilt direction changes
figure(24)
clf
plot(caldera_loc(:,1),caldera_loc(:,2),'k')
hold on
plot(AXCC1.x,AXCC1.y,'k^','markersize',15,'linewidth',2)
plot(AXID1.x,AXID1.y,'k^','markersize',15,'linewidth',2)
plot(AXEC2.x,AXEC2.y,'k^','markersize',15,'linewidth',2)
plot(ASHES.x,ASHES.y,'k^','markersize',15,'linewidth',2)
quiver(AXCC1.x,AXCC1.y,AXCC1.x_rot1*2,AXCC1.y_rot1*2,'b','linewidth',2)
quiver(AXEC2.x,AXEC2.y,AXEC2.x_rot1*2,AXEC2.y_rot1*2,'b','linewidth',2)
quiver(AXID1.x,AXID1.y,AXID1.x_rot1*2,AXID1.y_rot1*2,'b','linewidth',2)
h1=quiver(ASHES.x,ASHES.y,ASHES.x_rot1*2,ASHES.y_rot1*2,'b','linewidth',2);
quiver(AXCC1.x,AXCC1.y,AXCC1.x_rot2*2,AXCC1.y_rot2*2,'r','linewidth',2)
quiver(AXEC2.x,AXEC2.y,AXEC2.x_rot2*2,AXEC2.y_rot2*2,'r','linewidth',2)
quiver(AXID1.x,AXID1.y,AXID1.x_rot2*2,AXID1.y_rot2*2,'r','linewidth',2)
h2=quiver(ASHES.x,ASHES.y,ASHES.x_rot2*2,ASHES.y_rot2*2,'r','linewidth',2);
axis equal
legend([h1,h2],{'before','after'})
set(gca,'fontsize',18)
title(['Tilt change ' fldr],'fontsize',20)

orient tall
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print(['../tiltcompare/inflation_reversal/' fldr '/mapview_tilt'],'-dtiff')

%plot tilts
figure(23)
clf
subplot(411)
plot(AXCC1.time,AXCC1.LAX_rot-nanmean(AXCC1.LAX_rot),'linewidth',1)
hold on
plot(AXCC1.time,AXCC1.LAY_rot-nanmean(AXCC1.LAY_rot),'linewidth',1)
legend('x','y','location','southeast')
datetick('x',6)
set(gca,'xticklabel',[])
set(gca,'fontsize',18)
title('Central Caldera','fontsize',20)
subplot(412)
plot(AXEC2.time,AXEC2.LAX_rot-nanmean(AXEC2.LAX_rot),'linewidth',1)
hold on
plot(AXEC2.time,AXEC2.LAY_rot-nanmean(AXEC2.LAY_rot),'linewidth',1)
legend('x','y','location','southeast')
a=xlim;
xlim([a(1) a(2)])
datetick('x',6)
set(gca,'xticklabel',[])
set(gca,'fontsize',18)
ylabel('tilt (\murad)')
title('Eastern Caldera','fontsize',20)
subplot(413)
plot(ASHES.time,ASHES.LAX_rot-nanmean(ASHES.LAX_rot),'linewidth',1)
hold on
plot(ASHES.time,ASHES.LAY_rot-nanmean(ASHES.LAY_rot),'linewidth',1)
legend('x','y','location','southeast')
a=xlim;
xlim([a(1) a(2)])
datetick('x',6)
set(gca,'xticklabel',[])
set(gca,'fontsize',18)
title('ASHES','fontsize',20)
subplot(414)
plot(AXID1.time,AXID1.LAX_rot-nanmean(AXID1.LAX_rot),'linewidth',1)
hold on
plot(AXID1.time,AXID1.LAY_rot-nanmean(AXID1.LAY_rot),'linewidth',1)
legend('x','y','location','southeast')
datetick('x',6)
set(gca,'fontsize',18)
title('International District','fontsize',20)

%annotations
subplot(411)
b=ylim;
patch([AXCC1.time(int1a) AXCC1.time(int1b) AXCC1.time(int1b)...
    AXCC1.time(int1a)],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
    'linewidth',1)
patch([AXCC1.time(int2a) AXCC1.time(int2b) AXCC1.time(int2b)...
    AXCC1.time(int2a)],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
    'linewidth',1)
text(mean([AXCC1.time(int1a) AXCC1.time(int1b)]),b(1)+(b(2)-b(1))/4,...
    [num2str(round(sqrt(AXCC1.x_rot1^2+AXCC1.y_rot1^2),2)) ' \murad/day'],'fontsize',12)
text(mean([AXCC1.time(int2a) AXCC1.time(int2b)]),b(1)+(b(2)-b(1))/4,...
    [num2str(round(sqrt(AXCC1.x_rot2^2+AXCC1.y_rot2^2),2)) ' \murad/day'],'fontsize',12)
subplot(412)
b=ylim;
patch([AXEC2.time(int1a) AXEC2.time(int1b) AXEC2.time(int1b)...
    AXEC2.time(int1a)],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
    'linewidth',1)
patch([AXEC2.time(int2a) AXEC2.time(int2b) AXEC2.time(int2b)...
    AXEC2.time(int2a)],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
    'linewidth',1)
text(mean([AXEC2.time(int1a) AXEC2.time(int1b)]),b(1)+(b(2)-b(1))/4,...
    [num2str(round(sqrt(AXEC2.x_rot1^2+AXEC2.y_rot1^2),2)) ' \murad/day'],'fontsize',12)
text(mean([AXEC2.time(int2a) AXEC2.time(int2b)]),b(1)+(b(2)-b(1))/4,...
    [num2str(round(sqrt(AXEC2.x_rot2^2+AXEC2.y_rot2^2),2)) ' \murad/day'],'fontsize',12)
subplot(413)
b=ylim;
patch([ASHES.time(int1a) ASHES.time(int1b) ASHES.time(int1b)...
    ASHES.time(int1a)],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
    'linewidth',1)
patch([ASHES.time(int2a) ASHES.time(int2b) ASHES.time(int2b)...
    ASHES.time(int2a)],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
    'linewidth',1)
text(mean([ASHES.time(int1a) ASHES.time(int1b)]),b(1)+(b(2)-b(1))/4,...
    [num2str(round(sqrt(ASHES.x_rot1^2+ASHES.y_rot1^2),2)) ' \murad/day'],'fontsize',12)
text(mean([ASHES.time(int2a) ASHES.time(int2b)]),b(1)+(b(2)-b(1))/4,...
    [num2str(round(sqrt(ASHES.x_rot2^2+ASHES.y_rot2^2),2)) ' \murad/day'],'fontsize',12)
subplot(414)
b=ylim;
patch([AXID1.time(int1a) AXID1.time(int1b) AXID1.time(int1b)...
    AXID1.time(int1a)],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
    'linewidth',1)
patch([AXID1.time(int2a) AXID1.time(int2b) AXID1.time(int2b)...
    AXID1.time(int2a)],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
    'linewidth',1)
text(mean([AXID1.time(int1a) AXID1.time(int1b)]),b(1)+(b(2)-b(1))/4,...
    [num2str(round(sqrt(AXID1.x_rot1^2+AXID1.y_rot1^2),2)) ' \murad/day'],'fontsize',12)
text(mean([AXID1.time(int2a) AXID1.time(int2b)]),b(1)+(b(2)-b(1))/4,...
    [num2str(round(sqrt(AXID1.x_rot2^2+AXID1.y_rot2^2),2)) ' \murad/day'],'fontsize',12)

orient tall
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print(['../tiltcompare/inflation_reversal/' fldr '/tiltplot'],'-dtiff')

%plot pressures
figure(22)
clf
subplot(411)
plot(AXEC2.p_time,AXCC1.p_dif-nanmean(AXCC1.p_dif),'linewidth',1)
hold on
datetick('x',6)
set(gca,'xticklabel',[])
set(gca,'fontsize',18)
title('Central Caldera','fontsize',20)
subplot(412)
plot(AXEC2.p_time,AXEC2.p_dif-nanmean(AXEC2.p_dif),'linewidth',1)
hold on
a=xlim;
datetick('x',6)
set(gca,'xticklabel',[])
set(gca,'fontsize',18)
ylabel('Elevation (cm)')
title('Eastern Caldera','fontsize',20)
subplot(413)
plot(AXEC2.p_time,ASHES.p_dif-nanmean(ASHES.p_dif),'linewidth',1)
hold on
xlim([a(1) a(2)])
datetick('x',6,'keeplimits')
set(gca,'xticklabel',[])
set(gca,'fontsize',18)
title('ASHES','fontsize',20)
subplot(414)
plot(AXEC2.p_time,AXID1.p_dif-nanmean(AXID1.p_dif),'linewidth',1)
hold on
datetick('x',6)
set(gca,'fontsize',18)
title('International District','fontsize',20)

%annotations
subplot(411)
b=ylim;
patch([AXEC2.time(int1a) AXEC2.time(int1b) AXEC2.time(int1b)...
    AXEC2.time(int1a)],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
    'linewidth',1)
patch([AXEC2.time(int2a) AXEC2.time(int2b) AXEC2.time(int2b)...
    AXEC2.time(int2a)],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
    'linewidth',1)
text(mean([AXCC1.time(int1a) AXCC1.time(int1b)]),b(1)+(b(2)-b(1))/4,...
    [num2str(round(AXCC1.p_int1,2)) ' cm/day'],'fontsize',12)
text(mean([AXCC1.time(int2a) AXCC1.time(int2b)]),b(1)+(b(2)-b(1))/4,...
    [num2str(round(AXCC1.p_int2,2)) ' cm/day'],'fontsize',12)
subplot(412)
b=ylim;
patch([AXEC2.time(int1a) AXEC2.time(int1b) AXEC2.time(int1b)...
    AXEC2.time(int1a)],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
    'linewidth',1)
patch([AXEC2.time(int2a) AXEC2.time(int2b) AXEC2.time(int2b)...
    AXEC2.time(int2a)],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
    'linewidth',1)
text(mean([AXEC2.time(int1a) AXEC2.time(int1b)]),b(1)+(b(2)-b(1))/4,...
    [num2str(round(AXEC2.p_int1,2)) ' cm/day'],'fontsize',12)
text(mean([AXEC2.time(int2a) AXEC2.time(int2b)]),b(1)+(b(2)-b(1))/4,...
    [num2str(round(AXEC2.p_int2,2)) ' cm/day'],'fontsize',12)
subplot(413)
b=ylim;
patch([AXEC2.time(int1a) AXEC2.time(int1b) AXEC2.time(int1b)...
    AXEC2.time(int1a)],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
    'linewidth',1)
patch([AXEC2.time(int2a) AXEC2.time(int2b) AXEC2.time(int2b)...
    AXEC2.time(int2a)],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
    'linewidth',1)
text(mean([ASHES.time(int1a) ASHES.time(int1b)]),b(1)+(b(2)-b(1))/4,...
    [num2str(round(ASHES.p_int1,2)) ' cm/day'],'fontsize',12)
text(mean([ASHES.time(int2a) ASHES.time(int2b)]),b(1)+(b(2)-b(1))/4,...
    [num2str(round(ASHES.p_int2,2)) ' cm/day'],'fontsize',12)
subplot(414)
b=ylim;
patch([AXEC2.time(int1a) AXEC2.time(int1b) AXEC2.time(int1b)...
    AXEC2.time(int1a)],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
    'linewidth',1)
patch([AXEC2.time(int2a) AXEC2.time(int2b) AXEC2.time(int2b)...
    AXEC2.time(int2a)],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
    'linewidth',1)
text(mean([AXID1.time(int1a) AXID1.time(int1b)]),b(1)+(b(2)-b(1))/4,...
    [num2str(round(AXID1.p_int1,2)) ' cm/day'],'fontsize',12)
text(mean([AXID1.time(int2a) AXID1.time(int2b)]),b(1)+(b(2)-b(1))/4,...
    [num2str(round(AXID1.p_int2,2)) ' cm/day'],'fontsize',12)

orient tall
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print(['../tiltcompare/inflation_reversal/' fldr '/pressureplot'],'-dtiff')

%plot pressure changes
figure(21)
clf
plot(caldera_loc(:,1),caldera_loc(:,2),'k')
hold on
plot(AXCC1.x,AXCC1.y,'k^','markersize',15,'linewidth',2)
plot(AXID1.x,AXID1.y,'k^','markersize',15,'linewidth',2)
plot(AXEC2.x,AXEC2.y,'k^','markersize',15,'linewidth',2)
plot(ASHES.x,ASHES.y,'k^','markersize',15,'linewidth',2)
quiver(AXCC1.x,AXCC1.y,0,AXCC1.p_int1/abs(AXCC1.p_int1)*3,'b','linewidth',2)
quiver(AXEC2.x,AXEC2.y,0,AXEC2.p_int1/abs(AXCC1.p_int1)*3,'b','linewidth',2)
quiver(AXID1.x,AXID1.y,0,AXID1.p_int1/abs(AXCC1.p_int1)*3,'b','linewidth',2)
h1=quiver(ASHES.x,ASHES.y,0,ASHES.p_int1/abs(AXCC1.p_int1)*3,'b','linewidth',2);
quiver(AXCC1.x,AXCC1.y,0,AXCC1.p_int2/abs(AXCC1.p_int2)*3,'r','linewidth',2)
quiver(AXEC2.x,AXEC2.y,0,AXEC2.p_int2/abs(AXCC1.p_int2)*3,'r','linewidth',2)
quiver(AXID1.x,AXID1.y,0,AXID1.p_int2/abs(AXCC1.p_int2)*3,'r','linewidth',2)
h2=quiver(ASHES.x,ASHES.y,0,ASHES.p_int2/abs(AXCC1.p_int2)*3,'r','linewidth',2);
axis equal
legend([h1,h2],{'before','after'})
set(gca,'fontsize',18)
title(['Pressure change ' fldr],'fontsize',20)

orient tall
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print(['../tiltcompare/inflation_reversal/' fldr '/mapview_pressure'],'-dtiff')

if filtpres
    %plot non-differenced pressures
    figure(25)
    clf
    subplot(411)
    hold on
    plot(AXCC1.t_even,AXCC1.p_even_filt_filled-nanmean(AXCC1.p_even_filt),'r:','linewidth',1)
    plot(AXCC1.t_even,AXCC1.p_even_filt-nanmean(AXCC1.p_even_filt),'b','linewidth',1)
    ylim([-5 5])
    datetick('x',6)
    set(gca,'xticklabel',[])
    set(gca,'fontsize',15)
    title('Central Caldera')
    box on
    
    subplot(412)
    hold on
    plot(AXEC2.t_even,AXEC2.p_even_filt_filled-nanmean(AXEC2.p_even_filt),'r:','linewidth',1)
    plot(AXEC2.t_even,AXEC2.p_even_filt-nanmean(AXEC2.p_even_filt),'b','linewidth',1)
    ylim([-5 5])
    datetick('x',6)
    set(gca,'xticklabel',[])
    set(gca,'fontsize',15)
    title('Eastern Caldera')
    box on
    
    subplot(413)
    hold on
    plot(ASHES.t_even,ASHES.p_even_filt_filled-nanmean(ASHES.p_even_filt),'r:','linewidth',1)
    plot(ASHES.t_even,ASHES.p_even_filt-nanmean(ASHES.p_even_filt),'b','linewidth',1)
    ylim([-5 5])
    datetick('x',6)
    set(gca,'xticklabel',[])
    ylabel('Pressure (cm)');
    set(gca,'fontsize',15)
    title('Ashes Vent Field')
    box on
    
    subplot(414)
    hold on
    plot(AXID1.t_even,AXID1.p_even_filt_filled-nanmean(AXID1.p_even_filt),'r:','linewidth',1)
    plot(AXID1.t_even,AXID1.p_even_filt-nanmean(AXID1.p_even_filt),'b','linewidth',1)
    ylim([-5 5])
    datetick('x',6)
    set(gca,'fontsize',15)
    title('International District')
    box on
    
    orient tall
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print(['../tiltcompare/inflation_reversal/' fldr '/filtered_pressure'],'-dtiff')
    saveas(gcf,['../tiltcompare/inflation_reversal/' fldr '/filtered_pressure.fig'])
end