%compareSCTA_LILY.m
%
% Compares in time and frequency space the tilt measured by the SCTA and
% Lily tiltmeters
%

clear; close all

load('../compass_directions.mat')

t0=datenum(2018,12,9);
tf=datenum(2018,12,15);

%load data
sta={'AXCC1','AXCC2'};
cha1={'LAX','LAY'};
cha2={'BNE','BNN','BNZ','BKA'};

t1=t0;
AXCC1.time=[];AXCC1.LAX=[];AXCC1.LAY=[];
AXCC2.time=[];AXCC2.BNE=[];AXCC2.BNN=[];AXCC2.BNZ=[];AXCC2.BKA=[];
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
    for j=1:4
        file_string=['AXCC2_' cha2{j} '_' t1_s '.miniseed'];
        %attempt download if file not found
        if ~exist(['../tiltcompare/AXCC2/' file_string],'file')
            IRIS_data_pull('AXCC2',cha2{j},'--',t1,t1+1);
        end
        %some dates have no data (power failure, etc.)
        if exist(['../tiltcompare/AXCC2/' file_string],'file')
            temp=rdmseed(['../tiltcompare/AXCC2/' file_string]);
            dtemp=double(cat(1,temp.d));
            if length(dtemp)<3600 %short timeseries will throw errors
                t1=t1+1;
                continue
            end
            dtempd=decimate(dtemp,12,'fir');
            dtempd=decimate(dtempd,10,'fir');
            dtempd=decimate(dtempd,5,'fir');
            dtempd=decimate(dtempd,4,'fir');
            eval(['AXCC2.' cha2{j} '=[AXCC2.' cha2{j} ';dtempd];']);
            if j==1
                ttemp=cat(1,temp.t);
                ttempd=decimate(ttemp,12,'fir');
                ttempd=decimate(ttempd,10,'fir');
                ttempd=decimate(ttempd,5,'fir');
                ttempd=decimate(ttempd,4,'fir');
                AXCC2.time=[AXCC2.time;ttempd];
            end
        end
    end
    t1=t1+1;
end

%AXCC1 compass correction
crd_rot=[cosd(CCMP(4).plus_y) sind(CCMP(4).plus_y); -sind(CCMP(4).plus_y) cosd(CCMP(4).plus_y)];
temp=crd_rot*[AXCC1.LAX';AXCC1.LAY']; AXCC1.LAX_rot=temp(1,:)'; AXCC1.LAY_rot=temp(2,:)';

%AXCC2 temperature correction
AXCC2.BNE=AXCC2.BNE/10^7; AXCC2.BNN=AXCC2.BNN/10^7; AXCC2.BNZ=AXCC2.BNZ/10^7;...
    AXCC2.BKA=AXCC2.BKA/10^7;
me=get_Tdep_params(AXCC2.time,AXCC2.BNE,AXCC2.BKA); AXCC2.BNE=AXCC2.BNE-me(3)*AXCC2.BKA;
mn=get_Tdep_params(AXCC2.time,AXCC2.BNN,AXCC2.BKA); AXCC2.BNN=AXCC2.BNN-mn(3)*AXCC2.BKA;
mz=get_Tdep_params(AXCC2.time,AXCC2.BNZ,AXCC2.BKA); AXCC2.BNZ=AXCC2.BNZ-mz(3)*AXCC2.BKA;

%convert AXCC2 acceleration to tilt
AXCC2.LAX=asin(AXCC2.BNE/9.81)*10^6;
AXCC2.LAY=asin(AXCC2.BNN/9.81)*10^6;