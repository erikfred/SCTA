% longterm_tiltplots.m
%
% ***WORK IN PROGRESS***
%
% An improvement on 'monthly_tiltplots.m', discarding the separation of
% data into monthly units and providing data structures at different
% temporal resolution. 
%

clear; close all;

%%%%%%%%%%CONFIG%%%%%%%%%%
loaddata=true;
tf=datenum(date);

sta='AXCC2';
dec=[8 10 6 10 6];
%%%%%%%%END CONFIG%%%%%%%%

% Load pre-exisiting structure, if it exists
if loaddata && exist('../calibrations/Axial/axialdata_hr.mat','file')
    load('../calibrations/Axial/axialdata_hr.mat')
    load('../calibrations/Axial/axialdata_min.mat')
    t0=ceil(data_min.t(end));
else
    t0=datenum(2018,10,13);
end

% Determine datenums of calibrations
load('../calibrations/Axial/axialdata.mat','flipInfoAll')
[daylist,id,~]=unique(floor(flipInfoAll.t));

if t0==datenum(2018,10,13)
    data_min.t=[];data_min.MNE=[];data_min.MNN=[];data_min.MNZ=[];
    data_min.MKA=[];data_min.iflip=[];
    data_hr.t=[];data_hr.MNE=[];data_hr.MNN=[];data_hr.MNZ=[];
    data_hr.MKA=[];data_hr.iflip=[];
end

t1=t0;
while t1<tf
    t1_s=datestr(t1,31); t1_s=t1_s(1:10);
    MNN_string=[sta '_MNN_' t1_s '.miniseed'];
    MNE_string=[sta '_MNE_' t1_s '.miniseed'];
    MNZ_string=[sta '_MNZ_' t1_s '.miniseed'];
    MKA_string=[sta '_MKA_' t1_s '.miniseed'];

    %attempt download if file not found
    if ~exist(['../tiltcompare/' sta '/' MNN_string],'file') || ...
            ~exist(['../tiltcompare/' sta '/' MNE_string],'file') || ...
            ~exist(['../tiltcompare/' sta '/' MNZ_string],'file') || ...
            ~exist(['../tiltcompare/' sta '/' MKA_string],'file')
        IRIS_data_pull(sta,'MNN','--',t1,t1+1);
        IRIS_data_pull(sta,'MNE','--',t1,t1+1);
        IRIS_data_pull(sta,'MNZ','--',t1,t1+1);
        IRIS_data_pull(sta,'MKA','--',t1,t1+1);
    end
    %some dates have no data (power failure, etc.)
    if exist(['../tiltcompare/' sta '/' MNN_string],'file') && ...
            exist(['../tiltcompare/' sta '/' MNE_string],'file') && ...
            exist(['../tiltcompare/' sta '/' MNZ_string],'file') && ...
            exist(['../tiltcompare/' sta '/' MKA_string],'file')
        %MNN channel
        temp=rdmseed(['../tiltcompare/' sta '/' MNN_string]);
        ntemp=double(cat(1,temp.d));
        %decimate to 1 sample/min
        ntempd=decimate(ntemp,dec(1),'fir');
        ntempd=decimate(ntempd,dec(2),'fir');
        ntempd=decimate(ntempd,dec(3),'fir');
        %note calibration signal
        [checkcal,iflip]=max(abs(ntempd));
        if checkcal/10^7>9
            data_min.iflip=[data_min.iflip;length(data_min.t)+(iflip-10:iflip+10)'];
        end
        
        %MNE channel
        temp=rdmseed(['../tiltcompare/' sta '/' MNE_string]);
        etemp=double(cat(1,temp.d));
        %decimate to 1 sample/min
        etempd=decimate(etemp,dec(1),'fir');
        etempd=decimate(etempd,dec(2),'fir');
        etempd=decimate(etempd,dec(3),'fir');        
        
        %MNZ channel
        temp=rdmseed(['../tiltcompare/' sta '/' MNZ_string]);
        ztemp=double(cat(1,temp.d));
        %decimate to 1 sample/min
        ztempd=decimate(ztemp,dec(1),'fir');
        ztempd=decimate(ztempd,dec(2),'fir');
        ztempd=decimate(ztempd,dec(3),'fir');
        
        %MKA channel
        temp=rdmseed(['../tiltcompare/' sta '/' MKA_string]);
        Ttemp=double(cat(1,temp.d));
        %decimate to 1 sample/min
        Ttempd=decimate(Ttemp,dec(1),'fir');
        Ttempd=decimate(Ttempd,dec(2),'fir');
        Ttempd=decimate(Ttempd,dec(3),'fir');
        
        %time
        ttemp=cat(1,temp.t);
        ttempd=decimate(ttemp,dec(1),'fir');
        ttempd=decimate(ttempd,dec(2),'fir');
        ttempd=decimate(ttempd,dec(3),'fir');
        
        %append
        if length(ttempd)==length(etempd) && length(ttempd)==length(ntempd) ...
                && length(ttempd)==length(ztempd) && length(ttempd)==length(Ttempd)
            data_min.MNN=[data_min.MNN; ntempd/10^7];
            data_min.MNE=[data_min.MNE; etempd/10^7];
            data_min.MNZ=[data_min.MNZ; ztempd/10^7];
            data_min.MKA=[data_min.MKA; Ttempd/10^7];
            data_min.t=[data_min.t; ttempd];
        else
            warning('vectors have different lengths')
            keyboard
        end
        
        % hourly version
        if length(ntempd)>620 %minimum length for decimation filter
            %decimate all to a sample/hr
            ntempd=decimate(ntempd,dec(4),'fir'); ntempd=decimate(ntempd,dec(5),'fir');
            etempd=decimate(etempd,dec(4),'fir'); etempd=decimate(etempd,dec(5),'fir');
            ztempd=decimate(ztempd,dec(4),'fir'); ztempd=decimate(ztempd,dec(5),'fir');
            Ttempd=decimate(Ttempd,dec(4),'fir'); Ttempd=decimate(Ttempd,dec(5),'fir');
            ttempd=decimate(ttempd,dec(4),'fir'); ttempd=decimate(ttempd,dec(5),'fir');
            %append
            if length(ttempd)==length(etempd) && length(ttempd)==length(ntempd) ...
                    && length(ttempd)==length(ztempd) && length(ttempd)==length(Ttempd)
                data_hr.MNN=[data_hr.MNN; ntempd/10^7];
                data_hr.MNE=[data_hr.MNE; etempd/10^7];
                data_hr.MNZ=[data_hr.MNZ; ztempd/10^7];
                data_hr.MKA=[data_hr.MKA; Ttempd/10^7];
                data_hr.t=[data_hr.t; ttempd];
            else
                warning('vectors have different lengths')
                keyboard
            end
            
            %note calibration signal
            [checkcal,iflip]=max(abs(median(ntempd)-ntempd));
            if checkcal/10^7>0.01
                data_hr.iflip=[data_hr.iflip;length(data_hr.t)-length(ntempd)+(iflip-2:length(ntempd))'];
            end
        end
    end
    t1=t1+1;
end

save('../calibrations/Axial/axialdata_hr.mat','data_hr')
save('../calibrations/Axial/axialdata_min.mat','data_min')