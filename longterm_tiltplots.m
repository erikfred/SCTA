% longterm_tiltplots.m
%
% An improvement on 'monthly_tiltplots.m', discarding the separation of
% data into monthly units and providing data structures at different
% temporal resolution. 
%

clear; close all;

%%%%%%%%%%CONFIG%%%%%%%%%%
axial=true;
lily=false;
pf=false;

loaddata=0; % 0 - start from scratch, 1 - load from pre-move era, 2 - load from post-move era

%-----(un)comment as desired
flipfile='../calibrations/Axial/axialdata.mat'; % only used if loaddata==0
% flipfile='../calibrations/Axial/axialdata_newloc.mat'; % only used if loaddata==0
% flipfile='../calibrations/PF/PFdata.mat'; % only used if loaddata==0

startdate=datenum(2018,10,10); % [only used if loaddata==0]
% startdate=datenum(2020,09,11); % [only used if loaddata==0]
% startdate=datenum(2018,10,18); % [only used if loaddata==0]

tf=datenum('09/10/20'); % OOI SCTA moved 9/11/2020
% tf=datenum('08/27/21'); % OOI SCTA recovered 8/27/2021
% tf=datenum('03/05/20'); % PF SCTA recovered 3/5/2020
%%%%%%%%END CONFIG%%%%%%%%

if axial
    
    % Load pre-exisiting structure, if requested & it exists
    if loaddata==1
        flipfile='../calibrations/Axial/axialdata.mat';
        datafile_hr='../calibrations/Axial/axialdata_hr.mat';
        datafile_min='../calibrations/Axial/axialdata_min.mat';
        load(datafile_hr)
        load(datafile_min)
        t0=ceil(data_min.t(end));
    elseif loaddata==2
        flipfile='../calibrations/Axial/axialdata_newloc.mat';
        datafile_hr='../calibrations/Axial/axialdata_newloc_hr.mat';
        datafile_min='../calibrations/Axial/axialdata_newloc_min.mat';
        load(datafile_hr)
        load(datafile_min)
        t0=ceil(data_min.t(end));
    else
        t0=startdate;
    end

    sta='AXCC2';

    % Determine datenums of calibrations
    load(flipfile,'flipInfoAll')
    [daylist,id,~]=unique(floor(flipInfoAll.t));
    
    if t0==startdate
        fields={'t','MNE','MNN','MNZ','MXG','MKA','iflip'};
        for i=1:length(fields)
            eval(['data_min.' fields{i} '=[];'])
            eval(['data_hr.' fields{i} '=[];'])
        end
    end
    
    t1=t0;
    while t1<tf
        % skip dates when we don't expect data
        if t1>datenum(2020,05,09) && t1<datenum(2020,06,02) % OOI outage
            t1=t1+1;
            continue
        elseif t1>datenum(2021,01,13) && t1<datenum(2021,01,18) % OOI outage
            t1=t1+1;
            continue
        end
        
        % download and read in data from 8Hz channels
        data=[];
        t1_s=datestr(t1,31); t1_s=t1_s(1:10);
        cha={'MNE','MNN','MNZ','MXG','MKA'};
        chastr={'a(:,1)','a(:,2)','a(:,3)','as','T'};
        for m=1:length(cha)
            %   attempt download if file not found
            fstring=[sta '_' cha{m} '_' t1_s '.miniseed'];
            if ~exist(['../tiltcompare/' sta '/' fstring],'file')
                IRIS_data_pull(sta,cha{m},'--',t1,t1+1);
            end
            % this redundancy avoids crashing if files weren't found on IRIS
            if exist(['../tiltcompare/' sta '/' fstring],'file')
                % read file
                temp=rdmseed(['../tiltcompare/' sta '/' fstring]);
                if m==1 % only need to pull time once per data day
                    data.t=cat(1,temp.t);
                    eval(['data.' chastr{m} '=cat(1,temp.d)/10^7;']);
                else
                    atemp=cat(1,temp.d)/10^7;
                    if length(atemp)==length(data.a)
                        eval(['data.' chastr{m} '=atemp;']);
                    else
                        warning('Length discrepancy')
                        keyboard
                    end
                end
            end
        end
        
        if isempty(data)
            t1=t1+1;
            continue
        end
        
        % if day of move, remove first half day of data
        if t1==datenum(2020,09,11)
            cond=data.t>datenum(2020,09,11,12,0,0);
            data.t=data.t(cond);
            data.a=data.a(cond,:);
            data.as=data.as(cond);
            data.T=data.T(cond);
        end
        
        % decimate to 1 sample/min
        dataDecMin=decimate_SCTA(data,60);
        dataDecHr=decimate_SCTA(data,3600);
        
        % note calibration signal
        [checkcal,iflip]=max(abs(dataDecMin.a(:,1)));
        if checkcal>2
            i1=find(dataDecMin.t-floor(dataDecMin.t(1))>=datenum(0,0,0,20,58,0),1);
            i2=find(dataDecMin.t-floor(dataDecMin.t(1))>=datenum(0,0,0,21,11,0),1);
            dataDecMin.iflip=(i1:i2)'+length(data_min.t);
            i1=find(dataDecHr.t-floor(dataDecHr.t(1))>=datenum(0,0,0,20,0,0),1);
            i2=find(dataDecHr.t-floor(dataDecHr.t(1))>=datenum(0,0,0,22,0,0),1);
            dataDecHr.iflip=(i1:i2-1)'+length(data_hr.t);
        else
            dataDecMin.iflip=[];
            dataDecHr.iflip=[];
        end
        
        % append
        fields={'t','MNE','MNN','MNZ','MXG','MKA','iflip'};
        chastr={'t','a(:,1)','a(:,2)','a(:,3)','as','T','iflip'};
        for i=1:length(fields)
            eval(['data_min.' fields{i} '=[data_min.' fields{i} '; dataDecMin.' chastr{i} '];'])
            eval(['data_hr.' fields{i} '=[data_hr.' fields{i} '; dataDecHr.' chastr{i} '];'])
        end
        
%         if any(t1==floor(flipInfoAll.t))
%             keyboard
%         end
        
        % move on to next day
        t1=t1+1;
    end
    
    save(insertBefore(flipfile,'.mat','_hr'),'data_hr')
    save(insertBefore(flipfile,'.mat','_min'),'data_min')
end

if lily
    sta='AXCC1';
    
    % Load pre-exisiting structure, if it exists
    if loaddata>0 && exist('../tiltcompare/SCTA_Lily_comp/AXCC1data_hr.mat','file')
        load('../tiltcompare/SCTA_Lily_comp/AXCC1data_hr.mat')
        load('../tiltcompare/SCTA_Lily_comp/AXCC1data_min.mat')
        t0=ceil(lily_min.t(end));
    else
        t0=datenum(2018,10,13);
    end
    
    if t0==datenum(2018,10,13)
        lily_min.t=[];lily_min.LAX=[];lily_min.LAY=[];lily_min.BDO=[];
        lily_hr.t=[];lily_hr.LAX=[];lily_hr.LAY=[];lily_hr.BDO=[];
    end
    
    t1=t0;
    while t1<tf
        if t1>datenum(2020,05,09) && t1<datenum(2020,06,02) %OOI outage
            t1=t1+1;
            continue
        end
        
        t1_s=datestr(t1,31); t1_s=t1_s(1:10);
        LAY_string=[sta '_LAY_' t1_s '.miniseed'];
        LAX_string=[sta '_LAX_' t1_s '.miniseed'];
        BDO_string=[sta '_BDO_' t1_s '.miniseed'];
        
        %attempt download if file not found
        if ~exist(['../tiltcompare/' sta '/' LAY_string],'file') || ...
                ~exist(['../tiltcompare/' sta '/' LAX_string],'file') || ...
                ~exist(['../tiltcompare/' sta '/' BDO_string],'file')
            IRIS_data_pull(sta,'LAY','11',t1,t1+1);
            IRIS_data_pull(sta,'LAX','11',t1,t1+1);
            IRIS_data_pull(sta,'BDO','11',t1,t1+1);
        end
        %some dates have no data (power failure, etc.)
        if exist(['../tiltcompare/' sta '/' LAY_string],'file') && ...
                exist(['../tiltcompare/' sta '/' LAX_string],'file') && ...
                exist(['../tiltcompare/' sta '/' BDO_string],'file')
            
            %LAY channel
            temp=rdmseed(['../tiltcompare/' sta '/' LAY_string]);
            ttemp=cat(1,temp.t);
            ntemp=double(cat(1,temp.d));
            %downsample to 1 sample/min
            [ttempd,ntempd,~,~]=downsample_uneven5(ttemp,ntemp,1/(24*60));
            
            %LAX channel
            temp=rdmseed(['../tiltcompare/' sta '/' LAX_string]);
            ttemp=cat(1,temp.t);
            etemp=double(cat(1,temp.d));
            %downsample to 1 sample/min
            [~,etempd,~,~]=downsample_uneven5(ttemp,etemp,1/(24*60));
            
            %BDO channel
            temp=rdmseed(['../tiltcompare/' sta '/' BDO_string]);
            ttemp=cat(1,temp.t);
            ztemp=double(cat(1,temp.d));
            %downsample to 1 sample/min
            [ttest,ztempd,~,~]=downsample_uneven5(ttemp,ztemp,1/(24*60));

            %append
            if length(ttempd)==length(etempd) && length(ttempd)==length(ntempd) ...
                    && length(ttempd)==length(ztempd)
                lily_min.LAY=[lily_min.LAY; ntempd]; % [urad]
                lily_min.LAX=[lily_min.LAX; etempd]; % [urad]
                lily_min.BDO=[lily_min.BDO; ztempd/(1.45038e-4)]; % [Pa]
                lily_min.t=[lily_min.t; ttempd];
            else % auto-fix missing data at start or end
                tilttime1=datetime(ttempd(1),'convertfrom','datenum');
                prestime1=datetime(ttest(1),'convertfrom','datenum');
                tilttime2=datetime(ttempd(end),'convertfrom','datenum');
                prestime2=datetime(ttest(end),'convertfrom','datenum');
                if length(ztempd)-length(etempd)==1
                    ztempd(end)=[];
                elseif tilttime1~=prestime1
                    ldif=length(etempd)-length(ztempd);
                    ztempd=[NaN(ldif,1);ztempd];
                    ztempd(length(ttempd)+1:end)=[];
                elseif tilttime2~=prestime2
                    ldif=length(etempd)-length(ztempd);
                    ztempd=[ztempd;NaN(ldif,1)];
                    ztempd(length(ttempd)+1:end)=[];
                end
                % check if auto-fix worked
                if length(ttempd)==length(etempd) && length(ttempd)==length(ntempd) ...
                        && length(ttempd)==length(ztempd)
                    lily_min.LAY=[lily_min.LAY; ntempd]; % [urad]
                    lily_min.LAX=[lily_min.LAX; etempd]; % [urad]
                    lily_min.BDO=[lily_min.BDO; ztempd/(1.45038e-4)]; % [Pa]
                    lily_min.t=[lily_min.t; ttempd];
                else % keyboard control for more complicated issues
                    warning('vectors have different lengths')
                    keyboard
                end
            end
            
            % hourly version
            if length(ntempd)>620 %minimum length for decimation filter
                %downsample all to a sample/hr
                [~,ntempd,~,~]=downsample_uneven5(ttempd,ntempd,1/(24));
                [~,etempd,~,~]=downsample_uneven5(ttempd,etempd,1/(24));
                [ttempd,ztempd,~,~]=downsample_uneven5(ttempd,ztempd,1/(24));
                %append
                if length(ttempd)==length(etempd) && length(ttempd)==length(ntempd) ...
                        && length(ttempd)==length(ztempd)
                    lily_hr.LAY=[lily_hr.LAY; ntempd]; % [urad]
                    lily_hr.LAX=[lily_hr.LAX; etempd]; % [urad]
                    lily_hr.BDO=[lily_hr.BDO; ztempd/(1.45038e-4)]; % [Pa]
                    lily_hr.t=[lily_hr.t; ttempd];
                else
                    warning('vectors have different lengths')
                    keyboard
                end
            end
        end
        t1=t1+1;
    end
    
    save('../tiltcompare/SCTA_Lily_comp/AXCC1data_hr.mat','lily_hr')
    save('../tiltcompare/SCTA_Lily_comp/AXCC1data_min.mat','lily_min')
end

if pf
    
    % Load pre-exisiting structure, if requested & it exists
    if loaddata>0 && exist('../calibrations/PinonFlat/PFdata_hr.mat','file')
        load('../calibrations/PinonFlat/PFdata_hr.mat')
        load('../calibrations/PinonFlat/PFdata_min.mat')
        t0=ceil(data_min.t(end));
    else
        t0=datenum(2018,10,18);
    end
    
    % Determine datenums of calibrations
    load('../calibrations/PinonFlat/PFdata.mat','flipInfoAll')
    [daylist,id,~]=unique(floor(flipInfoAll.t));
    
    if t0==datenum(2018,10,18)
        data_min.t=[];data_min.MNE=[];data_min.MNN=[];data_min.MNZ=[];
        data_min.MKA=[];data_min.iflip=[];
        data_hr.t=[];data_hr.MNE=[];data_hr.MNN=[];data_hr.MNZ=[];
        data_hr.MKA=[];data_hr.iflip=[];
    end
    
    % import and decimate
    t1=t0;
    while t1<tf
        data=get_sctaDay('/Volumes/GoogleDrive/My Drive/Oceanography/SCTA-Share/OOI-PF/SCTA-PF/ParsedData',t1);
        if isempty(data.t)
            t1=t1+1;
            continue
        end
        
        % x
        xtemp=data.a(:,1);
        %decimate to 1 sample/min
        xtemp1=decimate(xtemp,4,'fir');
        xtemp2=decimate(xtemp1,10,'fir');
        xtemp3=decimate(xtemp2,6,'fir');
        xtempd=decimate(xtemp3,10,'fir');
        
        % y
        ytemp=data.a(:,2);
        %decimate to 1 sample/min
        ytemp1=decimate(ytemp,4,'fir');
        ytemp2=decimate(ytemp1,10,'fir');
        ytemp3=decimate(ytemp2,6,'fir');
        ytempd=decimate(ytemp3,10,'fir');
        
        % z
        ztemp=data.a(:,3);
        %decimate to 1 sample/min
        ztemp1=decimate(ztemp,4,'fir');
        ztemp2=decimate(ztemp1,10,'fir');
        ztemp3=decimate(ztemp2,6,'fir');
        ztempd=decimate(ztemp3,10,'fir');
        
        % T
        Ttemp=data.T;
        %decimate to 1 sample/min
        Ttemp1=decimate(Ttemp,4,'fir');
        Ttemp2=decimate(Ttemp1,10,'fir');
        Ttemp3=decimate(Ttemp2,6,'fir');
        Ttempd=decimate(Ttemp3,10,'fir');
        
        % time
        ttempd=linspace(data.t(1),data.t(1)+1439/1440,length(Ttempd))';
        
        %note calibration signal
        [checkcal,iflip]=max(abs(ytempd));
        if checkcal>9
            data_min.iflip=[data_min.iflip;length(data_min.t)+(iflip-10:iflip+10)'];
        end
        
        %append
        if length(ttempd)==length(xtempd) && length(ttempd)==length(ytempd) ...
                && length(ttempd)==length(ztempd) && length(ttempd)==length(Ttempd)
            data_min.MNN=[data_min.MNN; ytempd];
            data_min.MNE=[data_min.MNE; xtempd];
            data_min.MNZ=[data_min.MNZ; ztempd];
            data_min.MKA=[data_min.MKA; Ttempd];
            data_min.t=[data_min.t; ttempd];
        else
            warning('vectors have different lengths')
            keyboard
        end
        
        % hourly version
        if length(Ttempd)>620 %minimum length for decimation filter
            %decimate all to a sample/hr
            ytempd=decimate(ytempd,10,'fir'); ytempd=decimate(ytempd,6,'fir');
            xtempd=decimate(xtempd,10,'fir'); xtempd=decimate(xtempd,6,'fir');
            ztempd=decimate(ztempd,10,'fir'); ztempd=decimate(ztempd,6,'fir');
            Ttempd=decimate(Ttempd,10,'fir'); Ttempd=decimate(Ttempd,6,'fir');
            ttempd=linspace(data.t(1),data.t(1)+23/24,length(Ttempd))';
            %append
            if length(ttempd)==length(xtempd) && length(ttempd)==length(ytempd) ...
                    && length(ttempd)==length(ztempd) && length(ttempd)==length(Ttempd)
                data_hr.MNN=[data_hr.MNN; ytempd];
                data_hr.MNE=[data_hr.MNE; xtempd];
                data_hr.MNZ=[data_hr.MNZ; ztempd];
                data_hr.MKA=[data_hr.MKA; Ttempd];
                data_hr.t=[data_hr.t; ttempd];
            else
                warning('vectors have different lengths')
                keyboard
            end
            
            %note calibration signal
            [checkcal,iflip]=max(abs(median(ytempd)-ytempd));
            if checkcal>0.01
                data_hr.iflip=[data_hr.iflip;length(data_hr.t)-length(Ttempd)+(iflip-2:length(Ttempd))'];
            end
        end
        t1=t1+1;
    end
    
    save('../calibrations/PinonFlat/PFdata_hr.mat','data_hr')
    save('../calibrations/PinonFlat/PFdata_min.mat','data_min')
end