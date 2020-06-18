% tide_check.m
%
% Compares tidal signals between Axial SCTA and BOPT sensor at AXCC1. Also
% generates plots for the other BOPT locations. Will eventually include
% comparison to a tidal model. Should allow us to resolve
% orientation/polarity ambiguity between instruments.
%

clear; close all

load('../compass_directions.mat')
CCMP=CCMP([4 3 2 1]);

t0=datenum(2019,06,9);
tf=datenum(2019,06,19);

%load data
sta1={'AXCC2'};
cha1={'MNE','MNN'};
sta2={'AXCC1','AXEC2','AXID1'};
cha2={'LAX','LAY'};

t1=t0;
AXCC1.time=[];AXCC1.LAX=[];AXCC1.LAY=[];
AXEC2.time=[];AXEC2.LAX=[];AXEC2.LAY=[];
AXID1.time=[];AXID1.LAX=[];AXID1.LAY=[];
ASHES.time=[];ASHES.LAX=[];ASHES.LAY=[];
AXCC2.time=[];AXCC2.BNE=[];AXCC2.BNN=[];
while t1<=tf
    t1_s=datestr(t1,31); t1_s=t1_s(1:10);
    
    % grab all BOPT data
    for i=1:length(sta2)
        for j=1:length(cha2)
            IRIS_data_pull(sta2{i},cha2{j},'11',t1,t1+1)
            if ~exist(['../tiltcompare/' sta2{i} '/' sta2{i} '_' cha2{j} '_' datestr(t1,29) '.miniseed'],'file')
                if ~exist(['../tiltcompare/' sta2{i} '/' sta2{i} '_' cha2{j} '_' datestr(t1,29) '.miniseed'],'file')
                    IRIS_data_pull(sta2{i},cha2{j},'11',t1,t1+1);
                end
            end
            temp=rdmseed(['../tiltcompare/' sta2{i} '/' sta2{i} '_' cha2{j} '_' datestr(t1,29) '.miniseed']);
            data.t=cat(1,temp.t);
            data.a(:,j)=cat(1,temp.d);
        end
        % append to appropriate structure
        eval([sta2{i} '.time=[' sta2{i} '.time; data.t];']);
        eval([sta2{i} '.LAX=[' sta2{i} '.LAX; data.a(:,1)];']);
        eval([sta2{i} '.LAY=[' sta2{i} '.LAY; data.a(:,2)];']);
        clear data
    end
    
    % grab AXCC2 (SCTA) data
    for k=1:length(cha1)
        IRIS_data_pull('AXCC2',cha1{k},'--',t1,t1+1)
        if ~exist(['../tiltcompare/AXCC2/AXCC2_' cha1{k} '_' datestr(t1,29) '.miniseed'],'file')
            if ~exist(['../tiltcompare/AXCC2/AXCC2_' cha1{k} '_' datestr(t1,29) '.miniseed'],'file')
                IRIS_data_pull('AXCC2',cha1{k},'--',t1,t1+1);
            end
        end
        temp=rdmseed(['../tiltcompare/AXCC2/AXCC2_' cha1{k} '_' datestr(t1,29) '.miniseed']);
        data.t=cat(1,temp.t);
        data.a(:,k)=cat(1,temp.d)/10^7;
    end
    % append to AXCC2
    AXCC2.time=[AXCC2.time; data.t];
    AXCC2.BNE=[AXCC2.BNE; data.a(:,1)];
    AXCC2.BNN=[AXCC2.BNN; data.a(:,2)];
    clear data
    
    t1=t1+1;
end

% BOPT compass correction
for l=1:length(sta2)
    crd_rot=[cosd(CCMP(l).plus_y) sind(CCMP(l).plus_y); -sind(CCMP(l).plus_y) cosd(CCMP(l).plus_y)];
    eval(['temp=crd_rot*[' sta2{l} '.LAX'';' sta2{l} '.LAY''];']);
    eval([sta2{l} '.LAX=temp(1,:)'';']);
    eval([sta2{l} '.LAY=temp(2,:)'';']);
end

%convert AXCC2 acceleration to tilt
AXCC2.LAX=asin(AXCC2.BNE/9.81)*10^6;
AXCC2.LAY=asin(AXCC2.BNN/9.81)*10^6;
    
%     for i=1:2
%         file_string=['AXCC1_' cha1{i} '_' t1_s '.miniseed'];
%         %attempt download if file not found
%         if ~exist(['../tiltcompare/AXCC1/' file_string],'file')
%             IRIS_data_pull('AXCC1',cha1{i},'11',t1,t1+1);
%         end
%         %some dates have no data (power failure, etc.)
%         if exist(['../tiltcompare/AXCC1/' file_string],'file')
%             temp=rdmseed(['../tiltcompare/AXCC1/' file_string]);
%             dtemp=double(cat(1,temp.d));
%             if length(dtemp)<3600 %short timeseries will throw errors
%                 t1=t1+1;
%                 continue
%             end
%             dtempd=decimate(dtemp,6,'fir');
%             dtempd=decimate(dtempd,10,'fir');
%             eval(['AXCC1.' cha1{i} '=[AXCC1.' cha1{i} ';dtempd];']);
%             if i==1
%                 ttemp=cat(1,temp.t);
%                 ttempd=decimate(ttemp,6,'fir');
%                 ttempd=decimate(ttempd,10,'fir');
%                 AXCC1.time=[AXCC1.time;ttempd];
%             end
%         end
%     end
%     for j=1:4
%         file_string=['AXCC2_' cha2{j} '_' t1_s '.miniseed'];
%         %attempt download if file not found
%         if ~exist(['../tiltcompare/AXCC2/' file_string],'file')
%             IRIS_data_pull('AXCC2',cha2{j},'--',t1,t1+1);
%         end
%         %some dates have no data (power failure, etc.)
%         if exist(['../tiltcompare/AXCC2/' file_string],'file')
%             temp=rdmseed(['../tiltcompare/AXCC2/' file_string]);
%             dtemp=double(cat(1,temp.d));
%             if length(dtemp)<3600 %short timeseries will throw errors
%                 t1=t1+1;
%                 continue
%             end
%             dtempd=decimate(dtemp,12,'fir');
%             dtempd=decimate(dtempd,10,'fir');
%             dtempd=decimate(dtempd,5,'fir');
%             dtempd=decimate(dtempd,4,'fir');
%             eval(['AXCC2.' cha2{j} '=[AXCC2.' cha2{j} ';dtempd];']);
%             if j==1
%                 ttemp=cat(1,temp.t);
%                 ttempd=decimate(ttemp,12,'fir');
%                 ttempd=decimate(ttempd,10,'fir');
%                 ttempd=decimate(ttempd,5,'fir');
%                 ttempd=decimate(ttempd,4,'fir');
%                 AXCC2.time=[AXCC2.time;ttempd];
%             end
%         end
%     end
%     t1=t1+1;
% end
% 
% 
% %%
% dayn=daylist(i);
% datestr(dayn)
% data=[];
% cha={'MNE','MNN','MNZ'};
% for j=1:length(cha)
%     IRIS_data_pull('AXCC2',cha{j},'--',dayn,dayn+1)
%     if ~exist(['../tiltcompare/AXCC2/AXCC2_' cha{j} '_' datestr(dayn,29) '.miniseed'],'file')
%         cha={'BNE','BNN','BNZ'};
%         if ~exist(['../tiltcompare/AXCC2/AXCC2_' cha{j} '_' datestr(dayn,29) '.miniseed'],'file')
%             IRIS_data_pull('AXCC2',cha{j},'--',dayn,dayn+1);
%         end
%     end
%     temp=rdmseed(['../tiltcompare/AXCC2/AXCC2_' cha{j} '_' datestr(dayn,29) '.miniseed']);
%     data.t=cat(1,temp.t);
%     data.a(:,j)=cat(1,temp.d)/10^7;
% end
% temp=sqrt(data.a(:,1).^2+data.a(:,2).^2+data.a(:,3).^2);
% data.a=cat(2,temp,data.a);
% % probably worthwhile to add in something that trims data to just
% % the necessary 10 minute window
% 
% for j=1:5
%     iflips=find(floor(flipInfoSome.t)==dayn);
%     label_strings={'Total Accel (m/s^2)','x (m/s^2)','y (m/s^2)','z (m/s^2)'};
%     figure(j)
%     subplot(411)
%     title([orientation_strings{j} ' Calibration'])
%     for k=1:4
%         subplot(4,1,k)
%         hold on
%         plot(data.t-data.t(1),data.a(:,k),'c','linewidth',1);
%         xlim([flipInfoSome.t(iflips(j))-data.t(1)+15/60/60/24 flipInfoSome.t(iflips(j))-data.t(1)+(flipInfoSome.duration(iflips(j))-5)/60/60/24])
%         ylabel(label_strings{k})
%         set(gca,'xtick',[])
%         box on
%     end
%     datetick('x','keeplimits')
% end
% 
% %%
% %%%%%%%%%%%%CONFIG%%%%%%%%%%%%%%
% cha1='HNZ';
% cha2='BNZ';
% 
% fac1=314252;
% fac2=10^7;
% %%%%%%%%END CONFIG%%%%%%%%%%%%%%
% 
% for i=1:30
%     t0=datenum(2019,4,1)+i-1;
%     t0_s=datestr(t0,31); t0_s=t0_s(1:10);
%     
%     AXCC1=rdmseed(['../tiltcompare/AXCC1/AXCC1_' cha1 '_' t0_s '.miniseed']);
%     t1=cat(1,AXCC1.t);
%     z1=cat(1,AXCC1.d);
%     
%     AXCC2=rdmseed(['../tiltcompare/AXCC2/AXCC2_' cha2 '_' t0_s '.miniseed']);
%     t2=cat(1,AXCC2.t);
%     z2=cat(1,AXCC2.d);
%     
%     figure(77)
%     clf
%     plot(t1,(z1-mean(z1))/fac1,'linewidth',1)
%     hold on
%     plot(t2,(z2-mean(z2))/fac2,'linewidth',1)
%     set(gca,'fontsize',18)
%     ylabel('Acceleration (m/s^2)')
%     datetick('x')
%     title(t0_s)
%     
%     orient landscape
%     fh=gcf;
%     fh.PaperUnits='inches';
%     fh.PaperPosition=[0 0 11 8.5];
%     print(['../tiltcompare/dailycomps/' t0_s],'-dtiff')
% end
% 
% clear; close all;
% 
% %%
% %%%%%%%%%%CONFIG%%%%%%%%%%
% loaddata=true;
% tf=datenum(date);
% 
% sta='AXCC2';
% dec=[8 10 6 10 6];
% %%%%%%%%END CONFIG%%%%%%%%
% 
% % Load pre-exisiting structure, if it exists
% if loaddata && exist('../calibrations/Axial/axialdata_hr.mat','file')
%     load('../calibrations/Axial/axialdata_hr.mat')
%     load('../calibrations/Axial/axialdata_min.mat')
%     t0=ceil(data_min.t(end));
% else
%     t0=datenum(2018,10,13);
% end
% 
% % Determine datenums of calibrations
% load('../calibrations/Axial/axialdata.mat','flipInfoAll')
% [daylist,id,~]=unique(floor(flipInfoAll.t));
% 
% if t0==datenum(2018,10,13)
%     data_min.t=[];data_min.MNE=[];data_min.MNN=[];data_min.MNZ=[];
%     data_min.MKA=[];data_min.iflip=[];
%     data_hr.t=[];data_hr.MNE=[];data_hr.MNN=[];data_hr.MNZ=[];
%     data_hr.MKA=[];data_hr.iflip=[];
% end
% 
% t1=t0;
% while t1<tf
%     t1_s=datestr(t1,31); t1_s=t1_s(1:10);
%     MNN_string=[sta '_MNN_' t1_s '.miniseed'];
%     MNE_string=[sta '_MNE_' t1_s '.miniseed'];
%     MNZ_string=[sta '_MNZ_' t1_s '.miniseed'];
%     MKA_string=[sta '_MKA_' t1_s '.miniseed'];
% 
%     %attempt download if file not found
%     if ~exist(['../tiltcompare/' sta '/' MNN_string],'file') || ...
%             ~exist(['../tiltcompare/' sta '/' MNE_string],'file') || ...
%             ~exist(['../tiltcompare/' sta '/' MNZ_string],'file') || ...
%             ~exist(['../tiltcompare/' sta '/' MKA_string],'file')
%         IRIS_data_pull(sta,'MNN','--',t1,t1+1);
%         IRIS_data_pull(sta,'MNE','--',t1,t1+1);
%         IRIS_data_pull(sta,'MNZ','--',t1,t1+1);
%         IRIS_data_pull(sta,'MKA','--',t1,t1+1);
%     end
%     %some dates have no data (power failure, etc.)
%     if exist(['../tiltcompare/' sta '/' MNN_string],'file') && ...
%             exist(['../tiltcompare/' sta '/' MNE_string],'file') && ...
%             exist(['../tiltcompare/' sta '/' MNZ_string],'file') && ...
%             exist(['../tiltcompare/' sta '/' MKA_string],'file')
%         %MNN channel
%         temp=rdmseed(['../tiltcompare/' sta '/' MNN_string]);
%         ntemp=double(cat(1,temp.d));
%         %decimate to 1 sample/min
%         ntempd=decimate(ntemp,dec(1),'fir');
%         ntempd=decimate(ntempd,dec(2),'fir');
%         ntempd=decimate(ntempd,dec(3),'fir');
%         %note calibration signal
%         [checkcal,iflip]=max(abs(ntempd));
%         if checkcal/10^7>9
%             data_min.iflip=[data_min.iflip;length(data_min.t)+(iflip-10:iflip+10)'];
%         end
%     end
%     
% end
% 
