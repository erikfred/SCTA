% monthly_tiltplots.m
%
% Generates monthly plots of tilt, with calibrations marked, between
% specified dates. Uses 8Hz filtered data to avoid OOI gaps.
%

clear; close all;

%%%%%%%%%%CONFIG%%%%%%%%%%
t0=datenum(2020,02,1);
tf=datenum(date);

sta='AXCC2';
dec=[6 10 8];
%%%%%%%%END CONFIG%%%%%%%%

t1=t0;
monum=str2double(datestr(t0,5));
while t1<tf
    t1_s=datestr(t1,31); t1_s=t1_s(1:10);
    MNN_string=[sta '_MNN_' t1_s '.miniseed'];
    MNE_string=[sta '_MNE_' t1_s '.miniseed'];
    MNZ_string=[sta '_MNZ_' t1_s '.miniseed'];
    MKA_string=[sta '_MKA_' t1_s '.miniseed'];
    %establish new figure for each month
    if t1==t0
        AXCC2.time=[];
        AXCC2.MNE=[];
        AXCC2.MNN=[];
        AXCC2.MNZ=[];
        AXCC2.MKA=[];
        AXCC2.iflip=[];
    elseif monum~=str2double(t1_s(end-4:end-3))
        %make previous month's figure
        figure(101)
        clf
        plot(AXCC2.time,AXCC2.MNE-nanmean(AXCC2.MNE),'linewidth',1)
        hold on
        plot(AXCC2.time,AXCC2.MNN-nanmean(AXCC2.MNN),'linewidth',1)
        legend('East','North')
        datetick('x',6)
        set(gca,'fontsize',18)
        ylabel('Accel (m/s^2)')
        title([num2str(monum) '/' datestr(t1,10)])
        fh=gcf;
        fh.PaperUnits='inches';
        fh.PaperPosition=[0 0 11 8.5];

        figure(102)
        clf
        AXCC2.LAX=asin(AXCC2.MNE/9.81)*10^6;
        AXCC2.LAY=asin(AXCC2.MNN/9.81)*10^6;
        plot(AXCC2.time,AXCC2.LAX-nanmean(AXCC2.LAX),'linewidth',1)
        hold on
        plot(AXCC2.time,AXCC2.LAY-nanmean(AXCC2.LAY),'linewidth',1)
        legend('East','North')
        datetick('x',6)
        set(gca,'fontsize',18)
        ylabel('Tilt (\murad)')
        title([num2str(monum) '/' datestr(t1,10)])
        fh=gcf;
        fh.PaperUnits='inches';
        fh.PaperPosition=[0 0 11 8.5];
        
        if ~exist(['../monthly_plots/' datestr(t1-1,'mmm')],'dir')
            mkdir(['../monthly_plots/' datestr(t1-1,'mmm')])
        end
        %figure(101); print(['../monthly_plots/' datestr(t1-1,'mmm') '/accel'],'-dtiff')
        %figure(102); print(['../monthly_plots/' datestr(t1-1,'mmm') '/tilt'],'-dtiff')
        save(['../monthly_plots/' datestr(t1-1,'mmm') '/AXCC2'],'AXCC2')
        
        AXCC2.time=[];
        AXCC2.MNE=[];
        AXCC2.MNN=[];
        AXCC2.MNZ=[];
        AXCC2.MKA=[];
        AXCC2.iflip=[];
        
        monum=monum+1; if monum==13; monum=1; end
    end
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
        %if 1st of month, note calibration signal
        [checkcal,iflip]=max(abs(ntempd));
        if checkcal/10^7>9
            AXCC2.iflip=[AXCC2.iflip;(iflip-30:iflip+29)];
        end
        %append decimated and cleaned record to structure
        AXCC2.MNN=[AXCC2.MNN; ntempd/10^7];
        
        %MNE channel
        temp=rdmseed(['../tiltcompare/' sta '/' MNE_string]);
        etemp=double(cat(1,temp.d));
        %decimate to 1 sample/min
        etempd=decimate(etemp,dec(1),'fir');
        etempd=decimate(etempd,dec(2),'fir');
        etempd=decimate(etempd,dec(3),'fir');
        %append decimated record to structure
        AXCC2.MNE=[AXCC2.MNE; etempd/10^7];
        
        %MNZ channel
        temp=rdmseed(['../tiltcompare/' sta '/' MNZ_string]);
        ztemp=double(cat(1,temp.d));
        %decimate to 1 sample/min
        ztempd=decimate(ztemp,dec(1),'fir');
        ztempd=decimate(ztempd,dec(2),'fir');
        ztempd=decimate(ztempd,dec(3),'fir');
        %append decimated record to structure
        AXCC2.MNZ=[AXCC2.MNZ; ztempd/10^7];
        
        %MKA channel
        temp=rdmseed(['../tiltcompare/' sta '/' MKA_string]);
        Ttemp=double(cat(1,temp.d));
        %decimate to 1 sample/min
        Ttempd=decimate(Ttemp,dec(1),'fir');
        Ttempd=decimate(Ttempd,dec(2),'fir');
        Ttempd=decimate(Ttempd,dec(3),'fir');
        %append decimated record to structure
        AXCC2.MKA=[AXCC2.MKA; Ttempd/10^7];
        
        %time
        ttemp=cat(1,temp.t);
        ttempd=decimate(ttemp,dec(1),'fir');
        ttempd=decimate(ttempd,dec(2),'fir');
        ttempd=decimate(ttempd,dec(3),'fir');
        AXCC2.time=[AXCC2.time; ttempd];
    end
    t1=t1+1;
end