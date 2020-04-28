% all_deflation_events.m
%
% Examines all short-term deflation events at Axial Seamount from
% 2016-present. Calculates a number of variables associated with these
% events:
%   - Long-term inflation/tilt rate before event
%   - Inflation/tilt leading up to event
%   - Deflation/tilt of event
%
% These values are used in other scripts to validate forward models of magma
% chamber in-/deflation and caldera fault slip.
%

close all; clear

%%%%%CONFIG%%%%%%
eventnum=5; %1=08/2016 2=02/2017 3=12/2017 4=06/2018 5=05/2019

filtpres=true; %plots filtered pressures

strts=[datenum(2016,01,01),datenum(2016,09,20),datenum(2017,03,01),...
    datenum(2018,01,15),datenum(2018,08,15)];
enders=[datenum(2016,09,20),datenum(2017,03,01),datenum(2018,01,15),...
    datenum(2018,08,15),datenum(2019,09,01)];
fldrs={'aug2016','feb2017','dec2017','jun2018','may2019'};

ashlist(3).fl={'BOTPTA_20170815-20171226.nc';'BOTPTA_20171227-20180114.nc'};
ashlist(4).fl={'BOTPTA_20180115-20180526.nc';'BOTPTA_20180527-20180814.nc'};
ashlist(5).fl={'BOTPTA_20180815-20190113.nc';'BOTPTA_20190114-20190527.nc';'BOTPTA_20190528-20190831.nc'};
dr='/Users/erikfred/Documents/SCTA/tiltcompare/ASHES/tilt/';
%%%%%END CONFIG%%%%%

for ii=eventnum
    t0=strts(ii);
    tf=enders(ii);
    fldr=fldrs{ii};
    
    load('../compass_directions.mat')
    load('../calderacoordinates.mat')
    
    pfac=6894.745/10065*100; %6.89Pa/count / 10065Pa/m * 100cm/m
    
    %load data
    sta={'AXCC1','AXEC2','AXID1'};
    cha={'LAX','LAY'};
    dec=[6 6;5 5;2 2];
    
    % ASHES first
    if tf<=datenum(2017,08,17)
        exclude_ashes=true;
    else
        exclude_ashes=false;
        if ~exist(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/ASHES.mat'],'file')
            ASHES.t_time=[];
            ASHES.LAX=[];
            ASHES.LAY=[];
            ASHES.p_time=[];
            ASHES.BDO=[];
            % tilt
            fl=ashlist(ii).fl;
            for l=1:length(fl)
                ttemp=ncread([dr fl{l}],'time');
                ttemp1=decimate(ttemp,dec(1,1),'fir');
                ttemp2=decimate(ttemp1,dec(2,1),'fir');
                ttemp3=decimate(ttemp2,dec(3,1),'fir');
                ASHES.t_time=[ASHES.t_time; ttemp3/60/60/24+datenum(1900,1,1)];
                xtemp=ncread([dr fl{l}],'lily_x_tilt');
                xtemp1=decimate(xtemp,dec(1,1),'fir');
                xtemp2=decimate(xtemp1,dec(2,1),'fir');
                xtemp3=decimate(xtemp2,dec(3,1),'fir');
                ASHES.LAX=[ASHES.LAX; xtemp3];
                ytemp=ncread([dr fl{l}],'lily_y_tilt');
                ytemp1=decimate(ytemp,dec(1,1),'fir');
                ytemp2=decimate(ytemp1,dec(2,1),'fir');
                ytemp3=decimate(ytemp2,dec(3,1),'fir');
                ASHES.LAY=[ASHES.LAY; ytemp3];
            end
            % pressure
            yrstr1=datestr(t0,'yy');
            yrstr2=datestr(tf,'yy');
            if strcmp(yrstr1,yrstr2)
                load(['../Chadwick_files/data/matfiles/ASHES_' yrstr1])
                eval(['ASHES.p_time=ASHES_' yrstr1 '.time;'])
                eval(['ASHES.BDO=ASHES_' yrstr1 '.p_m*100;']) % [cm]
                eval(['clear ASHES_' yrstr1])
            else
                load(['../Chadwick_files/data/matfiles/ASHES_' yrstr1])
                load(['../Chadwick_files/data/matfiles/ASHES_' yrstr2])
                eval(['ASHES.p_time=ASHES_' yrstr1 '.time;'])
                eval(['ASHES.BDO=ASHES_' yrstr1 '.p_m*100;']) % [cm]
                eval(['ASHES.p_time=cat(1,ASHES.p_time,ASHES_' yrstr2 '.time);'])
                eval(['ASHES.BDO=cat(1,ASHES.BDO,ASHES_' yrstr2 '.p_m*100);']) % [cm]
                eval(['clear ASHES_' yrstr1])
                eval(['clear ASHES_' yrstr2])
            end
            ASHES.lon=-130.01368;ASHES.lat=45.93363;ASHES.depth=1542;
            it=find(ASHES.p_time>=t0 & ASHES.p_time<=tf);
            ASHES.p_time=ASHES.p_time(it);
            ASHES.BDO=ASHES.BDO(it);
            %decimate
            ASHES.p_time=decimate(ASHES.p_time,4,'fir');
            ASHES.BDO=decimate(ASHES.BDO,4,'fir');
        else
            load(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/ASHES.mat'])
        end
    end
    
    % now others
    if ~exist(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/AXEC2.mat'],'file')
        for i=1:length(sta)
            eval([sta{i} '.t_time=[];']);
            eval([sta{i} '.p_time=[];']);
            % tilt
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
                            eval([sta{i} '.t_time=[' sta{i} '.t_time;ttempd];']);
                        end
                    end
                    t1=t1+1;
                end
            end
            % pressure
            yrstr1=datestr(t0,'yy');
            yrstr2=datestr(tf,'yy');
            if strcmp(yrstr1,yrstr2)
                load(['../Chadwick_files/data/matfiles/' sta{i} '_' yrstr1])
                eval([sta{i} '.p_time=' sta{i} '_' yrstr1 '.time;'])
                eval([sta{i} '.BDO=' sta{i} '_' yrstr1 '.p_m*100;']) % [cm]
                eval(['clear ' sta{i} '_' yrstr1])
            else
                load(['../Chadwick_files/data/matfiles/' sta{i} '_' yrstr1])
                load(['../Chadwick_files/data/matfiles/' sta{i} '_' yrstr2])
                eval([sta{i} '.p_time=' sta{i} '_' yrstr1 '.time;'])
                eval([sta{i} '.BDO=' sta{i} '_' yrstr1 '.p_m*100;']) % [cm]
                eval([sta{i} '.p_time=cat(1,' sta{i} '.p_time,' sta{i} '_' yrstr2 '.time);'])
                eval([sta{i} '.BDO=cat(1,' sta{i} '.BDO,' sta{i} '_' yrstr2 '.p_m*100);']) % [cm]
                eval(['clear ' sta{i} '_' yrstr1])
                eval(['clear ' sta{i} '_' yrstr2])
            end
            eval(['it=find(' sta{i} '.p_time>=t0 & ' sta{i} '.p_time<=tf);'])
            eval([sta{i} '.p_time=' sta{i} '.p_time(it);'])
            eval([sta{i} '.BDO=' sta{i} '.BDO(it);'])
            %decimate
            eval([sta{i} '.p_time=decimate(' sta{i} '.p_time,4,''fir'');'])
            eval([sta{i} '.BDO=decimate(' sta{i} '.BDO,4,''fir'');'])
        end
        AXCC1.lat=45.954681;AXCC1.lon=-130.008896;AXCC1.depth=-1528;
        AXEC2.lat=45.939671;AXEC2.lon=-129.973801;AXEC2.depth=-1519;
        AXID1.lat=45.925732;AXID1.lon=-129.977997;AXID1.depth=-1527.5;
    else
        load(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/AXEC2.mat'])
        load(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/AXCC1.mat'])
        load(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/AXID1.mat'])
    end
    
    %save everything before changing
    save(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/AXEC2'],'AXEC2')
    save(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/AXCC1'],'AXCC1')
    save(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/AXID1'],'AXID1')
    if ~exclude_ashes
        save(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/ASHES'],'ASHES')
    end
    
    %truncate data
    if ~exist(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/AXEC2clean.mat'],'file')
        clean_axial_data
    else
        load(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/AXEC2clean.mat'])
        load(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/AXCC1clean.mat'])
        load(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/AXID1clean.mat'])
        if ~exclude_ashes
            load(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/ASHESclean.mat'])
        end
    end
    
    if ~exist(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/AXEC2qc.mat'],'file')
        %select data ranges for various observables
        xt_lt=[]; xt_u=[];xt_d=[]; m=0;
        while isempty(xt_d) || isempty(xt_u) || isempty(xt_lt)
            m=m+1; if m>3; m=1; end
            figure(43)
            subplot(311)
            eval(['plot(' sta{m} '.t_time,' sta{m} '.LAX,''linewidth'',1)'])
            datetick
            subplot(312)
            eval(['plot(' sta{m} '.t_time,' sta{m} '.LAY,''linewidth'',1)'])
            datetick
            lim_x=xlim;
            subplot(313)
            eval(['plot(' sta{m} '.p_time,' sta{m} '.BDO,''linewidth'',1)'])
            datetick
            xlim([lim_x(1) lim_x(2)])
            if isempty(xt_lt)
                disp('Select start and end of long-term uplift')
                [xt_lt,~]=ginput(2);
            end
            if isempty(xt_u)
                disp('Select start and end of uplift')
                [xt_u,~]=ginput(2);
            end
            if isempty(xt_d)
                disp('Select start and end of downdrop')
                [xt_d,~]=ginput(2);
            end
            clf
        end
        close
        
        % tilt and pressure magnitudes
        quantify_intervals
    else
        load(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/AXEC2qc.mat'])
        load(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/AXCC1qc.mat'])
        load(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/AXID1qc.mat'])
        if ~exclude_ashes
            load(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/ASHESqc.mat'])
        end
    end
    
    %% plotting
    
    ASHES.lon=-130.01368;ASHES.lat=45.93363;ASHES.depth=1542;
    
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
    h1=quiver(AXCC1.x,AXCC1.y,AXCC1.x_u_rot*2,AXCC1.y_u_rot*2,'b','linewidth',2,'autoscale','off');
    quiver(AXEC2.x,AXEC2.y,AXEC2.x_u_rot*2,AXEC2.y_u_rot*2,'b','linewidth',2,'autoscale','off')
    quiver(AXID1.x,AXID1.y,AXID1.x_u_rot*2,AXID1.y_u_rot*2,'b','linewidth',2,'autoscale','off')
    if ~exclude_ashes
        quiver(ASHES.x,ASHES.y,ASHES.x_u_rot*2,ASHES.y_u_rot*2,'b','linewidth',2,'autoscale','off')
    end
    h2=quiver(AXCC1.x,AXCC1.y,AXCC1.x_d_rot*2,AXCC1.y_d_rot*2,'r','linewidth',2,'autoscale','off');
    quiver(AXEC2.x,AXEC2.y,AXEC2.x_d_rot*2,AXEC2.y_d_rot*2,'r','linewidth',2,'autoscale','off')
    quiver(AXID1.x,AXID1.y,AXID1.x_d_rot*2,AXID1.y_d_rot*2,'r','linewidth',2,'autoscale','off')
    if ~exclude_ashes
        quiver(ASHES.x,ASHES.y,ASHES.x_d_rot*2,ASHES.y_d_rot*2,'r','linewidth',2,'autoscale','off')
    end
    axis equal
    legend([h1,h2],{'before','after'})
    set(gca,'fontsize',18)
    title(['Tilt change ' fldr],'fontsize',20)
    
    orient tall
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/mapview_tilt'],'-dtiff')
    
    %plot tilts
    figure(23)
    clf
    subplot(411)
    plot(AXCC1.t_time,AXCC1.LAX_rot-nanmean(AXCC1.LAX_rot),'linewidth',1)
    hold on
    plot(AXCC1.t_time,AXCC1.LAY_rot-nanmean(AXCC1.LAY_rot),'linewidth',1)
    datetick('x',6)
    set(gca,'xticklabel',[])
    set(gca,'fontsize',18)
    title('Central Caldera','fontsize',20)
    subplot(412)
    plot(AXEC2.t_time,AXEC2.LAX_rot-nanmean(AXEC2.LAX_rot),'linewidth',1)
    hold on
    plot(AXEC2.t_time,AXEC2.LAY_rot-nanmean(AXEC2.LAY_rot),'linewidth',1)
    datetick('x',6)
    a=xlim;
    set(gca,'xticklabel',[])
    set(gca,'fontsize',18)
    ylabel('tilt (\murad)')
    title('Eastern Caldera','fontsize',20)
    subplot(413)
    if ~exclude_ashes
        plot(ASHES.t_time,ASHES.LAX_rot-nanmean(ASHES.LAX_rot),'linewidth',1)
        hold on
        plot(ASHES.t_time,ASHES.LAY_rot-nanmean(ASHES.LAY_rot),'linewidth',1)
        xlim([a(1) a(2)])
        datetick('x',6,'keeplimits')
        set(gca,'xticklabel',[])
        set(gca,'fontsize',18)
    end
    title('ASHES','fontsize',20)
    subplot(414)
    plot(AXID1.t_time,AXID1.LAX_rot-nanmean(AXID1.LAX_rot),'linewidth',1)
    hold on
    plot(AXID1.t_time,AXID1.LAY_rot-nanmean(AXID1.LAY_rot),'linewidth',1)
    datetick('x',6)
    set(gca,'fontsize',18)
    title('International District','fontsize',20)
    
    %annotations
    subplot(411)
    b=ylim;
    patch([AXCC1.t_time(AXCC1.u_int(1)) AXCC1.t_time(AXCC1.u_int(2)) AXCC1.t_time(AXCC1.u_int(2))...
        AXCC1.t_time(AXCC1.u_int(1))],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
        'linewidth',1)
    patch([AXCC1.t_time(AXCC1.d_int(1)) AXCC1.t_time(AXCC1.d_int(2)) AXCC1.t_time(AXCC1.d_int(2))...
        AXCC1.t_time(AXCC1.d_int(1))],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
        'linewidth',1)
%     text(mean(AXCC1.t_time(AXCC1.u_int)),b(1)+(b(2)-b(1))/4,...
%         [num2str(round(sqrt(AXCC1.x_u_rot^2+AXCC1.y_u_rot^2),2)) ' \murad/day'],'fontsize',12)
%     text(mean(AXCC1.t_time(AXCC1.d_int)),b(1)+(b(2)-b(1))/4,...
%         [num2str(round(sqrt(AXCC1.x_d_rot^2+AXCC1.y_d_rot^2),2)) ' \murad/day'],'fontsize',12)
    legend('x','y','location','southeast')
    subplot(412)
    b=ylim;
    patch([AXEC2.t_time(AXEC2.u_int(1)) AXEC2.t_time(AXEC2.u_int(2)) AXEC2.t_time(AXEC2.u_int(2))...
        AXEC2.t_time(AXEC2.u_int(1))],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
        'linewidth',1)
    patch([AXEC2.t_time(AXEC2.d_int(1)) AXEC2.t_time(AXEC2.d_int(2)) AXEC2.t_time(AXEC2.d_int(2))...
        AXEC2.t_time(AXEC2.d_int(1))],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
        'linewidth',1)
%     text(mean(AXEC2.t_time(AXEC2.u_int)),b(1)+(b(2)-b(1))/4,...
%         [num2str(round(sqrt(AXEC2.x_u_rot^2+AXEC2.y_u_rot^2),2)) ' \murad/day'],'fontsize',12)
%     text(mean(AXEC2.t_time(AXEC2.d_int)),b(1)+(b(2)-b(1))/4,...
%         [num2str(round(sqrt(AXEC2.x_d_rot^2+AXEC2.y_d_rot^2),2)) ' \murad/day'],'fontsize',12)
    legend('x','y','location','southeast')
    subplot(413)
    if ~exclude_ashes
        b=ylim;
        patch([ASHES.t_time(ASHES.u_int(1)) ASHES.t_time(ASHES.u_int(2)) ASHES.t_time(ASHES.u_int(2))...
            ASHES.t_time(ASHES.u_int(1))],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
            'linewidth',1)
        patch([ASHES.t_time(ASHES.d_int(1)) ASHES.t_time(ASHES.d_int(2)) ASHES.t_time(ASHES.d_int(2))...
            ASHES.t_time(ASHES.d_int(1))],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
            'linewidth',1)
%         text(mean(ASHES.t_time(ASHES.u_int)),b(1)+(b(2)-b(1))/4,...
%             [num2str(round(sqrt(ASHES.x_u_rot^2+ASHES.y_u_rot^2),2)) ' \murad/day'],'fontsize',12)
%         text(mean(ASHES.t_time(ASHES.d_int)),b(1)+(b(2)-b(1))/4,...
%             [num2str(round(sqrt(ASHES.x_d_rot^2+ASHES.y_d_rot^2),2)) ' \murad/day'],'fontsize',12)
        legend('x','y','location','southeast')
    end
    subplot(414)
    b=ylim;
    patch([AXID1.t_time(AXID1.u_int(1)) AXID1.t_time(AXID1.u_int(2)) AXID1.t_time(AXID1.u_int(2))...
        AXID1.t_time(AXID1.u_int(1))],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
        'linewidth',1)
    patch([AXID1.t_time(AXID1.d_int(1)) AXID1.t_time(AXID1.d_int(2)) AXID1.t_time(AXID1.d_int(2))...
        AXID1.t_time(AXID1.d_int(1))],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
        'linewidth',1)
%     text(mean(AXID1.t_time(AXID1.u_int)),b(1)+(b(2)-b(1))/4,...
%         [num2str(round(sqrt(AXID1.x_u_rot^2+AXID1.y_u_rot^2),2)) ' \murad/day'],'fontsize',12)
%     text(mean(AXID1.t_time(AXID1.d_int)),b(1)+(b(2)-b(1))/4,...
%         [num2str(round(sqrt(AXID1.x_d_rot^2+AXID1.y_d_rot^2),2)) ' \murad/day'],'fontsize',12)
    legend('x','y','location','southeast')
    
    orient tall
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    keyboard % adjustments
    print(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/tiltplot'],'-dtiff')
    saveas(gcf,['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/tiltplot.fig'])
    
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
    if ~exclude_ashes
        plot(AXEC2.p_time,ASHES.p_dif-nanmean(ASHES.p_dif),'linewidth',1)
        hold on
        xlim([a(1) a(2)])
        datetick('x',6,'keeplimits')
        set(gca,'xticklabel',[])
        set(gca,'fontsize',18)
    end
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
    patch([AXEC2.p_time(AXCC1.pu_int(1)) AXEC2.p_time(AXCC1.pu_int(2)) AXEC2.p_time(AXCC1.pu_int(2))...
        AXEC2.p_time(AXCC1.pu_int(1))],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
        'linewidth',1)
    patch([AXEC2.p_time(AXCC1.pd_int(1)) AXEC2.p_time(AXCC1.pd_int(2)) AXEC2.p_time(AXCC1.pd_int(2))...
        AXEC2.p_time(AXCC1.pd_int(1))],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
        'linewidth',1)
    text(mean(AXEC2.p_time(AXCC1.pu_int)),b(1)+(b(2)-b(1))/4,...
        [num2str(round(AXCC1.p_u,2)) ' cm/day'],'fontsize',12)
    text(mean(AXCC1.p_time(AXCC1.pd_int)),b(1)+(b(2)-b(1))/4,...
        [num2str(round(AXCC1.p_d,2)) ' cm/day'],'fontsize',12)
    subplot(412)
    b=ylim;
    patch([AXEC2.p_time(AXEC2.pu_int(1)) AXEC2.p_time(AXEC2.pu_int(2)) AXEC2.p_time(AXEC2.pu_int(2))...
        AXEC2.p_time(AXEC2.pu_int(1))],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
        'linewidth',1)
    patch([AXEC2.p_time(AXEC2.pd_int(1)) AXEC2.p_time(AXEC2.pd_int(2)) AXEC2.p_time(AXEC2.pd_int(2))...
        AXEC2.p_time(AXEC2.pd_int(1))],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
        'linewidth',1)
    text(mean(AXEC2.p_time(AXEC2.pu_int)),b(1)+(b(2)-b(1))/4,...
        [num2str(round(AXEC2.p_u,2)) ' cm/day'],'fontsize',12)
    text(mean(AXEC2.p_time(AXEC2.pd_int)),b(1)+(b(2)-b(1))/4,...
        [num2str(round(AXEC2.p_d,2)) ' cm/day'],'fontsize',12)
    subplot(413)
    if ~exclude_ashes
        b=ylim;
        patch([AXEC2.p_time(ASHES.pu_int(1)) AXEC2.p_time(ASHES.pu_int(2)) AXEC2.p_time(ASHES.pu_int(2))...
            AXEC2.p_time(ASHES.pu_int(1))],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
            'linewidth',1)
        patch([AXEC2.p_time(ASHES.pd_int(1)) AXEC2.p_time(ASHES.pd_int(2)) AXEC2.p_time(ASHES.pd_int(2))...
            AXEC2.p_time(ASHES.pd_int(1))],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
            'linewidth',1)
        text(mean(AXEC2.p_time(ASHES.pu_int)),b(1)+(b(2)-b(1))/4,...
            [num2str(round(ASHES.p_u,2)) ' cm/day'],'fontsize',12)
        text(mean(AXEC2.p_time(ASHES.pd_int)),b(1)+(b(2)-b(1))/4,...
            [num2str(round(ASHES.p_d,2)) ' cm/day'],'fontsize',12)
    end
    subplot(414)
    b=ylim;
    patch([AXEC2.p_time(AXID1.pu_int(1)) AXEC2.p_time(AXID1.pu_int(2)) AXEC2.p_time(AXID1.pu_int(2))...
        AXEC2.p_time(AXID1.pu_int(1))],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
        'linewidth',1)
    patch([AXEC2.p_time(AXID1.pd_int(1)) AXEC2.p_time(AXID1.pd_int(2)) AXEC2.p_time(AXID1.pd_int(2))...
        AXEC2.p_time(AXID1.pd_int(1))],[b(1) b(1) b(2) b(2)],[-1 -1 -1 -1],[0.9 0.9 0.9],...
        'linewidth',1)
    text(mean(AXEC2.p_time(AXID1.pu_int)),b(1)+(b(2)-b(1))/4,...
        [num2str(round(AXID1.p_u,2)) ' cm/day'],'fontsize',12)
    text(mean(AXID1.p_time(AXID1.pd_int)),b(1)+(b(2)-b(1))/4,...
        [num2str(round(AXID1.p_d,2)) ' cm/day'],'fontsize',12)
    
    orient tall
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    keyboard % adjustments
    print(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/pressureplot'],'-dtiff')
    saveas(gcf,['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/pressureplot.fig'])
    
    %plot pressure changes
    figure(21)
    clf
    plot(caldera_loc(:,1),caldera_loc(:,2),'k')
    hold on
    plot(AXCC1.x,AXCC1.y,'k^','markersize',15,'linewidth',2)
    plot(AXID1.x,AXID1.y,'k^','markersize',15,'linewidth',2)
    plot(AXEC2.x,AXEC2.y,'k^','markersize',15,'linewidth',2)
    plot(ASHES.x,ASHES.y,'k^','markersize',15,'linewidth',2)
    h1=quiver(AXCC1.x,AXCC1.y,0,AXCC1.p_u*3,'b','linewidth',2,'autoscale','off');
    quiver(AXEC2.x,AXEC2.y,0,AXEC2.p_u*3,'b','linewidth',2,'autoscale','off')
    quiver(AXID1.x,AXID1.y,0,AXID1.p_u*3,'b','linewidth',2,'autoscale','off')
    if ~exclude_ashes
        quiver(ASHES.x,ASHES.y,0,ASHES.p_u*3,'b','linewidth',2,'autoscale','off')
    end
    h2=quiver(AXCC1.x,AXCC1.y,0,AXCC1.p_d*3,'r','linewidth',2,'autoscale','off');
    quiver(AXEC2.x,AXEC2.y,0,AXEC2.p_d*3,'r','linewidth',2,'autoscale','off')
    quiver(AXID1.x,AXID1.y,0,AXID1.p_d*3,'r','linewidth',2,'autoscale','off')
    if ~exclude_ashes
        quiver(ASHES.x,ASHES.y,0,ASHES.p_d*3,'r','linewidth',2,'autoscale','off')
    end
    axis equal
    legend([h1,h2],{'before','after'})
    set(gca,'fontsize',18)
    title(['Pressure change ' fldr],'fontsize',20)
    
    orient tall
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/mapview_pressure'],'-dtiff')
end