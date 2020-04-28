% check_eruption_tilts.m
%
% Calculates tilt directions and amplitudes 1) pre-eruption, 2) during dike
% intrusion, 3) during deflation, and 4) during re-inflation against the
% published results of Nooner and Chadwick.
%

close all; clear

%%%%%CONFIG%%%%%%
eventnum=1:4;

strts=[datenum(2014,9,18),datenum(2015,4,24,5,0,0),datenum(2015,4,24,9,0,0),datenum(2015,5,20)];
enders=[datenum(2015,4,24),datenum(2015,4,24,9,0,0),datenum(2015,5,19),datenum(2016,5,20)];
fldrs={'1-pre_eruption','2-dike_intrusion','3-deflation','4-re_inflation'};

cal_a=[87217 0 0 0];
cal_b=[87229 0 0 0];
%%%%%END CONFIG%%%%%

load('../calderacoordinates.mat');
load('../compass_directions.mat');

for m=1:length(eventnum)
    t0=strts(eventnum(m));
    tf=enders(eventnum(m));
    fldr=fldrs{eventnum(m)};
    sta={'AXCC1','AXEC2','AXID1'};
    cha={'LAX','LAY'};
    dec=[6 6 12;5 5 10;2 2 10];
    %load data
    if ~exist(['../tiltcompare/2015_eruption/' fldr '/AXEC2.mat'],'file')
        for i=1:length(sta)
            eval([sta{i} '.time=[];']);
            for j=1:length(cha)
                t1=floor(t0);
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
                        end
                    end
                    t1=t1+1;
                end
            end
            %trim to fit time range, where appropriate
            if rem(t0,1)~=0
                eval([sta{i} '.LAX(' sta{i} '.time<t0 | ' sta{i} '.time>tf)=[];']);
                eval([sta{i} '.LAY(' sta{i} '.time<t0 | ' sta{i} '.time>tf)=[];']);
                eval([sta{i} '.time(' sta{i} '.time<t0 | ' sta{i} '.time>tf)=[];']);
            end
        end
        AXCC1.lat=45.954681;AXCC1.lon=-130.008896;AXCC1.depth=-1528;
        AXEC2.lat=45.939671;AXEC2.lon=-129.973801;AXEC2.depth=-1519;
        AXID1.lat=45.925732;AXID1.lon=-129.977997;AXID1.depth=-1527.5;
        
        save(['../tiltcompare/2015_eruption/' fldr '/AXEC2'],'AXEC2')
        save(['../tiltcompare/2015_eruption/' fldr '/AXCC1'],'AXCC1')
        save(['../tiltcompare/2015_eruption/' fldr '/AXID1'],'AXID1')
    else
        load(['../tiltcompare/2015_eruption/' fldr '/AXEC2.mat'])
        load(['../tiltcompare/2015_eruption/' fldr '/AXCC1.mat'])
        load(['../tiltcompare/2015_eruption/' fldr '/AXID1.mat'])
    end
    
    %tilt magnitudes, with compass rotation
    AXCC1.x_mag=mean(AXCC1.LAX(end-9:end))-mean(AXCC1.LAX(1:10));
    AXCC1.y_mag=mean(AXCC1.LAY(end-9:end))-mean(AXCC1.LAY(1:10));
    crd_rot=[cosd(CCMP(4).plus_y) sind(CCMP(4).plus_y); -sind(CCMP(4).plus_y) cosd(CCMP(4).plus_y)];
    temp3=crd_rot*[AXCC1.x_mag;AXCC1.y_mag]; AXCC1.x_rot=temp3(1); AXCC1.y_rot=temp3(2);
    
    AXEC2.x_mag=mean(AXEC2.LAX(end-9:end))-mean(AXEC2.LAX(1:10));
    AXEC2.y_mag=mean(AXEC2.LAY(end-9:end))-mean(AXEC2.LAY(1:10));
    crd_rot=[cosd(CCMP(3).plus_y) sind(CCMP(3).plus_y); -sind(CCMP(3).plus_y) cosd(CCMP(3).plus_y)];
    temp3=crd_rot*[AXEC2.x_mag;AXEC2.y_mag]; AXEC2.x_rot=temp3(1); AXEC2.y_rot=temp3(2);
    
    if m~=2
        AXID1.x_mag=mean(AXID1.LAX(end-9:end))-mean(AXID1.LAX(1:10));
        AXID1.y_mag=mean(AXID1.LAY(end-9:end))-mean(AXID1.LAY(1:10));
        crd_rot=[cosd(CCMP(2).plus_y) sind(CCMP(2).plus_y); -sind(CCMP(2).plus_y) cosd(CCMP(2).plus_y)];
        temp3=crd_rot*[AXID1.x_mag;AXID1.y_mag]; AXID1.x_rot=temp3(1); AXID1.y_rot=temp3(2);
    end
    
    %prep for plots
    caldera_loc=llh2local(caldera_rim',[AXCC1.lon,AXCC1.lat])';
    crds=llh2local([AXCC1.lon,AXCC1.lat]',[AXCC1.lon,AXCC1.lat]);
    AXCC1.x=crds(1); AXCC1.y=crds(2);
    crds=llh2local([AXEC2.lon,AXEC2.lat]',[AXCC1.lon,AXCC1.lat]);
    AXEC2.x=crds(1); AXEC2.y=crds(2);
    crds=llh2local([AXID1.lon,AXID1.lat]',[AXCC1.lon,AXCC1.lat]);
    AXID1.x=crds(1); AXID1.y=crds(2);
    
    %plot tilt directions
    figure(42)
    clf
    plot(caldera_loc(:,1),caldera_loc(:,2),'k')
    hold on
    plot(AXCC1.x,AXCC1.y,'k^','markersize',15,'linewidth',2)
    plot(AXID1.x,AXID1.y,'k^','markersize',15,'linewidth',2)
    plot(AXEC2.x,AXEC2.y,'k^','markersize',15,'linewidth',2)
    quiver(AXCC1.x,AXCC1.y,AXCC1.x_rot/400,AXCC1.y_rot/400,'b','linewidth',2)
    quiver(AXEC2.x,AXEC2.y,AXEC2.x_rot/400,AXEC2.y_rot/400,'b','linewidth',2)
    if m~=2
        quiver(AXID1.x,AXID1.y,AXID1.x_rot/400,AXID1.y_rot/400,'b','linewidth',2)
    end
    quiver(-4,-3,1000/400,0,'r','linewidth',2)
    axis equal
    set(gca,'fontsize',15)
    undr=strfind(fldr,'_');
    title([fldr(3:undr-1) ' ' fldr(undr+1:end)],'fontsize',20)
    
    orient portrait
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print(['../tiltcompare/2015_eruption/' fldr '/mapview'],'-dtiff')
    
    %plotting tilts
    figure(23)
    clf
    subplot(311)
    plot(AXCC1.time,AXCC1.LAX-nanmean(AXCC1.LAX),'linewidth',1)
    hold on
    plot(AXCC1.time,AXCC1.LAY-nanmean(AXCC1.LAY),'linewidth',1)
    legend('x','y','location','southeast')
    datetick('x','keeplimits')
    set(gca,'xticklabel',[])
    set(gca,'fontsize',18)
    title('Central Caldera','fontsize',20)
    subplot(312)
    plot(AXEC2.time,AXEC2.LAX-nanmean(AXEC2.LAX),'linewidth',1)
    hold on
    plot(AXEC2.time,AXEC2.LAY-nanmean(AXEC2.LAY),'linewidth',1)
    legend('x','y','location','southeast')
    datetick('x','keeplimits')
    set(gca,'xticklabel',[])
    set(gca,'fontsize',18)
    ylabel('tilt (\murad)')
    title('Eastern Caldera','fontsize',20)
    subplot(313)
    plot(AXID1.time,AXID1.LAX-nanmean(AXID1.LAX),'linewidth',1)
    hold on
    plot(AXID1.time,AXID1.LAY-nanmean(AXID1.LAY),'linewidth',1)
    legend('x','y','location','southeast')
    datetick('x','keeplimits')
    set(gca,'fontsize',18)
    title('International District','fontsize',20)
    
    orient tall
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print(['../tiltcompare/2015_eruption/' fldr '/tiltplot'],'-dtiff')
    
end