% make_dailyspect.m
%
% Makes spectrograms from daily files at 40Hz. Options at beginning to
% specify range of dates, choose whether or not to make spectrograms from
% partial days. Warns if you try to make a spectrogram from a day with a
% calibration.
%

clear; close all

%%%%% CONFIG %%%%%

t0=datenum(2018,09,13); % beginning of interval
tf=datenum(2018,09,23); % end of interval (same as above if only 1 day)

partialday=false; % whether or not to make spectrograms from partial days

%%% END CONFIG %%%

load ../calibrations/Axial/detailed_flipInfo.mat

for dayn=t0:tf
    % get data
    daystr=datestr(dayn,29);
    flipday=false;
    partday=false;
    
    if sum(dayn==floor(flipInfoSome.t))
        warning([datestr(dayn) ' is a calibration day'])
        flipday=true;
    end
    data = get_sctaDay('/Volumes/GoogleDrive/My Drive/Oceanography/SCTA-Share/OOI-SCTA/ParsedData',dayn);
    
    if isempty(data)
        cha={'BNE','BNN','BNZ','BKA'};
        chastr={'a(:,1)','a(:,2)','a(:,3)','T'};
        for m=1:length(cha)
            IRIS_data_pull('AXCC2',cha{m},'--',dayn,dayn+1)
            temp=rdmseed(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayn,29) '.miniseed']);
            data.t=cat(1,temp.t);
            eval(['data40_pre.' chastr{m} '=cat(1,temp.d)/10^7;']);
        end
        data.as=sqrt(data.a(:,1).^2+data.a(:,2).^2+data.a(:,3).^2);
    end
    
    % check length of time series
    if isempty(data)
        warning([datestr(dayn) ' data unavailable'])
        continue
    elseif length(data.t)<576000
        warning([datestr(dayn) ' has < 4 hours of data available'])
        continue
    elseif length(data.t)>576000 && length(data.t)<3456000
        numhrs=floor(length(data.t)/40/60/60);
        warning([datestr(dayn) ' producing plot from available ' num2str(numhrs) ' hours of data'])
        partday=true;
    end
    
    % make spectrograms
    
    xdata=data.a(:,1)-mean(data.a(:,1));
    ydata=data.a(:,2)-mean(data.a(:,2));
    zdata=data.a(:,3)-mean(data.a(:,3));
    Tdata=data.T-mean(data.T);
    
    % X
    figure(31)
    clf
    spectrogram(xdata,400*3,200*3,2^10,40,'yaxis')
    title(['X channel ' datestr(dayn)])
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    if partday && flipday
        title(['X channel ' datestr(dayn) ' (calibration day)'])
        print(['../noise exploration/daily_spectrograms/' daystr '_' num2str(numhrs) 'hrs_calday_X'],'-dtiff','-r300')
        print(['../noise exploration/daily_spectrograms/' daystr '_' num2str(numhrs) 'hrs_calday_X'],'-djpeg')
    elseif partday
        title(['X channel ' datestr(dayn)])
        print(['../noise exploration/daily_spectrograms/' daystr '_' num2str(numhrs) 'hrs_X'],'-dtiff','-r300')
        print(['../noise exploration/daily_spectrograms/' daystr '_' num2str(numhrs) 'hrs_X'],'-djpeg')
    elseif flipday
        title(['X channel ' datestr(dayn) ' (calibration day)'])
        print(['../noise exploration/daily_spectrograms/' daystr '_calday_X'],'-dtiff','-r300')
        print(['../noise exploration/daily_spectrograms/' daystr '_calday_X'],'-djpeg')
    else
    print(['../noise exploration/daily_spectrograms/' daystr '_X'],'-dtiff','-r300')
    print(['../noise exploration/daily_spectrograms/' daystr '_X'],'-djpeg')
    end
    
    %Y
    figure(32)
    clf
    spectrogram(ydata,400*3,200*3,2^10,40,'yaxis')
    title(['Y channel ' datestr(dayn)])
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    if partday && flipday
        title(['Y channel ' datestr(dayn) ' (calibration day)'])
        print(['../noise exploration/daily_spectrograms/' daystr '_' num2str(numhrs) 'hrs_calday_Y'],'-dtiff','-r300')
        print(['../noise exploration/daily_spectrograms/' daystr '_' num2str(numhrs) 'hrs_calday_Y'],'-djpeg')
    elseif partday
        title(['Y channel ' datestr(dayn)])
        print(['../noise exploration/daily_spectrograms/' daystr '_' num2str(numhrs) 'hrs_Y'],'-dtiff','-r300')
        print(['../noise exploration/daily_spectrograms/' daystr '_' num2str(numhrs) 'hrs_Y'],'-djpeg')
    elseif flipday
        title(['Y channel ' datestr(dayn) ' (calibration day)'])
        print(['../noise exploration/daily_spectrograms/' daystr '_calday_Y'],'-dtiff','-r300')
        print(['../noise exploration/daily_spectrograms/' daystr '_calday_Y'],'-djpeg')
    else
    print(['../noise exploration/daily_spectrograms/' daystr '_Y'],'-dtiff','-r300')
    print(['../noise exploration/daily_spectrograms/' daystr '_Y'],'-djpeg')
    end
    
    figure(33)
    clf
    spectrogram(zdata,400*3,200*3,2^10,40,'yaxis')
    title(['Z channel ' datestr(dayn)])
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    if partday && flipday
        title(['Z channel ' datestr(dayn) ' (calibration day)'])
        print(['../noise exploration/daily_spectrograms/' daystr '_' num2str(numhrs) 'hrs_calday_Z'],'-dtiff','-r300')
        print(['../noise exploration/daily_spectrograms/' daystr '_' num2str(numhrs) 'hrs_calday_Z'],'-djpeg')
    elseif partday
        title(['Z channel ' datestr(dayn)])
        print(['../noise exploration/daily_spectrograms/' daystr '_' num2str(numhrs) 'hrs_Z'],'-dtiff','-r300')
        print(['../noise exploration/daily_spectrograms/' daystr '_' num2str(numhrs) 'hrs_Z'],'-djpeg')
    elseif flipday
        title(['Z channel ' datestr(dayn) ' (calibration day)'])
        print(['../noise exploration/daily_spectrograms/' daystr '_calday_Z'],'-dtiff','-r300')
        print(['../noise exploration/daily_spectrograms/' daystr '_calday_Z'],'-djpeg')
    else
    print(['../noise exploration/daily_spectrograms/' daystr '_Z'],'-dtiff','-r300')
    print(['../noise exploration/daily_spectrograms/' daystr '_Z'],'-djpeg')
    end
    
    figure(34)
    clf
    spectrogram(Tdata,400*3,200*3,2^10,40,'yaxis')
    title(['Temperature ' datestr(dayn)])
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    if partday && flipday
        title(['Temperature ' datestr(dayn) ' (calibration day)'])
        print(['../noise exploration/daily_spectrograms/' daystr '_' num2str(numhrs) 'hrs_calday_T'],'-dtiff','-r300')
        print(['../noise exploration/daily_spectrograms/' daystr '_' num2str(numhrs) 'hrs_calday_T'],'-djpeg')
    elseif partday
        title(['Temperature ' datestr(dayn)])
        print(['../noise exploration/daily_spectrograms/' daystr '_' num2str(numhrs) 'hrs_T'],'-dtiff','-r300')
        print(['../noise exploration/daily_spectrograms/' daystr '_' num2str(numhrs) 'hrs_T'],'-djpeg')
    elseif flipday
        title(['Temperature ' datestr(dayn) ' (calibration day)'])
        print(['../noise exploration/daily_spectrograms/' daystr '_calday_T'],'-dtiff','-r300')
        print(['../noise exploration/daily_spectrograms/' daystr '_calday_T'],'-djpeg')
    else
    print(['../noise exploration/daily_spectrograms/' daystr '_T'],'-dtiff','-r300')
    print(['../noise exploration/daily_spectrograms/' daystr '_T'],'-djpeg')
    end
        
end