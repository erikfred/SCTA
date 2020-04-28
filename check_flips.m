
% Simple script to look at all SCTA flips in x, y,
% and z to help identify anomalies that may have lead to bad calibrations.

load('../calibrations/Axial/axialdata','flipInfoAll');

multiplots=true;

[~,id,~]=unique(floor(flipInfoAll.t));
daylist=floor(flipInfoAll.t(id));
% daylist=daylist(end-4:end); % just a few of the most recent
% daylist=daylist(26:end); % eveything since start of 5-position calibration
keyboard

if ~multiplots
    legstr={};
    figure(22), hold on;
    figure(23), hold on;
    figure(24), hold on;
    figure(25), hold on;
    figure(26), hold on;
    figure(5), hold on;
    figure(6), hold on;
    figure(7), hold on;
    for i=1:length(daylist)
        dayn=daylist(i);
        datestr(dayn)
        data=[];
        cha={'MNE','MNN','MNZ','MKA'};
        chastr={'a(:,1)','a(:,2)','a(:,3)','T'};
        for m=1:length(cha)
            IRIS_data_pull('AXCC2',cha{m},'--',dayn,dayn+1)
            if ~exist(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayn,29) '.miniseed'],'file')
                cha={'BNE','BNN','BNZ','BKA'};
                if ~exist(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayn,29) '.miniseed'],'file')
                    IRIS_data_pull('AXCC2',cha{m},'--',dayn,dayn+1);
                end
            end
            temp=rdmseed(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayn,29) '.miniseed']);
            data.t=cat(1,temp.t);
            eval(['data.' chastr{m} '=cat(1,temp.d)/10^7;']);
        end
        data.as=sqrt(data.a(:,1).^2+data.a(:,2).^2+data.a(:,3).^2);
        
        legstr{i}=datestr(dayn,6);
        %time series
        figure(22)
        plot(data.t-floor(data.t(1)-daylist(1)),data.a(:,1),'linewidth',1)
        ylim([-0.15 0.25])
        title('X channel')
        figure(23)
        plot(data.t-floor(data.t(1)-daylist(1)),data.a(:,2),'linewidth',1)
        title('Y channel')
        figure(24)
        plot(data.t-floor(data.t(1)-daylist(1)),data.a(:,3),'linewidth',1)
        title('Z channel')
        figure(25)
        plot(data.t-floor(data.t(1)-daylist(1)),data.as,'linewidth',1)
        title('Total Acceleration')
        figure(26)
        plot(data.t-floor(data.t(1)-daylist(1)),data.T,'linewidth',1)
        title('Temperature')
        
        %ensure even sampling
        tsmooth=linspace(data.t(1),data.t(end),length(data.t))*24*60*60; %[s]
        asmooth=resample(data.as,data.t);
        %     xsmooth=interp1(data.t,data.a(:,1),tsmooth);
        %     ysmooth=interp1(data.t,data.a(:,2),tsmooth);
        %     zsmooth=interp1(data.t,data.a(:,3),tsmooth);
        xsmooth=resample(data.a(:,1),data.t);
        ysmooth=resample(data.a(:,2),data.t);
        zsmooth=resample(data.a(:,3),data.t);
        
        %spectra
        inti=1:5e5;
        intf=6.2e5:length(tsmooth);
        x1(1)=find(abs(asmooth-mean(asmooth(inti)))>1e-3,1)+110;
        x1(2)=x1(1)+800;
        intx1=x1(1):x1(2);
        x2(1)=find(abs(asmooth-mean(asmooth(intf)))>1e-3,1,'last')-95;
        x2(2)=x2(1)-795;
        intx2=x2(2):x2(1);
        
        figure(5)
        subplot(311,'xscale','log','yscale','log')
        hold on
        [f,a]=amplitudespectrum(tsmooth(inti),xsmooth(inti));
        loglog(f,a,'linewidth',1)
        ylabel({'X',''})
        title('Before Flip')
        subplot(312,'xscale','log','yscale','log')
        hold on
        [f,a]=amplitudespectrum(tsmooth(inti),ysmooth(inti));
        loglog(f,a,'linewidth',1)
        ylabel({'Y','Amplitude'})
        subplot(313,'xscale','log','yscale','log')
        hold on
        [f,a]=amplitudespectrum(tsmooth(inti),zsmooth(inti));
        loglog(f,a,'linewidth',1)
        ylabel({'Z',''})
        
        figure(6)
        subplot(311,'xscale','log','yscale','log')
        hold on
        [f,a]=amplitudespectrum(tsmooth(intx1),xsmooth(intx1));
        loglog(f,a,'linewidth',1)
        ylabel({'X',''})
        title('+X Orientation')
        subplot(312,'xscale','log','yscale','log')
        hold on
        [f,a]=amplitudespectrum(tsmooth(intx1),ysmooth(intx1));
        loglog(f,a,'linewidth',1)
        ylabel({'Y','Amplitude'})
        subplot(313,'xscale','log','yscale','log')
        hold on
        [f,a]=amplitudespectrum(tsmooth(intx1),zsmooth(intx1));
        loglog(f,a,'linewidth',1)
        ylabel({'Z',''})
        
        figure(7)
        subplot(311,'xscale','log','yscale','log')
        hold on
        [f,a]=amplitudespectrum(tsmooth(intx2),xsmooth(intx2));
        loglog(f,a,'linewidth',1)
        ylabel({'X',''})
        title('-X Orientation')
        subplot(312,'xscale','log','yscale','log')
        hold on
        [f,a]=amplitudespectrum(tsmooth(intx2),ysmooth(intx2));
        loglog(f,a,'linewidth',1)
        ylabel({'Y',''})
        subplot(313,'xscale','log','yscale','log')
        hold on
        [f,a]=amplitudespectrum(tsmooth(intx2),zsmooth(intx2));
        loglog(f,a,'linewidth',1)
        ylabel({'Z',''})
    end
    
    %% add legends, etc.
    
    figure(22)
    legend(legstr)
    datetick('x','keeplimits')
    ylabel('Acceleration (m/s^2)')
    figure(23)
    legend(legstr)
    datetick('x','keeplimits')
    ylabel('Acceleration (m/s^2)')
    figure(24)
    legend(legstr)
    datetick('x','keeplimits')
    ylabel('Acceleration (m/s^2)')
    figure(25)
    legend(legstr)
    datetick('x','keeplimits')
    ylabel('Acceleration (m/s^2)')
    figure(26)
    legend(legstr)
    datetick('x','keeplimits')
    ylabel([char(176) 'C'])
    
    figure(5)
    subplot(311),legend(legstr)
    subplot(312),legend(legstr)
    subplot(313),legend(legstr)
    figure(6)
    subplot(311),legend(legstr)
    subplot(312),legend(legstr)
    subplot(313),legend(legstr)
    figure(7)
    subplot(311),legend(legstr)
    subplot(312),legend(legstr)
    subplot(313),legend(legstr)
else
    for i=1:length(daylist)
        dayn=daylist(i);
        datestr(dayn)
        data=[];
        cha={'MNE','MNN','MNZ','MKA'};
        chastr={'a(:,1)','a(:,2)','a(:,3)','T'};
        for m=1:length(cha)
            IRIS_data_pull('AXCC2',cha{m},'--',dayn,dayn+1)
            if ~exist(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayn,29) '.miniseed'],'file')
                cha={'BNE','BNN','BNZ','BKA'};
                if ~exist(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayn,29) '.miniseed'],'file')
                    IRIS_data_pull('AXCC2',cha{m},'--',dayn,dayn+1);
                end
            end
            temp=rdmseed(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayn,29) '.miniseed']);
            data.t=cat(1,temp.t);
            eval(['data.' chastr{m} '=cat(1,temp.d)/10^7;']);
        end
        data.as=sqrt(data.a(:,1).^2+data.a(:,2).^2+data.a(:,3).^2);
                    
        %time series
        figure(22)
        hold on
        if floor(data.t(1))==datenum(2019,11,01) % -X bad
            h2=plot3(data.t-floor(data.t(1)-daylist(1)),data.a(:,1),2*ones(size(data.a(:,1))),'linewidth',1);
            set(h2,'color','b')
        elseif floor(data.t(1))==datenum(2019,12,01) || ...
                floor(data.t(1))==datenum(2020,01,08) || ...
                floor(data.t(1))==datenum(2020,01,15) % many bad
            h4=plot3(data.t-floor(data.t(1)-daylist(1)),data.a(:,1),2*ones(size(data.a(:,1))),'linewidth',1);
            set(h4,'color','k')
        elseif floor(data.t(1))==datenum(2020,01,22) || ...
                floor(data.t(1))==datenum(2020,03,18) % -Y bad
            h3=plot3(data.t-floor(data.t(1)-daylist(1)),data.a(:,1),2*ones(size(data.a(:,1))),'linewidth',1);
            set(h3,'color','m')
        else
            h1=plot3(data.t-floor(data.t(1)-daylist(1)),data.a(:,1),ones(size(data.a(:,1))),'linewidth',0.5);
            set(h1,'color',[255 160 122]/255)
        end
        title('X channel')
        
        figure(23)
        hold on
        if floor(data.t(1))==datenum(2019,11,01) % -X bad
            h2=plot3(data.t-floor(data.t(1)-daylist(1)),data.a(:,2),2*ones(size(data.a(:,2))),'linewidth',1);
            set(h2,'color','b')
        elseif floor(data.t(1))==datenum(2019,12,01) || ...
                floor(data.t(1))==datenum(2020,01,08) || ...
                floor(data.t(1))==datenum(2020,01,15) % many bad
            h4=plot3(data.t-floor(data.t(1)-daylist(1)),data.a(:,2),2*ones(size(data.a(:,2))),'linewidth',1);
            set(h4,'color','k')
        elseif floor(data.t(1))==datenum(2020,01,22) || ...
                floor(data.t(1))==datenum(2020,03,18) % -Y bad
            h3=plot3(data.t-floor(data.t(1)-daylist(1)),data.a(:,2),2*ones(size(data.a(:,2))),'linewidth',1);
            set(h3,'color','m')
        else
            h1=plot3(data.t-floor(data.t(1)-daylist(1)),data.a(:,2),ones(size(data.a(:,2))),'linewidth',0.5);
            set(h1,'color',[255 160 122]/255)
        end
        title('Y channel')
        
        figure(24)
        hold on
        if floor(data.t(1))==datenum(2019,11,01) % -X bad
            h2=plot3(data.t-floor(data.t(1)-daylist(1)),data.a(:,3),2*ones(size(data.a(:,3))),'linewidth',1);
            set(h2,'color','b')
        elseif floor(data.t(1))==datenum(2019,12,01) || ...
                floor(data.t(1))==datenum(2020,01,08) || ...
                floor(data.t(1))==datenum(2020,01,15) % many bad
            h4=plot3(data.t-floor(data.t(1)-daylist(1)),data.a(:,3),2*ones(size(data.a(:,3))),'linewidth',1);
            set(h4,'color','k')
        elseif floor(data.t(1))==datenum(2020,01,22) || ...
                floor(data.t(1))==datenum(2020,03,18) % -Y bad
            h3=plot3(data.t-floor(data.t(1)-daylist(1)),data.a(:,3),2*ones(size(data.a(:,3))),'linewidth',1);
            set(h3,'color','m')
        else
            h1=plot3(data.t-floor(data.t(1)-daylist(1)),data.a(:,3),ones(size(data.a(:,3))),'linewidth',0.5);
            set(h1,'color',[255 160 122]/255)
        end
        title('Z channel')
        
        figure(25)
        hold on
        if floor(data.t(1))==datenum(2019,11,01) % -X bad
            h2=plot3(data.t-floor(data.t(1)-daylist(1)),data.as,2*ones(size(data.as)),'linewidth',1);
            set(h2,'color','b')
        elseif floor(data.t(1))==datenum(2019,12,01) || ...
                floor(data.t(1))==datenum(2020,01,08) || ...
                floor(data.t(1))==datenum(2020,01,15) % many bad
            h4=plot3(data.t-floor(data.t(1)-daylist(1)),data.as,2*ones(size(data.as)),'linewidth',1);
            set(h4,'color','k')
        elseif floor(data.t(1))==datenum(2020,01,22) || ...
                floor(data.t(1))==datenum(2020,03,18) % -Y bad
            h3=plot3(data.t-floor(data.t(1)-daylist(1)),data.as,2*ones(size(data.as)),'linewidth',1);
            set(h3,'color','m')
        else
            h1=plot3(data.t-floor(data.t(1)-daylist(1)),data.as,ones(size(data.as)),'linewidth',0.5);
            set(h1,'color',[255 160 122]/255)
        end
        title('Total Acceleration')
        
        figure(26)
        hold on
        if floor(data.t(1))==datenum(2019,11,01) % -X bad
            h2=plot3(data.t-floor(data.t(1)-daylist(1)),data.T,2*ones(size(data.T)),'linewidth',1);
            set(h2,'color','b')
        elseif floor(data.t(1))==datenum(2019,12,01) || ...
                floor(data.t(1))==datenum(2020,01,08) || ...
                floor(data.t(1))==datenum(2020,01,15) % many bad
            h4=plot3(data.t-floor(data.t(1)-daylist(1)),data.T,2*ones(size(data.T)),'linewidth',1);
            set(h4,'color','k')
        elseif floor(data.t(1))==datenum(2020,01,22) || ...
                floor(data.t(1))==datenum(2020,03,18) % -Y bad
            h3=plot3(data.t-floor(data.t(1)-daylist(1)),data.T,2*ones(size(data.T)),'linewidth',1);
            set(h3,'color','m')
        else
            h1=plot3(data.t-floor(data.t(1)-daylist(1)),data.T,ones(size(data.T)),'linewidth',0.5);
            set(h1,'color',[255 160 122]/255)
        end
        title('Temperature')
    end
    
    figure(22)
    legend([h1,h2,h3,h4],'good','bad -X','bad -Y','many bad','location','northeast')
    box on
    
    figure(23)
    legend([h1,h2,h3,h4],'good','bad -X','bad -Y','many bad','location','northeast')
    box on
    
    figure(24)
    legend([h1,h2,h3,h4],'good','bad -X','bad -Y','many bad','location','northeast')
    box on
    
    figure(25)
    legend([h1,h2,h3,h4],'good','bad -X','bad -Y','many bad','location','northeast')
    box on
    
    figure(26)
    legend([h1,h2,h3,h4],'good','bad -X','bad -Y','many bad','location','northeast')
    box on
end
