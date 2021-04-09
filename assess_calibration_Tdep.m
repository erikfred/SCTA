% assess_calibration_Tdep.m
%
% Uses pre-deployment testing data to determine the temperature dependence
% of the x and y channels when in the vertical.
%

Axial=false; % Axial sensor, but test conducted at UW
PinonFlat=false; % PF sensor, but test conducted at UW
OldOcean=true; % simultaneous testing of both sensors

if Axial
    
elseif PinonFlat
    Tref=30;
    load('../temperature_dependence/PinonFlat/pre-deployment.mat')
    
    tfull=dataDec1.t;
    xfull=dataDec1.as;
    yfull=dataDec1.as;
    Tfull=dataDec1.T;
    lt=82500; % approx number of samples in a day
    
    % T dependence, trimming a day off the front with each loop
    for i=1:7
        ttemp=tfull(1053000+(i-1)*lt:1630000);
        xtemp=xfull(1053000+(i-1)*lt:1630000);
        xtemp=inpaint_nans(xtemp);
        Ttemp=Tfull(1053000+(i-1)*lt:1630000);
        Ttemp=inpaint_nans(Ttemp);
        
        m=get_Tdep_params(ttemp,xtemp,Ttemp,Tref);
        Tdep.x1(i)=m(3);
        
        figure(7); clf
        subplot(311)
        plot(ttemp,xtemp-mean(xtemp),'linewidth',1)
        hold on
        plot(ttemp,(m(3)*Ttemp)-mean(m(3)*Ttemp),'linewidth',1)
        datetick('x','keeplimits')
        ylabel('X1 (m/s^2)')
        title(['dx/dT = ' num2str(m(3))])
        legend('data','inv model','location','northeast')
        
        ttemp=tfull(1631000+(i-1)*lt:2230000-200000);
        ytemp=yfull(1631000+(i-1)*lt:2230000-200000);
        ytemp=inpaint_nans(ytemp);
        Ttemp=Tfull(1631000+(i-1)*lt:2230000-200000);
        Ttemp=inpaint_nans(Ttemp);
        
        m=get_Tdep_params(ttemp,ytemp,Ttemp,Tref);
        Tdep.y1(i)=m(3);
        
        subplot(312)
        plot(ttemp,ytemp-mean(ytemp),'linewidth',1)
        hold on
        plot(ttemp,(m(3)*Ttemp)-mean(m(3)*Ttemp),'linewidth',1)
        datetick('x','keeplimits')
        ylabel('Y (m/s^2)')
        title(['dy/dT = ' num2str(m(3))])
        legend('data','inv model','location','northeast')
        
        ttemp=tfull(2250000+(i-1)*lt:2830000-200000);
        xtemp=xfull(2250000+(i-1)*lt:2830000-200000);
        xtemp=inpaint_nans(xtemp);
        Ttemp=Tfull(2250000+(i-1)*lt:2830000-200000);
        Ttemp=inpaint_nans(Ttemp);
        
        m=get_Tdep_params(ttemp,xtemp,Ttemp,Tref);
        Tdep.x2(i)=m(3);
        
        subplot(313)
        plot(ttemp,xtemp-mean(xtemp),'linewidth',1)
        hold on
        plot(ttemp,(m(3)*Ttemp)-mean(m(3)*Ttemp),'linewidth',1)
        datetick('x','keeplimits')
        ylabel('X2 (m/s^2)')
        title(['dx/dT = ' num2str(m(3))])
        legend('data','inv model','location','northeast')
        keyboard
    end
    
elseif OldOcean
    % read in the data for both accelerometers
    t0=datenum(2017,09,19);
    tf=datenum(2017,10,02);
    accel=[]; % older sensor, housed in auto-flipper (for this time interval)
    accel2=[]; % newer sensor, in manual housing (for this time interval)
    for dayn=t0:tf
        data=get_sctaDay('/Volumes/GoogleDrive/My Drive/Oceanography/SCTA-Share/Data/ParsedData',dayn);
        [accel]=decimate_SCTA(data,60,accel);
        
        data=get_sctaDay2('/Volumes/GoogleDrive/My Drive/Oceanography/SCTA-Share/Data/ParsedData',dayn);
        [accel2]=decimate_SCTA(data,60,accel2);
    end
    
    lt=60*24; % number of samples in a day
    Tref=30; % reference temperature for inversion
    
    % T dependence, trimming a day off the front with each loop
    for i=1:6
        % X in vertical for both sensors
        % older sensor
        ttemp=accel.t(1+(i-1)*lt:7*lt-6*60);
        xtemp=accel.a(1+(i-1)*lt:7*lt-6*60,1); %accel.as(1+(i-1)*lt:7*lt-6*60);
        Ttemp=accel.T(1+(i-1)*lt:7*lt-6*60);
        
        m=get_Tdep_params(ttemp,xtemp,Ttemp,Tref);
        accel.Tdep_x(i)=m(3);
        
        figure(17); clf
        subplot(211)
        plot(ttemp,xtemp-mean(xtemp),'linewidth',1)
        hold on
        plot(ttemp,(m(1)+m(2)*(ttemp-ttemp(1))+m(3)*(Ttemp-Tref)),'linewidth',1)
        datetick('x','keeplimits')
        ylabel('X (m/s^2)')
        title({'OLDER SENSOR';['dx/dT = ' num2str(m(3))]})
        legend('data','inv model','location','northeast')
        set(gca,'fontsize',14)
        
        % newer sensor
        ttemp=accel2.t(1+(i-1)*lt:7*lt-6*60);
        xtemp=accel2.a(1+(i-1)*lt:7*lt-6*60,1); %accel2.as(1+(i-1)*lt:7*lt-6*60);
        Ttemp=accel2.T(1+(i-1)*lt:7*lt-6*60);
        
        m=get_Tdep_params(ttemp,xtemp,Ttemp,Tref);
        accel2.Tdep_x(i)=m(3);
        
        figure(27); clf
        subplot(211)
        plot(ttemp,xtemp-mean(xtemp),'linewidth',1)
        hold on
        plot(ttemp,(m(1)+m(2)*(ttemp-ttemp(1))+m(3)*(Ttemp-Tref)),'linewidth',1)
        datetick('x','keeplimits')
        ylabel('X (m/s^2)')
        title({'NEWER SENSOR';['dx/dT = ' num2str(m(3))]})
        legend('data','inv model','location','northeast')
        set(gca,'fontsize',14)
        
        % Y in vertical for both sensors
        % older sensor
        ttemp=accel.t(7*lt+(i-1)*lt:14*lt-6*60);
        ytemp=accel.a(7*lt+(i-1)*lt:14*lt-6*60,2); %accel.as(7*lt+(i-1)*lt:14*lt-6*60);
        Ttemp=accel.T(7*lt+(i-1)*lt:14*lt-6*60);
        
        m=get_Tdep_params(ttemp,ytemp,Ttemp,Tref);
        accel.Tdep_y(i)=m(3);
        
        figure(17)
        subplot(212)
        plot(ttemp,ytemp-mean(ytemp),'linewidth',1)
        hold on
        plot(ttemp,(m(1)+m(2)*(ttemp-ttemp(1))+m(3)*(Ttemp-Tref)),'linewidth',1)
        datetick('x','keeplimits')
        ylabel('Y (m/s^2)')
        title(['dy/dT = ' num2str(m(3))])
        legend('data','inv model','location','northeast')
        set(gca,'fontsize',14)
        
        fh=gcf;
        fh.PaperUnits='inches';
        fh.PaperPosition=[0 0 11 8.5];
        print(['../temperature_dependence/OldOceanTesting/' datestr(t0,29) '/old' num2str(8-i) 'day_Tdep'],'-djpeg','-r100')
        
        % newer sensor
        ttemp=accel2.t(7*lt+(i-1)*lt:14*lt-6*60);
        ytemp=accel2.a(7*lt+(i-1)*lt:14*lt-6*60,2); %accel2.as(7*lt+(i-1)*lt:14*lt-6*60);
        Ttemp=accel2.T(7*lt+(i-1)*lt:14*lt-6*60);
        
        m=get_Tdep_params(ttemp,ytemp,Ttemp,Tref);
        accel2.Tdep_y(i)=m(3);
        
        figure(27)
        subplot(212)
        plot(ttemp,ytemp-mean(ytemp),'linewidth',1)
        hold on
        plot(ttemp,(m(1)+m(2)*(ttemp-ttemp(1))+m(3)*(Ttemp-Tref)),'linewidth',1)
        datetick('x','keeplimits')
        ylabel('Y (m/s^2)')
        title(['dy/dT = ' num2str(m(3))])
        legend('data','inv model','location','northeast')
        set(gca,'fontsize',14)
        
        fh=gcf;
        fh.PaperUnits='inches';
        fh.PaperPosition=[0 0 11 8.5];
        print(['../temperature_dependence/OldOceanTesting/' datestr(t0,29) '/new_' num2str(8-i) 'day_Tdep'],'-djpeg','-r100')
    end
    
    % plot Tdep values to compare intervals
    % older sensor
    figure(18); clf
    subplot(211)
    plot(accel.Tdep_x,'kx','markersize',15,'linewidth',2)
    title('OLDER SENSOR')
    ylabel('dx/dT (m/s^2/C)')
    set(gca,'xtick',1:7)
    set(gca,'xticklabel',0:6)
    set(gca,'fontsize',14)
    subplot(212)
    plot(accel.Tdep_y,'k+','markersize',15,'linewidth',2)
    ylabel('dy/dT (m/s^2/C)')
    xlabel('Days removed from beginning')
    set(gca,'xtick',1:7)
    set(gca,'xticklabel',0:6)
    set(gca,'fontsize',14)
    
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print(['../temperature_dependence/OldOceanTesting/' datestr(t0,29) '/old_Tdep_summary'],'-djpeg','-r100')
    
    % newer sensor
    figure(28); clf
    subplot(211)
    plot(accel2.Tdep_x,'kx','markersize',15,'linewidth',2)
    title('NEWER SENSOR')
    ylabel('dx/dT (m/s^2/C)')
    set(gca,'xtick',1:7)
    set(gca,'xticklabel',0:6)
    set(gca,'fontsize',14)
    subplot(212)
    plot(accel2.Tdep_y,'k+','markersize',15,'linewidth',2)
    ylabel('dy/dT (m/s^2/C)')
    xlabel('Days removed from beginning')
    set(gca,'xtick',1:7)
    set(gca,'xticklabel',0:6)
    set(gca,'fontsize',14)
    
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print(['../temperature_dependence/OldOceanTesting/' datestr(t0,29) '/new_Tdep_summary'],'-djpeg','-r100')
end