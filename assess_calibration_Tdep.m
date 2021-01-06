% assess_calibration_Tdep.m
%
% Uses pre-deployment testing data to determine the temperature dependence
% of the x and y channels when in the vertical.
%

PinonFlat=true;
Axial=false;

if Axial
    
elseif PinonFlat
    Tref=30;
    load('../temperature_dependence/PinonFlat/pre-deployment.mat')
    
    tfull=dataDec1.t;
    xfull=dataDec1.as;
    yfull=dataDec1.as;
    Tfull=dataDec1.T;
    lt=82500; % approx number of samples in a day
    
    % T dependence, trimming a day off the front each loop
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
        legend('data','model','location','northeast')
        
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
        legend('data','model','location','northeast')
        
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
        legend('data','model','location','northeast')
        keyboard
    end
end