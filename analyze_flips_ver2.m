function flipInfoOut = analyze_flips_ver2(data,flipInfoIn,p,makePlots)
% Process flips to get calibration values
%
% Usage
%   flipInfoOut(iset) = analyze_flips(dataDec,flipInfoIn,p,makePlots)
%
% Input
%     data       - Data structure (normally decimated) with fields t, a, as, T
%     flipInfoIn - Structure with information about flips with following fields
%                   normalOrientation = The normal orientation of the sensor (usually 3)
%                   i0 - Vector with first sample of each flip
%                   i1 - Vector with last sample of each flip
%                   orientation - Vector with orientation of each flip
%                   t - Vector with serial data number of each flip
%                   aMed - Vector with median acceleration of flipped channel 
%                   range80 - Vector with 80th percentile range of accelerations during the flip
%                   complex - Indicates that the flip was complex (sensor moved while near vertical)
%      p         - Parameter structure with fields
%                   daStartCall - Calibration starts when successive samples change by less than this
%                   tLimCal     - Time limits (s) of calibrations samples to use relative to start time
%      makePlots - Controls what plots are made
%
% Outputs
%      flipInfoOut - As for flipInfoIn but with additional fields
%
% Need to add in temperature calibration

if nargin<4
  makePlots = false;
end

flipInfoOut = flipInfoIn;
for iset = 1:length(flipInfoIn)
  flipInfoOut(iset).i0g = NaN(size(flipInfoOut(iset).i0));
  flipInfoOut(iset).i1g = NaN(size(flipInfoOut(iset).i0));
  flipInfoOut(iset).xCal = NaN(size(flipInfoOut(iset).i0));
  flipInfoOut(iset).xCalTCor = NaN(size(flipInfoOut(iset).i0));
  flipInfoOut(iset).yCal = NaN(size(flipInfoOut(iset).i0));
  flipInfoOut(iset).yCalTCor = NaN(size(flipInfoOut(iset).i0));
  flipInfoOut(iset).zCal = NaN(size(flipInfoOut(iset).i0));
  flipInfoOut(iset).zCalTCor = NaN(size(flipInfoOut(iset).i0));
  flipInfoOut(iset).gCal = NaN(size(flipInfoOut(iset).i0));
  flipInfoOut(iset).gCalTCor = NaN(size(flipInfoOut(iset).i0));
  flipInfoOut(iset).gFrac = NaN(size(flipInfoOut(iset).i0));
  flipInfoOut(iset).length = NaN(size(flipInfoOut(iset).i0));
  flipInfoOut(iset).T = NaN(size(flipInfoOut(iset).i0));
  flipInfoOut(iset).Tstd = NaN(size(flipInfoOut(iset).i0));
  flipInfoOut(iset).duration = NaN(size(flipInfoOut(iset).i0));
end

%Loop through the sets of flips
for iset = 1:length(flipInfoOut)

  %Loop through flips
  for i = 1:length(flipInfoOut(iset).i0)

    if ~flipInfoOut(iset).complex(i)

      % Find the start and end of the calibration based on changes in acceleration
      index = find(abs(diff(data.a(flipInfoOut(iset).i0(i):flipInfoOut(iset).i1(i),abs(flipInfoOut(iset).orientation(i))))) < p.daMax);
      while (index(2) - index(1)) ~=1
        index = index(2:end);
      end
      while (index(end) - index(end-1)) ~=1
        index = index(1:end-1);
      end
      i0g = flipInfoOut(iset).i0(i)+index(1);
      i1g = flipInfoOut(iset).i0(i)+index(end);
%       hold on
%       yyaxis left
%       plot(data.t(i0g:i1g),data.a(i0g:i1g,abs(flipInfoOut(iset).orientation(i))),'g')

      tg = (data.t(i0g:i1g) - data.t(i0g))*86400;
      g = data.a(i0g:i1g,abs(flipInfoOut(iset).orientation(i)));
      Tg = data.T(i0g:i1g);
      x0 = data.a(i0g:i1g,1);
      y0 = data.a(i0g:i1g,2);
      z0 = data.a(i0g:i1g,3);
      g0 = data.as(i0g:i1g);
      x0_cor = x0 - (Tg-p.TRef)*p.dadT(1);
      y0_cor = y0 - (Tg-p.TRef)*p.dadT(2);
      z0_cor = z0; % UNKNOWN COEFFICIENT
      g0_cor = g0 - (Tg-p.TRef)*p.dadT(abs(flipInfoOut(iset).orientation(i)));
      gFrac = mean(g(tg>=p.tCalLim(1) &  tg<=p.tCalLim(2))) / mean(g0(tg>=p.tCalLim(1) &  tg<=p.tCalLim(2)));

      index = tg>=p.tCalLim(1) &  tg<=p.tCalLim(2);
      xCal = mean(x0(index));
      yCal = mean(y0(index));
      zCal = mean(z0(index));
      gCal = mean(g0(index));
      xCalTCor = mean(x0_cor(index));
      yCalTCor = mean(y0_cor(index));
      zCalTCor = mean(z0_cor(index));
      gCalTCor = mean(g0_cor(index));

      flipInfoOut(iset).i0g(i) = i0g;
      flipInfoOut(iset).i1g(i) = i1g;
      flipInfoOut(iset).xCal(i) = xCal;
      flipInfoOut(iset).yCal(i) = yCal;
      flipInfoOut(iset).zCal(i) = zCal;
      flipInfoOut(iset).gCal(i) = gCal;
      flipInfoOut(iset).xCalTCor(i) = xCalTCor;
      flipInfoOut(iset).yCalTCor(i) = yCalTCor;
      flipInfoOut(iset).zCalTCor(i) = zCalTCor;
      flipInfoOut(iset).gCalTCor(i) = gCalTCor;
      flipInfoOut(iset).gFrac(i) = gFrac;
      flipInfoOut(iset).T(i) = mean(Tg);
      flipInfoOut(iset).Tstd(i) = std(Tg);
      flipInfoOut(iset).duration(i) = tg(end)-tg(1)+1;
      
    end
    
    % Plot of acceleration on vertical channel and temperature.
    if makePlots
        orientation_strings={'x1','y','negy','x2','negx'};
        figure(101)
        clf
        
        subplot(341)
        hold on
        yyaxis left
        plot(tg,g0,'b');
        ylabel('g (m/s^2)')
        lim_y=ylim;
        hp=patch([60 90 90 60],[lim_y(1) lim_y(1) lim_y(2) lim_y(2)],[0.8 0.8 0.8]);
        hp.EdgeColor='none';
        hp.FaceVertexAlphaData=0.2; 
        hp.FaceAlpha='flat';
        ylim(lim_y)
        yyaxis right
        plot(tg,Tg,'r');
        title({datestr(floor(data.t(flipInfoOut(iset).i0(i)))), 'g'});
        set(gca,'yticklabel',[])
        box on
        
        subplot(345)
        hold on
        plot(tg(index),g0(index),'b','linewidth',1)
        lim_x=xlim;
        plot([lim_x(1) lim_x(2)],[gCal gCal]','k--','linewidth',1)
        ylabel('g (m/s^2)')
        yyaxis right
        plot(tg(index),(Tg(index)-p.TRef)*p.dadT(abs(flipInfoOut(iset).orientation(i))))
        title('Uncorrected')
        box on
        
        subplot(349)
        hold on
        plot(tg(index),g0_cor(index),'b','linewidth',1)
        lim_x=xlim;
        plot([lim_x(1) lim_x(2)],[gCalTCor gCalTCor]','k--','linewidth',1)
        xlabel('Time, s')
        ylabel('g (m/s^2)')
        title('Corrected')
        box on
        
        subplot(342)
        hold on
        yyaxis left
        plot(tg,x0,'b');
        lim_y=ylim;
        hp=patch([60 90 90 60],[lim_y(1) lim_y(1) lim_y(2) lim_y(2)],[0.8 0.8 0.8]);
        hp.EdgeColor='none';
        hp.FaceVertexAlphaData=0.2; 
        hp.FaceAlpha='flat';
        ylim(lim_y)
        yyaxis right
        plot(tg,Tg,'r');
        title('x');
        set(gca,'yticklabel',[])
        box on
        
        subplot(346)
        hold on
        plot(tg(index),x0(index),'b','linewidth',1)
        lim_x=xlim;
        plot([lim_x(1) lim_x(2)],[xCal xCal]','k--','linewidth',1)
        yyaxis right
        plot(tg(index),(Tg(index)-p.TRef)*p.dadT(1))
        title('Uncorrected')
        box on
        
        subplot(3,4,10)
        hold on
        plot(tg(index),x0_cor(index),'b','linewidth',1)
        lim_x=xlim;
        plot([lim_x(1) lim_x(2)],[xCalTCor xCalTCor]','k--','linewidth',1)
        xlabel('Time, s')
        title('Corrected')
        box on
        
        subplot(343)
        hold on
        yyaxis left
        plot(tg,y0,'b');
        lim_y=ylim;
        hp=patch([60 90 90 60],[lim_y(1) lim_y(1) lim_y(2) lim_y(2)],[0.8 0.8 0.8]);
        hp.EdgeColor='none';
        hp.FaceVertexAlphaData=0.2; 
        hp.FaceAlpha='flat';
        ylim(lim_y)
        yyaxis right
        plot(tg,Tg,'r');
        title('y');
        set(gca,'yticklabel',[])
        box on
        
        subplot(347)
        hold on
        plot(tg(index),y0(index),'b','linewidth',1)
        lim_x=xlim;
        plot([lim_x(1) lim_x(2)],[yCal yCal]','k--','linewidth',1)
        yyaxis right
        plot(tg(index),(Tg(index)-p.TRef)*p.dadT(2))
        title('Uncorrected')
        box on
        
        subplot(3,4,11)
        hold on
        plot(tg(index),y0_cor(index),'b','linewidth',1)
        lim_x=xlim;
        plot([lim_x(1) lim_x(2)],[yCalTCor yCalTCor]','k--','linewidth',1)
        xlabel('Time, s')
        title('Corrected')
        box on
        
        subplot(344)
        hold on
        yyaxis left
        plot(tg,z0,'b');
        lim_y=ylim;
        hp=patch([60 90 90 60],[lim_y(1) lim_y(1) lim_y(2) lim_y(2)],[0.8 0.8 0.8]);
        hp.EdgeColor='none';
        hp.FaceVertexAlphaData=0.2; 
        hp.FaceAlpha='flat';
        ylim(lim_y)
        yyaxis right
        plot(tg,Tg,'r');
        title('z');
        ylabel(['T (' char(176) 'C)'])
        box on
        
        subplot(348)
        hold on
        plot(tg(index),z0(index),'b','linewidth',1)
        lim_x=xlim;
        plot([lim_x(1) lim_x(2)],[zCal zCal]','k--','linewidth',1)
        yyaxis right
        plot(tg(index),(Tg(index)-p.TRef)*p.dadT(abs(flipInfoOut(iset).orientation(i)))) % needs to stay in for y labels
        ylabel('T correction term')
        title('Uncorrected')
        box on
        
        subplot(3,4,12)
        hold on
        plot(tg(index),z0_cor(index),'b','linewidth',1)
        lim_x=xlim;
        plot([lim_x(1) lim_x(2)],[zCalTCor zCalTCor]','k--','linewidth',1)
        xlabel('Time, s')
        title('Corrected')
        box on
        
        fh=gcf;
        fh.PaperUnits='inches';
        fh.PaperPosition=[0 0 11 8.5];
        print(['../calibrations/Axial/cal_plots/' datestr(floor(data.t(flipInfoOut(iset).i0(i))),29) ...
            '_' orientation_strings{i}],'-dtiff','-r300')
    end
    
  end
end
