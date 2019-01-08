function flipInfoOut = analyze_flips(data,flipInfoIn,p,makePlots)
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
%      flipInfoOut - As for flipInfoIn but with additional fiels
%
% Need to add in temperature calibration

if nargin<4
  makePlots = false;
end

flipInfoOut = flipInfoIn;
for iset = 1:length(flipInfoIn)
  flipInfoOut(iset).i0g = NaN(size(flipInfoOut(iset).i0));
  flipInfoOut(iset).i1g = NaN(size(flipInfoOut(iset).i0));
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

    % Plot of acceleration on vertical channel and temperature.
    if makePlots
      figure(101)
      clf
      subplot(311)
      yyaxis left
      plot(data.t(flipInfoOut(iset).i0(i):flipInfoOut(iset).i1(i)),data.a(flipInfoOut(iset).i0(i):flipInfoOut(iset).i1(i),abs(flipInfoOut(iset).orientation(i))),'b');
      yyaxis right
      plot(data.t(flipInfoOut(iset).i0(i):flipInfoOut(iset).i1(i)),data.T(flipInfoOut(iset).i0(i):flipInfoOut(iset).i1(i)),'r');
      datetick
      title({['Flip ' int2str(i) ';   Orientation ' int2str(flipInfoOut(iset).orientation(i)) ';   ' datestr(floor(data.t(flipInfoOut(iset).i0(i)))) ] , 'Complex Flip'});
      subplot(312)
      yyaxis left
      plot(data.t(flipInfoOut(iset).i0(i):flipInfoOut(iset).i1(i)),data.a(flipInfoOut(iset).i0(i):flipInfoOut(iset).i1(i),abs(flipInfoOut(iset).orientation(i))),'b');
      ylim([flipInfoOut(iset).aMed(i)-0.0001 flipInfoOut(iset).aMed(i)+0.0001])
      yyaxis right
      plot(data.t(flipInfoOut(iset).i0(i):flipInfoOut(iset).i1(i)),data.T(flipInfoOut(iset).i0(i):flipInfoOut(iset).i1(i)),'r');
      datetick
    end

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
      hold on
      yyaxis left
      plot(data.t(i0g:i1g),data.a(i0g:i1g,abs(flipInfoOut(iset).orientation(i))),'g')

      tg = (data.t(i0g:i1g) - data.t(i0g))*86400;
      g = data.a(i0g:i1g,abs(flipInfoOut(iset).orientation(i)));
      Tg = data.T(i0g:i1g);
      g0 = data.as(i0g:i1g);
      gFrac = mean(g(tg>=p.tCalLim(1) &  tg<=p.tCalLim(2))) / mean(g0(tg>=p.tCalLim(1) &  tg<=p.tCalLim(2)));

  %     %Temperature correction
  %     if abs(flipInfoOut(iset).orientation(i)) == 1
  %       g = g - dgdTX*(Tg-p.TRef);
  %       g0 = g0 - dgdTX*(Tg-p.TRef);
  %     end
  %     if abs(flipInfoOut(iset).orientation(i)) == 2
  %       g = g - dgdTY*(Tg-p.TRef);
  %       g0 = g0 - dgdTY*(Tg-p.TRef);
  %     end
      subplot(313)
      index = tg>=p.tCalLim(1) &  tg<=p.tCalLim(2);
      gCal = mean(g0(index));
      gCalTCor = mean(g0(index) - (Tg(index)-p.TRef)*p.dadT(flipInfoOut(iset).orientation(i)));
      plot(tg,g0,'g',[tg(1) tg(end)],[gCal gCal]','k') 
      xlabel('Time, s')

      flipInfoOut(iset).i0g(i) = i0g;
      flipInfoOut(iset).i1g(i) = i1g;
      flipInfoOut(iset).gCal(i) = gCal;
      flipInfoOut(iset).gCalTCor(i) = gCalTCor;
      flipInfoOut(iset).gFrac(i) = gFrac;
      flipInfoOut(iset).T(i) = mean(Tg);
      flipInfoOut(iset).Tstd(i) = std(Tg);
      flipInfoOut(iset).duration(i) = tg(end)-tg(1)+1;
      subplot(311)
      title({['Flip ' int2str(i) ';   Orientation ' int2str(flipInfoOut(iset).orientation(i)) ';   ' datestr(floor(data.t(flipInfoOut(iset).i0(i)))) ] , ...
             ['gCal = ' num2str(gCal,'%9.6f') '    gFrac = ' num2str(gFrac,'%9.6f')]});

    end
  end
%   %Save Plots
%   flipStr = int2str(i);
%   while length(flipStr)<4
%     flipStr = ['0' flipStr];
%   end 
%   orient tall
% %     eval(['print -djpeg ' p.dirPlot '/' int2str(flipInfoOut(iset).orientation(i)) '/flip' flipStr '.jpeg'])


end
