%get_span_PF.m
%
% Calculates the span of both x- and y-channels under the
% +Z->+X->+Y->-Y->+X->-X->+Z flipping schema.
%
%

clear; close all

% Temperature sensitivity parameters
p.dadT=[6.2023e-5 2.9562e-5 NaN];
p.TRef = 30;

% Find flip parameters
p.cosThreshVert = 0.99;         % 0.99 = 8 degrees from vertical - threshold for a flip
p.minTime4Flip = 100;           % Minimum duration in seconds for a flip to be counted
p.complexRange80 = 1e-2;        % If 90% value - 10% value is greater than this then not a simple flip into one orientation
p.cosThreshNorm = 0.9996;       % Threshold for normal (one channel vertical) orientation (2�)
p.tBufferNorm = 120;            % Make non-normal any sample within this of a non-normal orientation
p.nMadNorm = 6;                 % If total acceleration in normal position is this many MADs from the median, make it non-normal

% Process flips parameters
p.daMax = 1e-4;                 % During a calibration successive samples will not change by more than this
p.tCalLim = [60 90];            % Time limits for calibration in seconds since start of stable output

t0=datenum(2019,8,9);
tf=floor(now);

t1=t0;
count=0;
span.ty=[];span.yup=[];span.ydn=[];span.tx=[];span.xup=[];span.xdn=[];
span.yup_theta=[];span.ydn_theta=[];span.xup_theta=[];span.xdn_theta=[];
while t1<tf
    count=count+1;
    data=get_sctaDay('/Volumes/GoogleDrive/My Drive/Oceanography/SCTA-Share/OOI-PF/SCTA-PF (1)/ParsedData',t1);
    
    %plot to ensure flipped properly
    figure(6)
    clf
    subplot(211)
    plot(data.t,data.a(:,1),'linewidth',1)
    ylim([-10 10])
    xlim([t1+21/24 t1+21.25/24])
    datetick('x','keeplimits')
    set(gca,'fontsize',15)
    subplot(212)
    plot(data.t,data.a(:,2),'linewidth',1)
    ylim([-10 10])
    xlim([t1+21/24 t1+21.25/24])
    datetick('x','keeplimits')
    set(gca,'fontsize',15)
    
    %calculate spans
    a=data.a;
    as=data.as;
    t=data.t;
    % orientation for each sample is set to
    %    0 if no channel is vertical
    %    1, 2, or 3 if that channel is vertical up
    %    -1, -2, or -3 if that channel is vertical down
    [~,orientation] = max(abs(a),[],2);
    [~,temp] = max(a,[],2);
    orientation(orientation~=temp) = -orientation(orientation~=temp);
    orientation(max(abs(a),[],2)./as < p.cosThreshVert) = 0;
    orientation(~~sum(isnan(a),2)) = 0;
    
    % Normal (most common) orientation
    normalOrientation = mode(orientation(orientation~=0));
    
    % Start and End samples of intervals when the sensor is flipped
    % Flipped requires one sensor is vertical and it is not in the normal
    % orientation and it lasts long enough
    flipped = orientation~=normalOrientation & orientation~=0;
    i0 = find(diff(flipped)==1)+1;
    i1 = find(diff(flipped)==-1);
    
    if i0(1)>i1(1)
        i0 = [1 i0];
    end
    if i1(end)<i0(end)
        i1 = [i1 length(flipped)];
    end
    i = (t(i1)-t(i0))>=p.minTime4Flip/86400;
    i0 = i0(i);
    i1 = i1(i);
    
    % Now length of orientation is adjusted from lenght of time series to number of flips
    orientation = orientation(i0);
    
    % Set the flip Information structure
    flipInfo.normalOrientation = normalOrientation;
    flipInfo.i0 = i0;
    flipInfo.i1 = i1;
    flipInfo.orientation = orientation;
    flipInfo.t = t(i0);
    
    flipInfo.aMed = zeros(size(flipInfo.i0));
    flipInfo.range80 = zeros(size(flipInfo.i0));
    for i=1:length(flipInfo.i0)
        flipInfo.aMed(i) = median(a(i0(i):i1(i),abs(orientation(i))));
        temp = sort(a(i0(i):i1(i),abs(orientation(i))));
        flipInfo.range80(i) = temp(round(length(temp)*0.9)) - temp(round(length(temp)*0.1));
    end
    
    flipInfo.complex = flipInfo.range80>p.complexRange80;
    flipInfo.roughDuration = (t(i1) - t(i0)) * 86400;
    
    % calculate spans, etc.
    span.ty=[span.ty t(flipInfo.i0(2))];
    span.yup=[span.yup median(a(flipInfo.i0(2)+60*40:flipInfo.i0(2)+90*40,2))];
    span.ydn=[span.ydn median(a(flipInfo.i0(3)+60*40:flipInfo.i0(3)+90*40,2))];
    span.tx=[span.tx t(flipInfo.i0(4))];
    span.xup=[span.xup median(a(flipInfo.i0(4)+60*40:flipInfo.i0(4)+90*40,1))];
    span.xdn=[span.xdn median(a(flipInfo.i0(5)+60*40:flipInfo.i0(5)+90*40,1))];
    span.yup_theta=[span.yup_theta atand(median(a(flipInfo.i0(2)+60*40:flipInfo.i0(2)+90*40,1))/...
        median(a(flipInfo.i0(2)+60*40:flipInfo.i0(2)+90*40,2)))];
    span.ydn_theta=[span.ydn_theta atand(median(a(flipInfo.i0(3)+60*40:flipInfo.i0(3)+90*40,1))/...
        median(a(flipInfo.i0(3)+60*40:flipInfo.i0(3)+90*40,2)))];
    span.xup_theta=[span.xup_theta atand(median(a(flipInfo.i0(4)+60*40:flipInfo.i0(4)+90*40,2))/...
        median(a(flipInfo.i0(4)+60*40:flipInfo.i0(4)+90*40,1)))];
    span.xdn_theta=[span.xdn_theta atand(median(a(flipInfo.i0(5)+60*40:flipInfo.i0(5)+90*40,2))/...
        median(a(flipInfo.i0(5)+60*40:flipInfo.i0(5)+90*40,1)))];
    
    % separate fig to inpect flips closely
    figure(7)
    subplot(221)
    hold on
    h(count)=plot(data.t-t1,data.a(:,1),'linewidth',1);
    xlim([t(flipInfo.i0(4))-t1+10/60/60/24 t(flipInfo.i1(4))-t1-10/60/60/24])
    datetick('x','keeplimits')
    set(gca,'fontsize',12)
    title('+x')
    subplot(222)
    hold on
    plot(data.t-t1,data.a(:,1),'linewidth',1)
    xlim([t(flipInfo.i0(5))-t1+10/60/60/24 t(flipInfo.i1(5))-t1-10/60/60/24])
    datetick('x','keeplimits')
    set(gca,'fontsize',12)
    title('-x')
    subplot(223)
    hold on
    plot(data.t-t1,data.a(:,2),'linewidth',1)
    xlim([t(flipInfo.i0(2))-t1+10/60/60/24 t(flipInfo.i1(2))-t1-10/60/60/24])
    datetick('x','keeplimits')
    set(gca,'fontsize',12)
    title('+y')
    subplot(224)
    hold on
    plot(data.t-t1,data.a(:,2),'linewidth',1)
    xlim([t(flipInfo.i0(3))-t1+10/60/60/24 t(flipInfo.i1(3))-t1-10/60/60/24])
    datetick('x','keeplimits')
    set(gca,'fontsize',12)
    title('-y')
    
    subplot(221)
    plot(t(flipInfo.i0(4)+60*40:flipInfo.i0(4)+90*40)-t1,a(flipInfo.i0(4)+60*40:...
        flipInfo.i0(4)+90*40,1),'k','linewidth',3)
    subplot(222)
    plot(t(flipInfo.i0(5)+60*40:flipInfo.i0(5)+90*40)-t1,a(flipInfo.i0(5)+60*40:...
        flipInfo.i0(5)+90*40,1),'k','linewidth',3)
    subplot(223)
    plot(t(flipInfo.i0(2)+60*40:flipInfo.i0(2)+90*40)-t1,a(flipInfo.i0(2)+60*40:...
        flipInfo.i0(2)+90*40,2),'k','linewidth',3)
    subplot(224)
    plot(t(flipInfo.i0(3)+60*40:flipInfo.i0(3)+90*40)-t1,a(flipInfo.i0(3)+60*40:...
        flipInfo.i0(3)+90*40,2),'k','linewidth',3)
    
    keyboard
    t1=t1+7;
end

% add legend to inspection plot
figure(7)
subplot(221)
legend([h(1),h(2),h(3),h(4),h(5)],datestr(t1-35,'dd-mmm'),datestr(t1-28,'dd-mmm'),datestr(t1-21,'dd-mmm'),...
    datestr(t1-14,'dd-mmm'),datestr(t1-7,'dd-mmm'),'location','southeast')

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../calibrations/PinonFlat/span/inspect','-dtiff')
saveas(gcf,'../calibrations/PinonFlat/span/inspect.fig')

% Plot spans
figure(8)
clf
yyaxis left
plot(span.ty,span.yup-span.ydn,'o','markersize',15)
hold on
ylabel('Y (m/s^2)')
lim1=ylim;
yyaxis right
plot(span.tx,span.xup-span.xdn,'x','markersize',15)
hold on
datetick('x')
set(gca,'fontsize',15)
ylabel('X (m/s^2)')
lim2=ylim;
title('PinonFlat SCTA span')
%ensure equal axes
if (lim1(2)-lim1(1))>(lim2(2)-lim2(1))
    yyaxis right
    ylim([lim2(1) lim2(2)+(lim1(2)-lim1(1))-(lim2(2)-lim2(1))])
else
    yyaxis left
    ylim([lim1(1) lim1(2)+(lim2(2)-lim2(1))-(lim1(2)-lim1(1))])
end

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../calibrations/PinonFlat/span/spanVtime','-dtiff')
saveas(gcf,'../calibrations/PinonFlat/span/spanVtime.fig')

%save variables
save('../calibrations/PinonFlat/span/spanstruct','span')

%temperature correction?

%convert acceleration to tilt
AXCC2.LAX=asin(data.a(:,1)/9.81)*10^6;
AXCC2.LAY=asin(data.a(:,2)/9.81)*10^6;