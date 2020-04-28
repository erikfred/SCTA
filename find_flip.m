function [flipInfo,lNormOrt] = find_flip(t,a,as,p)
% Finds flips in an accelerometer time series
% 
% Usage
%   [flipInfo,lNormalOrtientation] = find_flip(t,a,as,p)
%
% Inputs
%   t  - time series
%   a  - acceleration (3 components)
%   as - scalar acceleration
%   p  - parameter structure with fields
%          cosThreshVert  - Accelerometer is in a vertical orientation when one acceleration is this fraction of total
%          minTime4Flip   - A flip requires at least this much time (s) in flipped orientation
%          complexRange80 - If 80% of the flipped time series is not within this range then the flip is complex with the 
%                             accelerometer in multiple near-vertical orientations
%          cosThreshNorm  - This is the fractional threshold for determining that the accelerometer is in its normal 
%                             orientation.  It is used for plotting and is normally set to a higher fraction than cosThreshVert  
%          tBufferNorm    - Samples within this time of an interval of  non-normal orientation are deemed to be non-normal
%          nMadNorm       - If total acceleration in normal position is this many MADs from the median, make it non-normal
%          
%
% Outputs
%  flipInfo - Structure with information about flips with following fields
%               normalOrientation = The normal orientation of the sensor (usually 3)
%               i0 - Vector with first sample of each flip
%               i1 - Vector with last sample of each flip
%               orientation - Vector with orientation of each flip
%               t - Vector with serial data number of each flip
%               aMed - Vector with median acceleration of flipped channel 
%               range80 - Vector with 80th percentile range of accelerations during the flip
%               complex - Indicates that the flip was complex (sensor moved while near vertical)
%               roughDuration - Rough duration of flip
%
%  lNormalOrtientation - Logical vector of length t set to true when accelerometer is in its normal orientation
%
% Modified 11/6 to flip break flips up into groups when not successive and to log duration

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

% No flips found
if isempty(i0) || isempty(i1)
  flipInfo.normalOrientation = normalOrientation;
  flipInfo.i0 = [];
  flipInfo.i1 = [];
  flipInfo.orientation =  [];
  flipInfo.t =  [];
  flipInfo.aMed =  [];
  flipInfo.range80 =  [];
  flipInfo.complex =  [];
  flipInfo.duration = [];
  lNormOrt = ones(size(t));
  return
end

if i0(1)>i1(1)
  i0 = [1 i0];
end
if i1(end)<i0(end)
  i1 = [i1 length(flipped)];
end
i = (t(i1)-t(i0))>=p.minTime4Flip/86400;
i0 = i0(i);
i1 = i1(i);

% Now length of orientation is adjusted from length of time series to number of flips
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

%% Breakup flips into sets
ibreak = find((diff(flipInfo.t)*86400 - flipInfo.roughDuration(1:end-1))>60);
if ~isempty(ibreak)
  ib0 = [1; ibreak+1];
  ib1 = [ibreak; length(i0)];
  for i=length(ibreak)+1:-1:1
    flipInfo(i).normalOrientation = flipInfo(1).normalOrientation;
    flipInfo(i).i0 = flipInfo(1).i0(ib0(i):ib1(i));
    flipInfo(i).i1 = flipInfo(1).i1(ib0(i):ib1(i));
    flipInfo(i).orientation = flipInfo(1).orientation(ib0(i):ib1(i));
    flipInfo(i).t = flipInfo(1).t(ib0(i):ib1(i));
    flipInfo(i).aMed = flipInfo(1).aMed(ib0(i):ib1(i));
    flipInfo(i).range80 = flipInfo(1).range80(ib0(i):ib1(i));
    flipInfo(i).complex = flipInfo(1).complex(ib0(i):ib1(i));
    flipInfo(i).roughDuration = flipInfo(1).roughDuration(ib0(i):ib1(i));   
  end
end

%% Create a logical that is true for samples that are in the normal orienation, set apart from flips and without extreme accelerations
medianAs = nanmedian(as);
madAs = mad(as(~isnan(as)));
% Normal orientation when near vertical and acceleration not to far from mean values
lNormOrt = a(:,normalOrientation)./as>p.cosThreshNorm & as>medianAs-p.nMadNorm*madAs & as<medianAs+p.nMadNorm*madAs & ~isnan(as);
% First samples when orientation is not normal
on  = find(diff(lNormOrt)==-1)+1;
% Last sample when orientation is not normal
off = find(diff(lNormOrt)==1);
if on(1)>off(1)
  on = [1; on];
end
if off(end)<on(end)
  off = [off; length(lNormOrt)];
end
% Make anything close to a non normal orientation, non-normal
for i=1:length(on)
  lNormOrt(t>=t(on(i))-p.tBufferNorm/86400 & t<=t(off(i))+p.tBufferNorm/86400) = false;
end

