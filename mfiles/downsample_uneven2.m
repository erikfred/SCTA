function [tds,xds,nds] = downsample_uneven2(t,x,dt)
% Downsamples an uneven time series by averaging within bins
% This second version does not specify output sample times but instead obtains 
% them from the input data.  It also returns the number of samples in each average
% 
% Usage
%   [tds,xds,nused] = downsample_uneven(t,x,dt)
%
% Inputs
%   t  - Times of samples in time series
%   x  - Time series values (matrix okay)
%   dt - Length of each averaging window (output sample interval)
%
% Outputs
%   tds - Times of samples in output time series (time of center of window)
%   xds - Output averaged time series (NaN's where no input samples in time interval)
%   nds - Number of samples averaged for each output in 

% Small values compared with sample interval (in days)
small = 1e-3/86400;

% Check inputs have consistent dimensions
xDim = size(x);
if ~any(length(t(:))==xDim)
  error('downsample_uneven2: Input dimensions of time and data inconsistent')
end

% Put time series into column orientation
t = t(:);
if xDim(1) == length(t)
  flip = false;
else
  flip = true;
  x = x';
  xDim = xDim([2 1]);
end

% Create time limits of averaged samples and output time series as NaNs 
t0 = floor(t(1)/dt)*dt;
t1 = ceil((t(end)+small)/dt)*dt;
tds = t0:dt:t1+dt/1e6;
tds = tds(:);
xds = NaN(length(tds)-1,xDim(2));
nds = zeros(length(tds)-1,1);

% For each input time determine the appropriate index in the output time
index = floor((t-t0)/dt+1);

% Find the indices corresponding to a change in output sample including the first and last input sample
iEnd = [true; ~~diff(index); true];
iEnd = find(iEnd);

% Cumulative sum with extra 0 at start
xc = [zeros(1,xDim(2)); cumsum(x)];

% Averaged values
n = iEnd(2:end)-iEnd(1:end-1);
xm = (xc(iEnd(2:end),:) - xc(iEnd(1:end-1),:)) ./ repmat(n,1,xDim(2));

% Determine the sample that each averaged values represents 
index = index(iEnd(1:end-1));

% Assign average values to output
for j = 1:xDim(2)
  xds(index,j) = xm(:,j);
end
nds(index) = n;

% Sample time is average of interval
tds = (tds(1:end-1) + tds(2:end))/2;

% Put outputs into same orientation as inputs
if flip
  xds = xds';
  tds = tds';
end