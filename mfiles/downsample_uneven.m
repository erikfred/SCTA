function [tds,xds] = downsample_uneven(t,x,t0,dt,n)
% Downsamples an uneven time series by averaging within bins
% 
% Usage
%   [tds,xds] = downsample_uneven(t,x,t0,dt,n)
%
% Inputs
%   t  - Times of samples in time series
%   x  - Time series values (matrix okay)
%   t0 - Start time of first averaging window
%   dt - Length of each averaging window (output sample interval)
%   n  - Number of samples in average time series
%
% Outputs
%   tds - Times of samples in output time series (time of center of window)
%   xds - Output averaged time series (NaN's where no input samples in time interval)

% Check inputs
xDim = size(x);
if ~any(length(t(:))==xDim)
  error('downsample_uneven: Input dimensions of time and data inconsistent')
end

% Put time series into columns
t = t(:);
if xDim(1) == length(t);
  flip = false;
else
  flip = true;
  x = x';
  xDim = xDim([2 1]);
end

% Create time limits of averaged samples and output time series as NaNs 
tds = t0+(0:n)*dt;
tds = tds(:);
xds = NaN(length(tds)-1,xDim(2));

% For each input time determine the appropriate index in the output time
index = floor((t-t0)/dt+1);

% Indices for samples outside time limits with an extra index at teh start
index(index<1) = 0;
index(index>n) = n+1;
index = [-1; index];

% Find the indices corresponding to a change in output sample including the last input sample
iEnd = [~~diff(index); true];
iEnd = find(iEnd);

% Cumulative sum with extra 0 at start
xc = [zeros(1,xDim(2)); cumsum(x)];

% Averaged values
xm = (xc(iEnd(2:end),:) - xc(iEnd(1:end-1),:)) ./ repmat(iEnd(2:end)-iEnd(1:end-1),1,xDim(2));

% Determine the sample that each averaged values represents and which are for samples between 1 and n
index = index(iEnd(2:end));
i = index>0 & index<=n;

% Assign average values to output
for j = 1:xDim(2)
  xds(index(i),j) = xm(i,j);
end

% Sample time is average of interval
tds = (tds(1:end-1) + tds(2:end))/2;

% Put outputs into same orientation as inputs
if flip
  xds = xds';
  tds = tds';
end