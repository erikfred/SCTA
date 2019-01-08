function [datafilt]=trace_filter(data,filtparam,samprate)
% Applies a butterworth filter to seismic data after detrending
%
% [datafilt]=trace_filter(data,filtparam,samprate)
% 
% Inputs
% DATA      - A column wise matrix of seismic traces
%             NaN's identify 'blank' values at the end of each trace
% FILTPARAM - Filter parameter structure with the following elements
%             TYPE  - String of filter type.  Valid options are:
%              'bandpass','low'(pass),'high'(pass) or (band)'stop'
%             CUT   - Filter limits in Hz
%             ORDER - Order of filter
%             PHASE - Phase of filter.  Valid options are:
%             'min'(imum phase) or 'zero'( phase)
% SAMPRATE  - Sample rate in Hz (Default is 1)
% NSAMP     - Number of samples in each trace.  Do not set if you want to 
%             determine this from the size of DATA and the presence of NaN's
%
% Outputs
% DATAFILT  - Filtered data
%
% Paramters set in function
% ERROR_CHECK - Determines whether inputs undergo error checking
% FTAPER      - Length of filter taper given by FTAPER * 1/min(FILTPARAM.CUT)
% FLEN_TAPER_MAX - Maximun allowed length of taper relative to trace length

ERROR_CHECK     = 1;
FTAPER          = 2;
FLEN_TAPER_MAX  = 0.2;

if ERROR_CHECK
  %Check number of arguments
  if nargin<2 | nargin>3
    error ('TRACE_FILTER.m requires 2-3 input arguments')
  end
  if nargout>1
    error('TRACE_FILTER.m requires <=1 output argument');
  end
  
  ntrace = size(data, 2);
  % Check filterparam structure
  if ~isstruct(filtparam)
    error('TRACE_FILTER.m - FILTPARAM must be a structure');
  end
  if ~isfield(filtparam,'type')
    error('TRACE_FILTER.m - FILTPARAM must have field ''type''');
  elseif ~strcmp(filtparam.type,'bandpass') & ~strcmp(filtparam.type,'low') & ...
     ~strcmp(filtparam.type,'high') & ~strcmp(filtparam.type,'stop') 
    error('TRACE_FILTER.m - FILTPARAM.TYPE - unrecognized option')
  end
  if ~isfield(filtparam,'cut')
    error('TRACE_FILTER.m - FILTPARAM must have field ''cut''');
  else
    if strcmp(filtparam.type,'bandpass') | strcmp(filtparam.type,'pass') 
      if length(filtparam.cut)~=2
        error('TRACE_FILTER.m - FILTPARAM.CUT must have 2 elements for bandpass/bandstop filters')
      end
    elseif length(filtparam.cut)~=1
      error('TRACE_FILTER.m - FILTPARAM.CUT must have 1 element for lowpass/highpass filters')
    end
  end
  if ~isfield(filtparam,'order')
    error('TRACE_FILTER.m - FILTPARAM must have field ''order''');
  else
    if length(filtparam.order)~=1
      error('TRACE_FILTER.m - FILTPARAM.ORDER must have 1 element')
    end
    if rem(filtparam.order,1) | filtparam.order<0
      error('TRACE_FILTER.m - FILTPARAM.ORDER must be a positive interger')
    end
  end
  if ~isfield(filtparam,'phase')
    error('TRACE_FILTER.m - FILTPARAM must have field ''phase''');
  elseif ~strcmp(filtparam.phase,'min') & ~strcmp(filtparam.phase,'zero')
    error('TRACE_FILTER.m - FILTPARAM.PHASE - unrecognized option')
  end

  % Check optional inputs
  if nargin<3;
    samprate = 2;
  else
    if length(samprate)~=1 & length(samprate) ~= ntrace
      error('TRACE_FILTER.m - SAMPRATE must be a scalar or have an element for each trace')
    end
  end
  if nargin<4;
    nsamp = [];
  else
    if length(nsamp)~=1 & length(nsamp) ~= ntrace
      error('TRACE_FILTER.m - NSAMP must be a scalar or have an element for each trace')
    end
  end
end

% Set values of NSAMP and SAMPRATE
logicalnan = ~~sum(sum(isnan(data)));
if length(nsamp)
  if length(samprate)~=1 | length(nsamp)~=1
    if length(samprate)==1
      samprate = repmat(samprate, 1, ntrace);
    end
    if length(nsamp)==1
      samprate = repmat(nsamp, 1, ntrace);
    end
  end
else
  if length(samprate) == 1
    if ~logicalnan
      nsamp = size(data,1);
    else
      samprate = repmat(samprate, 1, ntrace);
    end
  end
  if length(samprate)>1
    if logicalnan
      for i = 1:ntrace
        temp = max(find(~isnan(data(:,i))));
        if ~length(temp); temp = 0; end;
        nsamp(i) = temp;
      end
    else
      nsamp = repmat(size(data,1), 1, ntrace);
    end
  end
end

% Detrend data
if ~logicalnan
  datafilt = detrend(data);
else
  datafilt = zeros(size(data))+NaN;
  for i = 1:ntrace
    k = find(~isnan(data(:,i)));
    if length(k);
      datafilt(k,i) = detrend(data(k,i));
    end
  end
end

% Filter data
for i = 1:length(samprate);
  j = i;
  if length(samprate)==1; j = 1:ntrace; end
  wn = 2*filtparam.cut(1) / samprate(i);
  if ERROR_CHECK
    if wn>1
      error('TRACE_FILTER.m - FILTPARAM.CUT(1) > SAMPRATE/2')
    end
  end
  if ~strcmp(filtparam.type,'high') & ~strcmp(filtparam.type,'low')
    wn(2) = 2*filtparam.cut(2) / samprate(i);
    if ERROR_CHECK
      if wn(2)>1
        error('TRACE_FILTER.m - FILTPARAM.CUT(2) > SAMPRATE/2')
      end
    end
  end
  [b,a] = butter(filtparam.order,wn,filtparam.type);
%   ttaper = FTAPER / min(filtparam.cut);
%   ntaper = round(ttaper * samprate(i));
%   if ERROR_CHECK
%     if ntaper > FLEN_TAPER_MAX*nsamp(i) & ~sum(sum(isnan(data(:,j))));
%       error(['TRACE_FILTER.m - Taper extends too far into trace ' int2str(j(1))])
%     end
%   end
%   if ntaper > nsamp(i); ntaper = nsamp(i); end
%   taper = ones(size(data(:,j)));
%   temp = sin(linspace(0,pi/2,ntaper))'; 
%   temp = repmat(temp,1,length(j));
%   taper(1:ntaper,:) = taper(1:ntaper,:).*temp;
%   if strcmp(filtparam.type,'zero')
%     taper(end-ntaper+1:end,:) = taper(end-ntaper+1:end,:) .* flipud(temp);
%   end 
%   datafilt(:,j) = datafilt(:,j) .* taper;
  if nsamp(i)>0
    if ~logicalnan
      if strcmp(filtparam.phase,'min')
        datafilt(:,j) = filter(b,a,datafilt(:,j));
      elseif strcmp(filtparam.phase,'zero')
        datafilt(:,j) = filtfilt(b,a,datafilt(:,j));
      end
    else
      if j==200; keyboard; end
      k = find(~isnan(data(:,j)));
      if strcmp(filtparam.phase,'min')
        datafilt(k,j) = filter(b,a,datafilt(k,j));
      elseif strcmp(filtparam.phase,'zero')
        datafilt(k,j) = filtfile(b,a,datafilt(k,j));
      end
    end
  end
end
  
  
