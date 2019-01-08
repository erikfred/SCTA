function [dataDec] = decimate_SCTA(data,sampInt,dataDec)
% Function to decimate/downsample SCTA and tilt data by averaging 
%
% Usage
%   [dataDec] = decimate_SCTA(data,sampInt,dataDec)
%
% Inputs
%   data    - SCTA structure with fields t, a, as, and T
%               or
%             TILT data structure with fields t, a, T
%   sampInt - Sample interval(s) in seconds of decimated data
%   dataDec - Optionally provide decimated output for appending
%
% Outputs
%   dataDec - Decimated data structure with fields t, a, {as}, T, n
%             Will be a vector if sampInt is a vector
%
% Use Fast decimation of unevenly spaced time series by averaging

na = size(data.a,2);

% Automated flipping apparatus
for i=1:length(sampInt)
  if  ~isempty(data.t)
    [t,a,n] = downsample_uneven2(data.t,data.a,sampInt(i)/86400);
    [~,T]   = downsample_uneven2(data.t,data.T,sampInt(i)/86400);
    if isfield(data,'as')
      [~,as]   = downsample_uneven2(data.t,data.as,sampInt(i)/86400);
    end
    if nargin<3 || ~isfield(dataDec,'t') || isempty(dataDec(i).t)
      dataDec(i).t = t;
      dataDec(i).a = a;
      dataDec(i).T = T;
      if isfield(data,'as')
        dataDec(i).as = as;
      end
      dataDec(i).n = n;
    else
      % Add to end with one sample overlap (average)
      if abs(t(1)-dataDec(i).t(end))<0.5*sampInt(i)/86400  
        dataDec(i).a(end,:) = (dataDec(i).a(end,:)*dataDec(i).n(end) + a(1,:)*n(1)) ...
                                  ./repmat(dataDec(i).n(end)+n(1),1,na);
        if isfield(data,'as')
          dataDec(i).as(end) = (dataDec(i).as(end)*dataDec(i).n(end) + as(1)*n(1)) ...
                                  /(dataDec(i).n(end)+n(1));
        end
        dataDec(i).T(end) = (dataDec(i).T(end)*dataDec(i).n(end) + T(1)*n(1)) ...
                                  /(dataDec(i).n(end)+n(1));
        dataDec(i).n(end) = dataDec(i).n(end) + n(1);
        dataDec(i).t = [dataDec(i).t; t(2:end)];
        dataDec(i).a = [dataDec(i).a; a(2:end,:)];
        dataDec(i).n = [dataDec(i).n; n(2:end)];
        if isfield(data,'as')
          dataDec(i).as = [dataDec(i).as; as(2:end)];
        end
        dataDec(i).T = [dataDec(i).T; T(2:end)];
      % Add to end with no gap
      elseif abs(t(1)-dataDec(i).t(end)-sampInt(i)/86400)<0.5*sampInt(i)/86400
        dataDec(i).t = [dataDec(i).t; t];
        dataDec(i).a = [dataDec(i).a; a];
        dataDec(i).n = [dataDec(i).n; n];
        if isfield(data,'as')
          dataDec(i).as = [dataDec(i).as; as];
        end
        dataDec(i).T = [dataDec(i).T; T];
      % Add to end with overlap of more than one sample
      elseif abs(t(1)-dataDec(i).t(end))>0.5*sampInt(i)/86400 && t(1)>dataDec(i).t(1);
        [~,ilast] = min(abs(dataDec(i).t-t(1)));
        ilast = ilast-1;
        dataDec(i).t = [dataDec(i).t(1:ilast); t];
        dataDec(i).a = [dataDec(i).a(1:ilast,:); a];
        dataDec(i).n = [dataDec(i).n(1:ilast); n];
        if isfield(data,'as')
          dataDec(i).as = [dataDec(i).as(1:ilast); as];
        end
        dataDec(i).T = [dataDec(i).T(1:ilast); T];
      % Add to end with empty samples
      elseif t(1)>dataDec(i).t(end)
        tmiss = dataDec(i).t(end):sampInt(i):t(end)+0.5*sampInt(i);
        tmiss = tmiss(2:end-1)';
        nmiss = length(tmiss);
        dataDec(i).t = [dataDec(i).t; tmiss; t];
        dataDec(i).a = [dataDec(i).a; NaN(nmiss,na); a];
        dataDec(i).n = [dataDec(i).n; zeros(nmiss,1); n];
        if isfield(data,'as')
          dataDec(i).as = [dataDec(i).as; NaN(nmiss,1); as];
        end
        dataDec(i).T = [dataDec(i).T; NaN(nmiss,1); T];
      % Presently this is all that is allowed
      else
        disp('Problems Downsampling')
        keyboard
      end
    end
  end
end

