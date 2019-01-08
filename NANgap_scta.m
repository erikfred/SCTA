function dataOut = NANgap_scta(data,dt)
% Add gaps with NaNs into an SCTA time series
%
% Usage
%   dataOut = NANgap_scta(data)
%
% Inputs
%   data - Data structure with t, a, as, T and n fields
%   dt   - Time increment in seconds of samples (will calculate otherwise)
%
% Outputs
%   dataOut - Data structure with t, a, as, T and n fields where t is uniformly sampled 
%   through data gaps and other fields are filled in with Nans

if nargin<2
  dt = min(diff(data.t));
  n = round((data.t(end) - data.t(1))/dt);
  dt = (data.t(end) - data.t(1))/n;
end

dataOut.t = linspace(data.t(1),data.t(end),n+1)';
dataOut.a = NaN(n+1,3);
dataOut.as = NaN(n+1,1);
dataOut.T = NaN(n+1,1);
dataOut.n = NaN(n+1,1);

i = round((data.t-data.t(1))/dt) + 1;
dataOut.a(i,:) = data.a;
dataOut.as(i) = data.as;
dataOut.T(i) = data.T;
dataOut.n(i) = data.n;


