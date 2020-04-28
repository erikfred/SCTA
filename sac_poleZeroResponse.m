function trans = sac_poleZeroResponse(f,sacpz,type)
% Calculate the transfer ruction given a SAC response
%
% Usage
%   trans = sac_poleZeroResponse(f,sacpz)
%
% Inputs
%   F     - Vector of frequencies
%   SACPZ - Structure with SAC style displacement response with fields
%             constant - Scaling
%             zeros    - Vector of response zeroes
%             poles    - Vector of response poles
%           Note: SAC response is for DISPLACEMENT
%   TYPE  - Type of response
%             'Disp' - Displacement 
%             'Vel'  - Velocity (default)
%             'Acc'  - Acceleration
%
% Outputs
%   TRANS - Transfer function 
%
% July 18, 2017 - Fixed to take into account that SAC pole-zero response is
%                 for displacement
% July 19, 2017 - Added use of constantFix to apply correction to incorrect SP response

if nargin<3
  type = 'V';
end
if isempty(type)
  type = 'V';
end

s = sqrt(-1)*2*pi*f;

trans = ones(size(f))*sacpz.constant;
for i=1:length(sacpz.zeros)
  trans = trans.*(s-sacpz.zeros(i));
end
for i=1:length(sacpz.poles)
  trans = trans./(s-sacpz.poles(i));
end

if strcmpi(type(1),'v')
  trans = trans./s;
elseif strcmpi(type(1),'a')
  trans = trans./s./s;
end

if isfield(sacpz,'constantFix')
  trans = trans .* sacpz.constantFix;
end