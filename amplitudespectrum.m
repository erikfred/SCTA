function [f,a,z]=amplitudespectrum(x,y)
%Returns the amplidude spectrum of an evenly spaced series (Y(x))
%
%Usage
%  [f,a]=amplitudespectrum(x,y)
%Inputs
%  X - Evenly spaced x values of series
%  Y - Y values of the series at X
%Outputs
%  F - Frequencies
%  A - Amplitude Spectrum
%  Z - Angle

% Check inputs
if length(x)~=length(y)
  error('Lengths of ''x'' and ''y'' are not the same')
end
% if any((diff(diff(x)))>mean(diff(x))/1e6) %originally 1e6
%   error('amplitudespectrum - values of ''x'' are not evenly spaced')
% end

% Compute the FT (note scaling by dx to make it equivalent to the continuous FT 
%                  - see Section 4.4 in Bob Crosson's notes)
n = length(y);
dx = x(2)-x(1);
f = linspace(0,1/(2*(x(2)-x(1))),floor(n/2)+1);
df = f(2)-f(1);
a = dx*abs(fftshift(fft(y)));
z = angle(fftshift(fft(y)));
if ~rem(n,2)
  f = [-f(end:-1:2) f(1:end-1)];
else
  f = [-f(end:-1:2) f];
end
