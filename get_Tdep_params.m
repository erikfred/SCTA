function [m1]=get_Tdep_params(t,g,T,Tref)
% Solves for temperature dependence of vertical accelerometer signal
%
% Usage
%   [m1]=solve_Tdependence(t,a,T,tSkip)
%
% Let measured signal g(t) = a + bt + cT(t)
%   where a-c are unknown constants, 
%   b gives linear trend in underlying signal
%   d gives linear and temperature dependence
%
%  We can write g = A m where
%     A = [ 1  dT  T(1)
%           1 2dT  T(2)
%           1 3dT  T(3)
%           ......
%           1 NdT  T(N) ]
%      m = [a b c d e]'
%   
%   Solve as least squares following Menke section 3.4
%   Some normalization of A and m needed to avoid singular matrix (i.e., coefficients in 
%   linear equations must not be wildly different)

t = t-min(t);
  
N = length(t);

% Vertical acceleration demeaned
gmean = mean(g);
g = g-gmean;

T = T - Tref;

% All linear with time in days
A1 = [ones(N,1) t T];
m1 = inv(A1' * A1) * A1' * g;