function [m1]=plot_Tdep_params(t,g,T,Tref)
% Solves for temperature dependence of vertical accelerometer signal,
% produces plot of data and model for visual assessment
%
% Usage
%   [m1]=plot_Tdependence(t,a,T,Tref)
%
% Let measured signal g(t) = a + bt + c(T(t)-Tref)
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

%plotting
figure
plot(t,g,'linewidth',1)
hold on
plot(t,m1(1)+m1(2)*t+m1(3)*T,'linewidth',1)
legend('signal','model')