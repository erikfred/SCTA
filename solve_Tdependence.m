function [lGlT,qGlT]=solve_Tdependence(t,g,T,tSkip)
% Solves for temperature dependence of vertical accelerometer signal
%
% Usage
%   [lGlT qGlT lGqT qGqT]=solv_eTdependence(t,a,T,tSkip)
%
% Let measured signal g(t) = a + bt + [ct^2] + dT(t) + eT(t)^2
%   where a-e are unknown constants, 
%   b and c give linear and quadratic trend in underlying signal
%   d and e give linear and quadratic temperature dependence
%
%  We can write g = A m where
%     A = [ 1  dT   dT^2  T(1) T(1)^2
%           1 2dT (2dT)^2 T(2) T(2)^2
%           1 3dT (3dT)^2 T(3) T(3)^2
%           ......
%           ......
%           1 NdT (NdT)^2 T(N) T(N)^2 ]
%      m = [a b c d e]'
%   and the 3rd and 5th columns of A and rows of m may be dropped based on the choice of 
%   underlying trend.
%   
%   Solve as least squares following Menke section 3.4
%   Some normalization of A and m needed to avoid singular matrix (i.e., coefficients in 
%   linear equations must not be wildly different)

t = t-min(t);

if nargin>=4 && tSkip>0
  i = t>=tSkip;
  t = t(i);
  g = g(i);
  T = T(i);
end
  
N = length(t);

% Vertical acceleration demeaned
gmean = mean(g);
g = g-gmean;
% Temperature demeaned
Tmean = mean(T);
T = T - Tmean;

% All linear with time in days
A1 = [ones(N,1) t T];
m1 = inv(A1' * A1) * A1' * g;
g1 = A1 * m1;
lGlT=m1(3);
disp(['Linear g, Linear T: dg/dT = ' num2str(m1(3)) ' m/s^2/°C; misfit = ' num2str(mean(sqrt((g-g1).^2))) ' m/s^2'])

% Quadratic trend in g and linear T dependence
A2 = [ones(N,1) (0:N-1)'/(86400*4) (((0:N-1)/(86400*4)).^2)' T];
m2 = inv(A2' * A2) * A2' * g;
g2 = A2 * m2;
qGlT=m2(4);
disp(['Quadratic g, Linear T: dg/dT = ' num2str(m2(4)) ' m/s^2/°C; misfit = ' num2str(mean(sqrt((g-g2).^2))) ' m/s^2'])

% % Linear trend in g and quadratic T dependence
% A3 = [ones(N,1) (0:N-1)'/(86400*4) T T.^2];
% m3 = inv(A3' * A3) * A3' * g;
% g3 = A3 * m3;
% lGqT=m3(3);
% disp(['Linear g, Quadratic T: dg/dT = ' num2str(m3(3)) ' m/s^2/°C; d^2g/dT^2 = ' num2str(m3(4)) ...
%       ' m/s^2/°C^2; misfit = ' num2str(mean(sqrt((g-g3).^2))) ' m/s^2'])
% 
% % All quadratic
% A4 = [ones(N,1) (0:N-1)'/(86400*4) (((0:N-1)/(86400*4)).^2)' T T.^2];
% m4 = inv(A4' * A4) * A4' * g;
% g4 = A4 * m4;
% qGqT=m4(4);
% disp(['Quadratic g, Quadratic T: dg/dT = ' num2str(m4(4)) ' m/s^2/°C; d^2g/dT^2 = ' num2str(m4(5)) ...
%        ' m/s^2/°C^2; misfit = ' num2str(mean(sqrt((g-g4).^2))) ' m/s^2'])

%clf
plot(t,g)
hold all
plot(t,g1)
plot(t,g2)
% plot(t,g3)
% plot(t,g4)
legend('Data','Linear g, Linear T dep','Quadratic g, Linear T dep','Linear g, Quadratic T dep','Quadratic g, Quadratic T dep')
ylabel('Accel, m/s')
xlabel('Day')