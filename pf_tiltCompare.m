% Processes pf tilt
startDate = datenum('Nov 20, 2018');
endDate = datenum('Nov 27, 2018');
np = 3;


load ../calibrations/PinonFlat/pfdata

i = dataDec1.t>=startDate & dataDec1.t<endDate;
x1 = dataDec1.a(i,1);
x1 = detrend(x1);
y1 = dataDec1.a(i,2);
y1 = detrend(y1);
t1 = dataDec1.t(i);
n = length(t1);

i = tilt2Dec1.t>=startDate & tilt2Dec1.t<endDate;
x2 = tilt2Dec1.a(i,1);
x2 = detrend(x2);
y2 = tilt2Dec1.a(i,2);
y2 = detrend(y2);
t2 = tilt2Dec1.t(i);

x1 = x1 - polyval(polyfit(t1-t1(1),x1,np),t1-t1(1));
x2 = x2 - polyval(polyfit(t2,x2,np),t2);
y1 = y1 - polyval(polyfit(t1-t1(1),y1,np),t1-t1(1));
y2 = y2 - polyval(polyfit(t2,y2,np),t2);

figure
clf
subplot(211)
plot(t1,x1,t2,x2)
hold on
datetick
title('X - channel')
subplot(212)
plot(t1,y1,t2,y2)
datetick
title('Y - channel')
hold on


a = 1;
b = 1;
theta = -0.3;
damp = 0.7;
% x2pred = a*x1*cos(theta) + b*y1*sin(theta);
% y2pred = -a*x1*sin(theta) + b*y1*cos(theta);
% subplot(211)
% plot(t2,x2pred);
% subplot(212)
% plot(t2,y2pred);
% return

for i=1:10

  x2pred = a*x1*cos(theta) + b*y1*sin(theta);
  y2pred = -a*x1*sin(theta) + b*y1*cos(theta);
  if i==10
    subplot(211)
    plot(t2,x2pred,t2,x2-x2pred)
%     legend('x_{1}','x_{2}','x_{1}-x_{2}','x_{1,rot}','x_{1,rot}-x_{2}','location','northwest')
     legend('x_{1}','x_{2}','x_{1,rot}','x_{1,rot}-x_{2}','location','northwest')
     ylabel('Acceleration m/s^{2}')

    subplot(212)
    plot(t2,y2pred,t2,y2-y2pred)
%     legend('y_{1}','y_{2}','y_{1}-y_{2}','y_{1,rot}','y_{1,rot}-y_{2}','location','northwest')
     legend('y_{1}','y_{2}','y_{1,rot}','y_{1,rot}-y_{2}','location','northwest')
    xlabel(['Rotation applied = ' num2str(theta*180/pi) '°'])
     ylabel('Acceleration m/s^{2}')
  end
  
  d = [x2-x2pred; y2-y2pred];

  m = [a b theta]';

  G = [x1*cos(theta) y1*sin(theta) -a*x1*sin(theta)+b*y1*cos(theta);
       -x1*sin(theta) y1*cos(theta) -a*x1*cos(theta)-b*y1*cos(theta)];
   G = G(:,3);
     
  dm = inv(G'*G)*G'*d;

%   a = a + dm(1)*damp;
%   b = b + dm(2)*damp;
%   theta = theta + dm(3)*damp;
    theta = theta + dm(1)*damp;
  
end

save pf_tiltCompare theta
print -djpeg pf_tiltCompare.jpeg
print -dtiff pf_tiltCompare.tiff -r300