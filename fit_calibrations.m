load ../calibrations/Axial/axialdata

% range to include in fit
ix1=43:3:67; iy=44:3:68; ix2=45:3:69;

% first +X
t_x1=flipInfoAll.t(ix1);
a_x1=flipInfoAll.gCal(ix1);
px1=polyfit(t_x1-t_x1(1),a_x1,1);
x1_star=(t_x1-t_x1(1))*px1(1)+px1(2);
rmse_x1=sqrt(mean((a_x1-x1_star).^2))

figure(11)
clf
subplot(311)
plot(t_x1,a_x1,'bo')
hold on
plot(t_x1,x1_star,'k*')

% second +X
t_x2=flipInfoAll.t(ix2);
a_x2=flipInfoAll.gCal(ix2);
px2=polyfit(t_x2-t_x2(1),a_x2,1);
x2_star=(t_x2-t_x2(1))*px2(1)+px2(2);
rmse_x2=sqrt(mean((a_x2-x2_star).^2))

subplot(312)
plot(t_x2,a_x2,'bo')
hold on
plot(t_x2,x2_star,'k*')

% +Y
t_y=flipInfoAll.t(iy);
a_y=flipInfoAll.gCal(iy);
py=polyfit(t_y-t_y(1),a_y,1);
y_star=(t_y-t_y(1))*py(1)+py(2);
rmse_y=sqrt(mean((a_y-y_star).^2))

subplot(313)
plot(t_y,a_y,'bo')
hold on
plot(t_y,y_star,'k*')