% Axial_map_simple.m

load('calderacoordinates.mat')
calxy=llh2local(caldera_rim',[-130.1,45.9]);

% names [FI1,FO1,FO2,FI2]
dike_centers=[[7.184 10.179];[8.838 3.995];[6.319 5.123];[6.92 2.6]]; % km
dike_angles=[-70,113,-67,112]; % CCW from East
dike_lengths=[[0 4.5];[0.95 6];[4,3];[0,6]]; % in km, to the left and to the right
dike_tops=[0.0,0.2,0.5,0.0]; % km
dike_bots=[1.0,2.0,1.5,0.3];
dike_deltas=[80,58,64,45];

%----- Make a map
figure(94); clf; hold on
% caldera
plot(calxy(1,:)*1000,calxy(2,:)*1000,'k','linewidth',1)
axis('equal')
% AXCC1 location
plot(7140,5980,'xk','markersize',10,'linewidth',1)
% patch([-38 38 38 -38]+7140,[-26 -26 26 26]+5980,'k') % too small to display
% proposed inflation centers
hefxy=llh2local([-130.01;45.94],[-130.1,45.9]); % Hefner model
noonxy=llh2local([-129.99;45.95],[-130.1,45.9]); % Nooner model
plot(hefxy(1)*1000,hefxy(2)*1000,'ok','markerfacecolor','b','markersize',10,'linewidth',1)
plot(noonxy(1)*1000,noonxy(2)*1000,'sk','markerfacecolor','r','markersize',10,'linewidth',1)
%fault locations
for j=1:length(dike_centers)
    if j==1 || j==4
        sty='k:'; wd=1;
    else
        sty='k--'; wd=0.5;
    end
    quiver(dike_centers(j,1)*1000,dike_centers(j,2)*1000,...
        dike_lengths(j,2)*cosd(dike_angles(j))*1000,...
        dike_lengths(j,2)*sind(dike_angles(j))*1000,sty,'linewidth',wd,...
        'autoscale','off','showarrowhead','off')
    quiver(dike_centers(j,1)*1000,dike_centers(j,2)*1000,...
        -dike_lengths(j,1)*cosd(dike_angles(j))*1000,...
        -dike_lengths(j,1)*sind(dike_angles(j))*1000,sty,'linewidth',wd,...
        'autoscale','off','showarrowhead','off')
end

xlim([0 10000])
ylim([2000 12000])
box on; grid on
set(gca,'fontsize',14)
legend('Caldera rim','Inset center','Inflation center A','Inflation center B',...
    'Inward-dipping fault','Outward-dipping faults','location','northwest')

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print('map_test','-dtiff','-r300')