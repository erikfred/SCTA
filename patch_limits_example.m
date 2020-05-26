% patch_limits_example.m
%
% simplified example of 'patch' function causing ylim to change
%

figure(90)
clf
subplot(411)
hold on
plot(1:100,rand(1,100),'o')
box on
lim_y=ylim;
hp=patch([60 80 80 60],[lim_y(1) lim_y(1) lim_y(2) lim_y(2)],[0.8 0.8 0.8]);
hp.EdgeColor='none';
hp.FaceVertexAlphaData=0.2;
hp.FaceAlpha='flat';

subplot(412)
hold on
plot(1:100,rand(1)*rand(1,100),'o')
box on
lim_y=ylim;
hp=patch([60 80 80 60],[lim_y(1) lim_y(1) lim_y(2) lim_y(2)],[0.8 0.8 0.8]);
hp.EdgeColor='none';
hp.FaceVertexAlphaData=0.2;
hp.FaceAlpha='flat';

subplot(413)
hold on
plot(1:100,rand(1)*10^6*rand(1,100),'o')
box on
lim_y=ylim;
hp=patch([60 80 80 60],[lim_y(1) lim_y(1) lim_y(2) lim_y(2)],[0.8 0.8 0.8]);
hp.EdgeColor='none';
hp.FaceVertexAlphaData=0.2;
hp.FaceAlpha='flat';

subplot(414)
hold on
plot(1:100,rand(1,100)/rand(1),'o')
box on
lim_y=ylim;
hp=patch([60 80 80 60],[lim_y(1) lim_y(1) lim_y(2) lim_y(2)],[0.8 0.8 0.8]);
hp.EdgeColor='none';
hp.FaceVertexAlphaData=0.2;
hp.FaceAlpha='flat';

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
% orient tall
print ../test -djpeg
!open ../test.jpg