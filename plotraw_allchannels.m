% plotraw_allchannels.m
%
% For each deployment, plot X, Y, Z, A at 1 sample/100 s without any kind
% of alignment, drift correction, etc.
%

clear; close all

% Axial pre-move
load('../calibrations/Axial/axialdata.mat','dataDec100')

figure(1); clf
subplot(411); hold on
plot(dataDec100.t,(dataDec100.a(:,1)+0.035)*10^5,'linewidth',1)
datetick('x',3)
ylim(([-0.04 -0.03]+0.035)*10^5)
ylabel('X (\mug)')
box on; grid on
subplot(412)
plot(dataDec100.t,(dataDec100.a(:,2)+0.024)*10^5,'linewidth',1)
datetick('x',3)
ylim(([-0.029 -0.019]+0.024)*10^5)
ylabel('Y (\mug)')
box on; grid on
subplot(413)
plot(dataDec100.t,(dataDec100.a(:,3)-9.815)*10^5,'linewidth',1)
datetick('x',3)
ylim(([9.81 9.82]-9.815)*10^5)
ylabel('Z-g (\mug)')
box on; grid on
subplot(414)
plot(dataDec100.t,(dataDec100.as-9.815)*10^5,'linewidth',1)
datetick('x',3)
ylim(([9.81 9.82]-9.815)*10^5)
ylabel('A_t_o_t-g (\mug)')
box on; grid on
xline(datenum(2019,08,13),'k','linewidth',2);
xline(datenum(2019,11,01),'k','linewidth',0.5);
xline(datenum(2019,12,01),'k','linewidth',0.5);
xline(datenum(2020,01,08),'k','linewidth',0.5);
xline(datenum(2020,01,15),'k','linewidth',0.5);
xline(datenum(2020,01,22),'k','linewidth',0.5);
xline(datenum(2020,01,29),'k','linewidth',0.5);
xline(datenum(2020,02,19),'k','linewidth',0.5);
xline(datenum(2020,03,18),'k','linewidth',0.5);
xline(datenum(2020,04,29),'k','linewidth',0.5);
xline(datenum(2020,06,03),'k','linewidth',0.5);
xline(datenum(2020,01,17),'k','linewidth',0.5);
xline(datenum(2020,09,02),'k','linewidth',0.5);
xline(datenum(2020,09,09),'k','linewidth',0.5);

% Axial new location
load('../calibrations/Axial/axialdata_newloc.mat','dataDec100')

figure(2); clf
subplot(411); hold on
plot(dataDec100.t,(dataDec100.a(:,1)-0.0275)*10^5,'linewidth',1)
datetick('x',3)
ylim(([0.02 0.035]-0.0275)*10^5)
ylabel('X (\mug)')
box on; grid on
subplot(412)
plot(dataDec100.t,(dataDec100.a(:,2)+0.143)*10^5,'linewidth',1)
datetick('x',3)
ylim(([-0.1401 -0.1395]+0.14)*10^5)
ylabel('Y (\mug)')
box on; grid on
subplot(413)
plot(dataDec100.t,(dataDec100.a(:,3)-9.8125)*10^5,'linewidth',1)
datetick('x',3)
ylim(([9.811 9.812]-9.8125)*10^5)
ylabel('Z-g (\mug)')
box on; grid on
subplot(414)
plot(dataDec100.t,(dataDec100.as-9.8125)*10^5,'linewidth',1)
datetick('x',3)
ylim(([9.812 9.813]-9.8125)*10^5)
ylabel('A_t_o_t-g (\mug)')
box on; grid on