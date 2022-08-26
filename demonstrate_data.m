% demonstrate_data.m
%
% For each of Axial and PF, generates example time series plots of 1) a full
% calibration sequence and 2) a non-calibration day.
%

%% Pinon Flat 3-flip

load('../calibrations/PinonFlat/PFdata.mat','flipInfoAll')

%--- calibration examples
dayPF=floor(flipInfoAll.t(20));
dataPF=get_sctaDay('/Users/erikfred/Google Drive/My Drive/Oceanography/SCTA-Share/OOI-PF/SCTA-PF/ParsedData',dayPF);

cal_plots(dataPF,3,[datenum(0000,00,00,20,59,00),datenum(0000,00,00,21,08,00)],...
    [datenum(0000,00,00,21,00,20),datenum(0000,00,00,21,02,05)],...
    [datenum(0000,00,00,21,01,20),datenum(0000,00,00,21,01,50)])

% full calibration sequence
figure(10);
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../calibrations/example_sequences/cal_sequence_3_PF','-dtiff','-r300')
print('../calibrations/example_sequences/cal_sequence_3_PF','-depsc','-painters')

% one-orientation interval
figure(11);
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../calibrations/example_sequences/x_cal_3_PF','-dtiff','-r300')
print('../calibrations/example_sequences/x_cal_3_PF','-depsc','-painters')

%% Pinon Flat 5-flip

%--- calibration examples
dayPF=floor(flipInfoAll.t(91));
dataPF=get_sctaDay('/Users/erikfred/Google Drive/My Drive/Oceanography/SCTA-Share/OOI-PF/SCTA-PF/ParsedData',dayPF);

cal_plots(dataPF,5,[datenum(0000,00,00,21,00,00),datenum(0000,00,00,21,12,00)],...
    [datenum(0000,00,00,21,00,50) datenum(0000,00,00,21,02,35)],...
    [datenum(0000,00,00,21,01,50) datenum(0000,00,00,21,02,20)])

% full calibration sequence
figure(10);
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../calibrations/example_sequences/cal_sequence_5_PF','-dtiff','-r300')
print('../calibrations/example_sequences/cal_sequence_5_PF','-depsc','-painters')

% one-orientation interval
figure(11);
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../calibrations/example_sequences/x_cal_5_PF','-dtiff','-r300')
print('../calibrations/example_sequences/x_cal_5_PF','-depsc','-painters')

%% Pinon Flat non-calibration intervals

%--- before data gap
% get a week of data
dayPF=datenum(2018,12,09);
dataPF=get_sctaDay('/Users/erikfred/Google Drive/My Drive/Oceanography/SCTA-Share/OOI-PF/SCTA-PF/ParsedData',dayPF);
for i=1:6
    dataday=get_sctaDay('/Users/erikfred/Google Drive/My Drive/Oceanography/SCTA-Share/OOI-PF/SCTA-PF/ParsedData',dayPF+i);
    dataPF=cell2struct(cellfun(@vertcat,struct2cell(dataPF),struct2cell(dataday),'uni',0),fieldnames(dataPF),1);
end

noncal_plots(dataPF,40)

% time series of week of data
figure(12);
ylim([-2 10]*10^-6)
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../calibrations/example_sequences/example_week_time_PF','-dtiff','-r300')
print('../calibrations/example_sequences/example_week_time_PF','-depsc','-painters')

% spectra of data week
figure(13);
ylim([10^-14 10^-11])
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../calibrations/example_sequences/example_week_spectral_PF','-dtiff','-r300')
print('../calibrations/example_sequences/example_week_spectral_PF','-depsc','-painters')

%--- after data gap
% get a week of data
dayPF=datenum(2019,05,12);
dataPF=get_sctaDay('/Users/erikfred/Google Drive/My Drive/Oceanography/SCTA-Share/OOI-PF/SCTA-PF/ParsedData',dayPF);
for i=1:6
    dataday=get_sctaDay('/Users/erikfred/Google Drive/My Drive/Oceanography/SCTA-Share/OOI-PF/SCTA-PF/ParsedData',dayPF+i);
    dataPF=cell2struct(cellfun(@vertcat,struct2cell(dataPF),struct2cell(dataday),'uni',0),fieldnames(dataPF),1);
end

noncal_plots(dataPF,40)

% time series of week of data
figure(12);
ylim([-2 10]*10^-6)
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../calibrations/example_sequences/example_week2_time_PF','-dtiff','-r300')
print('../calibrations/example_sequences/example_week2_time_PF','-depsc','-painters')

% spectra of data week
figure(13);
ylim([10^-14 10^-11])
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../calibrations/example_sequences/example_week2_spectral_PF','-dtiff','-r300')
print('../calibrations/example_sequences/example_week2_spectral_PF','-depsc','-painters')

%% Axial 3-flip

load('../calibrations/Axial/axialdata.mat','flipInfoAll')

%--- calibration examples
dayA=floor(flipInfoAll.t(20));
dataA=get_sctaDay('/Users/erikfred/Google Drive/My Drive/Oceanography/SCTA-Share/OOI-SCTA/ParsedData',dayA);

cal_plots(dataA,3,[datenum(0000,00,00,20,59,00),datenum(0000,00,00,21,08,00)],...
    [datenum(0000,00,00,21,00,20),datenum(0000,00,00,21,02,05)],...
    [datenum(0000,00,00,21,01,20),datenum(0000,00,00,21,01,50)])

% full calibration sequence
figure(10);
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../calibrations/example_sequences/cal_sequence_3_Axial','-dtiff','-r300')
print('../calibrations/example_sequences/cal_sequence_3_Axial','-depsc','-painters')

% one-orientation interval
figure(11);
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../calibrations/example_sequences/x_cal_3_Axial','-dtiff','-r300')
print('../calibrations/example_sequences/x_cal_3_Axial','-depsc','-painters')

%% Axial 5-flip

%--- calibration examples
dayA=floor(flipInfoAll.t(101));
dataA=get_sctaDay('/Users/erikfred/Google Drive/My Drive/Oceanography/SCTA-Share/OOI-SCTA/ParsedData',dayA);

cal_plots(dataA,5,[datenum(0000,00,00,20,59,30),datenum(0000,00,00,21,11,30)],...
    [datenum(0000,00,00,21,00,10) datenum(0000,00,00,21,01,55)],...
    [datenum(0000,00,00,21,01,10) datenum(0000,00,00,21,01,40)])

% full calibration sequence
figure(10);
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../calibrations/example_sequences/cal_sequence_5_Axial','-dtiff','-r300')
print('../calibrations/example_sequences/cal_sequence_5_Axial','-depsc','-painters')

% one-orientation interval
figure(11);
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../calibrations/example_sequences/x_cal_5_Axial','-dtiff','-r300')
print('../calibrations/example_sequences/x_cal_5_Axial','-depsc','-painters')

%% Axial non-calibration intervals

%--- original location (pre September 09, 2020)
% get a week of data
dayA=datenum(2018,12,11);
dataA=pull_data_Axial(dayA);
for i=1:6
    dataday=pull_data_Axial(dayA+i);
    dataA=cell2struct(cellfun(@vertcat,struct2cell(dataA),struct2cell(dataday),'uni',0),fieldnames(dataA),1);
end

dataA.a=[dataA.MNE,dataA.MNN,dataA.MNZ];

noncal_plots(dataA,8)

% time series of week of data
figure(12);
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../calibrations/example_sequences/example_week_time_Axial','-dtiff','-r300')
print('../calibrations/example_sequences/example_week_time_Axial','-depsc','-painters')

% spectra of data week
figure(13);
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../calibrations/example_sequences/example_week_spectral_Axial','-dtiff','-r300')
print('../calibrations/example_sequences/example_week_spectral_Axial','-depsc','-painters')

%--- relocation (post September 09, 2020)
load('../calibrations/Axial/axialdata_newloc.mat','flipInfoAll')
% get a week of data
dayA=datenum(2021,07,11);
dataA=pull_data_Axial(dayA);
for i=1:6
    dataday=pull_data_Axial(dayA+i);
    dataA=cell2struct(cellfun(@vertcat,struct2cell(dataA),struct2cell(dataday),'uni',0),fieldnames(dataA),1);
end

dataA.a=[dataA.MNE,dataA.MNN,dataA.MNZ];

noncal_plots(dataA,8)

% time series of week of data
figure(12);
ylim([-0.06e-3 0.14e-3])
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../calibrations/example_sequences/example_week2_time_Axial','-dtiff','-r300')
print('../calibrations/example_sequences/example_week2_time_Axial','-depsc','-painters')

% spectra of data week
figure(13);
ylim([10^-15 10^-8])
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../calibrations/example_sequences/example_week2_spectral_Axial','-dtiff','-r300')
print('../calibrations/example_sequences/example_week2_spectral_Axial','-depsc','-painters')

%%%%%%%%%%%%%%%%%%%%%

function cal_plots(data_in,n,box0_lim,box1_lim,box2_lim)

data_day=floor(data_in.t(1));

% full calibration sequence
figure(10); clf
plot(data_in.t,data_in.a(:,1),'b','linewidth',1)
hold on
plot(data_in.t,data_in.a(:,2),'r','linewidth',1)
xlim(data_day+[box0_lim(1) box0_lim(2)])
if n==3
    ylim([-0.1 10]); lim_y=ylim;
elseif n==5
    ylim([-10 10]); lim_y=ylim;
end
hp1=patch(data_day+[box1_lim(1) box1_lim(2) box1_lim(2) box1_lim(1)],...
    [lim_y(1) lim_y(1) lim_y(2) lim_y(2)],[0.9 0.9 0.9]);
hp1.EdgeColor='none';
hp1.ZData=-[1;1;1;1];
datetick('x','keeplimits')
xtickangle(45)
legend('X','Y','location','northeast')
ylabel('Acceleration (m/s^2)')
% title(['Example Calibration Sequence at Pinon Flat on ' datestr(data_day)])
set(gca,'fontsize',14)
box on; grid on

% one-orientation interval
figure(11); clf
plot(data_in.t,data_in.as,'k','linewidth',1)
xlim(data_day+[box1_lim(1) box1_lim(2)])
datetick('x','keeplimits')
xtickangle(45)
lim_y=ylim;
hold on
hp2=patch(data_day+[box2_lim(1) box2_lim(2) box2_lim(2) box2_lim(1)],...
    [lim_y(1) lim_y(1) lim_y(2) lim_y(2)],[0.9 0.9 0.9]);
hp2.EdgeColor='none';
hp2.ZData=-[1;1;1;1];
aint=data_in.as(data_in.t>=data_day+box2_lim(1) & data_in.t<data_day+box2_lim(2));
plot(mean([box2_lim(1), box2_lim(2)])+data_day,mean(aint),...
    'x','color',[0.95 0.95 0.95],'linewidth',2,'markersize',10)
ylim(lim_y)
legend('A','location','northeast')
ylabel('Acceleration (m/s^2)')
% title(['Example Calibration +X calibration at Pinon Flat on ' datestr(data_day)])
set(gca,'fontsize',14)
box on; grid on

end

function noncal_plots(data_in,fs)

data_day=floor(data_in.t(1));

% remove EQs, anomalies
data_in.x=filloutliers(detrend(data_in.a(:,1)),'linear','thresholdfactor',6);
data_in.y=filloutliers(detrend(data_in.a(:,2)),'linear','thresholdfactor',6);
data_in.z=filloutliers(detrend(data_in.a(:,3)),'linear','thresholdfactor',6);
% decimate
if fs==40
    d=[4,10,6,10];
elseif fs==8
    d=[4,2,6,10];
end
data_in.xd=decimate(data_in.x,d(1),'fir'); data_in.xd=decimate(data_in.xd,d(2),'fir');
data_in.xd=decimate(data_in.xd,d(3),'fir'); data_in.xd=decimate(data_in.xd,d(4),'fir');
data_in.yd=decimate(data_in.y,d(1),'fir'); data_in.yd=decimate(data_in.yd,d(2),'fir');
data_in.yd=decimate(data_in.yd,d(3),'fir'); data_in.yd=decimate(data_in.yd,d(4),'fir');
data_in.zd=decimate(data_in.z,d(1),'fir'); data_in.zd=decimate(data_in.zd,d(2),'fir');
data_in.zd=decimate(data_in.zd,d(3),'fir'); data_in.zd=decimate(data_in.zd,d(4),'fir');
data_in.td=linspace(data_in.t(1),data_in.t(end),length(data_in.xd));
% box average to 1 sample/min
data_in.xm=movmean(data_in.x,fs*60);
data_in.ym=movmean(data_in.y,fs*60);
data_in.zm=movmean(data_in.z,fs*60);

% time series of day of data
if fs==40
    c=10^-6;
elseif fs==8
    c=2*10^-5;
end
figure(12); clf
plot(data_in.td,data_in.xd+6*c,'b','linewidth',1)
hold on
plot(data_in.td,data_in.yd+3*c,'r','linewidth',1)
plot(data_in.td,data_in.zd,'k','linewidth',1)
xlim([data_day data_day+6])
datetick('x','keeplimits')
xtickangle(45)
legend('X','Y','Z-9.81','location','northeast')
ylabel('Acceleration (m/s^2)')
% title(['Example Day at Pinon Flat on ' datestr(data_day)])
set(gca,'fontsize',14)
box on; grid on

% spectra
if fs==40
    [data_in.pxx,data_in.fx]=pwelch(data_in.x,1024*32,1024*16,1024*32,fs);
    [data_in.pyy,data_in.fy]=pwelch(data_in.y,1024*32,1024*16,1024*32,fs);
    [data_in.pyz,data_in.fz]=pwelch(data_in.z,1024*32,1024*16,1024*32,fs);
elseif fs==8
    [data_in.pxx,data_in.fx]=pwelch(data_in.x,1024*8,1024*4,1024*8,fs);
    [data_in.pyy,data_in.fy]=pwelch(data_in.y,1024*8,1024*4,1024*8,fs);
    [data_in.pyz,data_in.fz]=pwelch(data_in.z,1024*8,1024*4,1024*8,fs);
end

figure(13); clf
loglog(data_in.fx,data_in.pxx,'b','linewidth',1)
hold on
loglog(data_in.fy,data_in.pyy,'r','linewidth',1)
loglog(data_in.fz,data_in.pyz,'k','linewidth',1)
xlim([10^-2 10^1])
legend('X','Y','Z')
ylabel('PSD ((m/s^2)^2/Hz)')
xlabel('Frequency (Hz)')
% title(['Example Day at Pinon Flat on ' datestr(data_day)])
set(gca,'fontsize',14)
box on; grid on

end

function [data_out]=pull_data_Axial(dayn)

t_s=datestr(dayn,31); t_s=t_s(1:10);
MNN_string=['AXCC2_MNN_' t_s '.miniseed'];
MNE_string=['AXCC2_MNE_' t_s '.miniseed'];
MNZ_string=['AXCC2_MNZ_' t_s '.miniseed'];
MKA_string=['AXCC2_MKA_' t_s '.miniseed'];

data_out.t=[]; data_out.MNN=[]; data_out.MNE=[]; data_out.MNZ=[]; data_out.MKA=[];

%attempt download if file not found
if ~exist(['../tiltcompare/AXCC2/' MNN_string],'file') || ...
        ~exist(['../tiltcompare/AXCC2/' MNE_string],'file') || ...
        ~exist(['../tiltcompare/AXCC2/' MNZ_string],'file') || ...
        ~exist(['../tiltcompare/AXCC2/' MKA_string],'file')
    IRIS_data_pull('AXCC2','MNN','--',dayn,dayn+1);
    IRIS_data_pull('AXCC2','MNE','--',dayn,dayn+1);
    IRIS_data_pull('AXCC2','MNZ','--',dayn,dayn+1);
    IRIS_data_pull('AXCC2','MKA','--',dayn,dayn+1);
end

%some dates have no data (power failure, etc.)
if exist(['../tiltcompare/AXCC2/' MNN_string],'file') && ...
        exist(['../tiltcompare/AXCC2/' MNE_string],'file') && ...
        exist(['../tiltcompare/AXCC2/' MNZ_string],'file') && ...
        exist(['../tiltcompare/AXCC2/' MKA_string],'file')
    
    %MNN channel
    temp=rdmseed(['../tiltcompare/AXCC2/' MNN_string]);
    ntemp=double(cat(1,temp.d));
    data_out.MNN=ntemp/10^7;
    
    %MNE channel
    temp=rdmseed(['../tiltcompare/AXCC2/' MNE_string]);
    etemp=double(cat(1,temp.d));
    data_out.MNE=etemp/10^7;
    
    %MNZ channel
    temp=rdmseed(['../tiltcompare/AXCC2/' MNZ_string]);
    ztemp=double(cat(1,temp.d));
    data_out.MNZ=ztemp/10^7;
    
    %MKA channel
    temp=rdmseed(['../tiltcompare/AXCC2/' MKA_string]);
    Ttemp=double(cat(1,temp.d));
    data_out.MKA=Ttemp/10^7;
    
    %time
    ttemp=cat(1,temp.t);
    data_out.t=ttemp;
end
end