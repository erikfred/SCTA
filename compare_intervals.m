% compare_intervals.m
%
% Plot user-selected data intervals and compare in time and frequency
% space. Used to highlight anomalous noise characteristics found in OOI
% SCTA data
%

l=0.5; % length of interval, in days
fs=40; % desired sample rate
if fs==40
    cha={'BNE','BNN','BNZ','BKA'};
elseif fs==8
    cha={'MNE','MNN','MNZ','MKA'};
else
    keyboard
end

t0=[datenum(2019,09,06,12,00,00); datenum(2019,09,22,12,00,00); ...
    datenum(2019,10,15); datenum(2019,11,02,01,00,00); ...
    datenum(2019,12,03,01,00,00); datenum(2020,01,10,03,00,00); ...
    datenum(2020,09,04); datenum(2020,09,12)];

data=[];
for i=1:length(t0)
    % get the data
    chastr={'a(:,1)','a(:,2)','a(:,3)','T'};
    for j=1:length(cha)
        if ~exist(['../tiltcompare/AXCC2/AXCC2_' cha{j} '_' datestr(floor(t0(i)),29) '.miniseed'],'file')
            IRIS_data_pull('AXCC2',cha{j},'--',floor(t0(i)),floor(t0(i))+1)
        end
        temp=rdmseed(['../tiltcompare/AXCC2/AXCC2_' cha{j} '_' datestr(floor(t0(i)),29) '.miniseed']);
        tempt=cat(1,temp.t);
        tempd=cat(1,temp.d)/10^7;
        % isolate window of interest
        eval(['data(i).' chastr{j} '=interp1(tempt,tempd,linspace(t0(i),t0(i)+l,fs*60*60*24*l)'');']);
    end
    data(i).t=linspace(t0(i),t0(i)+l,fs*60*60*24*l)';
    data(i).as=sqrt(data(i).a(:,1).^2+data(i).a(:,2).^2+data(i).a(:,3).^2);    
    
    [data(i).pwr,data(i).fr]=pwelch(data(i).as-mean(data(i).as),40*60*60,40*60*30,2^18,40);
end

%% plotting

% time space
figure(43)
clf; hold on
plot(data(1).t-data(1).t(1),data(1).as-mean(data(1).as),'-k','linewidth',1)
plot(data(4).t-data(4).t(1),data(4).as-mean(data(4).as),'-c','linewidth',1)
plot(data(7).t-data(7).t(1),data(7).as-mean(data(7).as),'-r','linewidth',1)
plot(data(8).t-data(8).t(1),data(8).as-mean(data(8).as),'-b','linewidth',1)
legend('Sep 06 2019','Dec 03 2019','Sep 04 2020','Sep 12 2020','location','northwest')
ylabel('Acceleration - mean (m/s^2)')
datetick('x')
set(gca,'fontsize',14)
box on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../noise exploration/Axial_noise_evolution/timeseries','-dtiff','-r300')

figure(44)
clf; hold on
plot(data(1).t-data(1).t(1),data(1).T,':k','linewidth',2)
plot(data(4).t-data(4).t(1),data(4).T,':c','linewidth',2)
plot(data(7).t-data(7).t(1),data(7).T,':r','linewidth',2)
plot(data(8).t-data(8).t(1),data(8).T,':b','linewidth',2)
legend('Sep 06 2019','Dec 03 2019','Sep 04 2020','Sep 12 2020','location','northwest')
ylabel('Temperature (C)')
datetick('x')
set(gca,'fontsize',14)
box on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../noise exploration/Axial_noise_evolution/temperature','-dtiff','-r300')

%frequency space
figure(45)
clf
loglog(data(1).fr,data(1).pwr,'k','linewidth',1)
hold on
loglog(data(4).fr,data(4).pwr,'c','linewidth',1)
loglog(data(7).fr,data(7).pwr,'r','linewidth',1)
loglog(data(8).fr,data(8).pwr,'b','linewidth',1)
legend('Sep 06 2019','Dec 03 2019','Sep 04 2020','Sep 12 2020','location','northeast')
ylabel('Amplitude')
xlabel('Frequency (Hz)')
set(gca,'fontsize',14)
box on

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../noise exploration/Axial_noise_evolution/spectra','-dtiff','-r300')