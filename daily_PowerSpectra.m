% daily_PowerSpectra.m
%
% Computes power spectra for single day of 40Hz SCTA data, with option to
% compare multiple days.
%


Axial=true;
PF=false;
fs=40;

%% Axial
if Axial
    startDate = datenum('January 10, 2020');
    endDate = startDate+1;
    p1 = [];    
    
    if fs==40
        scta = get_sctaDay('/Volumes/GoogleDrive/My Drive/Oceanography/SCTA-Share/OOI-SCTA/ParsedData',startDate);
        
        if length(scta.t)<3456000
            disp('Incomplete day of data')
            return
        end
        
        for i=1:3
            scta.a(:,i) = detrend(scta.a(:,i));
        end
        [p1,f1] = pwelch(scta.a,1024*32,1024*16,1024*32,40);
        p1 = p1/(endDate-startDate);
    elseif fs==8
        scta=[];
        oneHz=[];
        cha={'MNE','MNN','MNZ','MKA'};
        chastr={'a(:,1)','a(:,2)','a(:,3)','T'};
        for m=1:length(cha)
            IRIS_data_pull('AXCC2',cha{m},'--',startDate,endDate)
            temp=rdmseed(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(startDate,29) '.miniseed']);
            scta.t=cat(1,temp.t);
            eval(['scta.' chastr{m} '=cat(1,temp.d)/10^7;']);
            % decimate to 1 Hz (due to odd filtering at 8Hz)
            eval(['oneHz(:,m)=decimate(scta.' chastr{m} ',8,''fir'');']);
        end
        scta.a=oneHz(:,1:3);
        scta.T=oneHz(:,4);
        scta.t=linspace(scta.t(1),scta.t(end),length(scta.a));
        scta.as=sqrt(scta.a(:,1).^2+scta.a(:,2).^2+scta.a(:,3).^2);
        for i=1:3
            scta.a(:,i) = detrend(scta.a(:,i));
        end
        [p1,f1] = pwelch(scta.a,round(1024*32/40),round(1024*16/40),round(1024*32/40),1);
        p1 = p1/(endDate-startDate);
    else
        disp('Invalid frequency selection')
        return
    end
    
    figure(102)
    clf
    subplot(311)
    loglog(f1,p1(:,3),'linewidth',1)
    hold on
    legend([datestr(startDate,6) ' - Z'],'location','northeast')
    xlabel('Frequency, Hz')
    ylabel('Acceleration Spectra (m/s^{2})^{2}/Hz')
%     title(['Axial Seamount ' datestr(startDate,'mmm dd, yyyy') ' to ' datestr(endDate,'mmm dd, yyyy')])
    xlim([0.0001 10])
    ylim([10^-12 10^-10])
    set(gca,'fontsize',12)
    subplot(312)
    loglog(f1,p1(:,1),'linewidth',1)
    hold on
    legend([datestr(startDate,6) ' - X'],'location','northeast')
    xlabel('Frequency, Hz')
    ylabel('Acceleration Spectra (m/s^{2})^{2}/Hz')
    xlim([0.0001 10])
    ylim([10^-12 10^-10])
    set(gca,'fontsize',12)
    subplot(313)
    loglog(f1,p1(:,2),'linewidth',1)
    hold on
    legend([datestr(startDate,6) ' - Y'],'location','northeast')
    xlabel('Frequency, Hz')
    ylabel('Acceleration Spectra (m/s^{2})^{2}/Hz')
    xlim([0.0001 10])
    ylim([10^-12 10^-10])
    set(gca,'fontsize',12)
    
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print -djpeg ../calibrations/Axial/weird_cal/spectra/noisecomp.jpg
    print -dtiff ../calibrations/Axial/weird_cal/spectra/noisecomp.tiff -r300
end
