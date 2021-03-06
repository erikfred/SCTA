
% Simple script to look at all SCTA flips in x, y,
% and z to help identify anomalies that may have lead to bad calibrations.

clear; close all

load('../calibrations/Axial/axialdata','flipInfoAll');
load('../calibrations/Axial/detailed_flipInfo');

multiplots=true;

if ~multiplots
    [~,id,~]=unique(floor(flipInfoAll.t));
    daylist=floor(flipInfoAll.t(id));
    % daylist=daylist(end-4:end); % just a few of the most recent
    daylist=daylist([31,32]); % eveything since start of 5-position calibration

    legstr={};
    figure(22), hold on;
    figure(23), hold on;
    figure(24), hold on;
    figure(25), hold on;
    figure(26), hold on;
    figure(5), hold on;
    figure(6), hold on;
    figure(7), hold on;
    for i=1:length(daylist)
        dayn=daylist(i);
        datestr(dayn)
        data=[];
        cha={'MNE','MNN','MNZ','MKA'};
        chastr={'a(:,1)','a(:,2)','a(:,3)','T'};
        for m=1:length(cha)
            IRIS_data_pull('AXCC2',cha{m},'--',dayn,dayn+1)
            if ~exist(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayn,29) '.miniseed'],'file')
                cha={'BNE','BNN','BNZ','BKA'};
                if ~exist(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayn,29) '.miniseed'],'file')
                    IRIS_data_pull('AXCC2',cha{m},'--',dayn,dayn+1);
                end
            end
            temp=rdmseed(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayn,29) '.miniseed']);
            data.t=cat(1,temp.t);
            eval(['data.' chastr{m} '=cat(1,temp.d)/10^7;']);
        end
        data.as=sqrt(data.a(:,1).^2+data.a(:,2).^2+data.a(:,3).^2);
        
        legstr{i}=datestr(dayn,6);
        %time series
        figure(22)
        plot(data.t-floor(data.t(1)-daylist(1)),data.a(:,1),'linewidth',1)
        ylim([-0.15 0.25])
        title('X channel')
        figure(23)
        plot(data.t-floor(data.t(1)-daylist(1)),data.a(:,2),'linewidth',1)
        title('Y channel')
        figure(24)
        plot(data.t-floor(data.t(1)-daylist(1)),data.a(:,3),'linewidth',1)
        title('Z channel')
        figure(25)
        plot(data.t-floor(data.t(1)-daylist(1)),data.as,'linewidth',1)
        title('Total Acceleration')
        figure(26)
        plot(data.t-floor(data.t(1)-daylist(1)),data.T,'linewidth',1)
        title('Temperature')
        
        %ensure even sampling
        tsmooth=linspace(data.t(1),data.t(end),length(data.t))*24*60*60; %[s]
        asmooth=resample(data.as,data.t);
        %     xsmooth=interp1(data.t,data.a(:,1),tsmooth);
        %     ysmooth=interp1(data.t,data.a(:,2),tsmooth);
        %     zsmooth=interp1(data.t,data.a(:,3),tsmooth);
        xsmooth=resample(data.a(:,1),data.t);
        ysmooth=resample(data.a(:,2),data.t);
        zsmooth=resample(data.a(:,3),data.t);
        
        %spectra
        inti=1:5e5;
        intf=6.2e5:length(tsmooth);
        x1(1)=find(abs(asmooth-mean(asmooth(inti)))>1e-3,1)+110;
        x1(2)=x1(1)+800;
        intx1=x1(1):x1(2);
        x2(1)=find(abs(asmooth-mean(asmooth(intf)))>1e-3,1,'last')-95;
        x2(2)=x2(1)-795;
        intx2=x2(2):x2(1);
        
        figure(5)
        subplot(311,'xscale','log','yscale','log')
        hold on
        [f,a]=amplitudespectrum(tsmooth(inti),xsmooth(inti));
        loglog(f,a,'linewidth',1)
        ylabel({'X',''})
        title('Before Flip')
        subplot(312,'xscale','log','yscale','log')
        hold on
        [f,a]=amplitudespectrum(tsmooth(inti),ysmooth(inti));
        loglog(f,a,'linewidth',1)
        ylabel({'Y','Amplitude'})
        subplot(313,'xscale','log','yscale','log')
        hold on
        [f,a]=amplitudespectrum(tsmooth(inti),zsmooth(inti));
        loglog(f,a,'linewidth',1)
        ylabel({'Z',''})
        
        figure(6)
        subplot(311,'xscale','log','yscale','log')
        hold on
        [f,a]=amplitudespectrum(tsmooth(intx1),xsmooth(intx1));
        loglog(f,a,'linewidth',1)
        ylabel({'X',''})
        title('+X Orientation')
        subplot(312,'xscale','log','yscale','log')
        hold on
        [f,a]=amplitudespectrum(tsmooth(intx1),ysmooth(intx1));
        loglog(f,a,'linewidth',1)
        ylabel({'Y','Amplitude'})
        subplot(313,'xscale','log','yscale','log')
        hold on
        [f,a]=amplitudespectrum(tsmooth(intx1),zsmooth(intx1));
        loglog(f,a,'linewidth',1)
        ylabel({'Z',''})
        
        figure(7)
        subplot(311,'xscale','log','yscale','log')
        hold on
        [f,a]=amplitudespectrum(tsmooth(intx2),xsmooth(intx2));
        loglog(f,a,'linewidth',1)
        ylabel({'X',''})
        title('-X Orientation')
        subplot(312,'xscale','log','yscale','log')
        hold on
        [f,a]=amplitudespectrum(tsmooth(intx2),ysmooth(intx2));
        loglog(f,a,'linewidth',1)
        ylabel({'Y',''})
        subplot(313,'xscale','log','yscale','log')
        hold on
        [f,a]=amplitudespectrum(tsmooth(intx2),zsmooth(intx2));
        loglog(f,a,'linewidth',1)
        ylabel({'Z',''})
    end
    
    %% add legends, etc.
    
    figure(22)
    legend(legstr)
    datetick('x','keeplimits')
    ylabel('Acceleration (m/s^2)')
    figure(23)
    legend(legstr)
    datetick('x','keeplimits')
    ylabel('Acceleration (m/s^2)')
    figure(24)
    legend(legstr)
    datetick('x','keeplimits')
    ylabel('Acceleration (m/s^2)')
    figure(25)
    legend(legstr)
    datetick('x','keeplimits')
    ylabel('Acceleration (m/s^2)')
    figure(26)
    legend(legstr)
    datetick('x','keeplimits')
    ylabel([char(176) 'C'])
    
    figure(5)
    subplot(311),legend(legstr)
    subplot(312),legend(legstr)
    subplot(313),legend(legstr)
    figure(6)
    subplot(311),legend(legstr)
    subplot(312),legend(legstr)
    subplot(313),legend(legstr)
    figure(7)
    subplot(311),legend(legstr)
    subplot(312),legend(legstr)
    subplot(313),legend(legstr)
else
    daylist=unique(floor(flipInfoSome.t)); % every calibration day since 5-orientation scheme
    
    load('../calibrations/Axial/badcaldates') % bad calibrations determined/stored in 'drift_sapndrift.m'

%     % identify bad calibrations (or time series of interest) manually
%     bad_x1=[datenum(2019,11,01) datenum(2019,12,01)];
%     bad_y=[datenum(2019,11,01) datenum(2019,12,01)];
%     bad_negy=[datenum(2019,11,01) datenum(2019,12,01)];
%     bad_x2=[datenum(2019,11,01) datenum(2019,12,01)];
%     bad_negx=[datenum(2019,11,01) datenum(2019,12,01)];
%     bads={bad_x1, bad_y, bad_negy, bad_x2, bad_negx};
%     ubads=cat(2,bad_x1,bad_y,bad_negy,bad_x2,bad_negx); ubads=unique(ubads);
    
    orientation_strings={'+X1','+Y','-Y','+X2','-X'};
    save_strings={'x1','y','negy','x2','negx'};
    
    % find linear best fits to only good calibrations
    [~,ix1]=setdiff(floor(flipInfoSome.t),bad_x1);
    px1=polyfit(flipInfoSome.t(ix1)-flipInfoSome.t(ix1(1)),flipInfoSome.gCal(ix1),1);
    x1_m=px1(1)*(flipInfoSome.t(1:5:end)-flipInfoSome.t(1))+px1(2);
    
    [~,iy]=setdiff(floor(flipInfoSome.t),bad_y); iy=iy+1;
    py=polyfit(flipInfoSome.t(iy)-flipInfoSome.t(iy(1)),flipInfoSome.gCal(iy),1);
    y_m=py(1)*(flipInfoSome.t(2:5:end)-flipInfoSome.t(2))+py(2);
    
    [~,inegy]=setdiff(floor(flipInfoSome.t),bad_negy); inegy=inegy+2;
    pnegy=polyfit(flipInfoSome.t(inegy)-flipInfoSome.t(inegy(1)),flipInfoSome.gCal(inegy),1);
    negy_m=pnegy(1)*(flipInfoSome.t(3:5:end)-flipInfoSome.t(3))+pnegy(2);
    
    [~,ix2]=setdiff(floor(flipInfoSome.t),bad_x2); ix2=ix2+3;
    px2=polyfit(flipInfoSome.t(ix2)-flipInfoSome.t(ix2(1)),flipInfoSome.gCal(ix2),1);
    x2_m=px2(1)*(flipInfoSome.t(4:5:end)-flipInfoSome.t(4))+px2(2);
    
    [~,inegx]=setdiff(floor(flipInfoSome.t),bad_negx); inegx=inegx+4;
    pnegx=polyfit(flipInfoSome.t(inegx)-flipInfoSome.t(inegx(1)),flipInfoSome.gCal(inegx),1);
    negx_m=pnegx(1)*(flipInfoSome.t(5:5:end)-flipInfoSome.t(5))+pnegx(2);
    
    % combine models into single matrix
    mods=[x1_m y_m negy_m x2_m negx_m];
    
    if exist(['../calibrations/Axial/weird_cals/multiplots/good_cals_' save_strings{1} '.fig'],'file')
        for v=1:5
            openfig(['../calibrations/Axial/weird_cals/multiplots/goodcals_' save_strings{v} '.fig']);
        end
    else
        cmap=cmocean('deep'); % smoothly-varying color scheme to indicate time
        cmap=cmap(51:end-50,:); % take off first and last ~20% for better visibility
        for i=1:length(daylist)
            dayn=daylist(i);
            if sum(ubads==dayn) % skip the bad cals for now
                continue
            end
            datestr(dayn)
            icmap=round((dayn-floor(flipInfoSome.t(1)))/(flipInfoSome.t(end)-flipInfoSome.t(1))*length(cmap));
            disp(icmap)
            if icmap==0
                icmap=1;
            end
            
            data=[];
            cha={'MNE','MNN','MNZ'};
            for m=1:length(cha)
                IRIS_data_pull('AXCC2',cha{m},'--',dayn,dayn+1)
                if ~exist(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayn,29) '.miniseed'],'file')
                    cha={'BNE','BNN','BNZ'};
                    if ~exist(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayn,29) '.miniseed'],'file')
                        IRIS_data_pull('AXCC2',cha{m},'--',dayn,dayn+1);
                    end
                end
                temp=rdmseed(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayn,29) '.miniseed']);
                data.t=cat(1,temp.t);
                data.a(:,m)=cat(1,temp.d)/10^7;
            end
            % trim down to only the calibration window
            icut=data.t>dayn+datenum(0,0,0,20,59,00) & data.t<dayn+datenum(0,0,0,21,11,00);
            data.t=data.t(icut);
            data.a=data.a(icut,:);
            temp=sqrt(data.a(:,1).^2+data.a(:,2).^2+data.a(:,3).^2);
            data.a=cat(2,temp,data.a);
            
            iflips=find(floor(flipInfoSome.t)==dayn); % set relative time scale
            label_strings={'Total Accel (m/s^2)','x (m/s^2)','y (m/s^2)','z (m/s^2)'};
            for j=1:5
                figure(j)
                subplot(411)
                title([orientation_strings{j} ' Calibration'])
                for k=1:4
                    subplot(4,1,k)
                    hold on
                    plot(data.t-flipInfoSome.t(iflips(1)),data.a(:,k),'color',cmap(icmap,:),'linewidth',1);
                    xlim([flipInfoSome.t(iflips(j))-flipInfoSome.t(iflips(1))+15/60/60/24 flipInfoSome.t(iflips(j))-flipInfoSome.t(iflips(1))+(flipInfoSome.duration(iflips(j))-5)/60/60/24])
                    ylabel(label_strings{k})
                    set(gca,'xtick',[])
                    box on
                end
                datetick('x','keeplimits')
            end
        end
        
        for l=1:5
            figure(l)
            fh=gcf;
            fh.PaperUnits='inches';
            fh.PaperPosition=[0 0 8.5 11];
            print(['../calibrations/Axial/weird_cals/multiplots/goodcals_' save_strings{l}],'-dtiff','-r300')
            saveas(gcf,['../calibrations/Axial/weird_cals/multiplots/goodcals_' save_strings{l} '.fig'])
        end
    end
    
    %% relative time scale has been changed in the above section; now I must change it below!
    
    % add in curves from bad calibrations, saving individually
    for k=1:length(ubads)
        dayn=ubads(k);
        data=[];
        cha={'MNE','MNN','MNZ'};
        for m=1:length(cha)
            IRIS_data_pull('AXCC2',cha{m},'--',dayn,dayn+1)
            if ~exist(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayn,29) '.miniseed'],'file')
                cha={'BNE','BNN','BNZ'};
                if ~exist(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayn,29) '.miniseed'],'file')
                    IRIS_data_pull('AXCC2',cha{m},'--',dayn,dayn+1);
                end
            end
            temp=rdmseed(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayn,29) '.miniseed']);
            data.t=cat(1,temp.t);
            data.a(:,m)=cat(1,temp.d)/10^7;
        end
        temp=sqrt(data.a(:,1).^2+data.a(:,2).^2+data.a(:,3).^2);
        data.a=cat(2,temp,data.a);
        
        ib=find(dayn==floor(flipInfoSome.t(1:5:end))); % determine which day we're looking at
        iflips=find(floor(flipInfoSome.t)==dayn); % set relative time scale
        for n=1:5
            badcheck=bads{n};
            if sum(dayn==badcheck)
                % calculate theoretical z for good cal
                zg=sqrt(mods(ib,n)^2-flipInfoSome.xCal(n+(ib-1)*5)^2-flipInfoSome.yCal(n+(ib-1)*5)^2);
                if flipInfoSome.zCal(n+(ib-1)*5)<0
                    zg=-zg;
                end
                figure(n)
                for q=1:4
                    subplot(4,1,q)
                    h(q)=plot(data.t-flipInfoSome.t(iflips(1)),data.a(:,q),'k','linewidth',1);
                    if q==1
                        h(9)=plot([data.t(1) data.t(end)]-flipInfoSome.t(iflips(1)),[mods(ib,n) mods(ib,n)],'k--','linewidth',1);
                        title([datestr(flipInfoSome.t(n+(ib-1)*5)) ' ' orientation_strings{n} ' Calibration'])
                    elseif q==4
                        if isreal(zg)
                            h(10)=plot([data.t(1) data.t(end)]-flipInfoSome.t(iflips(1)),[zg zg],'k--','linewidth',1);
                            title('')
                        else
                            title('Corrected z is imaginary')
                        end
                    end
                    ylim auto
                    lim_y=ylim;
                    h(4+q)=patch([flipInfoSome.t(iflips(n))+60/60/60/24 flipInfoSome.t(iflips(n))+90/60/60/24 flipInfoSome.t(iflips(n))+90/60/60/24 ...
                        flipInfoSome.t(iflips(n))+60/60/60/24]-flipInfoSome.t(iflips(1)),[lim_y(1) lim_y(1) lim_y(2) lim_y(2)],...
                        [0.8 0.8 0.8]);
                    h(4+q).EdgeColor='none';
                    h(4+q).FaceVertexAlphaData=0.2;
                    h(4+q).FaceAlpha='flat';
                    ylim([lim_y(1) lim_y(2)])
                    set(gca,'fontsize',14)
                end
                fh=gcf;
                fh.PaperUnits='inches';
                fh.PaperPosition=[0 0 8.5 11];
                print(['../calibrations/Axial/weird_cals/multiplots/' datestr(dayn,29) '_' save_strings{n}],'-dtiff','-r300')
                
                delete(h)
            end
        end
    end
end
