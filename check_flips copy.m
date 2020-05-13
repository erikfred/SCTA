
% Simple script to look at all SCTA flips in x, y,
% and z to help identify anomalies that may have lead to bad calibrations.

clear; close all

load('../calibrations/Axial/axialdata','flipInfoAll');
load('../calibrations/Axial/detailed_flipInfo');

multiplots=true;

[~,id,~]=unique(floor(flipInfoAll.t));
daylist=floor(flipInfoAll.t(id));
% daylist=daylist(end-4:end); % just a few of the most recent
% daylist=daylist(26:end); % eveything since start of 5-position calibration
keyboard

if ~multiplots
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
    % identify bad calibrations
    bad_x1=[datenum(2019,12,01) datenum(2020,01,08) datenum(2020,01,15) datenum(2020,01,22) datenum(2020,01,29)];
    bad_y=[datenum(2019,12,01) datenum(2020,01,08)];
    bad_negy=[datenum(2019,12,01) datenum(2020,01,08) datenum(2020,01,15) datenum(2020,01,22) datenum(2020,02,19) datenum(2020,03,18)];
    bad_x2=[datenum(2019,12,01) datenum(2020,01,08) datenum(2020,01,15) datenum(2020,01,22) datenum(2020,04,29)];
    bad_negx=[datenum(2019,10,15) datenum(2019,11,01) datenum(2019,12,01) datenum(2020,01,08) datenum(2020,01,22)];
    bads=cat(2,bad_x1,bad_y,bad_negy,bad_x2,bad_negx); bads=unique(bads);
    
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
    x2_m=px2(1)*(flipInfoSome.t(3:5:end)-flipInfoSome.t(3))+px2(2);
    
    [~,inegx]=setdiff(floor(flipInfoSome.t),bad_negx); inegx=inegx+4;
    pnegx=polyfit(flipInfoSome.t(inegx)-flipInfoSome.t(inegx(1)),flipInfoSome.gCal(inegx),1);
    negx_m=pnegx(1)*(flipInfoSome.t(4:5:end)-flipInfoSome.t(4))+pnegx(2);
    
    for i=1:length(daylist)
        dayn=daylist(i);
        if sum(bads==dayn) % skip the bad cals for now
            continue
        end
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
        % probably worthwhile to add in something that trims data to just
        % the necessary 10 minute window
        
        for j=1:5
            iflips=find(floor(flipInfoSome.t)==dayn);
            
            figure(j)
            subplot(411)
            hold on
            plot(data.t-data.t(1),data.as,'c','linewidth',1);
            xlim([flipInfoSome.t(iflips(j))-data.t(1)+15/60/60/24 flipInfoSome.t(iflips(j))-data.t(1)+(flipInfoSome.duration(iflips(j))-5)/60/60/24])
            ylabel('g_t_o_t (m/s^2)')
            set(gca,'xtick',[])
            title([orientation_strings{j} ' Calibration'])
            box on
            
            subplot(412)
            hold on
            plot(data.t-data.t(1),data.a(:,1),'c','linewidth',1);
            xlim([flipInfoSome.t(iflips(j))-data.t(1)+15/60/60/24 flipInfoSome.t(iflips(j))-data.t(1)+(flipInfoSome.duration(iflips(j))-5)/60/60/24])
            ylabel('x (m/s^2)')
            set(gca,'xtick',[])
            box on
            
            subplot(413)
            hold on
            plot(data.t-data.t(1),data.a(:,2),'c','linewidth',1);
            xlim([flipInfoSome.t(iflips(j))-data.t(1)+15/60/60/24 flipInfoSome.t(iflips(j))-data.t(1)+(flipInfoSome.duration(iflips(j))-5)/60/60/24])
            ylabel('y (m/s^2)')
            set(gca,'xtick',[])
            box on
            
            subplot(414)
            hold on
            plot(data.t-data.t(1),data.a(:,3),'c','linewidth',1);
            xlim([flipInfoSome.t(iflips(j))-data.t(1)+15/60/60/24 flipInfoSome.t(iflips(j))-data.t(1)+(flipInfoSome.duration(iflips(j))-5)/60/60/24])
            ylabel('z (m/s^2)')
            datetick('x','keeplimits')
            box on
        end
    end
    
%     for l=1:5
%         figure(l)
%         fh=gcf;
%         fh.PaperUnits='inches';
%         fh.PaperPosition=[0 0 8.5 11];
%         print(['../calibrations/Axial/weird_cals/multiplots/goodcals_' save_strings{l}],'-dtiff','-r300')
%         saveas(gcf,['../calibrations/Axial/weird_cals/multiplots/goodcals_' save_strings{l} '.fig'])
%     end
    
    for k=1:length(bads)
        dayn=bads(k);
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
        
        ib=find(dayn==floor(flipInfoSome.t(1:5:end)));
        if sum(dayn==bad_x1)
            n=1;
            % calculate theoretical z for good cal
            zg=sqrt(x1_m(ib)^2-flipInfoSome.xCal(n+(ib-1)*5)^2-flipInfoSome.yCal(n+(ib-1)*5)^2);
        elseif sum(dayn==bad_y)
            n=2;
            % calculate theoretical z for good cal
            zg=sqrt(y_m(ib)^2-flipInfoSome.xCal(n+(ib-1)*5)^2-flipInfoSome.yCal(n+(ib-1)*5)^2);
        elseif sum(dayn==bad_negy)
            n=3;
            % calculate theoretical z for good cal
            zg=sqrt(negy_m(ib)^2-flipInfoSome.xCal(n+(ib-1)*5)^2-flipInfoSome.yCal(n+(ib-1)*5)^2);
        elseif sum(dayn==bad_x2)
            n=4;
            % calculate theoretical z for good cal
            zg=sqrt(x2_m(ib)^2-flipInfoSome.xCal(n+(ib-1)*5)^2-flipInfoSome.yCal(n+(ib-1)*5)^2);
        elseif sum(dayn==bad_negx)
            n=5;
            % calculate theoretical z for good cal
            zg=sqrt(negx_m(ib)^2-flipInfoSome.xCal(n+(ib-1)*5)^2-flipInfoSome.yCal(n+(ib-1)*5)^2);
        end
            
        figure(n)
        subplot(411)
        h1=plot(data.t-data.t(1),data.as,'k','linewidth',1);
        h1b=plot([data.t(1) data.t(end)]-data.t(1),[negx_m(ib) negx_m(ib)],'k--','linewidth',1);
        title([datestr(dayn) ' -X Calibration'])
        lim_y=ylim;
        hp1=patch([flipInfoSome.t(n)+60/60/60/24 flipInfoSome.t(n)+90/60/60/24 flipInfoSome.t(n)+90/60/60/24 ...
            flipInfoSome.t(n)+60/60/60/24]-floor(flipInfoSome.t(n)),[lim_y(1) lim_y(1) lim_y(2) lim_y(2)],...
            [0.8 0.8 0.8]);
        hp1.EdgeColor='none';
        hp1.FaceVertexAlphaData=0.2;
        hp1.FaceAlpha='flat';
        ylim(lim_y)
        
        subplot(412)
        h2=plot(data.t-data.t(1),data.a(:,1),'k','linewidth',1);
        lim_y=ylim;
        hp2=patch([flipInfoSome.t(n)+60/60/60/24 flipInfoSome.t(n)+90/60/60/24 flipInfoSome.t(n)+90/60/60/24 ...
            flipInfoSome.t(n)+60/60/60/24]-floor(flipInfoSome.t(n)),[lim_y(1) lim_y(1) lim_y(2) lim_y(2)],...
            [0.8 0.8 0.8]);
        hp2.EdgeColor='none';
        hp2.FaceVertexAlphaData=0.2;
        hp2.FaceAlpha='flat';
        ylim(lim_y)
        
        subplot(413)
        h3=plot(data.t-data.t(1),data.a(:,2),'k','linewidth',1);
        lim_y=ylim;
        hp3=patch([flipInfoSome.t(n)+60/60/60/24 flipInfoSome.t(n)+90/60/60/24 flipInfoSome.t(n)+90/60/60/24 ...
            flipInfoSome.t(n)+60/60/60/24]-floor(flipInfoSome.t(n)),[lim_y(1) lim_y(1) lim_y(2) lim_y(2)],...
            [0.8 0.8 0.8]);
        hp3.EdgeColor='none';
        hp3.FaceVertexAlphaData=0.2;
        hp3.FaceAlpha='flat';
        ylim(lim_y)
        
        subplot(414)
        h4=plot(data.t-data.t(1),data.a(:,3),'k','linewidth',1);
        h4b=plot([data.t(1) data.t(end)]-data.t(1),[zg zg],'k--','linewidth',1);
        lim_y=ylim;
        hp4=patch([flipInfoSome.t(n)+60/60/60/24 flipInfoSome.t(n)+90/60/60/24 flipInfoSome.t(n)+90/60/60/24 ...
            flipInfoSome.t(n)+60/60/60/24]-floor(flipInfoSome.t(n)),[lim_y(1) lim_y(1) lim_y(2) lim_y(2)],...
            [0.8 0.8 0.8]);
        hp4.EdgeColor='none';
        hp4.FaceVertexAlphaData=0.2;
        hp4.FaceAlpha='flat';
        ylim(lim_y)
        
        fh=gcf;
        fh.PaperUnits='inches';
        fh.PaperPosition=[0 0 8.5 11];
        print(['../calibrations/Axial/weird_cals/multiplots/' datestr(dayn,29) '_' save_strings{n}],'-dtiff','-r300')
        
        delete([h1 h1b h2 h3 h4 h4b hp1 hp2 hp3 hp4])
    end
end
