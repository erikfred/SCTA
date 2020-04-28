doWhat = [false true false];
getData = [false true false];


%% Pinon Flat Acceleration
if doWhat(1)
    startDate = datenum('Nov 20, 2018');
    endDate = datenum('Nov 27, 2018');
    p1 = [];
    p2a = [];
    p2b = [];
    
    
    if getData(1)
        for testDate = startDate:endDate-1
            scta = get_sctaDay('/Volumes/GoogleDrive/My Drive/Oceanography/SCTA-Share/OOI-PF/SCTA-PF (1)/ParsedData',testDate);
            
            for i=1:3
                scta.a(:,i) = detrend(scta.a(:,i));
            end
            
            [ptemp,f1] = pwelch(scta.a,1024*32,1024*16,1024*32,40);
            if isempty(p1)
                p1 = ptemp;
            else
                p1 = p1+ptemp;
            end
            
            % Station PFO , network II, locatiion 00, BH1, BH2, BHZ - Streckeisen STS-1 Seismometer 20 Hz
            % trace = irisFetch.Traces(p.network.network,stationString,p.network.location,channelString,str1,str2);
            str1 = datestr(testDate, 'yyyy-mm-dd HH:MM:SS.FFF');
            str2 = datestr(testDate+1-1e-10, 'yyyy-mm-dd HH:MM:SS.FFF');
            %   trace = irisFetch.Traces('II','PFO','00','BHZ,BH1,BH2',str1,str2,'includePZ'); %'verbose');
            tracea = irisFetch.Traces('II','PFO','00','BHZ,BH1,BH2',str1,str2,'includePZ'); %'verbose');
            pfo = [];
            for i=1:3
                pfo = [pfo tracea(i).data];
            end
            [ptemp,f2a] = pwelch(pfo,1024*16,1024*8,1024*16,20);
            i = (f2a>0 & f2a<8);
            f2a = f2a(i);
            ptemp = ptemp(i,:);
            if isempty(p2a)
                p2a = ptemp;
            else
                p2a = p2a+ptemp;
            end
            
            % Station PFO , network II, locatiion 10, BH1, BH2, BHZ - Nanometrics Trillium 240 Seismometer 40 Hz
            %   trace = irisFetch.Traces('II','PFO','10','BHZ,BH1,BH2',str1,str2,'includePZ'); %'verbose');
            traceb = irisFetch.Traces('II','PFO','10','BHZ,BH1,BH2',str1,str2,'includePZ'); %'verbose');
            pfo = [];
            for i=1:3
                pfo = [pfo traceb(i).data];
            end
            [ptemp,f2b] = pwelch(pfo,1024*32,1024*16,1024*32,40);
            i = (f2b>0 & f2b<18);
            f2b = f2b(i);
            ptemp = ptemp(i,:);
            if isempty(p2b)
                p2b = ptemp;
            else
                p2b = p2b+ptemp;
            end
            
        end
        
        p1 = p1/(endDate-startDate);
        p2a = p2a/(endDate-startDate);
        p2b = p2b/(endDate-startDate);
        
        trans1 = sac_poleZeroResponse(f2a,tracea(1).sacpz,'acc');
        trans2 = sac_poleZeroResponse(f2a,tracea(2).sacpz,'acc');
        trans3 = sac_poleZeroResponse(f2a,tracea(3).sacpz,'acc');
        p2a(:,1) = p2a(:,1)./abs(trans1).^2;
        p2a(:,2) = p2a(:,2)./abs(trans2).^2;
        p2a(:,3) = p2a(:,3)./abs(trans3).^2;
        trans1 = sac_poleZeroResponse(f2b,traceb(1).sacpz,'acc');
        trans2 = sac_poleZeroResponse(f2b,traceb(2).sacpz,'acc');
        trans3 = sac_poleZeroResponse(f2b,traceb(3).sacpz,'acc');
        p2b(:,1) = p2b(:,1)./abs(trans1).^2;
        p2b(:,2) = p2b(:,2)./abs(trans2).^2;
        p2b(:,3) = p2b(:,3)./abs(trans3).^2;
        save ../powerspectra/pf_accel_PowerSpectra f1 f2a f2b p1 p2a p2b
        
    else
        load ../powerspectra/pf_accel_PowerSpectra
        
    end
    
    figure(101)
    clf
    subplot(211)
    h1 = loglog(f1,p1(:,3),'-r');
    hold on
    h2 = loglog(f2a,p2a(:,3),'-b',f2b,p2b(:,3),'-g');
    loglog(f1,p1(:,3),'-r')
    loglog(f2a,p2a(:,3),'-b');
    legend([h1 h2'],'SCTA - Z','PFO STS-1 - Z','PFO T240 - Z')
    xlabel('Frequency, Hz')
    ylabel('Acceleration Spectra (m/s^{2})^{2}/Hz')
    title(['Piñon Flat Observatory ' datestr(startDate,'mmm dd, yyyy') ' to ' datestr(endDate,'mmm dd, yyyy')])
    xlim([0.01 10])
    ylim([1e-17 1e-11])
    subplot(212)
    h1 = loglog(f1,p1(:,1),'-r');
    hold on
    h2 = loglog(f2a,p2a(:,1),'-b',f2b,p2b(:,1),'-g');
    loglog(f2b,p2b(:,2),'-g',f2a,p2a(:,2),'-b');
    loglog(f1,p1(:,1),'-r');
    legend([h1 h2'],'SCTA - X/Y','PFO STS-1 - 1/2','PFO T240 - 1/2')
    xlabel('Frequency, Hz')
    ylabel('Acceleration Spectra (m/s^{2})^{2}/Hz')
    xlim([0.01 10])
    ylim([1e-17 1e-11])
    print -djpeg ../powerspectra/pf_accelPowerSpectra.jpg
    print -dtiff ../powerspectra/pf_accelPowerSpectra.tiff -r300
end

%% Axial
if doWhat(2)
    startDate = datenum('September 04, 2019');
    endDate = datenum('September 05, 2019');
    p1 = [];
    p2 = [];
    
    
    if getData(2)
        for testDate = startDate:endDate-1
            scta = get_sctaDay('/Volumes/GoogleDrive/My Drive/Oceanography/SCTA-Share/OOI-SCTA/ParsedData',testDate);
            
            for i=1:3
                scta.a(:,i) = detrend(scta.a(:,i));
            end
            
            [ptemp,f1] = pwelch(scta.a,1024*32,1024*16,1024*32,40);
            if isempty(p1)
                p1 = ptemp;
            else
                p1 = p1+ptemp;
            end
            
            % Station AXCC1 , network OO, locatiion --, BHE, BHN, BHZ - Guralp CMG-1T
            str1 = datestr(testDate, 'yyyy-mm-dd HH:MM:SS.FFF');
            str2 = datestr(testDate+1-1e-10, 'yyyy-mm-dd HH:MM:SS.FFF');
            trace = irisFetch.Traces('OO','AXCC1','--','BHE,BHN,BHZ',str1,str2,'includePZ'); %'verbose');
            c1 = [];
            c2 = [];
            c3 = [];
            for i=1:length(trace)
                if strcmp(trace(i).channel,'BHE')
                    c1 = [c1; trace(i).data];
                elseif strcmp(trace(i).channel,'BHN')
                    c2 = [c2; trace(i).data];
                elseif strcmp(trace(i).channel,'BHZ')
                    c3 = [c3; trace(i).data];
                end
            end
            %       axcc1 = [c1 c2 c3];
            [ptemp1,f2] = pwelch(c1,1024*32,1024*16,1024*32,40);
            [ptemp2,f2] = pwelch(c2,1024*32,1024*16,1024*32,40);
            [ptemp3,f2] = pwelch(c3,1024*32,1024*16,1024*32,40);
            ptemp = [ptemp1 ptemp2 ptemp3];
            i = (f2>0 & f2<18);
            f2 = f2(i);
            ptemp = ptemp(i,:);
            if isempty(p2)
                p2 = ptemp;
            else
                p2 = p2+ptemp;
            end
        end
        
        p1 = p1/(endDate-startDate);
        p2 = p2/(endDate-startDate);
        
        trans1 = sac_poleZeroResponse(f2,trace(1).sacpz,'acc');
        trans2 = sac_poleZeroResponse(f2,trace(2).sacpz,'acc');
        trans3 = sac_poleZeroResponse(f2,trace(3).sacpz,'acc');
        p2(:,1) = p2(:,1)./abs(trans1).^2;
        p2(:,2) = p2(:,2)./abs(trans2).^2;
        p2(:,3) = p2(:,3)./abs(trans3).^2;
        save ../powerspectra/axial_accel_PowerSpectra f1 f2 p1 p2
        
    else
        load ../powerspectra/axial_accel_PowerSpectra
        
    end
    
    figure(102)
    clf
    subplot(211)
    loglog(f1,p1(:,3),'-b')
    hold on
    loglog(f2,p2(:,3),'-r')
    legend('SCTA - Z','AXCC1 CMG-1 - Z','location','south')
    xlabel('Frequency, Hz')
    ylabel('Acceleration Spectra (m/s^{2})^{2}/Hz')
    title(['Axial Seamount ' datestr(startDate,'mmm dd, yyyy') ' to ' datestr(endDate,'mmm dd, yyyy')])
    xlim([0.01 10])
    ylim([1e-15 1e-10])
    subplot(212)
    loglog(f1,p1(:,1:2),'-b')
    hold on
    loglog(f2,p2(:,1:2),'-r')
    legend('SCTA - X','SCTA - Y','AXCC1 CMG-1 - N','AXCC1 CMG-1 - E','location','south')
    xlabel('Frequency, Hz')
    ylabel('Acceleration Spectra (m/s^{2})^{2}/Hz')
    xlim([0.01 10])
    ylim([1e-15 1e-10])
    
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 11 8.5];
    print -djpeg ../powerspectra/axial_accelPowerSpectra.jpg
    print -dtiff ../powerspectra/axial_accelPowerSpectra.tiff -r300
end

%% Pinon Flat Tilt
if doWhat(3)
    %   startDate = datenum('Nov 18, 2018');
    %   endDate = datenum('Nov 24, 2018');
    startDate = datenum('Nov 20, 2018');
    endDate = datenum('Nov 24, 2018');
    p1_1 = [];
    p1_2 = [];
    p1Rot = [];
    p1Diff = [];
    p2 = [];
    load pf_tiltCompare
    
    if getData(3)
        for testDate = startDate:endDate-1
            tilt1 = get_tiltDay('/Users/wilcock/Mydrive/APL/SCTA-Share/OOI-PF/SCTA-PF (1)/ParsedData','SCTA-Tilt_20_19_',testDate);
            tilt2 = get_tiltDay('/Users/wilcock/Mydrive/APL/SCTA-Share/OOI-PF/SCTA-PF (1)/ParsedData','SCTA-Tilt_22_25_',testDate);
            
            % Rotate tiltmeter 1 into direction of 2
            tiltRot = tilt1;
            tiltRot.a(:,1) = tilt1.a(:,1)*cos(theta) + tilt1.a(:,2)*sin(theta);
            tiltRot.a(:,2) = -tilt1.a(:,1)*sin(theta) + tilt1.a(:,2)*cos(theta);
            % Difference in tilts
            tiltDiff = tilt1;
            t1 = tiltRot.t(1);
            n1 = length(tiltRot.t);
            t2 = tilt2.t(1);
            n2 = length(tilt2.t);
            if abs(t1-t2)<1e-6
                tiltDiff.a = tilt2.a(1:min(n1,n2),:)-tiltRot.a(1:min(n1,n2),:);
                tiltDiff.t = tiltDiff.t(1:min(n1,n2));
            elseif t1>t2
                i1 = round((t1-t2)*86400);
                tiltDiff.a = tilt2.a(1:min(n1,n2)-i1,:)-tiltRot.a(1+i1:min(n1,n2),:);
                tiltDiff.t = tilt2.t(1:min(n1,n2)-i1);
            elseif t2>t1
                i1 = round((t2-t1)*86400);
                tiltDiff.a = tilt2.a(1+i1:min(n1,n2),:)-tiltRot.a(1:min(n1,n2)-i1,:);
                tiltDiff.t = tilt2.t(1+i1:min(n1,n2));
            end
            
            for i=1:2
                tilt1.a(:,i) = detrend(tilt1.a(:,i));
                tilt2.a(:,i) = detrend(tilt2.a(:,i));
                tiltRot.a(:,i) = detrend(tiltRot.a(:,i));
                tiltDiff.a(:,i) = detrend(tiltDiff.a(:,i));
            end
            
            [ptemp,f1_1] = pwelch(tilt1.a,1024*4,512*4,1024*4,1);
            if isempty(p1_1)
                p1_1 = ptemp;
            else
                p1_1 = p1_1+ptemp;
            end
            [ptemp,f1_2] = pwelch(tilt2.a,1024*4,512*4,1024*4,1);
            if isempty(p1_2)
                p1_2 = ptemp;
            else
                p1_2 = p1_2+ptemp;
            end
            [ptemp,f1Rot] = pwelch(tiltRot.a,1024*4,512*4,1024*4,1);
            if isempty(p1Rot)
                p1Rot = ptemp;
            else
                p1Rot = p1Rot+ptemp;
            end
            [ptemp,f1Diff] = pwelch(tiltDiff.a,1024*4,512*4,1024*4,1);
            if isempty(p1Diff)
                p1Diff = ptemp;
            else
                p1Diff = p1Diff+ptemp;
            end
            
            % Station PFO , network II, locatiion 00, BH1, BH2, BHZ - Streckeisen STS-1 Seismometer 20 Hz
            % trace = irisFetch.Traces(p.network.network,stationString,p.network.location,channelString,str1,str2);
            str1 = datestr(testDate, 'yyyy-mm-dd HH:MM:SS.FFF');
            str2 = datestr(testDate+1-1e-10, 'yyyy-mm-dd HH:MM:SS.FFF');
            %   trace = irisFetch.Traces('II','PFO','00','BH1,BH2',str1,str2,'includePZ'); %'verbose');
            trace = irisFetch.Traces('II','PFO','00','BHZ,BH1,BH2',str1,str2,'includePZ'); %'verbose');
            pfo = [];
            for i=1:2
                pfo = [pfo trace(i).data];
            end
            [ptemp,f2] = pwelch(pfo,1024*16*4,1024*8*4,1024*16*4,20);
            i = (f2>0 & f2<8);
            f2 = f2(i);
            ptemp = ptemp(i,:);
            if isempty(p2)
                p2 = ptemp;
            else
                p2 = p2+ptemp;
            end
            
        end
        
        p1_1 = p1_1/(endDate-startDate);
        p1_2 = p1_2/(endDate-startDate);
        p1Rot = p1Rot/(endDate-startDate);
        p1Diff = p1Diff/(endDate-startDate);
        p2 = p2/(endDate-startDate);
        
        trans1 = sac_poleZeroResponse(f2,trace(1).sacpz,'acc');
        trans2 = sac_poleZeroResponse(f2,trace(2).sacpz,'acc');
        p2(:,1) = p2(:,1)./abs(trans1).^2;
        p2(:,2) = p2(:,2)./abs(trans2).^2;
        save pf_tilt_PowerSpectra f1_1 f1_2 f1Rot f1Diff f2 p1_1 p1_2 p1Rot p1Diff p2
        
    else
        load pf_tilt_PowerSpectra
        
    end
    
    figure(103)
    clf
    subplot(211)
    h1 = loglog(f1_1,p1_1(:,1),'-r',f1_2,p1_2(:,1),'-m');
    hold on
    %   h3 = loglog(f1Diff,p1Diff(:,1),'-g');
    h2 = loglog(f2,p2(:,1),'-b');
    loglog(f1_1,p1_1(:,1),'-r',f1_2,p1_2(:,1),'-m');
    legend([h1' h2],'Tilt 20-19','Tilt 22-25','PFO STS-1','location','north')
    xlabel('Frequency, Hz')
    ylabel('Acceleration Spectra (m/s^{2})^{2}/Hz')
    title(['X - Piñon Flat Observatory ' datestr(startDate,'mmm dd, yyyy') ' to ' datestr(endDate,'mmm dd, yyyy')])
    xlim([0.0005 0.3])
    ylim([0.3e-16 3e-12])
    subplot(212)
    h1 = loglog(f1_1,p1_1(:,2),'-r',f1_2,p1_2(:,2),'-m');
    hold on
    %   h3 = loglog(f1Diff,p1Diff(:,2),'-g');
    h2 = loglog(f2,p2(:,2),'-b');
    loglog(f1_1,p1_1(:,2),'-r',f1_2,p1_2(:,2),'-m');
    legend([h1' h2],'Tilt 20-19','Tilt 22-25','PFO STS-1','location','north')
    xlabel('Frequency, Hz')
    ylabel('Acceleration Spectra (m/s^{2})^{2}/Hz')
    title(['Y - Piñon Flat Observatory ' datestr(startDate,'mmm dd, yyyy') ' to ' datestr(endDate,'mmm dd, yyyy')])
    xlim([0.0005 0.3])
    ylim([0.3e-16 3e-12])
    print -djpeg pf_tiltPowerSpectra.jpg
    print -dtiff pf_tiltPowerSpectra.tiff -r300
    
    
    
end

