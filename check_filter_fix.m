% check_filter_fix.m
%
% Analysis to confirm that corrected 8Hz filter is doing what we expect it
% to. Grabs a day of data before and after the fix (on Apr-25-2020) and
% makes spectrograms and amplitude spectra. Also includes a day from Pinon
% Flat for comparison.
%

% get the data
dayn=datenum(2020,04,23);
data40_pre=[];
cha={'BNE','BNN','BNZ','BKA'};
chastr={'a(:,1)','a(:,2)','a(:,3)','T'};
for m=1:length(cha)
    IRIS_data_pull('AXCC2',cha{m},'--',dayn,dayn+1)
    temp=rdmseed(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayn,29) '.miniseed']);
    data40_pre.t=cat(1,temp.t);
    eval(['data40_pre.' chastr{m} '=cat(1,temp.d)/10^7;']);
end
data40_pre.as=sqrt(data40_pre.a(:,1).^2+data40_pre.a(:,2).^2+data40_pre.a(:,3).^2);

dayn=datenum(2020,04,23);
data8_pre=[];
cha={'MNE','MNN','MNZ','MKA'};
chastr={'a(:,1)','a(:,2)','a(:,3)','T'};
for m=1:length(cha)
    IRIS_data_pull('AXCC2',cha{m},'--',dayn,dayn+1)
    temp=rdmseed(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayn,29) '.miniseed']);
    data8_pre.t=cat(1,temp.t);
    eval(['data8_pre.' chastr{m} '=cat(1,temp.d)/10^7;']);
end
data8_pre.as=sqrt(data8_pre.a(:,1).^2+data8_pre.a(:,2).^2+data8_pre.a(:,3).^2);

dayn=datenum(2020,04,26);
data40_post=[];
cha={'BNE','BNN','BNZ','BKA'};
chastr={'a(:,1)','a(:,2)','a(:,3)','T'};
for m=1:length(cha)
    IRIS_data_pull('AXCC2',cha{m},'--',dayn,dayn+1)
    temp=rdmseed(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayn,29) '.miniseed']);
    data40_post.t=cat(1,temp.t);
    eval(['data40_post.' chastr{m} '=cat(1,temp.d)/10^7;']);
end
data40_post.as=sqrt(data40_post.a(:,1).^2+data40_post.a(:,2).^2+data40_post.a(:,3).^2);

dayn=datenum(2020,04,26);
data8_post=[];
cha={'MNE','MNN','MNZ','MKA'};
chastr={'a(:,1)','a(:,2)','a(:,3)','T'};
for m=1:length(cha)
    IRIS_data_pull('AXCC2',cha{m},'--',dayn,dayn+1)
    temp=rdmseed(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayn,29) '.miniseed']);
    data8_post.t=cat(1,temp.t);
    eval(['data8_post.' chastr{m} '=cat(1,temp.d)/10^7;']);
end
data8_post.as=sqrt(data8_post.a(:,1).^2+data8_post.a(:,2).^2+data8_post.a(:,3).^2);

dayn=datenum(2019,12,29);
dataPF = get_sctaDay('/Volumes/GoogleDrive/My Drive/Oceanography/SCTA-Share/OOI-PF/SCTA-PF/ParsedData',dayn);

% grab 4hr of E channels and make spectrograms
e40_pre=data40_pre.a(1:576001,1); e40_pre=e40_pre-mean(e40_pre);
e8_pre=data8_pre.a(1:115201,1); e8_pre=e8_pre-mean(e8_pre);
e40_post=data40_post.a(1:576001,1); e40_post=e40_post-mean(e40_post);
e8_post=data8_post.a(1:115201,1); e8_post=e8_post-mean(e8_post);
e_PF=dataPF.a(1:576001,1); e_PF=e_PF-mean(e_PF);

figure(1)
spectrogram(e40_pre,400*3,200*3,2^10,40,'yaxis')
title('BNE (40Hz) 23-April-2020')
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../noise exploration/filter_fix/pre40','-dtiff','-r300')
print('../noise exploration/filter_fix/pre40','-djpeg')

figure(2)
spectrogram(e8_pre,80*3,40*3,2^10,8,'yaxis')
title('MNE (8Hz) 23-April-2020')
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../noise exploration/filter_fix/pre8','-dtiff','-r300')
print('../noise exploration/filter_fix/pre8','-djpeg')

figure(3)
spectrogram(e40_post,400*3,200*3,2^10,40,'yaxis')
title('BNE (40Hz) 26-April-2020')
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../noise exploration/filter_fix/post40','-dtiff','-r300')
print('../noise exploration/filter_fix/post40','-djpeg')

figure(4)
spectrogram(e8_post,80*3,40*3,2^10,8,'yaxis')
title('MNE (8Hz) 26-April-2020')
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../noise exploration/filter_fix/post8','-dtiff','-r300')
print('../noise exploration/filter_fix/post8','-djpeg')

figure(5)
spectrogram(e_PF,400*3,200*3,2^10,40,'yaxis')
title('Pinon Flat (40Hz) 29-Dec-2019')
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../noise exploration/filter_fix/PF','-dtiff','-r300')
print('../noise exploration/filter_fix/PF','-djpeg')