% assess_calibrations.m
%
% Allows user to browse through calibration intervals and assess the
% accuracy of the automated data selection used in the calibration.
% Currently only works for calibrations under the 5-flip scheme.
%

load('../calibrations/Axial/axialdata','flipInfoAll');

[~,id,~]=unique(floor(flipInfoAll.t));
daylist=floor(flipInfoAll.t(id));
daylist=daylist(end-11:end);

for i=1:length(daylist)
    dayn=daylist(i);
    data=[];
    flipn=find(floor(flipInfoAll.t)==dayn);
    
    cha={'MNE','MNN','MNZ','MKA'};
    chastr={'a(:,1)','a(:,2)','a(:,3)','T'};
    for m=1:length(cha)
        IRIS_data_pull('AXCC2',cha{m},'--',dayn,dayn+1)
        temp=rdmseed(['../tiltcompare/AXCC2/AXCC2_' cha{m} '_' datestr(dayn,29) '.miniseed']);
        data.t=cat(1,temp.t);
        eval(['data.' chastr{m} '=cat(1,temp.d)/10^7;']);
    end
    data.as=sqrt(data.a(:,1).^2+data.a(:,2).^2+data.a(:,3).^2);
    
    [data1DayDec] = decimate_SCTA(data,1);
    
    figure(i)
    clf
    lbls={'+X','+Y','-Y','+X','-X'};
    for k=1:length(flipn)
        subplot(3,2,k)
        plot(data1DayDec.t(flipInfoAll.i0(flipn(k)):flipInfoAll.i1(flipn(k))),...
            data1DayDec.as(flipInfoAll.i0(flipn(k)):flipInfoAll.i1(flipn(k))),'k')
        hold on
        plot(data1DayDec.t(flipInfoAll.i0g(flipn(k)):flipInfoAll.i1g(flipn(k))),...
            data1DayDec.as(flipInfoAll.i0g(flipn(k)):flipInfoAll.i1g(flipn(k))),...
            'r','linewidth',2)
        ylim([flipInfoAll.gCal(flipn(k))-10e-5 flipInfoAll.gCal(flipn(k))+10e-5])
        a=xlim;
        plot([a(1) a(2)],[flipInfoAll.gCal(flipn(k)) flipInfoAll.gCal(flipn(k))],'b')
        datetick('x','keeplimits')
        title(lbls{k})
        set(gca,'fontsize',12)
    end
    keyboard
end