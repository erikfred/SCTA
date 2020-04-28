file_str=['AXEC2_BDO_'];
t1=datenum(2015,09,15);
tf=datenum(2015,11,01);
AXEC2.time=[];
AXEC2.LAX=[];
AXEC2.LAY=[];
AXEC2.BDO=[];
dec=[12;10;10];
while t1<tf
    t1_s=datestr(t1,31); t1_s=t1_s(1:10);
    file_string=[file_str t1_s '.miniseed'];
    %some dates have no data (power failure, etc.)
    if exist(['../tiltcompare/AXEC2/' file_string],'file')
        temp=rdmseed(['../tiltcompare/AXEC2/' file_string]);
        dtemp=double(cat(1,temp.d));
        dtempd=decimate(dtemp,dec(1),'fir');
        dtempd=decimate(dtempd,dec(2),'fir');
        dtempd=decimate(dtempd,dec(3),'fir');
        AXEC2.BDO=[AXEC2.BDO;dtempd];
        
%         ttemp=cat(1,temp.t);
%         ttempd=decimate(ttemp,dec(1),'fir');
%         ttempd=decimate(ttempd,dec(2),'fir');
%         ttempd=decimate(ttempd,dec(3),'fir');
%         AXEC2.time=[AXEC2.time;ttempd];
    end
    t1=t1+1;
end