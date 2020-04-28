% getaweekofdata.m
%
% Reads miniseed files to get a week's worth of data
%

sta='AXCC2';
cha='MNZ';
t0=datenum(2019,11,01);
t0_s=datestr(t0,31); t0_s=t0_s(1:10);

Z=[];
T=[];
for i=1:7
    temp=rdmseed(['../tiltcompare/' sta '/' sta '_' cha '_' t0_s '.miniseed']);
    t1=cat(1,temp.t);
    z1=cat(1,temp.d);
    
    if strcmp(cha,'HNZ')
        z2=decimate(z1,5,'fir');
        t2=decimate(t1,5,'fir');
    else
        z2=z1;
        t2=t1;
    end
    
    Z=cat(1,Z,z2);
    T=cat(1,T,t2);
    
    t0=t0+1;
    t0_s=datestr(t0,31); t0_s=t0_s(1:10);
end
