function IRIS_data_pull(sta,cha,loc,t_0,t_f)
%
% Queries IRIS database for data in specified range. Saves as daily files.
%   sta = station name as string
%   cha = channel ID as string
%   loc = location code as string
%   t_0 = start time as datenum
%   t_f = end time as datenum
%

inrange=true;
t1=t_0;
t2=t_0+1;
while inrange
    t1_s=datestr(t1,31); t1_s(11)='T';
    t2_s=datestr(t2,31); t2_s(11)='T';
    url=['http://service.iris.edu/irisws/timeseries/1/query?net=OO&sta='...
        sta '&cha=' cha '&start=' t1_s '&end=' t2_s '&format=miniseed&loc=' loc];
    
    savestr=[sta '_' cha '_' t1_s(1:10)];
    if ~exist(['../tiltcompare/' sta '/' savestr '.miniseed'],'file')
        [~,status]=urlread(url);
        if status
            websave(['../tiltcompare/' sta '/' savestr '.miniseed'],url);
        else
            warning(['Unsuccessful download for ' savestr])
        end
    end
        
    t1=t1+1;
    t2=t2+1;
    if t2>t_f
        inrange=false;
    end
end
end