% quantify_intervals.m
%
% Script called by 'all_deflation_events.m' to quantify the tilt and
% pressure changes associated with the deflation events at Axial
%

%%%%% tilts and long-term pressure %%%%%
%% AXCC1
% long-term tilt
[~,AXCC1.lt_int(1)]=min(abs(xt_lt(1)-AXCC1.t_time));[~,AXCC1.lt_int(2)]=min(abs(xt_lt(2)-AXCC1.t_time));
temp=polyfit([AXCC1.lt_int(1):AXCC1.lt_int(2)]',AXCC1.LAX(AXCC1.lt_int(1):AXCC1.lt_int(2)),1);temp2=polyval(temp,AXCC1.lt_int(1):AXCC1.lt_int(2));
AXCC1.x_lt=(temp2(end)-temp2(1))/(AXCC1.t_time(AXCC1.lt_int(2))-AXCC1.t_time(AXCC1.lt_int(1)));
temp=polyfit([AXCC1.lt_int(1):AXCC1.lt_int(2)]',AXCC1.LAY(AXCC1.lt_int(1):AXCC1.lt_int(2)),1);temp2=polyval(temp,AXCC1.lt_int(1):AXCC1.lt_int(2));
AXCC1.y_lt=(temp2(end)-temp2(1))/(AXCC1.t_time(AXCC1.lt_int(2))-AXCC1.t_time(AXCC1.lt_int(1)));

% event inflation tilt
[~,AXCC1.u_int(1)]=min(abs(xt_u(1)-AXCC1.t_time));[~,AXCC1.u_int(2)]=min(abs(xt_u(2)-AXCC1.t_time));
temp=polyfit([AXCC1.u_int(1):AXCC1.u_int(2)]',AXCC1.LAX(AXCC1.u_int(1):AXCC1.u_int(2)),1);temp2=polyval(temp,AXCC1.u_int(1):AXCC1.u_int(2));
AXCC1.x_u=(temp2(end)-temp2(1))/(AXCC1.t_time(AXCC1.u_int(2))-AXCC1.t_time(AXCC1.u_int(1)));
temp=polyfit([AXCC1.u_int(1):AXCC1.u_int(2)]',AXCC1.LAY(AXCC1.u_int(1):AXCC1.u_int(2)),1);temp2=polyval(temp,AXCC1.u_int(1):AXCC1.u_int(2));
AXCC1.y_u=(temp2(end)-temp2(1))/(AXCC1.t_time(AXCC1.u_int(2))-AXCC1.t_time(AXCC1.u_int(1)));

% event deflation tilt
[~,AXCC1.d_int(1)]=min(abs(xt_d(1)-AXCC1.t_time));[~,AXCC1.d_int(2)]=min(abs(xt_d(2)-AXCC1.t_time));
temp=polyfit([AXCC1.d_int(1):AXCC1.d_int(2)]',AXCC1.LAX(AXCC1.d_int(1):AXCC1.d_int(2)),1);temp2=polyval(temp,AXCC1.d_int(1):AXCC1.d_int(2));
AXCC1.x_d=(temp2(end)-temp2(1))/(AXCC1.t_time(AXCC1.d_int(2))-AXCC1.t_time(AXCC1.d_int(1)));
temp=polyfit([AXCC1.d_int(1):AXCC1.d_int(2)]',AXCC1.LAY(AXCC1.d_int(1):AXCC1.d_int(2)),1);temp2=polyval(temp,AXCC1.d_int(1):AXCC1.d_int(2));
AXCC1.y_d=(temp2(end)-temp2(1))/(AXCC1.t_time(AXCC1.d_int(2))-AXCC1.t_time(AXCC1.d_int(1)));

% long-term elevation
[~,AXCC1.plt_int(1)]=min(abs(xt_lt(1)-AXCC1.p_time));[~,AXCC1.plt_int(2)]=min(abs(xt_lt(2)-AXCC1.p_time));
temp=polyfit([AXCC1.plt_int(1):AXCC1.plt_int(2)]',AXCC1.BDO(AXCC1.plt_int(1):AXCC1.plt_int(2)),1);temp2=polyval(temp,AXCC1.plt_int(1):AXCC1.plt_int(2));
AXCC1.p_lt=-(temp2(end)-temp2(1))/(AXCC1.p_time(AXCC1.plt_int(2))-AXCC1.p_time(AXCC1.plt_int(1))); % sign change to make elevation

%% AXEC2
% long-term tilt
[~,AXEC2.lt_int(1)]=min(abs(xt_lt(1)-AXEC2.t_time));[~,AXEC2.lt_int(2)]=min(abs(xt_lt(2)-AXEC2.t_time));
temp=polyfit([AXEC2.lt_int(1):AXEC2.lt_int(2)]',AXEC2.LAX(AXEC2.lt_int(1):AXEC2.lt_int(2)),1);temp2=polyval(temp,AXEC2.lt_int(1):AXEC2.lt_int(2));
AXEC2.x_lt=(temp2(end)-temp2(1))/(AXEC2.t_time(AXEC2.lt_int(2))-AXEC2.t_time(AXEC2.lt_int(1)));
temp=polyfit([AXEC2.lt_int(1):AXEC2.lt_int(2)]',AXEC2.LAY(AXEC2.lt_int(1):AXEC2.lt_int(2)),1);temp2=polyval(temp,AXEC2.lt_int(1):AXEC2.lt_int(2));
AXEC2.y_lt=(temp2(end)-temp2(1))/(AXEC2.t_time(AXEC2.lt_int(2))-AXEC2.t_time(AXEC2.lt_int(1)));

% event inflation tilt
[~,AXEC2.u_int(1)]=min(abs(xt_u(1)-AXEC2.t_time));[~,AXEC2.u_int(2)]=min(abs(xt_u(2)-AXEC2.t_time));
temp=polyfit([AXEC2.u_int(1):AXEC2.u_int(2)]',AXEC2.LAX(AXEC2.u_int(1):AXEC2.u_int(2)),1);temp2=polyval(temp,AXEC2.u_int(1):AXEC2.u_int(2));
AXEC2.x_u=(temp2(end)-temp2(1))/(AXEC2.t_time(AXEC2.u_int(2))-AXEC2.t_time(AXEC2.u_int(1)));
temp=polyfit([AXEC2.u_int(1):AXEC2.u_int(2)]',AXEC2.LAY(AXEC2.u_int(1):AXEC2.u_int(2)),1);temp2=polyval(temp,AXEC2.u_int(1):AXEC2.u_int(2));
AXEC2.y_u=(temp2(end)-temp2(1))/(AXEC2.t_time(AXEC2.u_int(2))-AXEC2.t_time(AXEC2.u_int(1)));

% event deflation tilt
[~,AXEC2.d_int(1)]=min(abs(xt_d(1)-AXEC2.t_time));[~,AXEC2.d_int(2)]=min(abs(xt_d(2)-AXEC2.t_time));
temp=polyfit([AXEC2.d_int(1):AXEC2.d_int(2)]',AXEC2.LAX(AXEC2.d_int(1):AXEC2.d_int(2)),1);temp2=polyval(temp,AXEC2.d_int(1):AXEC2.d_int(2));
AXEC2.x_d=(temp2(end)-temp2(1))/(AXEC2.t_time(AXEC2.d_int(2))-AXEC2.t_time(AXEC2.d_int(1)));
temp=polyfit([AXEC2.d_int(1):AXEC2.d_int(2)]',AXEC2.LAY(AXEC2.d_int(1):AXEC2.d_int(2)),1);temp2=polyval(temp,AXEC2.d_int(1):AXEC2.d_int(2));
AXEC2.y_d=(temp2(end)-temp2(1))/(AXEC2.t_time(AXEC2.d_int(2))-AXEC2.t_time(AXEC2.d_int(1)));

% long-term elevation
[~,AXEC2.plt_int(1)]=min(abs(xt_lt(1)-AXEC2.p_time));[~,AXEC2.plt_int(2)]=min(abs(xt_lt(2)-AXEC2.p_time));
temp=polyfit([AXEC2.plt_int(1):AXEC2.plt_int(2)]',AXEC2.BDO(AXEC2.plt_int(1):AXEC2.plt_int(2)),1);temp2=polyval(temp,AXEC2.plt_int(1):AXEC2.plt_int(2));
AXEC2.p_lt=-(temp2(end)-temp2(1))/(AXEC2.p_time(AXEC2.plt_int(2))-AXEC2.p_time(AXEC2.plt_int(1))); % sign change to make elevation

%% AXID1
% long-term tilt
[~,AXID1.lt_int(1)]=min(abs(xt_lt(1)-AXID1.t_time));[~,AXID1.lt_int(2)]=min(abs(xt_lt(2)-AXID1.t_time));
temp=polyfit([AXID1.lt_int(1):AXID1.lt_int(2)]',AXID1.LAX(AXID1.lt_int(1):AXID1.lt_int(2)),1);temp2=polyval(temp,AXID1.lt_int(1):AXID1.lt_int(2));
AXID1.x_lt=(temp2(end)-temp2(1))/(AXID1.t_time(AXID1.lt_int(2))-AXID1.t_time(AXID1.lt_int(1)));
temp=polyfit([AXID1.lt_int(1):AXID1.lt_int(2)]',AXID1.LAY(AXID1.lt_int(1):AXID1.lt_int(2)),1);temp2=polyval(temp,AXID1.lt_int(1):AXID1.lt_int(2));
AXID1.y_lt=(temp2(end)-temp2(1))/(AXID1.t_time(AXID1.lt_int(2))-AXID1.t_time(AXID1.lt_int(1)));

%event inflation tilt
[~,AXID1.u_int(1)]=min(abs(xt_u(1)-AXID1.t_time));[~,AXID1.u_int(2)]=min(abs(xt_u(2)-AXID1.t_time));
temp=polyfit([AXID1.u_int(1):AXID1.u_int(2)]',AXID1.LAX(AXID1.u_int(1):AXID1.u_int(2)),1);temp2=polyval(temp,AXID1.u_int(1):AXID1.u_int(2));
AXID1.x_u=(temp2(end)-temp2(1))/(AXID1.t_time(AXID1.u_int(2))-AXID1.t_time(AXID1.u_int(1)));
temp=polyfit([AXID1.u_int(1):AXID1.u_int(2)]',AXID1.LAY(AXID1.u_int(1):AXID1.u_int(2)),1);temp2=polyval(temp,AXID1.u_int(1):AXID1.u_int(2));
AXID1.y_u=(temp2(end)-temp2(1))/(AXID1.t_time(AXID1.u_int(2))-AXID1.t_time(AXID1.u_int(1)));

% event deflation tilt
[~,AXID1.d_int(1)]=min(abs(xt_d(1)-AXID1.t_time));[~,AXID1.d_int(2)]=min(abs(xt_d(2)-AXID1.t_time));
temp=polyfit([AXID1.d_int(1):AXID1.d_int(2)]',AXID1.LAX(AXID1.d_int(1):AXID1.d_int(2)),1);temp2=polyval(temp,AXID1.d_int(1):AXID1.d_int(2));
AXID1.x_d=(temp2(end)-temp2(1))/(AXID1.t_time(AXID1.d_int(2))-AXID1.t_time(AXID1.d_int(1)));
temp=polyfit([AXID1.d_int(1):AXID1.d_int(2)]',AXID1.LAY(AXID1.d_int(1):AXID1.d_int(2)),1);temp2=polyval(temp,AXID1.d_int(1):AXID1.d_int(2));
AXID1.y_d=(temp2(end)-temp2(1))/(AXID1.t_time(AXID1.d_int(2))-AXID1.t_time(AXID1.d_int(1)));

% long-term elevation
[~,AXID1.plt_int(1)]=min(abs(xt_lt(1)-AXID1.p_time));[~,AXID1.plt_int(2)]=min(abs(xt_lt(2)-AXID1.p_time));
temp=polyfit([AXID1.plt_int(1):AXID1.plt_int(2)]',AXID1.BDO(AXID1.plt_int(1):AXID1.plt_int(2)),1);temp2=polyval(temp,AXID1.plt_int(1):AXID1.plt_int(2));
AXID1.p_lt=-(temp2(end)-temp2(1))/(AXID1.p_time(AXID1.plt_int(2))-AXID1.p_time(AXID1.plt_int(1))); % sign change to make elevation

%% ASHES
if ~exclude_ashes
    % long-term tilt
    [~,ASHES.lt_int(1)]=min(abs(xt_lt(1)-ASHES.t_time));[~,ASHES.lt_int(2)]=min(abs(xt_lt(2)-ASHES.t_time));
    temp=polyfit([ASHES.lt_int(1):ASHES.lt_int(2)]',ASHES.LAX(ASHES.lt_int(1):ASHES.lt_int(2)),1);temp2=polyval(temp,ASHES.lt_int(1):ASHES.lt_int(2));
    ASHES.x_lt=(temp2(end)-temp2(1))/(ASHES.t_time(ASHES.lt_int(2))-ASHES.t_time(ASHES.lt_int(1)));
    temp=polyfit([ASHES.lt_int(1):ASHES.lt_int(2)]',ASHES.LAY(ASHES.lt_int(1):ASHES.lt_int(2)),1);temp2=polyval(temp,ASHES.lt_int(1):ASHES.lt_int(2));
    ASHES.y_lt=(temp2(end)-temp2(1))/(ASHES.t_time(ASHES.lt_int(2))-ASHES.t_time(ASHES.lt_int(1)));
    
    % event inflation tilt
    [~,ASHES.u_int(1)]=min(abs(xt_u(1)-ASHES.t_time));[~,ASHES.u_int(2)]=min(abs(xt_u(2)-ASHES.t_time));
    temp=polyfit([ASHES.u_int(1):ASHES.u_int(2)]',ASHES.LAX(ASHES.u_int(1):ASHES.u_int(2)),1);temp2=polyval(temp,ASHES.u_int(1):ASHES.u_int(2));
    ASHES.x_u=(temp2(end)-temp2(1))/(ASHES.t_time(ASHES.u_int(2))-ASHES.t_time(ASHES.u_int(1)));
    temp=polyfit([ASHES.u_int(1):ASHES.u_int(2)]',ASHES.LAY(ASHES.u_int(1):ASHES.u_int(2)),1);temp2=polyval(temp,ASHES.u_int(1):ASHES.u_int(2));
    ASHES.y_u=(temp2(end)-temp2(1))/(ASHES.t_time(ASHES.u_int(2))-ASHES.t_time(ASHES.u_int(1)));
    
    % event deflation tilt
    [~,ASHES.d_int(1)]=min(abs(xt_d(1)-ASHES.t_time));[~,ASHES.d_int(2)]=min(abs(xt_d(2)-ASHES.t_time));
    temp=polyfit([ASHES.d_int(1):ASHES.d_int(2)]',ASHES.LAX(ASHES.d_int(1):ASHES.d_int(2)),1);temp2=polyval(temp,ASHES.d_int(1):ASHES.d_int(2));
    ASHES.x_d=(temp2(end)-temp2(1))/(ASHES.t_time(ASHES.d_int(2))-ASHES.t_time(ASHES.d_int(1)));
    temp=polyfit([ASHES.d_int(1):ASHES.d_int(2)]',ASHES.LAY(ASHES.d_int(1):ASHES.d_int(2)),1);temp2=polyval(temp,ASHES.d_int(1):ASHES.d_int(2));
    ASHES.y_d=(temp2(end)-temp2(1))/(ASHES.t_time(ASHES.d_int(2))-ASHES.t_time(ASHES.d_int(1)));
    
    % long-term elevation
    [~,ASHES.plt_int(1)]=min(abs(xt_lt(1)-ASHES.p_time));[~,ASHES.plt_int(2)]=min(abs(xt_lt(2)-ASHES.p_time));
    temp=polyfit([ASHES.plt_int(1):ASHES.plt_int(2)]',ASHES.BDO(ASHES.plt_int(1):ASHES.plt_int(2)),1);temp2=polyval(temp,ASHES.plt_int(1):ASHES.plt_int(2));
    ASHES.p_lt=-(temp2(end)-temp2(1))/(ASHES.p_time(ASHES.plt_int(2))-ASHES.p_time(ASHES.plt_int(1))); % sign change to make elevation
end
%% 
%%%%% pressures relative to AXEC2 %%%%%

% difference inflation relative to AXEC2, add sign reversal (decreased pressure = uplift)
AXCC1.p_dif=AXEC2.BDO-interp1(AXCC1.p_time,AXCC1.BDO,AXEC2.p_time);
AXEC2.p_dif=AXEC2.BDO*0;
AXID1.p_dif=AXEC2.BDO-interp1(AXID1.p_time,AXID1.BDO,AXEC2.p_time);
if ~exclude_ashes
    ASHES.p_dif=AXEC2.BDO-interp1(ASHES.p_time,ASHES.BDO,AXEC2.p_time);
end

% % cut out high-amplitude signals/blips
% AXCC1.p_dif(abs(AXCC1.p_dif-nanmedian(AXCC1.p_dif))>4)=nan;
% AXEC2.p_dif(abs(AXEC2.p_dif-median(AXEC2.p_dif))>4)=nan;
% AXID1.p_dif(abs(AXID1.p_dif-nanmedian(AXID1.p_dif))>4)=nan;
% if ~exclude_ashes
%     ASHES.p_dif(abs(ASHES.p_dif-nanmedian(ASHES.p_dif))>4)=nan;
% end

% 12-hr median filter to smooth remaining blips
AXCC1.p_dif=medfilt1(AXCC1.p_dif,60*12,'truncate');
AXEC2.p_dif=medfilt1(AXEC2.p_dif,60*12,'truncate');
AXID1.p_dif=medfilt1(AXID1.p_dif,60*12,'truncate');
if ~exclude_ashes
    ASHES.p_dif=medfilt1(ASHES.p_dif,60*12,'truncate');
end

% relative elevation changes
% AXCC1
[~,AXCC1.pu_int(1)]=min(abs(xt_u(1)-AXEC2.p_time));[~,AXCC1.pu_int(2)]=min(abs(xt_u(2)-AXEC2.p_time));
int=AXCC1.pu_int(1):AXCC1.pu_int(2); inan=~isnan(AXCC1.p_dif(int)); int=int(inan);
temp=polyfit(int',AXCC1.p_dif(int),1);temp2=polyval(temp,int);
AXCC1.p_u=(temp2(end)-temp2(1))/(AXEC2.p_time(AXCC1.pu_int(2))-AXEC2.p_time(AXCC1.pu_int(1)));

[~,AXCC1.pd_int(1)]=min(abs(xt_d(1)-AXEC2.p_time));[~,AXCC1.pd_int(2)]=min(abs(xt_d(2)-AXEC2.p_time));
int=AXCC1.pd_int(1):AXCC1.pd_int(2); inan=~isnan(AXCC1.p_dif(int)); int=int(inan);
temp=polyfit(int',AXCC1.p_dif(int),1);temp2=polyval(temp,int);
AXCC1.p_d=(temp2(end)-temp2(1))/(AXEC2.p_time(AXCC1.pd_int(2))-AXEC2.p_time(AXCC1.pd_int(1)));

% AXEC2
[~,AXEC2.pu_int(1)]=min(abs(xt_u(1)-AXEC2.p_time));[~,AXEC2.pu_int(2)]=min(abs(xt_u(2)-AXEC2.p_time));
int=AXEC2.pu_int(1):AXEC2.pu_int(2); inan=~isnan(AXEC2.p_dif(int)); int=int(inan);
temp=polyfit(int',AXEC2.p_dif(int),1);temp2=polyval(temp,int);
AXEC2.p_u=(temp2(end)-temp2(1))/(AXEC2.p_time(AXEC2.pu_int(2))-AXEC2.p_time(AXEC2.pu_int(1)));

[~,AXEC2.pd_int(1)]=min(abs(xt_d(1)-AXEC2.p_time));[~,AXEC2.pd_int(2)]=min(abs(xt_d(2)-AXEC2.p_time));
int=AXEC2.pd_int(1):AXEC2.pd_int(2); inan=~isnan(AXEC2.p_dif(int)); int=int(inan);
temp=polyfit(int',AXEC2.p_dif(int),1);temp2=polyval(temp,int);
AXEC2.p_d=(temp2(end)-temp2(1))/(AXEC2.p_time(AXEC2.pd_int(2))-AXEC2.p_time(AXEC2.pd_int(1)));

% AXID1
[~,AXID1.pu_int(1)]=min(abs(xt_u(1)-AXEC2.p_time));[~,AXID1.pu_int(2)]=min(abs(xt_u(2)-AXEC2.p_time));
int=AXID1.pu_int(1):AXID1.pu_int(2); inan=~isnan(AXID1.p_dif(int)); int=int(inan);
temp=polyfit(int',AXID1.p_dif(int),1);temp2=polyval(temp,int);
AXID1.p_u=(temp2(end)-temp2(1))/(AXEC2.p_time(AXID1.pu_int(2))-AXEC2.p_time(AXID1.pu_int(1)));

[~,AXID1.pd_int(1)]=min(abs(xt_d(1)-AXEC2.p_time));[~,AXID1.pd_int(2)]=min(abs(xt_d(2)-AXEC2.p_time));
int=AXID1.pd_int(1):AXID1.pd_int(2); inan=~isnan(AXID1.p_dif(int)); int=int(inan);
temp=polyfit(int',AXID1.p_dif(int),1);temp2=polyval(temp,int);
AXID1.p_d=(temp2(end)-temp2(1))/(AXEC2.p_time(AXID1.pd_int(2))-AXEC2.p_time(AXID1.pd_int(1)));

if ~exclude_ashes
    % ASHES
    [~,ASHES.pu_int(1)]=min(abs(xt_u(1)-AXEC2.p_time));[~,ASHES.pu_int(2)]=min(abs(xt_u(2)-AXEC2.p_time));
    int=ASHES.pu_int(1):ASHES.pu_int(2); inan=~isnan(ASHES.p_dif(int)); int=int(inan);
    temp=polyfit(int',ASHES.p_dif(int),1);temp2=polyval(temp,int);
    ASHES.p_u=(temp2(end)-temp2(1))/(AXEC2.p_time(ASHES.pu_int(2))-AXEC2.p_time(ASHES.pu_int(1)));
    
    [~,ASHES.pd_int(1)]=min(abs(xt_d(1)-AXEC2.p_time));[~,ASHES.pd_int(2)]=min(abs(xt_d(2)-AXEC2.p_time));
    int=ASHES.pd_int(1):ASHES.pd_int(2); inan=~isnan(ASHES.p_dif(int)); int=int(inan);
    temp=polyfit(int',ASHES.p_dif(int),1);temp2=polyval(temp,int);
    ASHES.p_d=(temp2(end)-temp2(1))/(AXEC2.p_time(ASHES.pd_int(2))-AXEC2.p_time(ASHES.pd_int(1)));
end

%% compass corrections
%compass corrections
crd_rot=[cosd(CCMP(4).plus_y) sind(CCMP(4).plus_y); -sind(CCMP(4).plus_y) cosd(CCMP(4).plus_y)];
temp3=crd_rot*[AXCC1.x_lt;AXCC1.y_lt]; AXCC1.x_lt_rot=temp3(1); AXCC1.y_lt_rot=temp3(2);
temp4=crd_rot*[AXCC1.x_u;AXCC1.y_u]; AXCC1.x_u_rot=temp4(1); AXCC1.y_u_rot=temp4(2);
temp5=crd_rot*[AXCC1.x_d;AXCC1.y_d]; AXCC1.x_d_rot=temp5(1); AXCC1.y_d_rot=temp5(2);
temp6=crd_rot*[AXCC1.LAX';AXCC1.LAY']; AXCC1.LAX_rot=temp6(1,:); AXCC1.LAY_rot=temp6(2,:);

crd_rot=[cosd(CCMP(3).plus_y) sind(CCMP(3).plus_y); -sind(CCMP(3).plus_y) cosd(CCMP(3).plus_y)];
temp3=crd_rot*[AXEC2.x_lt;AXEC2.y_lt]; AXEC2.x_lt_rot=temp3(1); AXEC2.y_lt_rot=temp3(2);
temp4=crd_rot*[AXEC2.x_u;AXEC2.y_u]; AXEC2.x_u_rot=temp4(1); AXEC2.y_u_rot=temp4(2);
temp5=crd_rot*[AXEC2.x_d;AXEC2.y_d]; AXEC2.x_d_rot=temp5(1); AXEC2.y_d_rot=temp5(2);
temp6=crd_rot*[AXEC2.LAX';AXEC2.LAY']; AXEC2.LAX_rot=temp6(1,:); AXEC2.LAY_rot=temp6(2,:);

crd_rot=[cosd(CCMP(2).plus_y) sind(CCMP(2).plus_y); -sind(CCMP(2).plus_y) cosd(CCMP(2).plus_y)];
temp3=crd_rot*[AXID1.x_lt;AXID1.y_lt]; AXID1.x_lt_rot=temp3(1); AXID1.y_lt_rot=temp3(2);
temp4=crd_rot*[AXID1.x_u;AXID1.y_u]; AXID1.x_u_rot=temp4(1); AXID1.y_u_rot=temp4(2);
temp5=crd_rot*[AXID1.x_d;AXID1.y_d]; AXID1.x_d_rot=temp5(1); AXID1.y_d_rot=temp5(2);
temp6=crd_rot*[AXID1.LAX';AXID1.LAY']; AXID1.LAX_rot=temp6(1,:); AXID1.LAY_rot=temp6(2,:);

if ~exclude_ashes
    crd_rot=[cosd(CCMP(1).plus_y) sind(CCMP(1).plus_y); -sind(CCMP(1).plus_y) cosd(CCMP(1).plus_y)];
    temp3=crd_rot*[ASHES.x_lt;ASHES.y_lt]; ASHES.x_lt_rot=temp3(1); ASHES.y_lt_rot=temp3(2);
    temp4=crd_rot*[ASHES.x_u;ASHES.y_u]; ASHES.x_u_rot=temp4(1); ASHES.y_u_rot=temp4(2);
    temp5=crd_rot*[ASHES.x_d;ASHES.y_d]; ASHES.x_d_rot=temp5(1); ASHES.y_d_rot=temp5(2);
    temp6=crd_rot*[ASHES.LAX';ASHES.LAY']; ASHES.LAX_rot=temp6(1,:); ASHES.LAY_rot=temp6(2,:);
end

%% save
save(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/AXEC2qc'],'AXEC2')
save(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/AXCC1qc'],'AXCC1')
save(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/AXID1qc'],'AXID1')
if ~exclude_ashes
    save(['../tiltcompare/inflation_reversal/in_progress_matfiles/' fldr '/ASHESqc'],'ASHES')
end