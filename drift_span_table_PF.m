function M = drift_span_table2
%
% Calculate drift rates, spans, and misfits for all locations and channels.
% Enforces a simplistic Temperature-dependence constraint within each
% channel.
%

%% Pinon Flat

load('../calibrations/PinonFlat/PFdata','flipInfoAll')

% identify and separate each of the calibrations
i_x=find(flipInfoAll.orientation==1);
i_x1=i_x(1:2:end); i_x2=i_x(2:2:end);
i_y1=find(flipInfoAll.orientation==2);
i_xneg=find(flipInfoAll.orientation==-1);
i_yneg=find(flipInfoAll.orientation==-2);

i_x1b=[]; i_x2b=[]; i_xnegb=[];
i_y1b=[]; i_ynegb=[];

%----- X CHANNEL
% fit +X1 to model y = a*t + b*T + f*(dT/dt) + c1(t<tgap) + c2(t>tgap) + d*exp(-t/f)
M.PF.X1=fit_cals(flipInfoAll,i_x1,i_x1b,2,5.3e-5,0.003);

% fit +X2 to model y = a*t + b*T + f*(dT/dt) + c1(t<tgap) + c2(t0>t>tgap) c3(t>t0) + d*exp(-t/f)
M.PF.X2=fit_cals(flipInfoAll,i_x2,i_x2b,3,5.3e-5,0.003);

% fit -X to model y = a*t + b*T + f*(dT/dt) + c [+ d*exp(-t/f)]
M.PF.negX=fit_cals(flipInfoAll,i_xneg,i_xnegb,4,-5.3e-5,-0.003);

% X1 span
tspan1=flipInfoAll.t(intersect(i_x1,i_xneg-4));
xspan1=(flipInfoAll.gCal(intersect(i_x1,i_xneg-4))+flipInfoAll.gCal(intersect(i_x1+4,i_xneg)));
p_xspan1=polyfit(tspan1-tspan1(1),xspan1,1);
m_xspan1=polyval(p_xspan1,tspan1-tspan1(1));
disp(['X1 Span Slope = ' num2str(p_xspan1(1)*365*10^5) ' ug/yr'])
disp(['X1 Span Misfit = ' num2str(std(xspan1-m_xspan1)*10^5) ' ug'])

% X2 span
tspan2=flipInfoAll.t(intersect(i_x2,i_xneg-1));
xspan2=(flipInfoAll.gCal(intersect(i_x2,i_xneg-1))+flipInfoAll.gCal(intersect(i_x2+1,i_xneg)));
p_xspan2=polyfit(tspan2-tspan2(1),xspan2,1);
m_xspan2=polyval(p_xspan2,tspan2-tspan2(1));
disp(['X2 Span Slope = ' num2str(p_xspan2(1)*365*10^5) ' ug/yr'])
disp(['X2 Span Misfit = ' num2str(std(xspan2-m_xspan2)*10^5) ' ug'])

%----- Y CHANNEL
% fit +Y to model y = a*t + b*T + f*(dT/dt) + c1(t<tgap) + c2(t>tgap) + d*exp(-t/f)
M.PF.Y=fit_cals(flipInfoAll,i_y1,i_y1b,2,8.44e-5,0.0026);

% fit -Y to model y = a*t + b*T + f*(dT/dt) + c [+ d*exp(-t/f)]
M.PF.negY=fit_cals(flipInfoAll,i_yneg,i_ynegb,4,-8.44e-5,-0.0026);

% Y span
tspan=flipInfoAll.t(intersect(i_y1,i_yneg-1));
yspan=flipInfoAll.gCal(intersect(i_y1,i_yneg-1))+flipInfoAll.gCal(intersect(i_y1+1,i_yneg));
p_yspan=polyfit(tspan-tspan(1),yspan,1);
m_yspan=polyval(p_yspan,tspan-tspan(1));
disp(['Y Span Slope = ' num2str(p_yspan(1)*365*10^5) ' ug/yr'])
disp(['Y Span Misfit = ' num2str(std(yspan-m_yspan)*10^5) ' ug'])

end

%% Additional functions

function M = fit_cals(flipInfoAll,i_good,i_bad,model,varargin)
% fits calibrations to model of chosen form:
%   0:  y = a*t + c + d*exp(-t/f)
%   1:  y = a*t + c1(t<t0) + c2(t>t0) + d*exp(-t/f)
%   2:  y = a*t + b*T + f*(dT/dt) + c1(t<tgap) + c2(t>tgap) + d*exp(-t/f)
%   3:  y = a*t + b*T + f*(dT/dt) + c1(t<tgap) + c2(t0>t>tgap) + c3(t>t0) + d*exp(-t/f)
%   4:  y = a*t + b*T + f*(dT/dt) + c + d*exp(-t/f)
%

t_good1=flipInfoAll.t(i_good)-flipInfoAll.t(i_good(1));
t_good=t_good1/max(t_good1);
a_good=flipInfoAll.gCal(i_good);
T_good=flipInfoAll.T(i_good);

if model==2
    % separate either side of gap
    a1=a_good(1:15);
    t1=t_good(1:15);
    a2=a_good(16:end);
    t2=t_good(16:end);
    
    % calculate dT/dt from daily values
    T_day=interp1(t_good1,T_good,(t_good1(1):t_good1(end)+0.1)','pchip');
    dT_day=diff(T_day); dT_day=[dT_day(1);dT_day];
    [~,ig,~]=intersect(round((t_good1(1):t_good1(end)+0.1)'),round(t_good1));
    dT_good=dT_day(ig);
    
    % temperature-independent calibrations
    a_noT=a_good-varargin{1}*T_good-varargin{2}*dT_good;
    
    % first guess at parameters, based off residual to non-exponential model
    G_lin(:,1)=t_good;
    G_lin(:,2)=[ones(size(t1));zeros(size(t2))];
    G_lin(:,3)=[zeros(size(t1));ones(size(t2))];
    m_lin=(G_lin'*G_lin)\G_lin'*a_noT;
    a_lin=G_lin*m_lin;
    % use as correction, solve for exponent
    a_temp=a_noT-a_lin;
    fit_exp=fit(t_good,a_temp-min(a_temp),'exp1','Lower',[0 -1000],'Upper',[1 0],...
        'StartPoint',[(max(a_temp)-min(a_temp))*2 -1000/100]);
    d_guess=fit_exp.a;
    e_guess=-1/fit_exp.b;
    % use above model parameters as guess model
    m_star=[m_lin;d_guess;e_guess];
    a_star=a_lin+m_star(4)*exp(-t_good/m_star(5));
    % solve for final model
    m_star2=nonlin_inv({t1;t2},a_noT,m_star,a_star);
    a_star2=m_star2(1)*t_good+...
        m_star2(2)*[ones(size(t1));zeros(size(t2))]+m_star2(3)*[zeros(size(t1));ones(size(t2))]+...
        m_star2(4)*exp(-t_good/m_star2(5));
    % rescale model to inputs
    m_star2(1)=m_star2(1)/max(t_good1);
    
    M.m=[m_star2(1);varargin{1};varargin{2};m_star2(2:end)];
    M.a=a_star2+varargin{1}*T_good+varargin{2}*dT_good;
    M.drift=[num2str(M.m(1)*365*10^5) ' ug/yr'];
    M.misfit=[num2str(std(a_good-M.a)*10^5) ' ug'];
elseif model==3 % y = a*t + b*T + f*(dT/dt) + c1(t<tgap) + c2(t0>t>tgap) c3(t>t0) + d*exp(-t/f)
    % find indices of gap and offset
    a1=a_good(1:15);
    t1=t_good(1:15);
    i_5cal=find(flipInfoAll.orientation==-2,1);
    i2=find(i_good>=i_5cal,1);
    a2=a_good(16:i2-1);
    t2=t_good(16:i2-1);
    a3=a_good(i2:end);
    t3=t_good(i2:end);
    
    % calculate dT/dt from daily values
    T_day=interp1(t_good1,T_good,(t_good1(1):t_good1(end)+0.1)','pchip');
    dT_day=diff(T_day); dT_day=[dT_day(1);dT_day];
    [~,ig,~]=intersect(round((t_good1(1):t_good1(end)+0.1)'),round(t_good1));
    dT_good=dT_day(ig);
    
    % temperature-independent calibrations
    a_noT=a_good-varargin{1}*T_good-varargin{2}*dT_good;
    
    % first guess at parameters, based off residual to non-exponential model
    G_lin(:,1)=t_good;
    G_lin(:,2)=[ones(size(t1));zeros(size(t2));zeros(size(t3))];
    G_lin(:,3)=[zeros(size(t1));ones(size(t2));zeros(size(t3))];
    G_lin(:,4)=[zeros(size(t1));zeros(size(t2));ones(size(t3))];
    m_lin=(G_lin'*G_lin)\G_lin'*a_noT;
    a_lin=G_lin*m_lin;
    % use as correction, solve for exponent
    a_temp=a_noT-a_lin;
    fit_exp=fit(t_good,a_temp-min(a_temp),'exp1','Lower',[0 -1000],'Upper',[1 0],...
        'StartPoint',[(max(a_temp)-min(a_temp))*2 -2000/100]);
    d_guess=fit_exp.a;
    e_guess=-1/fit_exp.b;
    % use above model parameters as guess model
    m_star=[m_lin;d_guess;e_guess];
    a_star=a_lin+m_star(5)*exp(-t_good/m_star(6));
    % solve for final model
    m_star2=nonlin_inv({t1;t2;t3},a_noT,m_star,a_star);
    a_star2=m_star2(1)*t_good+...
        m_star2(2)*[ones(size(t1));zeros(size(t2));zeros(size(t3))]+...
        m_star2(3)*[zeros(size(t1));ones(size(t2));zeros(size(t3))]+...
        m_star2(4)*[zeros(size(t1));zeros(size(t2));ones(size(t3))]+...
        m_star2(5)*exp(-t_good/m_star2(6));
    % rescale model to inputs
    m_star2(1)=m_star2(1)/max(t_good1);
    
    M.m=[m_star2(1);varargin{1};varargin{2};m_star2(2:end)];
    M.a=a_star2+varargin{1}*T_good+varargin{2}*dT_good;
    M.drift=[num2str(M.m(1)*365*10^5) ' ug/yr'];
    M.misfit=[num2str(std(a_good-M.a)*10^5) ' ug'];
elseif model==4 % y = a*t + b*T + f*(dT/dt) + ct [+ d*exp(-t/e)]
    % calculate dT/dt from daily values
    T_day=interp1(t_good1,T_good,(t_good1(1):t_good1(end)+0.1)','pchip');
    dT_day=diff(T_day); dT_day=[dT_day(1);dT_day];
    [~,ig,~]=intersect(round((t_good1(1):t_good1(end)+0.1)'),round(t_good1));
    dT_good=dT_day(ig);
    
    % temperature-independent calibrations
    a_noT=a_good-varargin{1}*T_good-varargin{2}*dT_good;
    
    % because t>>e, no need to include exponential term
    G_lin(:,1)=t_good;
    G_lin(:,2)=ones(size(t_good));
    m_lin=(G_lin'*G_lin)\G_lin'*a_noT;
    a_lin=G_lin*m_lin;
    
    m_star2=[m_lin;0;0];
    a_star2=a_lin;
    % rescale model to inputs
    m_star2(1)=m_star2(1)/max(t_good1);
    
    M.m=[m_star2(1);varargin{1};varargin{2};m_star2(2:end)];
    M.a=a_star2+varargin{1}*T_good+varargin{2}*dT_good;
    M.drift=[num2str(M.m(1)*365*10^5) ' ug/yr'];
    M.misfit=[num2str(std(a_good-M.a)*10^5) ' ug'];
end

figure(1); clf; hold on
plot(t_good,a_good,'o')
plot(t_good,a_star2+varargin{1}*T_good+varargin{2}*dT_good,'s')
keyboard

end

function m = nonlin_inv(t_obs,a_obs,m_star,a_star,varargin)
% use standard nonlinear methods to solve inversions

iter=true;
redo=false;
count=0;
count2=0;
disp('starting misfit:')
disp(norm(a_obs-a_star))
while iter
    if ~redo
        norm_check=norm(a_obs-a_star); %condition check
        
        % matrix of partial derivatives
        if length(m_star)==5
            % y = a*(t') + c1(t'<t0) + c2(t'>t0) + d*exp(-(t')/e)
            if iscell(t_obs)
                t1=t_obs{1}; t2=t_obs{2};
                t_obs=cat(1,t_obs{:});
                dyda=t_obs;
                dydc1=[ones(size(t1));zeros(size(t2))];
                dydc2=[zeros(size(t1));ones(size(t2))];
                G_lin=cat(2,dyda,dydc1,dydc2);
            end
            dyde=exp(-t_obs/m_star(5));
            dydf=m_star(4)*t_obs/m_star(5)^2.*exp(-t_obs/m_star(5));
            G_prime=cat(2,G_lin,dyde,dydf*10^4);
        elseif length(m_star)==6
            % y = a*(t') + c1(t'<tgap) + c2(t0>t'>tgap) + c3(t'>t0) + d*exp(-(t')/e)
            if iscell(t_obs)
                t1=t_obs{1}; t2=t_obs{2}; t3=t_obs{3};
                t_obs=cat(1,t_obs{:});
                dyda=t_obs;
                dydc1=[ones(size(t1));zeros(size(t2));zeros(size(t3))];
                dydc2=[zeros(size(t1));ones(size(t2));zeros(size(t3))];
                dydc3=[zeros(size(t1));zeros(size(t2));ones(size(t3))];
                G_lin=cat(2,dyda,dydc1,dydc2,dydc3);
            end
            dyde=exp(-t_obs/m_star(6));
            dydf=m_star(5)*t_obs/m_star(6)^2.*exp(-t_obs/m_star(6));
            G_prime=cat(2,G_lin,dyde,dydf*10^4);
        end
        
        dm_star=inv(G_prime'*G_prime)*G_prime'*(a_obs-a_star);

    else
        dm_star=dm_star/2; %checking for overstep
    end
    
    m_star2=m_star+[dm_star(1:end-1);dm_star(end)*10^4];
    while m_star2(end)<0 %ensures exponential decay
        dm_star(end)=dm_star(end)/2;
        m_star2(end)=m_star(end)+dm_star(end)*10^4;
    end
    
    a_star2=G_lin*m_star2(1:end-2)+m_star2(end-1)*exp(-t_obs/m_star2(end));
    
    new_norm=norm(a_obs-a_star2);
    disp('new misfit:')
    disp(new_norm)
    
    if new_norm>=norm_check
        if count<25
            redo=true;
            count=count+1;
        else
            iter=false;
        end
    elseif abs(norm_check-new_norm)/norm_check<0.001
        if count2>5
            m_star=m_star2;
            a_star=a_star2;
            iter=false;
        else
            redo=true;
            count2=count2+1;
        end
    else
        m_star=m_star2;
        a_star=a_star2;
        
        redo=false;
        count=0;
    end
end
disp('final misfit:')
disp(norm(a_obs-a_star))
m=m_star;

end