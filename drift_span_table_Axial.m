function M = drift_span_table_Axial
%
% Calculate drift rates, spans, and misfits for all locations and channels
%

%% Axial Location 1

load('../calibrations/Axial/axialdata','flipInfoAll')

% identify and separate each of the calibrations
i_x=find(flipInfoAll.orientation==1);
i_y=find(flipInfoAll.orientation==2);
i_negx=find(flipInfoAll.orientation==-1);
i_negy=find(flipInfoAll.orientation==-2);

%----- X CHANNEL
% separate good calibrations from bad
i_x1=i_x([1:2:61,73:2:127]); i_x1b=i_x([63,65,67,69,71,129]);
i_x2=i_x([2:2:62,70:2:96,100:2:126,130]); i_x2b=i_x([64,66,68,98,128]);
i_xneg=i_negx([1:6,9:39,41]); i_xnegb=i_negx([7,8,40]);

% fit +X1 to model y = a*t + c + d*exp(-t/f)
M.Axial1.X1=fit_cals(flipInfoAll,i_x1,i_x1b,0);

% fit +X2 to model y = a*t + c1(t<t0) + c2(t>t0) + d*exp(-t/f)
M.Axial1.X2=fit_cals(flipInfoAll,i_x2,i_x2b,1);

% fit -X to model y = a*t + c + d*exp(-t/f)
M.Axial1.negX=fit_cals(flipInfoAll,i_xneg,i_xnegb,0);

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
% separate good calibrations from bad
i_y1=i_y([1:31,34:50,52:end]); i_y1b=i_y([32,33,51]);
i_yneg=i_negy([1:7,12:14,16:18,20:26,28,30:40]); i_ynegb=i_negy([8,9,10,11,15,19,27,29,41]);

% fit +Y to model y = a*t + c + d*exp(-t/f)
M.Axial1.Y=fit_cals(flipInfoAll,i_y1,i_y1b,0);

% fit -Y to model y = a*t + c + d*exp(-t/f)
M.Axial1.negY=fit_cals(flipInfoAll,i_yneg,i_ynegb,0);

% Y span
tspan=flipInfoAll.t(intersect(i_y1,i_yneg-1));
yspan=flipInfoAll.gCal(intersect(i_y1,i_yneg-1))+flipInfoAll.gCal(intersect(i_y1+1,i_yneg));
p_yspan=polyfit(tspan-tspan(1),yspan,1);
m_yspan=polyval(p_yspan,tspan-tspan(1));
disp(['Y Span Slope = ' num2str(p_yspan(1)*365*10^5) ' ug/yr'])
disp(['Y Span Misfit = ' num2str(std(yspan-m_yspan)*10^5) ' ug'])

%% Axial Location 2

load('../calibrations/Axial/axialdata_newloc','flipInfoAll')

% identify and separate each of the calibrations
i_x=find(flipInfoAll.orientation==1);
i_y=find(flipInfoAll.orientation==2);
i_negx=find(flipInfoAll.orientation==-1);
i_negy=find(flipInfoAll.orientation==-2);

%----- X CHANNEL
% separate good calibrations from bad
i_x1=i_x([1:2:127,131:2:end]); i_x1b=i_x(129);
i_x2=i_x(2:2:end);
i_xneg=i_negx([2:38,44:50,52:63,65:end]); i_xnegb=i_negx([1,39:43,51,64]);

% fit +X1 to model y = a*t + c + d*exp(-t/f)
M.Axial2.X1=fit_cals(flipInfoAll,i_x1,i_x1b,0);

% fit +X2 to model y = a*t + c + d*exp(-t/f)
M.Axial2.X2=fit_cals(flipInfoAll,i_x2,i_x2b,0);

% fit -X to model y = a*t + c + d*exp(-t/f)
M.Axial2.negX=fit_cals(flipInfoAll,i_xneg,i_xnegb,0);

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
% separate good calibrations from bad
i_y1=i_y([2:17,19:39,41:end]); i_y1b=i_y([1,18,40]);
i_yneg=i_negy([1:38,40:63,65:end]); i_ynegb=i_negy([39,64]);

% fit +Y to model y = a*t + c + d*exp(-t/f)
M.Axial2.Y=fit_cals(flipInfoAll,i_y1,i_y1b,0);

% fit -Y to model y = a*t + c + d*exp(-t/f)
M.Axial2.negY=fit_cals(flipInfoAll,i_yneg,i_ynegb,0);

% Y span
tspan=flipInfoAll.t(intersect(i_y1,i_yneg-1));
yspan=flipInfoAll.gCal(intersect(i_y1,i_yneg-1))+flipInfoAll.gCal(intersect(i_y1+1,i_yneg));
p_yspan=polyfit(tspan-tspan(1),yspan,1);
m_yspan=polyval(p_yspan,tspan-tspan(1));
disp(['Y Span Slope = ' num2str(p_yspan(1)*365*10^5) ' ug/yr'])
disp(['Y Span Misfit = ' num2str(std(yspan-m_yspan)*10^5) ' ug'])
end

%% Additional functions

function M = fit_cals(flipInfoAll,i_good,i_bad,model)
% fits calibrations to model of chosen form:
%   0:  y = a*t + c + d*exp(-t/f)
%   1:  y = a*t + c1(t<t0) + c2(t>t0) + d*exp(-t/f)
%

t_good1=flipInfoAll.t(i_good)-flipInfoAll.t(i_good(1));
t_good=t_good1/max(t_good1);
a_good=flipInfoAll.gCal(i_good);
T_good=flipInfoAll.T(i_good);

if model == 0
    % get simple linear fit
    p=polyfit(t_good,a_good,1);
    a_lin=polyval(p,t_good);
    % use as correction, solve for exponent
    a_temp=a_good-a_lin;
    fit_exp=fit(t_good,a_temp-max(a_temp),'exp1','Lower',[-1 -1000],'Upper',[0 0],...
        'StartPoint',[min(a_temp)-max(a_temp)*2 -800/100]);
    d_guess=fit_exp.a;
    e_guess=-1/fit_exp.b;
    % use above model parameters as guess model
    m_star=[p';d_guess;e_guess];
    a_star=a_lin+m_star(3)*exp(-t_good/m_star(4));
    % solve for final model
    m_star2=nonlin_inv(t_good,a_good,m_star,a_star);
    a_star2=m_star2(1)*t_good+m_star2(2)+m_star2(3)*exp(-t_good/m_star2(4));
    % rescale model to inputs
    m_star2(1)=m_star2(1)/max(t_good1);
    
    M.m=m_star2;
    M.a=a_star2;
    M.drift=[num2str(M.m(1)*365*10^5) ' ug/yr'];
    M.misfit=[num2str(std(a_good-M.a)*10^5) ' ug'];
elseif model==1
    % find index of offset
    i_5cal=find(flipInfoAll.orientation==-2,1);
    a1=a_good(i_good<i_5cal);
    t1=t_good(i_good<i_5cal);
    a2=a_good(i_good>=i_5cal);
    t2=t_good(i_good>=i_5cal);
    % get simple linear fit
    p=polyfit(t_good,a_good,1);
    a_lin=polyval(p,t_good);
    % use as correction, solve for exponent
    a_temp=a_good-a_lin;
    fit_exp=fit(t_good,a_temp-max(a_temp),'exp1','Lower',[-1 -1000],'Upper',[0 0],...
        'StartPoint',[min(a_temp)-max(a_temp)*2 -800/100]);
    d_guess=fit_exp.a;
    e_guess=-1/fit_exp.b;
    % use above model parameters as guess model
    m_star=[p';p(2);d_guess;e_guess];
    a_star=a_lin+m_star(4)*exp(-t_good/m_star(5));
    % solve for final model
    m_star2=nonlin_inv({t1;t2},a_good,m_star,a_star);
    a_star2=m_star2(1)*t_good+m_star2(2)*[ones(size(t1));zeros(size(t2))]+...
        m_star2(3)*[zeros(size(t1));ones(size(t2))]+m_star2(4)*exp(-t_good/m_star2(5));
    % rescale model to inputs
    m_star2(1)=m_star2(1)/max(t_good1);
    
    M.m=m_star2;
    M.a=a_star2;
    M.drift=[num2str(M.m(1)*365*10^5) ' ug/yr'];
    M.misfit=[num2str(std(a_good-M.a)*10^5) ' ug'];
end

figure(1); clf; hold on
plot(t_good,a_good,'o')
plot(t_good,a_star2,'s')
% keyboard

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
        if length(m_star)==4
            % y = a*(t') + c + d*exp(-(t')/e)
            dyda=t_obs;
            dydc=ones(size(t_obs));
            G_lin=cat(2,dyda,dydc);
            dyde=exp(-t_obs/m_star(4));
            dydf=m_star(3)*t_obs/m_star(4)^2.*exp(-t_obs/m_star(4));
            G_prime=cat(2,G_lin,dyde,dydf*10^4);
        elseif length(m_star)==5
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