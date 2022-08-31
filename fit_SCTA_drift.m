function [M] = fit_SCTA_drift(time,accel,Temp,dTemp)
%
% FIT_DRIFT.m takes a pressure series and associated time as inputs and
% returns coefficients for a combined exponential+linear fit by least
% squares, assumed to represent the instrument drift.
%
% Inputs:
%   t - time in datenum format
%   a - acceleration
%   T - temperature
%   dT - dT/dt, as daily values
%
% Outputs:
%   	M = [a b c d e f], from
%       y = a*(t') + b*(T') + [c*(dT')] + d + e*exp(-(t')/f)
%

%% inversion without dT/dt dependence
if isempty(dTemp)
    
    % first guess at parameters, based off residual to non-exponential model
    
    % build starting matrix
    num_ind=length(time); % number of independent time series
    g_build=[];
    for i=1:num_ind
        g_build{i}=zeros(length(time{i}),1+2*num_ind);
        g_build{i}(:,i)=time{i};
        g_build{i}(:,1+num_ind)=Temp{i};
        g_build{i}(:,1+num_ind+i)=ones(size(time{i}));
    end
    G_lin=cat(1,g_build{:});
    
    m_lin=(G_lin'*G_lin)\G_lin'*cat(1,accel{:});
    a_lin_star=G_lin*m_lin;
    
    % assume no significant exponential in X2b or negX, solve only for X1 and X2a
    
    % fit X1 rediduals with exponential curve to get first guess
    res=accel{1}(1:15)-a_lin_star(1:15);
    fit_exp=fit(time{1}(1:15),res-min(res),'exp1');
    e1_guess=fit_exp.a;
    f1_guess=-1/fit_exp.b;
    
    if num_ind>1
        % fit X2a rediduals with exponential curve to get first guess
        res=accel{2}(1:15)-a_lin_star(length(accel{1})+1:length(accel{1})+15);
        fit_exp=fit(time{2}(1:15),res-min(res),'exp1');
        e2_guess=fit_exp.a;
        f2_guess=-1/fit_exp.b;
    end
    
    % generate "guess" model
    if num_ind==1
        m_star=[m_lin;e1_guess;f1_guess];
        a_star=a_lin_star+m_star(4)*exp(-time{1}/m_star(5));
    else
        m_star=[m_lin;e1_guess;e2_guess;f1_guess;f2_guess];
        a_star=a_lin_star+[m_star(10)*exp(-time{1}/m_star(12));...
            m_star(11)*exp(-time{2}/m_star(13));zeros(length(time{3}),1);...
            zeros(length(time{4}),1)];
    end
    
    %% plots to confirm things are working
    % compare guess against data and against initial model
    
    if num_ind==4
        figure(91); clf
        subplot(221)
        plot(time{1},accel{1}-mean(accel{1}),'k','linewidth',1)
        hold on
        plot(time{1},a_lin_star(1:length(time{1}))-mean(accel{1}),'r','linewidth',1)
        plot(time{1},a_star(1:length(time{1}))-mean(accel{1}),'b','linewidth',1)
        plot(time{1},accel{1}-a_lin_star(1:length(time{1})),'rx','linewidth',2)
        plot(time{1},accel{1}-a_star(1:length(time{1})),'bx','linewidth',2)
        xlabel('Days since first flip')
        ylabel('Acceleration (m/s^2)')
        title('X1')
        yyaxis right
        plot(time{1},Temp{1},'c:','linewidth',1)
        ylabel('Temperature (C)')
        legend('data','linear model','guess model','linear residual','guess residual','location','northwest')
        
        subplot(222)
        plot(time{2},accel{2}-mean(accel{2}),'k','linewidth',1)
        hold on
        plot(time{2},a_lin_star(length(time{1})+1:length([time{1};time{2}]))-mean(accel{2}),'r','linewidth',1)
        plot(time{2},a_star(length(time{1})+1:length([time{1};time{2}]))-mean(accel{2}),'b','linewidth',1)
        plot(time{2},accel{2}-a_lin_star(length(time{1})+1:length([time{1};time{2}])),'rx','linewidth',2)
        plot(time{2},accel{2}-a_star(length(time{1})+1:length([time{1};time{2}])),'bx','linewidth',2)
        xlabel('Days since first flip')
        ylabel('Acceleration (m/s^2)')
        title('X2a')
        yyaxis right
        plot(time{1},Temp{1},'c:','linewidth',1)
        ylabel('Temperature (C)')
        legend('data','linear model','guess model','linear residual','guess residual','location','northwest')
        
        subplot(223)
        plot(time{3},accel{3}-mean(accel{3}),'k','linewidth',1)
        hold on
        plot(time{3},a_lin_star(length([time{1};time{2}])+1:length([time{1};time{2};time{3}]))-mean(accel{3}),'r','linewidth',1)
        plot(time{3},a_star(length([time{1};time{2}])+1:length([time{1};time{2};time{3}]))-mean(accel{3}),'b','linewidth',1)
        plot(time{3},accel{3}-a_lin_star(length([time{1};time{2}])+1:length([time{1};time{2};time{3}])),'rx','linewidth',2)
        plot(time{3},accel{3}-a_star(length([time{1};time{2}])+1:length([time{1};time{2};time{3}])),'bx','linewidth',2)
        xlabel('Days since first flip')
        ylabel('Acceleration (m/s^2)')
        title('X2b')
        yyaxis right
        plot(time{1},Temp{1},'c:','linewidth',1)
        ylabel('Temperature (C)')
        legend('data','linear model','guess model','linear residual','guess residual','location','northwest')
        
        subplot(224)
        plot(time{4},accel{4}-mean(accel{4}),'k','linewidth',1)
        hold on
        plot(time{4},a_lin_star(length([time{1};time{2};time{3}])+1:length([time{1};time{2};time{3};time{4}]))-mean(accel{4}),'r','linewidth',1)
        plot(time{4},a_star(length([time{1};time{2};time{3}])+1:length([time{1};time{2};time{3};time{4}]))-mean(accel{4}),'b','linewidth',1)
        plot(time{4},accel{4}-a_lin_star(length([time{1};time{2};time{3}])+1:length([time{1};time{2};time{3};time{4}])),'rx','linewidth',2)
        plot(time{4},accel{4}-a_star(length([time{1};time{2};time{3}])+1:length([time{1};time{2};time{3};time{4}])),'bx','linewidth',2)
        xlabel('Days since first flip')
        ylabel('Acceleration (m/s^2)')
        title('-X')
        yyaxis right
        plot(time{1},Temp{1},'c:','linewidth',1)
        ylabel('Temperature (C)')
        legend('data','linear model','guess model','linear residual','guess residual','location','northwest')
    end
    
    %% iterate to solve for exponential terms
    % define useful variables
    a_obs=cat(1,accel{:});
    e_ind=2+2*num_ind;
    
    iter=true;
    redo=false;
    count=0;
    count2=0;
    disp('starting misfit:')
    disp(norm(a_obs-a_star))
    while iter
        if ~redo
            norm_check=norm(a_obs-a_star); %condition check
            
            % matrix of partial derivatives from equation of the form:
            % y = a*(t') + b*(T') + [c*(dT')] + d + e*exp(-(t')/f)
            % need to scale for well conditioned inversion
            G_prime=G_lin;
            if num_ind==1
                dAde1=exp(-time{1}/m_star(e_ind+1));
                dAdf1=m_star(e_ind)*time{1}/m_star(e_ind+1)^2.*exp(-time{1}/m_star(e_ind+1));
                
                G_prime=cat(2,G_prime,dAde1,dAdf1);
            else
                dAde1=[exp(-time{1}/m_star(e_ind+2));zeros(length(time{2}),1);zeros(length(time{3}),1);...
                    zeros(length(time{4}),1)];
                dAde2=[zeros(length(time{1}),1);exp(-time{2}/m_star(e_ind+3));zeros(length(time{3}),1);...
                    zeros(length(time{4}),1)];
                dAdf1=(10^4*m_star(e_ind)*[time{1};zeros(length(time{2}),1);zeros(length(time{3}),1);...
                    zeros(length(time{4}),1)]/m_star(e_ind+2)^2).*[exp(-time{1}/m_star(e_ind+2));...
                    zeros(length(time{2}),1);zeros(length(time{3}),1);zeros(length(time{4}),1)];
                dAdf2=(10^4*m_star(e_ind+1)*[zeros(length(time{1}),1);time{2};zeros(length(time{3}),1);...
                    zeros(length(time{4}),1)]/m_star(e_ind+3)^2).*[zeros(length(time{1}),1);...
                    exp(-time{2}/m_star(e_ind+3));zeros(length(time{3}),1);zeros(length(time{4}),1)];
                G_prime=cat(2,G_prime,dAde1,dAde2,dAdf1,dAdf2);
            end
            
            %setup inverse problem
            dm_star=inv(G_prime'*G_prime)*G_prime'*(a_obs-a_star);
%             dm_star=-dm_star/10;
        else
            dm_star=dm_star/2; %checking for overstep
        end
        
        m_star2=m_star+dm_star;
        if m_star2(e_ind+1)<0 %ensures exponential decay
            m_star2(e_ind+1)=m_star(e_ind+1)/100;
        end
        if num_ind>1
            if m_star2(e_ind+3)<0 %ensures exponential decay
                m_star2(e_ind+3)=m_star(e_ind+3)/100;
            end
            
            a_star2=[m_star2(1)*time{1};m_star2(2)*time{2};m_star2(3)*time{3};m_star2(4)*time{4}]...
                +m_star2(5)*[Temp{1};Temp{2};Temp{3};Temp{4}]...
                +[m_star2(6)*ones(length(time{1}),1);m_star2(7)*ones(length(time{2}),1);...
                m_star2(8)*ones(length(time{3}),1);m_star2(9)*ones(length(time{4}),1)]...
                +[m_star2(10)*exp(-time{1}/m_star2(12));m_star2(11)*exp(-time{2}/m_star2(13));...
                zeros(length(time{3}),1);zeros(length(time{4}),1)];
    
        else
            a_star2=m_star2(1)*time{1}+m_star2(2)*Temp{1}+m_star2(3)*ones(length(time{1}),1)+...
                m_star2(e_ind)*exp(-time{1}/m_star2(e_ind+1));
        end
        
        new_norm=norm(a_obs-a_star2);
        disp('new misfit:')
        disp(new_norm)
        
        if new_norm>norm_check
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
    keyboard
    disp('final misfit:')
    disp(norm(a_obs-a_star))
    M=a_star;

%% inversion with  dT/dt dependence
else
    
    % first guess at parameters, based off residual to non-exponential model
    
    % build starting matrix
    num_ind=length(time); % number of independent time series
    g_build=[];
    for i=1:num_ind
        g_build{i}=zeros(length(time{i}),2+2*num_ind);
        g_build{i}(:,i)=time{i};
        g_build{i}(:,1+num_ind)=Temp{i};
        g_build{i}(:,2+num_ind)=dTemp{i};
        g_build{i}(:,2+num_ind+i)=ones(size(time{i}));
    end
    G_lin=cat(1,g_build{:});
    
    m_lin=(G_lin'*G_lin)\G_lin'*cat(1,accel{:});
    a_lin_star=G_lin*m_lin;
    
    % assume no significant exponential in X2b or negX, solve only for X1 and X2a
    
    % fit X1 or Y rediduals with exponential curve to get first guess
    res=accel{1}(1:10)-a_lin_star(1:10);
    fit_exp=fit(time{1}(1:10),res-min(res),'exp1');
    e1_guess=fit_exp.a;
    f1_guess=-1/fit_exp.b;
    
    if num_ind==4
        % fit X2a rediduals with exponential curve to get first guess
        res=accel{2}(1:10)-a_lin_star(length(accel{1})+1:length(accel{1})+10);
        fit_exp=fit(time{2}(1:10),res-min(res),'exp1');
        e2_guess=fit_exp.a;
        f2_guess=-1/fit_exp.b;
    end
    
    % generate "guess" model
    if num_ind==4
        m_star=[m_lin;e1_guess;e2_guess;f1_guess;f2_guess];
        a_star=a_lin_star+[m_star(11)*exp(-time{1}/m_star(13));...
            m_star(12)*exp(-time{2}/m_star(14));zeros(length(time{3}),1);...
            zeros(length(time{4}),1)];
    elseif num_ind==2
        m_star=[m_lin;e1_guess;f1_guess];
        a_star=a_lin_star+[m_star(7)*exp(-time{1}/m_star(8));zeros(length(time{2}),1)];
    else
        m_star=[m_lin;e1_guess;f1_guess];
        a_star=a_lin_star+m_star(5)*exp(-time{1}/m_star(6));
    end
    
    %% plots to confirm things are working
    % compare guess against data and against initial model
    
    if num_ind==4
        
        figure(91); clf
        subplot(221)
        plot(time{1},accel{1}-mean(accel{1}),'k','linewidth',1)
        hold on
        plot(time{1},a_lin_star(1:length(time{1}))-mean(accel{1}),'r','linewidth',1)
        plot(time{1},a_star(1:length(time{1}))-mean(accel{1}),'b','linewidth',1)
        plot(time{1},accel{1}-a_lin_star(1:length(time{1})),'rx','linewidth',2)
        plot(time{1},accel{1}-a_star(1:length(time{1})),'bx','linewidth',2)
        xlabel('Days since first flip')
        ylabel('Acceleration (m/s^2)')
        title('X1')
        yyaxis right
        plot(time{1},Temp{1},'c:','linewidth',1)
        ylabel('Temperature (C)')
        legend('data','linear model','guess model','linear residual','guess residual','location','northwest')
        
        subplot(222)
        plot(time{2},accel{2}-mean(accel{2}),'k','linewidth',1)
        hold on
        plot(time{2},a_lin_star(length(time{1})+1:length([time{1};time{2}]))-mean(accel{2}),'r','linewidth',1)
        plot(time{2},a_star(length(time{1})+1:length([time{1};time{2}]))-mean(accel{2}),'b','linewidth',1)
        plot(time{2},accel{2}-a_lin_star(length(time{1})+1:length([time{1};time{2}])),'rx','linewidth',2)
        plot(time{2},accel{2}-a_star(length(time{1})+1:length([time{1};time{2}])),'bx','linewidth',2)
        xlabel('Days since first flip')
        ylabel('Acceleration (m/s^2)')
        title('X2a')
        yyaxis right
        plot(time{1},Temp{1},'c:','linewidth',1)
        ylabel('Temperature (C)')
        legend('data','linear model','guess model','linear residual','guess residual','location','northwest')
        
        subplot(223)
        plot(time{3},accel{3}-mean(accel{3}),'k','linewidth',1)
        hold on
        plot(time{3},a_lin_star(length([time{1};time{2}])+1:length([time{1};time{2};time{3}]))-mean(accel{3}),'r','linewidth',1)
        plot(time{3},a_star(length([time{1};time{2}])+1:length([time{1};time{2};time{3}]))-mean(accel{3}),'b','linewidth',1)
        plot(time{3},accel{3}-a_lin_star(length([time{1};time{2}])+1:length([time{1};time{2};time{3}])),'rx','linewidth',2)
        plot(time{3},accel{3}-a_star(length([time{1};time{2}])+1:length([time{1};time{2};time{3}])),'bx','linewidth',2)
        xlabel('Days since first flip')
        ylabel('Acceleration (m/s^2)')
        title('X2b')
        yyaxis right
        plot(time{1},Temp{1},'c:','linewidth',1)
        ylabel('Temperature (C)')
        legend('data','linear model','guess model','linear residual','guess residual','location','northwest')
        
        subplot(224)
        plot(time{4},accel{4}-mean(accel{4}),'k','linewidth',1)
        hold on
        plot(time{4},a_lin_star(length([time{1};time{2};time{3}])+1:length([time{1};time{2};time{3};time{4}]))-mean(accel{4}),'r','linewidth',1)
        plot(time{4},a_star(length([time{1};time{2};time{3}])+1:length([time{1};time{2};time{3};time{4}]))-mean(accel{4}),'b','linewidth',1)
        plot(time{4},accel{4}-a_lin_star(length([time{1};time{2};time{3}])+1:length([time{1};time{2};time{3};time{4}])),'rx','linewidth',2)
        plot(time{4},accel{4}-a_star(length([time{1};time{2};time{3}])+1:length([time{1};time{2};time{3};time{4}])),'bx','linewidth',2)
        xlabel('Days since first flip')
        ylabel('Acceleration (m/s^2)')
        title('-X')
        yyaxis right
        plot(time{1},Temp{1},'c:','linewidth',1)
        ylabel('Temperature (C)')
        legend('data','linear model','guess model','linear residual','guess residual','location','northwest')
    end
    
    %% iterate to solve for exponential terms
    % define useful variables
    a_obs=cat(1,accel{:});
    e_ind=3+2*num_ind;
    
    iter=true;
    redo=false;
    count=0;
    count2=0;
    disp('starting misfit:')
    disp(norm(a_obs-a_star))
    while iter
        if ~redo
            norm_check=norm(a_obs-a_star); %condition check
            
            % matrix of partial derivatives from equation of the form:
            % y = a*(t') + b*(T') + [c*(dT')] + d + e*exp(-(t')/f)
            % need to scale for well conditioned inversion
            G_prime=G_lin;
            if num_ind==4
                dAde1=[exp(-time{1}/m_star(e_ind+2));zeros(length(time{2}),1);zeros(length(time{3}),1);...
                    zeros(length(time{4}),1)];
                dAde2=[zeros(length(time{1}),1);exp(-time{2}/m_star(e_ind+3));zeros(length(time{3}),1);...
                    zeros(length(time{4}),1)];
                dAdf1=(10^4*m_star(e_ind)*[time{1};zeros(length(time{2}),1);zeros(length(time{3}),1);...
                    zeros(length(time{4}),1)]/m_star(e_ind+2)^2).*[exp(-time{1}/m_star(e_ind+2));...
                    zeros(length(time{2}),1);zeros(length(time{3}),1);zeros(length(time{4}),1)];
                dAdf2=(10^4*m_star(e_ind+1)*[zeros(length(time{1}),1);time{2};zeros(length(time{3}),1);...
                    zeros(length(time{4}),1)]/m_star(e_ind+3)^2).*[zeros(length(time{1}),1);...
                    exp(-time{2}/m_star(e_ind+3));zeros(length(time{3}),1);zeros(length(time{4}),1)];
                
                G_prime=cat(2,G_prime,dAde1,dAde2,dAdf1,dAdf2);
            elseif num_ind==2
                dAde1=[exp(-time{1}/m_star(e_ind+1));zeros(length(time{2}),1)];
                dAdf1=10^4*m_star(e_ind)*[time{1};zeros(length(time{2}),1)]/m_star(e_ind+1)^2.*...
                    [exp(-time{1}/m_star(e_ind+1));zeros(length(time{2}),1)];
                
                G_prime=cat(2,G_prime,dAde1,dAdf1);
            else
                dAde1=exp(-time{1}/m_star(e_ind+1));
                dAdf1=m_star(e_ind)*time{1}/m_star(e_ind+1)^2.*exp(-time{1}/m_star(e_ind+1));
                
                G_prime=cat(2,G_prime,dAde1,dAdf1*10^7);
            end
            
            %setup inverse problem
            dm_star=inv(G_prime'*G_prime)*G_prime'*(a_obs-a_star);
%             dm_star=-dm_star/10;
        else
            dm_star=dm_star/2; %checking for overstep
        end
        
        m_star2=m_star+[dm_star(1:end-1); dm_star(end)*10^7];
        if num_ind==4
            while m_star2(e_ind+3)<0 %ensures exponential decay
                dm_star(e_ind+3)=dm_star(e_ind+3)/2;
                m_star2(e_ind+3)=m_star(e_ind+3)+dm_star(e_ind+3);
            end
            
            a_star2=[m_star2(1)*time{1};m_star2(2)*time{2};m_star2(3)*time{3};m_star2(4)*time{4}]...
                +m_star2(5)*[Temp{1};Temp{2};Temp{3};Temp{4}]...
                +m_star(6)*[dTemp{1};dTemp{2};dTemp{3};dTemp{4}]...
                +[m_star2(7)*ones(length(time{1}),1);m_star2(8)*ones(length(time{2}),1);...
                m_star2(9)*ones(length(time{3}),1);m_star2(10)*ones(length(time{4}),1)]...
                +[m_star2(11)*exp(-time{1}/m_star2(13));m_star2(12)*exp(-time{2}/m_star2(14));...
                zeros(length(time{3}),1);zeros(length(time{4}),1)];
    
        elseif num_ind==2
            while m_star2(e_ind+1)<0 %ensures exponential decay
                dm_star(e_ind+1)=dm_star(e_ind+1)/2;
                m_star2(e_ind+1)=m_star(e_ind+1)+dm_star(e_ind+1);
            end
            a_star2=[m_star2(1)*time{1};m_star2(2)*time{2}]+m_star2(3)*[Temp{1};Temp{2}]+...
                m_star2(4)*[dTemp{1};dTemp{2}]+...
                [m_star2(5)*ones(length(time{1}),1);m_star2(6)*ones(length(time{2}),1)]+...
                [m_star2(e_ind)*exp(-time{1}/m_star2(e_ind+1));zeros(length(time{2}),1)];
        else
            a_star2=m_star2(1)*time{1}+m_star2(2)*Temp{1}+m_star2(3)*dTemp{1}+...
                m_star2(4)*ones(length(time{1}),1)+m_star2(e_ind)*exp(-time{1}/m_star2(e_ind+1));
        end
        
        new_norm=norm(a_obs-a_star2);
        disp('new misfit:')
        disp(new_norm)
        
        if new_norm>norm_check
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
    M=a_star;
    
    % plots to demonstrate fit improvement
    if num_ind>1
        figure(89); clf
        y_it=accel{1}; y_it=y_it-y_it(1); y_its=y_it*10^5;
        plot(time{1},y_its,'ok','markersize',15)
        hold on
        y_it=y_it-m_star(5)*Temp{1}; y_it=y_it-y_it(1); y_its=y_it*10^5;
        plot(time{1},y_its,'rs','markersize',15)
        y_it=y_it-m_star(6)*dTemp{1}; y_it=y_it-y_it(1); y_its=y_it*10^5;
        plot(time{1},y_its,'b^','markersize',15)
        y_it=y_it-m_star(11)*exp(-time{1}/m_star(13)); y_it=y_it-y_it(1); y_its=y_it*10^5;
        plot(time{1},y_its,'kd','markersize',15)
        y_it=y_it-m_star(1)*time{1}; y_it=y_it-y_it(1); y_its=y_it*10^5;
        plot(time{1},y_its,'gp','markersize',15)
        xlabel('Days since first flip')
        ylabel('Acceleration (\mug)')
        title('Calibrations')
        yyaxis right
        plot(time{1},Temp{1},':c','linewidth',2)
        ylabel('Temperature (C)')
        legend('calibration','T-corrected','dT/dt-corrected','exp-corrected','lin-corrected','location','northwest')
        set(gca,'fontsize',16)
        box on; grid on; grid minor
        
        figure(90); clf
        plot(time{1},y_its,'gp','markersize',15)
        xlabel('Days since first flip')
        ylabel('Acceleration (\mug)')
        title('Fully-corrected calibrations')
        set(gca,'fontsize',16)
        box on; grid on; grid minor
        
        figure(91); clf; hold on
        plot(time{1},accel{1}-a_star(1:length(time{1})),'ko','linewidth',2)
        xlabel('Days since first flip')
        ylabel('Acceleration (m/s^2)')
        title(['+X1 \sigma = ' num2str(std(accel{1}-a_star(1:length(time{1})))) ' m/s^2'])
        yyaxis right
        plot(time{1},Temp{1},'c:','linewidth',2)
        ylabel('Temperature (C)')
        legend('data','linear model','final model','linear residual','final residual','location','northwest')
        
        subplot(222)
        plot(time{2},accel{2}-mean(accel{2}),'-ok','linewidth',1)
        hold on
        plot(time{2},a_lin_star(length(time{1})+1:length([time{1};time{2}]))-mean(accel{2}),'-or','linewidth',1)
        plot(time{2},a_star(length(time{1})+1:length([time{1};time{2}]))-mean(accel{2}),'-ob','linewidth',1)
        plot(time{2},accel{2}-a_lin_star(length(time{1})+1:length([time{1};time{2}])),'rx','linewidth',2)
        plot(time{2},accel{2}-a_star(length(time{1})+1:length([time{1};time{2}])),'bx','linewidth',2)
        xlabel('Days since first flip')
        ylabel('Acceleration (m/s^2)')
        title(['+X2a \sigma = ' num2str(std(accel{2}-a_star(length(time{1})+1:length([time{1};time{2}])))) ' m/s^2'])
        yyaxis right
        plot(time{1},Temp{1},'c:','linewidth',2)
        ylabel('Temperature (C)')
        legend('data','linear model','final model','linear residual','final residual','location','northwest')
        
        subplot(223)
        plot(time{3},accel{3}-mean(accel{3}),'-ok','linewidth',1)
        hold on
        plot(time{3},a_lin_star(length([time{1};time{2}])+1:length([time{1};time{2};time{3}]))-mean(accel{3}),'-or','linewidth',1)
        plot(time{3},a_star(length([time{1};time{2}])+1:length([time{1};time{2};time{3}]))-mean(accel{3}),'-ob','linewidth',1)
        plot(time{3},accel{3}-a_lin_star(length([time{1};time{2}])+1:length([time{1};time{2};time{3}])),'rx','linewidth',2)
        plot(time{3},accel{3}-a_star(length([time{1};time{2}])+1:length([time{1};time{2};time{3}])),'bx','linewidth',2)
        xlabel('Days since first flip')
        ylabel('Acceleration (m/s^2)')
        title(['+X2b \sigma = ' num2str(std(accel{3}-a_star(length([time{1};time{2}])+1:length([time{1};time{2};time{3}])))) ' m/s^2'])
        yyaxis right
        plot(time{1},Temp{1},'c:','linewidth',2)
        ylabel('Temperature (C)')
        legend('data','linear model','final model','linear residual','final residual','location','northwest')
        
        subplot(224)
        plot(time{4},accel{4}-mean(accel{4}),'-ok','linewidth',1)
        hold on
        plot(time{4},a_lin_star(length([time{1};time{2};time{3}])+1:length([time{1};time{2};time{3};time{4}]))-mean(accel{4}),'-or','linewidth',1)
        plot(time{4},a_star(length([time{1};time{2};time{3}])+1:length([time{1};time{2};time{3};time{4}]))-mean(accel{4}),'-ob','linewidth',1)
        plot(time{4},accel{4}-a_lin_star(length([time{1};time{2};time{3}])+1:length([time{1};time{2};time{3};time{4}])),'rx','linewidth',2)
        plot(time{4},accel{4}-a_star(length([time{1};time{2};time{3}])+1:length([time{1};time{2};time{3};time{4}])),'bx','linewidth',2)
        xlabel('Days since first flip')
        ylabel('Acceleration (m/s^2)')
        title(['-X \sigma = ' num2str(std(accel{4}-a_star(length([time{1};time{2};time{3}])+1:length([time{1};time{2};time{3};time{4}])))) ' m/s^2'])
        yyaxis right
        plot(time{1},Temp{1},'c:','linewidth',2)
        ylabel('Temperature (C)')
        legend('data','linear model','final model','linear residual','final residual','location','northwest')
    end
    
    figure(70); clf; hold on
    cal_plots(time,accel,Temp,dTemp,a_star,m_star)
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print('../calibrations/PinonFlat/process_PinonFlat_Tmodel_y','-dtiff','-r300')
    print('../calibrations/PinonFlat/process_PinonFlat_Tmodel_y','-depsc','-painters')
    
    figure(71); clf; hold on
    calcor_plots(time,accel,Temp,dTemp,a_star,m_star)
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 8.5 11];
    print('../calibrations/PinonFlat/process_PinonFlat_Tcorrected_y','-dtiff','-r300')
    print('../calibrations/PinonFlat/process_PinonFlat_Tcorrected_y','-depsc','-painters')
    
end

figure(2); clf
plot(time{1},a_obs-9.793,'bo','markersize',10,'linewidth',1)
hold on
plot(time{1},a_lin_star-9.793,'r^','markersize',10,'linewidth',1)
plot(time{1},a_star-9.793,'s','markersize',10,'linewidth',1)
plot(time{1},a_obs-a_star,'kx','markersize',10,'linewidth',2)
ylabel('Acceleration - 9.793 (m/s^2)')
xlabel('Days since first Calibration')
title(['+Y Only \sigma = ' num2str(std(a_obs-a_star)) ' m/s^2'])
set(gca,'FontSize',16)
legend('Calibrations','Linear Model','Exponential Model','Misfit','location','northwest')

temp={['A = ' num2str(m_star(1)) ' (m/s^2)/d'];['B = ' num2str(m_star(2)) ' (m/s^2)/C'];...
    ['C = ' num2str(m_star(3)) ' (m/s^2)/(C/d)'];['D = ' num2str(m_star(4)) ' m/s^2'];...
    ['E = ' num2str(m_star(5)) ' m/s^2'];['F = ' num2str(m_star(6)) ' d^-^1']};
text(50,4e-4,temp,'fontsize',16)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 11 8.5];
print('../temperature_dependence/PinonFlat/by_orientation/PF_expInversion_dT_Y','-dtiff','-r300')

end

%%%%%%%%%%

function cal_plots(time,accel,Temp,dTemp,a_star,m_star)

if length(time)==4
        
    % Plot calibrations on separate axes
    subplot(311); hold on
    plot(time{1},accel{1},'ko','markersize',12)
    plot(time{1},a_star(1:length(time{1})),'k--','linewidth',1)
    ylabel('+X1 (m/s^2)')
    ylim([9.79250 9.79400])
    % title({'Pinon Flat SCTA X Calibrations',[datestr(startDate,'mmm dd, yyyy') ' - ' datestr(endDate,'mmm dd, yyyy')]})
    yyaxis right
    plot([time{2};time{3}],[accel{2};accel{3}],'ro','markersize',12)
    plot([time{2};time{3}],a_star(length(time{1})+1:length([time{1};time{2};time{3}])),'r--','linewidth',1)
    datetick('x',3)
    xtickangle(45)
    ylabel('+X2 (m/s^2)')
    set(gca,'fontsize',12)
    set(gca,'Ycolor','r')
    ylim([9.79275 9.79425])    
    % xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.0005 -0.0005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'100 \mug','fontsize',12)
    xl=xlim;
    legend('+X1','+X1 model','+X2','+X2 model','location','northwest')
    box on; grid on
    subplot(312); hold on
    plot(time{4},accel{4},'ks','markersize',12)
    plot(time{4},a_star(length([time{1};time{2};time{3}])+1:length([time{1};time{2};time{3};time{4}])),'k--','linewidth',1)
    xlim(xl)
    datetick('x',3,'keeplimits')
    xtickangle(45)
    ylabel('-X (m/s^2)')
    set(gca,'fontsize',12)
    ylim([-9.79300 -9.79150])
    % xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.0005 -0.0005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'100 \mug','fontsize',12)
    legend('-X','-X model','location','northwest')
    box on; grid on
    subplot(313); hold on
    xspan1=accel{1}(end-length(accel{4})+1:end)-accel{4};
    plot(time{4},xspan1-xspan1(1),'k^','markersize',12);
    ylabel('\Delta X1 span (m/s^2)')
    ylim([-0.00004 0.00006])
    yyaxis right
    xspan2=accel{3}-accel{4};
    plot(time{4},xspan2-xspan2(1),'r^','markersize',12);
    xlim(xl)
    datetick('x',3,'keeplimits')
    xtickangle(45)
    ylabel('\Delta X2 span (m/s^2)')
    set(gca,'fontsize',12)
    ylim([-0.00002 0.00008])
    % xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.00005 -0.00005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'10 \mug','fontsize',12)
    legend('X1','X2','location','northwest')    
    box on; grid on
    
elseif length(time)==2
    
    % Plot calibrations on separate axes
    subplot(311); hold on
    plot(time{1},accel{1},'ko','markersize',12)
    plot(time{1},a_star(1:length(time{1})),'k--','linewidth',1)
    datetick('x',3)
    xtickangle(45)
    ylabel('+Y (m/s^2)')
    % title({'Pinon Flat SCTA X Calibrations',[datestr(startDate,'mmm dd, yyyy') ' - ' datestr(endDate,'mmm dd, yyyy')]})
    set(gca,'fontsize',12)
    ylim([9.79375 9.79525])
    % xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.0005 -0.0005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'100 \mug','fontsize',12)
    xl=xlim;
    legend('+Y','+Y model','location','northwest')
    box on; grid on
    subplot(312); hold on
    plot(time{2},accel{2},'ks','markersize',12)
    plot(time{2},a_star(length(time{1})+1:length(a_star)),'k--','linewidth',1)
    xlim(xl)
    datetick('x',3,'keeplimits')
    xtickangle(45)
    ylabel('-Y (m/s^2)')
    set(gca,'fontsize',12)
    ylim([-9.79200 -9.79050])
    % xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.0005 -0.0005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'100 \mug','fontsize',12)
    legend('-Y','-Y model','location','northwest')
    box on; grid on
    subplot(313); hold on
    yspan=accel{1}(end-length(accel{2})+1:end)-accel{2};
    plot(time{2},yspan-yspan(1),'k^','markersize',12);
    xlim(xl)
    datetick('x',3,'keeplimits')
    xtickangle(45)
    ylabel('\Delta Y span (m/s^2)')
    set(gca,'fontsize',12)
    ylim([-0.00002 0.00008])
    % xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.00005 -0.00005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'10 \mug','fontsize',12)
    box on; grid on

end
end

function calcor_plots(time,accel,Temp,dTemp,a_star,m_star)

if length(time)==4
    
    x1=accel{1}-a_star(1:length(time{1}))+(time{1}-time{1}(1))*m_star(1)+m_star(7);
    x2a=accel{2}-a_star(length(time{1})+1:length([time{1};time{2}]))+(time{2}-time{2}(1))*m_star(2)+m_star(8);
    x2b=accel{3}-a_star(length([time{1};time{2}])+1:length([time{1};time{2};time{3}]))+(time{3}-time{3}(1))*m_star(3)+m_star(9)+(m_star(9)-m_star(8));
    negx=accel{4}-a_star(length([time{1};time{2};time{3}])+1:length([time{1};time{2};time{3};time{4}]))+(time{4}-time{4}(1))*m_star(4)+m_star(10);
    
    disp(['+X1 Slope = ' num2str(m_star(1)*365*10^5) ' ug/yr'])
    disp(['+X1  Misfit = ' num2str(std(accel{1}-a_star(1:length(time{1})))*10^5) ' ug'])
    
    p_x2=polyfit([time{2};time{3}]-time{2}(1),[x2a;x2b],1);
    m_x2=polyval(p_x2,[time{2};time{3}]-time{2}(1));
    disp(['+X2 Slope = ' num2str(p_x2(1)*365*10^5) ' ug/yr'])
    disp(['+X2 Misfit = ' num2str(std([x2a;x2b]-m_x2)*10^5) ' ug'])
    
    disp(['-X Slope = ' num2str(m_star(4)*365*10^5) ' ug/yr'])
    disp(['-X  Misfit = ' num2str(std(accel{4}-a_star(length([time{1};time{2};time{3}])+1:length([time{1};time{2};time{3};time{4}])))*10^5) ' ug'])
    
    % Plot calibrations on separate axes
    subplot(311); hold on
    plot(time{1},x1,'ko','markersize',12)
    ylabel('+X1 (m/s^2)')
    ylim([9.79075 9.79225])
    % title({'Pinon Flat SCTA X Calibrations',[datestr(startDate,'mmm dd, yyyy') ' - ' datestr(endDate,'mmm dd, yyyy')]})
    yyaxis right
    plot([time{2};time{3}],[x2a;x2b],'ro','markersize',12)
    datetick('x',3)
    xtickangle(45)
    ylabel('+X2 (m/s^2)')
    set(gca,'fontsize',12)
    ylim([9.79075 9.79225])
    % xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.0005 -0.0005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'100 \mug','fontsize',12)
    xl=xlim;
    legend('X1','X2','location','northwest')
    box on; grid on
    subplot(312); hold on
    plot(time{4},negx,'ks','markersize',12)
    xlim(xl)
    datetick('x',3,'keeplimits')
    xtickangle(45)
    ylabel('-X (m/s^2)')
    set(gca,'fontsize',12)
    ylim([-9.79500 -9.79350])
    % xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.0005 -0.0005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'100 \mug','fontsize',12)
    box on; grid on
    subplot(313); hold on
    xspan1=x1(end-length(negx)+1:end)-negx;
    plot(time{4},xspan1-xspan1(1),'k^','markersize',12);
    ylabel('\Delta X1 span (m/s^2)')
    ylim([-0.00005 0.00005])
    yyaxis right
    xspan2=x2b-negx;
    plot(time{4},xspan2-xspan2(1),'r^','markersize',12);
    xlim(xl)
    datetick('x',3,'keeplimits')
    xtickangle(45)
    ylabel('\Delta X2 span (m/s^2)')
    set(gca,'fontsize',12)
    ylim([-0.00002 0.00008])
    % xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.00005 -0.00005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'10 \mug','fontsize',12)
    legend('X1','X2','location','northwest')    
    box on; grid on
    
    p_xspan1=polyfit(time{4}-time{4}(1),xspan1,1);
    m_xspan1=polyval(p_xspan1,time{4}-time{4}(1));
    disp(['X1 Span Slope = ' num2str(p_xspan1(1)*365*10^5) ' ug/yr'])
    disp(['X1 Span Misfit = ' num2str(std(xspan1-m_xspan1)*10^5) ' ug'])
    
    p_xspan2=polyfit(time{4}-time{4}(1),xspan2,1);
    m_xspan2=polyval(p_xspan2,time{4}-time{4}(1));
    disp(['X2 Span Slope = ' num2str(p_xspan2(1)*365*10^5) ' ug/yr'])
    disp(['X2 Span Misfit = ' num2str(std(xspan2-m_xspan2)*10^5) ' ug'])
elseif length(time)==2
    
    y1=accel{1}-a_star(1:length(time{1}))+(time{1}-time{1}(1))*m_star(1)+m_star(5);
    negy=accel{2}-a_star(length(time{1})+1:length([time{1};time{2}]))+(time{2}-time{2}(1))*m_star(2)+m_star(6);
    
    disp(['+Y Slope = ' num2str(m_star(1)*365*10^5) ' ug/yr'])
    disp(['+Y  Misfit = ' num2str(std(accel{1}-a_star(1:length(time{1})))*10^5) ' ug'])
    
    disp(['-Y Slope = ' num2str(m_star(2)*365*10^5) ' ug/yr'])
    disp(['-Y  Misfit = ' num2str(std(accel{2}-a_star(length(time{1})+1:length([time{1};time{2}])))*10^5) ' ug'])
    
    % Plot calibrations on separate axes
    subplot(311); hold on
    plot(time{1},y1,'ko','markersize',12)
    datetick('x',3)
    xtickangle(45)
    ylabel('+Y (m/s^2)')
    % title({'Pinon Flat SCTA X Calibrations',[datestr(startDate,'mmm dd, yyyy') ' - ' datestr(endDate,'mmm dd, yyyy')]})
    set(gca,'fontsize',12)
    ylim([9.79150 9.79300])
    % xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.0005 -0.0005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'100 \mug','fontsize',12)
    xl=xlim;
    box on; grid on
    subplot(312); hold on
    plot(time{2},negy,'ks','markersize',12)
    xlim(xl)
    datetick('x',3,'keeplimits')
    xtickangle(45)
    ylabel('-Y (m/s^2)')
    set(gca,'fontsize',12)
    ylim([-9.79500 -9.79350])
    % xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.0005 -0.0005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'100 \mug','fontsize',12)
    box on; grid on
    subplot(313); hold on
    yspan=y1(end-length(negy)+1:end)-negy;
    plot(time{2},yspan-yspan(1),'k^','markersize',12);
    xlim(xl)
    datetick('x',3,'keeplimits')
    xtickangle(45)
    ylabel('\Delta Y span (m/s^2)')
    set(gca,'fontsize',12)
    ylim([-0.00005 0.00005])
    % xl=xlim; yl=ylim; plot([0 0]+xl(1)+diff(xl)/10,mean(yl)-[0.00005 -0.00005],'-k','linewidth',2); text(xl(1)+diff(xl)/9,mean(yl),'10 \mug','fontsize',12)
    box on; grid on
    
    p_yspan=polyfit(time{2}-time{2}(1),yspan,1);
    m_yspan=polyval(p_yspan,time{2}-time{2}(1));
    disp(['Y Span Slope = ' num2str(p_yspan(1)*365*10^5) ' ug/yr'])
    disp(['Y Span Misfit = ' num2str(std(yspan-m_yspan)*10^5) ' ug'])
end
end