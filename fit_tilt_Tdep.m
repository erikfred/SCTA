function [M] = fit_tilt_Tdep(time,accel,Temp,dTemp)
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

model=1; % 0 - unique offsets allowed, 1 - unique slopes allowed, 2 - dT/dt dependence

% scale datenums better for inversions
tstart=time{1}(1);
time{1}=time{1}-tstart;
time{2}=time{2}-tstart;

%% inversion without unique offsets on either side of gap
if model==0
    
    %----- assume form y = a*(t') + b*(T') + c1(t'<tgap) + c2(t'>tgap) + e*exp(-(t')/f)
    % starting matrix using only data after gap (avoids exponential and offset) to get T-dependence
    g2_lin=time{2};
    g2_T=Temp{2};
    g2_c=ones(size(time{2}));
    
    G2_lin=[g2_lin, g2_T, g2_c];
    m2_lin=(G2_lin'*G2_lin)\G2_lin'*accel{2};
    a2_lin_star=G2_lin*m2_lin;
    
    % determine exponential dependence by fitting first segment of temperature-corrected data 
    res=accel{1}-m2_lin(2)*Temp{1};
    if res(1)<res(end)
        fit_exp=fit(time{1},res-max(res),'exp1','Lower',[-1 -1],'Upper',[0 0],'StartPoint',[-1e-4 -1/5]);
    else
        fit_exp=fit(time{1},res-min(res),'exp1','Lower',[0 -1],'Upper',[1 0],'StartPoint',[-1e-4 -1/5]);
    end
    e1_guess=fit_exp.a;
    f1_guess=-1/fit_exp.b;
    
    % now apply T- and exponent- corrections to entire time series
    a_obs=cat(1,accel{:})-m2_lin(2)*cat(1,Temp{:})-e1_guess*exp(-cat(1,time{:})/f1_guess);
    
    % new matrix to get linear and offset terms
    g_lin=cat(1,time{:});
    g_c1=[ones(size(time{1})); zeros(size(time{2}))];
    g_c2=[zeros(size(time{1})); ones(size(time{2}))];

    G_lin=[g_lin, g_c1, g_c2];
    
    m_lin=(G_lin'*G_lin)\G_lin'*a_obs;
    a_lin_star=G_lin*m_lin;
    
    % compile into "guess" model
    m_star=[m_lin(1);m2_lin(2);m_lin(2:3);e1_guess;f1_guess];
    a_star=m_star(1)*cat(1,time{:})+m_star(2)*cat(1,Temp{:})+m_star(3)*g_c1+...
        m_star(4)*g_c2+m_star(5)*exp(-cat(1,time{:})/m_star(6));
    
    %% iterate to solve for exponential terms
    % define useful variables
    a_obs=cat(1,accel{:});
    t_obs=cat(1,time{:});
    t_obs1=[time{1}; zeros(size(time{2}))];
    t_obs2=[zeros(size(time{1})); time{2}];
    T_obs=cat(1,Temp{:});
    
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
            % y = a*(t') + b*(T') + c1(t'<tgap) + c2(t'>tgap) + e*exp(-(t')/f)
            dyda=t_obs;
            dydb=T_obs;
            dydc1=[ones(size(time{1})); zeros(size(time{2}))];
            dydc2=[zeros(size(time{1})); ones(size(time{2}))];
            G_lin=cat(2,dyda,dydb,dydc1,dydc2);
            dyde=exp(-t_obs/m_star(6));
            dydf=m_star(5)*t_obs/m_star(6)^2.*exp(-t_obs/m_star(6));
            
            G_prime=cat(2,G_lin,dyde,dydf*10^7);
            
            % invert for model perturbation
            if rcond(G_prime'*G_prime)<10^-16
%                 keyboard
            end
            dm_star=inv(G_prime'*G_prime)*G_prime'*(a_obs-a_star);
%             dm_star=-dm_star/10;
        else
            dm_star=dm_star/2; %checking for overstep
        end
        
        m_star2=m_star+[dm_star(1:5);dm_star(6)*10^7];
        while m_star2(6)<0 %ensures exponential decay
            dm_star(6)=dm_star(6)/2;
            m_star2(6)=m_star(6)+dm_star(6)*10^7;
        end
        
%         % plots to see how scaling works
%         figure(3); clf
%         subplot(321); hold on
%         plot(t_obs,a_obs-a_star)
%         plot(t_obs,a_obs-(m_star2(1)*dyda+m_star(2)*dydb+m_star(3)*dydc1+...
%             m_star(4)*dydc2+m_star(5)*exp(-t_obs/m_star(6))))
%         subplot(322); hold on
%         plot(t_obs,a_obs-a_star)
%         plot(t_obs,a_obs-(m_star(1)*dyda+m_star2(2)*dydb+m_star(3)*dydc1+...
%             m_star(4)*dydc2+m_star(5)*exp(-t_obs/m_star(6))))
%         subplot(323); hold on
%         plot(t_obs,a_obs-a_star)
%         plot(t_obs,a_obs-(m_star(1)*dyda+m_star(2)*dydb+m_star2(3)*dydc1+...
%             m_star(4)*dydc2+m_star(5)*exp(-t_obs/m_star(6))))
%         subplot(324); hold on
%         plot(t_obs,a_obs-a_star)
%         plot(t_obs,a_obs-(m_star(1)*dyda+m_star(2)*dydb+m_star(3)*dydc1+...
%             m_star2(4)*dydc2+m_star(5)*exp(-t_obs/m_star(6))))
%         subplot(325); hold on
%         plot(t_obs,a_obs-a_star)
%         plot(t_obs,a_obs-(m_star(1)*dyda+m_star(2)*dydb+m_star(3)*dydc1+...
%             m_star(4)*dydc2+m_star2(5)*exp(-t_obs/m_star(6))))
%         subplot(326); hold on
%         plot(t_obs,a_obs-a_star)
%         plot(t_obs,a_obs-(m_star(1)*dyda+m_star(2)*dydb+m_star(3)*dydc1+...
%             m_star(4)*dydc2+m_star(5)*exp(-t_obs/m_star2(6))))
%         % viewed the other way around
%         figure(4); clf
%         subplot(321); hold on
%         plot(t_obs,a_obs-a_star)
%         plot(t_obs,a_obs-(m_star(1)*dyda+m_star2(2)*dydb+m_star2(3)*dydc1+...
%             m_star2(4)*dydc2+m_star2(5)*exp(-t_obs/m_star2(6))))
%         subplot(322); hold on
%         plot(t_obs,a_obs-a_star)
%         plot(t_obs,a_obs-(m_star2(1)*dyda+m_star(2)*dydb+m_star2(3)*dydc1+...
%             m_star2(4)*dydc2+m_star2(5)*exp(-t_obs/m_star2(6))))
%         subplot(323); hold on
%         plot(t_obs,a_obs-a_star)
%         plot(t_obs,a_obs-(m_star2(1)*dyda+m_star2(2)*dydb+m_star(3)*dydc1+...
%             m_star2(4)*dydc2+m_star2(5)*exp(-t_obs/m_star2(6))))
%         subplot(324); hold on
%         plot(t_obs,a_obs-a_star)
%         plot(t_obs,a_obs-(m_star2(1)*dyda+m_star2(2)*dydb+m_star2(3)*dydc1+...
%             m_star(4)*dydc2+m_star2(5)*exp(-t_obs/m_star2(6))))
%         subplot(325); hold on
%         plot(t_obs,a_obs-a_star)
%         plot(t_obs,a_obs-(m_star2(1)*dyda+m_star2(2)*dydb+m_star2(3)*dydc1+...
%             m_star2(4)*dydc2+m_star(5)*exp(-t_obs/m_star2(6))))
%         subplot(326); hold on
%         plot(t_obs,a_obs-a_star)
%         plot(t_obs,a_obs-(m_star2(1)*dyda+m_star2(2)*dydb+m_star2(3)*dydc1+...
%             m_star2(4)*dydc2+m_star2(5)*exp(-t_obs/m_star(6))))
%         keyboard
        
        a_star2=G_lin*m_star2(1:4)+m_star2(5)*exp(-t_obs/m_star2(6));
        
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
    M=a_star;
    
    %----- PLOTTING
    % show goodness of fit and exponential-/temperature-corrected tilt
    figure(2); clf
    subplot(121); hold on
    plot(t_obs+tstart,a_obs,'k','linewidth',1)
    plot(t_obs+tstart,a_star,'r:','linewidth',2)
    datetick('x',3,'keeplimits')
    xtickangle(45)
    ylabel('Acceleration (m/s^2)')
    set(gca,'fontsize',12)
    plot(t_obs+tstart,a_obs-a_star,'b','linewidth',1)
    title(['X tilt (fit \sigma = ' num2str(std(a_obs-a_star)) ' m/s^2)'])
    legend('Input','Exp + T Model','Residual','location','northwest')
    grid on; box on
    subplot(122);
    plot(t_obs+tstart,a_obs-(m_star(2)*G_lin(:,2)+m_star(3)*G_lin(:,3)+...
        m_star2(4)*G_lin(:,4)+m_star2(5)*exp(-t_obs/m_star2(6))),'k','linewidth',1)
    datetick('x',3,'keeplimits')
    xtickangle(45)
    ylabel('Acceleration (m/s^2)')
    title('Drift-, Temperature-, and Exponent-corrected x tilt')
    grid on; box on
    set(gca,'fontsize',12)
    
    % show modeled alignment of segments (?)
    figure(3); clf
    plot(time{1}+tstart,accel{1}-(m_star2(1)*time{1}+m_star2(5)*exp(-time{1}/m_star2(6))),'k','linewidth',1)
    hold on
    plot(time{2}+tstart,accel{2}-m_star2(1)*time{2}+(m_star2(3)-m_star2(4)),'b','linewidth',1)
    yyaxis right
    plot(t_obs+tstart,-T_obs,'r.','linewidth',1)
    set(gca,'YColor','r')
    
    % temp={['A = ' num2str(m_star(1)) ' (m/s^2)/d'];['B = ' num2str(m_star(2)) ' (m/s^2)/C'];...
    %     ['C = ' num2str(m_star(3)) ' (m/s^2)/(C/d)'];['D = ' num2str(m_star(4)) ' m/s^2'];...
    %     ['E = ' num2str(m_star(5)) ' m/s^2'];['F = ' num2str(m_star(6)) ' d^-^1']};
    % text(50,4e-4,temp,'fontsize',12)
    %
    % fh=gcf;
    % fh.PaperUnits='inches';
    % fh.PaperPosition=[0 0 11 8.5];
    % print('../temperature_dependence/PinonFlat/by_orientation/PF_expInversion_dT_Y','-dtiff','-r300')

%% inversion with unique offsets AND unique linear terms
elseif model==1
    
    %----- assume form y = a1*(t'<tgap) + a2*(t'>tgap) b*(T') + c1(t'<tgap) + c2(t'>tgap) + e*exp(-(t')/f)
    % starting matrix using only data after gap (avoids exponential and offset) to get T-dependence
    g2_lin=time{2};
    g2_T=Temp{2};
    g2_c=ones(size(time{2}));
    
    G2_lin=[g2_lin, g2_T, g2_c];
    m2_lin=(G2_lin'*G2_lin)\G2_lin'*accel{2};
    a2_lin_star=G2_lin*m2_lin;
    
    % determine exponential dependence by fitting first segment of temperature-corrected data 
    res=accel{1}-m2_lin(2)*Temp{1};
    if res(1)<res(end)
        fit_exp=fit(time{1},res-max(res),'exp1','Lower',[-1 -1],'Upper',[0 0],'StartPoint',[-1e-4 -1/5]);
    else
        fit_exp=fit(time{1},res-min(res),'exp1','Lower',[0 -1],'Upper',[1 0],'StartPoint',[-1e-4 -1/5]);
    end
    e1_guess=fit_exp.a;
    f1_guess=-1/fit_exp.b;
    
    % now apply T- and exponent- corrections to entire time series
%     a_obs=cat(1,accel{:})-m2_lin(2)*cat(1,Temp{:})-e1_guess*exp(-cat(1,time{:})/f1_guess);
    a_obs=cat(1,accel{:})-e1_guess*exp(-cat(1,time{:})/f1_guess);
    
    % new matrix to get linear and offset terms
    g_lin1=[time{1}; zeros(size(time{2}))];
    g_lin2=[zeros(size(time{1})); time{2}];
    g_b=cat(1,Temp{:});
    g_c1=[ones(size(time{1})); zeros(size(time{2}))];
    g_c2=[zeros(size(time{1})); ones(size(time{2}))];

    G_lin=cat(2,g_lin1,g_lin2,g_b,g_c1,g_c2);
    
    m_lin=(G_lin'*G_lin)\G_lin'*a_obs;
    a_lin_star=G_lin*m_lin;
    
    % compile into "guess" model
%     m_star=[m_lin(1:2);m2_lin(2);m_lin(3:4);e1_guess;f1_guess];
    m_star=[m_lin;e1_guess;f1_guess];
    a_star=m_star(1)*g_lin1+m_star(2)*g_lin2+m_star(3)*cat(1,Temp{:})+m_star(4)*g_c1+...
        m_star(5)*g_c2+m_star(6)*exp(-cat(1,time{:})/m_star(7));
    
    %% iterate to solve for exponential terms
    % define useful variables
    a_obs=cat(1,accel{:});
    t_obs=cat(1,time{:});
    t_obs1=[time{1}; zeros(size(time{2}))];
    t_obs2=[zeros(size(time{1})); time{2}];
    T_obs=cat(1,Temp{:});
    T_obs1=[Temp{1}; zeros(size(Temp{2}))];
    T_obs2=[zeros(size(Temp{1})); Temp{2}];
    
    m_star=[m_star(1:3);m_star(3:end)];
    
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
            % y = a1*(t'<tgap) + a2*(t'>tgap) + b*(T') + c1(t'<tgap) + c2(t'>tgap) + e*exp(-(t')/f)
            dyda1=t_obs1;
            dyda2=t_obs2;
            dydb1=T_obs1;
            dydb2=T_obs2;
            dydc1=[ones(size(time{1})); zeros(size(time{2}))];
            dydc2=[zeros(size(time{1})); ones(size(time{2}))];
            G_lin=cat(2,dyda1,dyda2,dydb1,dydb2,dydc1,dydc2);
            dyde=exp(-t_obs/m_star(8));
            dydf=m_star(7)*t_obs/m_star(8)^2.*exp(-t_obs/m_star(8));
            
            G_prime=cat(2,G_lin,dyde,dydf*10^7);
            
            % invert for model perturbation
            if rcond(G_prime'*G_prime)<10^-16
                keyboard
            end
            dm_star=inv(G_prime'*G_prime)*G_prime'*(a_obs-a_star);
%             dm_star=-dm_star/10;
        else
            dm_star=dm_star/2; %checking for overstep
        end
        
        m_star2=m_star+[dm_star(1:7);dm_star(8)*10^7];
        while m_star2(8)<0 %ensures exponential decay
            dm_star(8)=dm_star(8)/2;
            m_star2(8)=m_star(8)+dm_star(8)*10^7;
        end
        
        a_star2=G_lin*m_star2(1:6)+m_star2(7)*exp(-t_obs/m_star2(8));
        
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
    disp(std(a_obs-a_star))
    M=a_star;
    
    %----- PLOTTING
    crd='East';
    if strcmp(crd,'North')
        tiltfac=-10^5;
        lg='southeast';
        svstr='N_accel';
    else
        tiltfac=10^5;
        lg='northeast';
        svstr='E_accel';
    end
    % show goodness of fit and exponential-/temperature-corrected tilt
    figure(2); clf; hold on
    plot(t_obs+tstart,a_obs*tiltfac,'k','linewidth',1)
    plot(t_obs+tstart,a_star*tiltfac,'r:','linewidth',2)
    datetick('x',3,'keeplimits')
    xtickangle(45)
    ylabel([crd ' Tilt (\murad)'])
    set(gca,'fontsize',12)
    plot(t_obs+tstart,(a_obs-a_star)*tiltfac,'b','linewidth',1)
    yyaxis right
    plot(t_obs+tstart,cat(1,Temp{:}),'m','linewidth',1)
    ylabel(['Temperature (' char(176) 'C)'])
    disp(['Full model (x tilt fit \sigma = ' num2str(std(a_obs-a_star)) ' m/s^2)'])
    legend('Input','Model','Residual','location',lg)
    grid on; box on
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 5.5 4.25];
    print(['../longterm_tilt/PinonFlat/manuscript/' svstr '_model'],'-dtiff','-r300')
    print(['../longterm_tilt/PinonFlat/manuscript/' svstr '_model'],'-depsc','-painters')

    figure(3); clf; hold on
    a_cor=((a_obs-a_star)+m_star(1)*G_lin(:,1)+m_star(2)*G_lin(:,2))*tiltfac;
    plot(t_obs+tstart,a_cor,'k','linewidth',1)
    datetick('x',3,'keeplimits')
    xtickangle(45)
    ylabel([crd ' Tilt (\murad)'])
%     title('Drift-, Temperature-, and Exponent-corrected x tilt')
    grid on; box on
    set(gca,'fontsize',12)
    fh=gcf;
    fh.PaperUnits='inches';
    fh.PaperPosition=[0 0 5.5 4.25];
    print(['../longterm_tilt/PinonFlat/manuscript/' svstr],'-dtiff','-r300')
    print(['../longterm_tilt/PinonFlat/manuscript/' svstr],'-depsc','-painters')
    
    % show modeled alignment of segments without other corrections
    figure(4); clf
    plot(time{1}+tstart,accel{1}-(m_star(4)+m_star(6)*exp(-time{1}/m_star(7))),'k','linewidth',1)
    hold on
    plot(time{2}+tstart,accel{2}-(m_star(5)),'b','linewidth',1)
    yyaxis right
    plot(t_obs+tstart,-T_obs,'r.','linewidth',1)
    set(gca,'YColor','r')
    
    % temp={['A = ' num2str(m_star(1)) ' (m/s^2)/d'];['B = ' num2str(m_star(2)) ' (m/s^2)/C'];...
    %     ['C = ' num2str(m_star(3)) ' (m/s^2)/(C/d)'];['D = ' num2str(m_star(4)) ' m/s^2'];...
    %     ['E = ' num2str(m_star(5)) ' m/s^2'];['F = ' num2str(m_star(6)) ' d^-^1']};
    % text(50,4e-4,temp,'fontsize',12)
    
elseif model==2
        
    %----- assume form y = a1*(t'<tgap) + a2*(t'>tgap) b*(T') + c1(t'<tgap) + c2(t'>tgap) + d*(dT'/dt) + e*exp(-(t')/f)
    % use only data after gap to get T and dT/dt dependence
    g2_lin=time{2};
    g2_T=Temp{2};
    g2_c=ones(size(time{2}));
    g2_d=dTemp{2};
    
    G2_lin=[g2_lin, g2_T, g2_c, g2_d];
    m2_lin=(G2_lin'*G2_lin)\G2_lin'*accel{2};
    a2_lin_star=G2_lin*m2_lin;
    
    % determine exponential dependence by fitting first segment of T- and dT- corrected data 
    res=accel{1}-m2_lin(2)*Temp{1}-m2_lin(4)*dTemp{1};
    if res(1)<res(end)
        fit_exp=fit(time{1},res-max(res),'exp1','Lower',[-1 -1],'Upper',[0 0],'StartPoint',[-1e-4 -1/5]);
    else
        fit_exp=fit(time{1},res-min(res),'exp1','Lower',[0 -1],'Upper',[1 0],'StartPoint',[-1e-4 -1/5]);
    end
    e1_guess=fit_exp.a;
    f1_guess=-1/fit_exp.b;
    
    % now apply exponent corrections to entire time series and redo fit to other terms
%     a_obs=cat(1,accel{:})-m2_lin(2)*cat(1,Temp{:})-e1_guess*exp(-cat(1,time{:})/f1_guess);
    a_obs=cat(1,accel{:})-e1_guess*exp(-cat(1,time{:})/f1_guess);
    
    % new matrix to get linear and offset terms
    g_lin1=[time{1}; zeros(size(time{2}))];
    g_lin2=[zeros(size(time{1})); time{2}];
    g_b=cat(1,Temp{:});
    g_c1=[ones(size(time{1})); zeros(size(time{2}))];
    g_c2=[zeros(size(time{1})); ones(size(time{2}))];
    g_d=cat(1,dTemp{:});

    G_lin=cat(2,g_lin1,g_lin2,g_b,g_c1,g_c2,g_d);
    
    m_lin=(G_lin'*G_lin)\G_lin'*a_obs;
    a_lin_star=G_lin*m_lin;
    
    % compile into "guess" model
%     m_star=[m_lin(1:2);m2_lin(2);m_lin(3:4);e1_guess;f1_guess];
    m_star=[m_lin;e1_guess;f1_guess];
    a_star=m_star(1)*g_lin1+m_star(2)*g_lin2+m_star(3)*g_b+m_star(4)*g_c1+...
        m_star(5)*g_c2+m_star(6)*g_d+m_star(7)*exp(-cat(1,time{:})/m_star(8));
    
    %% iterate to solve for exponential terms
    % define useful variables
    a_obs=cat(1,accel{:});
    t_obs=cat(1,time{:});
    t_obs1=[time{1}; zeros(size(time{2}))];
    t_obs2=[zeros(size(time{1})); time{2}];
    T_obs=cat(1,Temp{:});
    dT_obs=cat(1,dTemp{:});
    
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
            % y = a1*(t'<tgap) + a2*(t'>tgap) + b*(T') + c1(t'<tgap) + c2(t'>tgap) + d*(dT*/dt) + e*exp(-(t')/f)
            dyda1=t_obs1;
            dyda2=t_obs2;
            dydb=T_obs;
            dydc1=[ones(size(time{1})); zeros(size(time{2}))];
            dydc2=[zeros(size(time{1})); ones(size(time{2}))];
            dydd=dT_obs;
            G_lin=cat(2,dyda1,dyda2,dydb,dydc1,dydc2,dydd);
            dyde=exp(-t_obs/m_star(8));
            dydf=m_star(7)*t_obs/m_star(8)^2.*exp(-t_obs/m_star(8));
            
            G_prime=cat(2,G_lin,dyde,dydf*10^7);
            
            % invert for model perturbation
            if rcond(G_prime'*G_prime)<10^-16
                keyboard
            end
            dm_star=inv(G_prime'*G_prime)*G_prime'*(a_obs-a_star);
%             dm_star=-dm_star/10;
        else
            dm_star=dm_star/2; %checking for overstep
        end
        
        m_star2=m_star+[dm_star(1:7);dm_star(8)*10^7];
        if m_star2(8)<0 %ensures exponential decay
            dm_star(8)=dm_star(8)/2;
            m_star2(8)=m_star(8)+dm_star(8)*10^7;
        end
        
        a_star2=G_lin*m_star2(1:6)+m_star2(7)*exp(-t_obs/m_star2(8));
        
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
    
    %----- PLOTTING
    % show goodness of fit and exponential-/temperature-corrected tilt
    figure(2); clf
    subplot(121); hold on
    plot(t_obs+tstart,a_obs,'k','linewidth',1)
    plot(t_obs+tstart,a_star,'r:','linewidth',2)
    datetick('x',3,'keeplimits')
    xtickangle(45)
    ylabel('Acceleration (m/s^2)')
    set(gca,'fontsize',12)
    plot(t_obs+tstart,a_obs-a_star,'b','linewidth',1)
%     yyaxis right
%     plot(t_obs+tstart,cat(1,Temp{:}),'m','linewidth',1)
%     ylabel(['Temperature (' char(176) 'C)'])
    title(['Full model (x tilt fit \sigma = ' num2str(std(a_obs-a_star)) ' m/s^2)'])
    legend('Input','Exp + T + dT Model','Residual','location','northwest')
    grid on; box on
    subplot(122);
    a_cor=(a_obs-a_star)+m_star(1)*G_lin(:,1)+m_star(2)*G_lin(:,2);
    plot(t_obs+tstart,a_cor,'k','linewidth',1)
    datetick('x',3,'keeplimits')
    xtickangle(45)
    ylabel('Acceleration (m/s^2)')
    title('Drift-, T-, dT- and Exponent-corrected x tilt')
    grid on; box on
    set(gca,'fontsize',12)
    
    % show modeled alignment of segments without other corrections
    figure(3); clf
    plot(time{1}+tstart,accel{1}-(m_star(4)+m_star(7)*exp(-time{1}/m_star(8))),'k','linewidth',1)
    hold on
    plot(time{2}+tstart,accel{2}-(m_star(5)),'b','linewidth',1)
    yyaxis right
    plot(t_obs+tstart,-T_obs,'r.','linewidth',1)
    set(gca,'YColor','r')
    
    % temp={['A = ' num2str(m_star(1)) ' (m/s^2)/d'];['B = ' num2str(m_star(2)) ' (m/s^2)/C'];...
    %     ['C = ' num2str(m_star(3)) ' (m/s^2)/(C/d)'];['D = ' num2str(m_star(4)) ' m/s^2'];...
    %     ['E = ' num2str(m_star(5)) ' m/s^2'];['F = ' num2str(m_star(6)) ' d^-^1']};
    % text(50,4e-4,temp,'fontsize',12)
    %
    % fh=gcf;
    % fh.PaperUnits='inches';
    % fh.PaperPosition=[0 0 11 8.5];
    % print('../temperature_dependence/PinonFlat/by_orientation/PF_expInversion_dT_Y','-dtiff','-r300')
    
end