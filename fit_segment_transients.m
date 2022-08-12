% fit_segment_transients.m
%   Model of the form
%       am = a*(t') + c + b*exp(-(t')/d)
%

clear; close all

load('../longterm_tilt/PinonFlat/WW_method/postflip_transients.mat')

% try one exponential term
for i=1:33
    if i==1 || i==17
        continue
    elseif i<=22
        x1=transients(i).x(17:end);
    else
        x1=transients(i).x(21:end);
    end
    x=x1(5:200);
    x1=x1(5:end);
    
    % test a variety of time constants
    lamda=1:0.1:15;
    tj=(1:length(x))';
    norm_list=[];
    for j=1:length(lamda)
        gexp1=exp(-tj/lamda(j));
            %construct matrix for inversion
            Gj=[tj,ones(size(tj)),gexp1];
            mj{j}=inv(Gj'*Gj)*Gj'*x;
            xm{j}=Gj*mj{j};
            norm_list(j)=norm(x-xm{j});
    end
    
    [~,im]=min(norm_list);
    lamda_best(:,i)=lamda(im);
    norm_best(:,i)=norm_list(im);
    m_best(:,i)=mj{im};
    x_best(:,i)=xm{im};
    exp_best(:,i)=m_best(3,i)*exp(-tj/lamda(im));
    
%     % plots to compare between segments
%     figure(3)
%     subplot(311); hold on
%     plot(x,'o')
%     subplot(312); hold on
%     plot(x_best(:,i)-m_best(2,i),'linewidth',1)
%     subplot(313); hold on
%     pcolor(lamda,lamda,norm_list)
    
%     % plots for visually inspecting fit quality
%     figure(4); clf
%     subplot(211); hold on
%     plot(x,'linewidth',1)
%     plot(x_best(:,i),'linewidth',1)
%     plot(x1,'--')
%     plot(m_best(1,i)*(1:length(x1))'+m_best(2,i)+m_best(3,i)*exp(-(1:length(x1))'/lamda_best(:,i)),'--')
% %     xlim([0 50])
%     subplot(212); hold on
%     plot(detrend(x-x_best(:,i)+m_best(1,i)*tj),'linewidth',1)
%     plot(detrend(x1-(m_best(2,i)+m_best(3,i)*exp(-(1:length(x1))'/lamda_best(:,i)))),':')
% %     xlim([0 50])

    % plots for visually inspecting fit quality
    figure(i); clf; hold on
    x_all=[x-m_best(3,i)*exp(-tj/lamda_best(:,i));x1(196:end)];
    plot(x_all,'linewidth',1)
    ly=ylim;
    plot(x1,':')
    ylim([ly(1) ly(1)+2*(ly(2)-ly(1))])
%     xlim([0 50])
end

clear
load('../longterm_tilt/PinonFlat/WW_method/postflip_transients.mat')

% try two exponential terms
for i=1:33
    if i==1 || i==17
        continue
    elseif i<=22
        x1=transients(i).x(17:end);
    else
        x1=transients(i).x(21:end);
    end
    x=x1(5:200);
    x1=x1(5:end);
    
    % test a variety of time constants
    lamda=1:0.1:15;
    tj=(1:length(x))';
    norm_list=[];
    for j=1:length(lamda)
        gexp1=exp(-tj/lamda(j));
        for k=1:length(lamda)
            if k==j
                norm_list(j,k)=1;
                continue
            end
            gexp2=exp(-tj/lamda(k));
            %construct matrix for inversion
            Gj=[tj,ones(size(tj)),gexp1,gexp2];
            %     rcond(Gj'*Gj)
            mj{j,k}=inv(Gj'*Gj)*Gj'*x;
            xm{j,k}=Gj*mj{j,k};
            norm_list(j,k)=norm(x-xm{j,k});
        end
    end
    
    [~,im1]=min(min(norm_list));
    [~,im2]=min(norm_list(im1,:));
    lamda1_best(:,i)=lamda(im1);
    lamda2_best(:,i)=lamda(im2);
    norm_best(:,i)=norm_list(im1,im2);
    m_best(:,i)=mj{im1,im2};
    x_best(:,i)=xm{im1,im2};
    exp1_best(:,i)=m_best(3,i)*exp(-tj/lamda(im1));
    exp2_best(:,i)=m_best(4,i)*exp(-tj/lamda(im2));
    
%     % plots to compare between segments
%     figure(3)
%     subplot(311); hold on
%     plot(x,'o')
%     subplot(312); hold on
%     plot(x_best(:,i)-m_best(2,i),'linewidth',1)
%     subplot(313); hold on
%     pcolor(lamda,lamda,norm_list)
    
%     % plots for visually inspecting fit quality
%     figure(4); clf
%     subplot(211); hold on
%     plot(x,'linewidth',1)
%     plot(x_best(:,i),'linewidth',1)
%     plot(x1,'--')
%     plot(m_best(1,i)*(1:length(x1))'+m_best(2,i)+m_best(3,i)*exp(-(1:length(x1))'/lamda1_best(:,i))+m_best(4,i)*exp(-(1:length(x1))'/lamda2_best(:,i)),'--')
% %     xlim([0 50])
%     subplot(212); hold on
%     plot(detrend(x-x_best(:,i)+m_best(1,i)*tj),'linewidth',1)
%     plot(detrend(x1-(m_best(2,i)+m_best(3,i)*exp(-(1:length(x1))'/lamda1_best(:,i))+m_best(4,i)*exp(-(1:length(x1))'/lamda2_best(:,i)))),':')
% %     xlim([0 50])

    % plots for visually inspecting fit quality
    figure(i);
    x_all=[x-(m_best(3,i)*exp(-tj/lamda1_best(:,i))+m_best(4,i)*exp(-tj/lamda2_best(:,i)));x1(196:end)];
    plot(x_all,'linewidth',1)
    xline(96)
end

% manually determine lags that best align segment transients
exp_best2(:,15)=interp1(1:length(exp_best),exp_best(:,15),1:0.1:length(exp_best));
for i=1:33
    if i==1 || i==17
        continue
    end
    
    exp_best2(:,i)=interp1(1:length(exp_best),exp_best(:,i),1:0.1:length(exp_best));
    
    figure(5)
    [~,im2(:,i)]=min(abs(exp_best2(1,15)-exp_best2(:,i)));
    i2a=1;
    i2b=length(exp_best2(:,15))-im2(:,i)+1;
    iia=im2(:,i);
    iib=length(exp_best2(:,i));
    t=(iia:iib)';
    y1=exp_best2(i2a:i2b,15);
    y2=exp_best2(iia:iib,i);
    subplot(211); hold on
    plot(y2,'o')
    subplot(212); hold on
    plot(t,y1-y2,'linewidth',1)
end

