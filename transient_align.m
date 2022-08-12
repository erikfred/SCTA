% transient_align.m
%
% finds optimal alignment via cross correlation
%

clear; close all

load('../longterm_tilt/PinonFlat/WW_method/postflip_transients.mat')

% visual inspection of initial similarity
figure(3); hold on
for i=1:33
    if i~=1 && i~=17
        plot(transients(i).x_cor-transients(i).x_cor(end),'o')
    end
end

% interpolate truncated transients for higher res
figure(4); hold on
for i=1:33
    if i~=1 && i~=17
        temp(i).x_cor=interp1(1:length(transients(i).x_cor(21:end)),transients(i).x_cor(21:end),1:1/60:length(transients(i).x_cor(21:end)));
        temp(i).x_cor=temp(i).x_cor-temp(i).x_cor(end);
        plot(temp(i).x_cor,'o')
    end
end

% % cross correlate and save lags
% figure(5)
% for i=1:33
%     if i~=1 && i~=2 && i~=17
%         [c,lags]=xcorr(temp(2).x_cor,temp(i).x_cor,'coeff');
%         subplot(211)
%         plot(1:length(temp(2).x_cor),temp(2).x_cor,'o',1:length(temp(i).x_cor),temp(i).x_cor,'s')
%         subplot(212)
% %         plot(lags,c,'o');
% %         plot(1:length(temp(2).x_cor),temp(2).x_cor-temp(i).x_cor,'.')
%         keyboard
%         [~,im]=max(c); temp(i).lag=lags(im);
%     end
% end

% manually align in time and save lag info
figure(5)
for i=1:33
    if i~=1 && i~=2 && i~=17
        if temp(2).x_cor(1)>temp(i).x_cor(1)
            [~,im]=min(abs(temp(2).x_cor-temp(i).x_cor(1)));
            i2a=im;
            i2b=length(temp(2).x_cor);
            iia=1;
            iib=length(temp(i).x_cor)-im+1;
            t=iia:iib;
            y1=temp(2).x_cor(i2a:i2b);
            y2=temp(i).x_cor(iia:iib);
        else
            [~,im]=min(abs(temp(2).x_cor(1)-temp(i).x_cor));
            i2a=1;
            i2b=length(temp(2).x_cor)-im+1;
            iia=im;
            iib=length(temp(i).x_cor);
            t=iia:iib;
            y1=temp(2).x_cor(i2a:i2b);
            y2=temp(i).x_cor(iia:iib);
        end
        subplot(211)
        plot(1:length(temp(2).x_cor),temp(2).x_cor,'o',1:length(temp(i).x_cor),temp(i).x_cor,'s'); hold on
        plot(t,y1,'o',t,y2,'s'); hold off
        subplot(212)
        plot(1:length(temp(2).x_cor),detrend(temp(2).x_cor-temp(i).x_cor),'.'); hold on
        plot(t,detrend(y1-y2),'.')
        text(25000,median(detrend(y1-y2)),num2str(std(detrend(temp(2).x_cor-temp(i).x_cor))),'color','b');
        text(25000,median(detrend(y1-y2))-1e-6,num2str(std(detrend(y1-y2))),'color','r'); hold off
        keyboard
%         [~,im]=max(c); temp(i).lag=lags(im);
    end
end