% inv_residual_signal.m
%
% Sets up an inversion problem for each deflation event to determine the
% component that can't be explained purely by a  reversal of  the long-term
% inflation
%

clear; close all

%%%%%%%%%%%%
%  CONFIG  %
%%%%%%%%%%%%

% load_dir='dec2017';
% load_dir='jun2018';
load_dir='may2019';

%%%%%%%%%%%%

load(['../tiltcompare/inflation_reversal/in_progress_matfiles/' load_dir '/AXCC1qc.mat'])
load(['../tiltcompare/inflation_reversal/in_progress_matfiles/' load_dir '/AXEC2qc.mat'])
load(['../tiltcompare/inflation_reversal/in_progress_matfiles/' load_dir '/ASHESqc.mat'])
load(['../tiltcompare/inflation_reversal/in_progress_matfiles/' load_dir '/AXID1qc.mat'])

%% Version 1a
% NOT WELL POSED, DO NOT USE

% set up variables

U=[0.4;0.1;0.2;0.2]; %need to properly calculate these for each reversal [m/yr]

Ru=[AXCC1.p_int1;AXEC2.p_int1;ASHES.p_int1;AXID1.p_int1]*365*10^-2; % [m/yr]

Rd=[AXCC1.p_int2;AXEC2.p_int2;ASHES.p_int2;AXID1.p_int2]*365*10^-2; % [m/yr]

%applying sign reversal to AXCC1 x-channel
eu=[-AXCC1.x_rot1;AXEC2.x_rot1;ASHES.x_rot1;AXID1.x_rot1]; % [urad/d]
nu=[AXCC1.y_rot1;AXEC2.y_rot1;ASHES.y_rot1;AXID1.y_rot1]; % [urad/d]

Eu=eu*10^-6*3000*365; % [m/yr]
Nu=nu*10^-6*3000*365; % [m/yr]

%applying sign reversal to AXCC1 x-channel
ed=[-AXCC1.x_rot2;AXEC2.x_rot2;ASHES.x_rot2;AXID1.x_rot2]; % [urad/d]
nd=[AXCC1.y_rot2;AXEC2.y_rot2;ASHES.y_rot2;AXID1.y_rot2]; % [urad/d]

Ed=ed*10^-6*3000*365; % [m/yr]
Nd=nd*10^-6*3000*365; % [m/yr]

% set up inverse problem

up=[Ru;Eu;Nu];
dwn=[Rd;Ed;Nd];
multipliers=[[1;1;1;1;0;0;0;0;0;0;0;0],[-1;-1;-1;-1;0;0;0;0;0;0;0;0]];

G=[up,multipliers]; % dwn = G * m; m = [b b*Cu Cd]'

% % solve
% 
% m=(G'*G)\G'*dwn;
% 
% % determine residual signal
% 
% non_sym=dwn-G*m;

%% Version 1b
% assumes that long-term uplift is equal to short-term uplift of event
% solves (Pdwn_rel; Edwn; Ndwn) = e*(Pup; Eup; Nup) - Cd

% set up inverse problem

Gb=[[U;Eu;Nu],[-1;-1;-1;-1;0;0;0;0;0;0;0;0]]; % dwn = Gb * mb; mb = [e Cd]'

% solve

mb=(Gb'*Gb)\Gb'*dwn;

% determine residual signal

non_symb=dwn-Gb*mb;

non_symb(1:4)=non_symb(1:4)/(365*10^-2); % [cm/d]
non_symb(5:end)=non_symb(5:end)/(10^-6*3000*365); % [urad/d]

%% Version 1c
% DO NOT USE, can't be solved for absolute subsidentce because Cu/Cd can't
% be separated out

% % set up inverse problem
% 
% Gc=[[Ru;Eu;Nu],[1;1;1;1;0;0;0;0;0;0;0;0]]; % dwn = Gc * mc; mc = [f D]'; D = b*Cu - Cd
% 
% % solve
% 
% mc=(Gc'*Gc)\Gc'*dwn;
% 
% % determine residual signal
% 
% non_symc=dwn-Gc*mc;
% 
% non_symc(1:4)=non_symc(1:4)/(365*10^-2); % [cm/d]
% non_symc(5:end)=non_symc(5:end)/(10^-6*3000*365); % [urad/d]

%% Version 2
% solves Pup_rel = a*P_longterm - Cu to get Pup,
% then solves (Pdwn_rel; Edwn; Ndwn) = b*(Pup; Eup; Nup) - Cd to get a
% model of Pdwn_rel, Edwn, and Ndwn

% set up inverse problem

G2=[U,[-1;-1;-1;-1]]; % Ru = G2 * m2; m2 = [a Cu]'

% solve

m2=(G2'*G2)\G2'*Ru;
a=m2(1);
Cu=m2(2);

Vu=m2(1)*U;

% second part

G3=[[Vu;Eu;Nu],[-1;-1;-1;-1;0;0;0;0;0;0;0;0]]; % dwn = G3 * m3; m3 = [c Cd]'

m3=(G3'*G3)\G3'*dwn;
dwn_mod=G3*m3;
c=m3(2);
Cd=m3(2);

% determine residual signal

non_sym2=dwn-dwn_mod;

% convert back to meaningful units
non_sym2(1:4)=non_sym2(1:4)/(365*10^-2); % [cm/d]
non_sym2(5:end)=non_sym2(5:end)/(10^-6*3000*365); % [urad/d]
up(1:4)=up(1:4)/(365*10^-2); % [cm/d]
up(5:end)=up(5:end)/(10^-6*3000*365); % [urad/d]
dwn(1:4)=dwn(1:4)/(365*10^-2); % [cm/d]
dwn(5:end)=dwn(5:end)/(10^-6*3000*365); % [urad/d]
dwn_mod(1:4)=dwn_mod(1:4)/(365*10^-2); % [cm/d]
dwn_mod(5:end)=dwn_mod(5:end)/(10^-6*3000*365); % [urad/d]
U(1:4)=U(1:4)/(365*10^-2); % [cm/d]

%% plotting

load('calderacoordinates.mat')
calxy=llh2local(caldera_rim',[-130.1,45.9]);
stalocs=llh2local([[AXCC1.lon;AXCC1.lat],[AXEC2.lon;AXEC2.lat],[ASHES.lon;ASHES.lat],...
    [AXID1.lon;AXID1.lat]],[-130.1,45.9]);

% pre
figure(11)
clf
hold on
plot(stalocs(1,:),stalocs(2,:),'^k','markerfacecolor','k','markersize',10)
plot(calxy(1,:),calxy(2,:),'k')

quiver(stalocs(1,:),stalocs(2,:),up(5:8)'*2,up(9:end)'*2,'k','linewidth',2,'autoscale','off')
quiver(stalocs(1,:),stalocs(2,:),zeros(1,4)*20^1,up(1:4)'*20^1,'b','linewidth',2,'autoscale','off')

lbl1={[num2str(abs(round(up(1),3))) ' cm/d'];[num2str(abs(round(up(2),3))) ' cm/d'];...
    [num2str(abs(round(up(3),3))) ' cm/d'];[num2str(abs(round(up(4),3))) ' cm/d']};
lbl2={[num2str(abs(round(sqrt(up(5)^2+up(9)^2),3))) ' urad/d'];...
    [num2str(abs(round(sqrt(up(6)^2+up(10)^2),3))) ' urad/d'];...
    [num2str(abs(round(sqrt(up(7)^2+up(11)^2),3))) ' urad/d'];...
    [num2str(abs(round(sqrt(up(8)^2+up(12)^2),3))) ' urad/d']};

text(stalocs(1,:)-1.5,stalocs(2,:)+0.6,lbl1,'color','b','fontsize',18)
text(stalocs(1,:)-1.5,stalocs(2,:)+0.3,lbl2,'color','k','fontsize',18)

title(['"Pre" ' load_dir])
set(gca,'fontsize',18)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print(['../tiltcompare/inflation_reversal/in_progress_matfiles/' load_dir '/pre_signal'],'-dtiff')

% post observed
figure(12)
clf
hold on
plot(stalocs(1,:),stalocs(2,:),'^k','markerfacecolor','k','markersize',10)
plot(calxy(1,:),calxy(2,:),'k')

quiver(stalocs(1,:),stalocs(2,:),dwn(5:8)'*2,dwn(9:end)'*2,'k','linewidth',2,'autoscale','off')
quiver(stalocs(1,:),stalocs(2,:),zeros(1,4)*20^1,dwn(1:4)'*20^1,'b','linewidth',2,'autoscale','off')

lbl1={[num2str(abs(round(dwn(1),3))) ' cm/d'];[num2str(abs(round(dwn(2),3))) ' cm/d'];...
    [num2str(abs(round(dwn(3),3))) ' cm/d'];[num2str(abs(round(dwn(4),3))) ' cm/d']};
lbl2={[num2str(abs(round(sqrt(dwn(5)^2+dwn(9)^2),3))) ' urad/d'];...
    [num2str(abs(round(sqrt(dwn(6)^2+dwn(10)^2),3))) ' urad/d'];...
    [num2str(abs(round(sqrt(dwn(7)^2+dwn(11)^2),3))) ' urad/d'];...
    [num2str(abs(round(sqrt(dwn(8)^2+dwn(12)^2),3))) ' urad/d']};

text(stalocs(1,:)-1.5,stalocs(2,:)+0.6,lbl1,'color','b','fontsize',18)
text(stalocs(1,:)-1.5,stalocs(2,:)+0.3,lbl2,'color','k','fontsize',18)

title(['"Post" (Observed) ' load_dir])
set(gca,'fontsize',18)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print(['../tiltcompare/inflation_reversal/in_progress_matfiles/' load_dir '/post_modeled_signal'],'-dtiff')

% post modeled
figure(13)
clf
hold on
plot(stalocs(1,:),stalocs(2,:),'^k','markerfacecolor','k','markersize',10)
plot(calxy(1,:),calxy(2,:),'k')

quiver(stalocs(1,:),stalocs(2,:),dwn_mod(5:8)'*2,dwn_mod(9:end)'*2,'k','linewidth',2,'autoscale','off')
quiver(stalocs(1,:),stalocs(2,:),zeros(1,4)*20^1,dwn_mod(1:4)'*20^1,'b','linewidth',2,'autoscale','off')

lbl1={[num2str(abs(round(dwn_mod(1),3))) ' cm/d'];[num2str(abs(round(dwn_mod(2),3))) ' cm/d'];...
    [num2str(abs(round(dwn_mod(3),3))) ' cm/d'];[num2str(abs(round(dwn_mod(4),3))) ' cm/d']};
lbl2={[num2str(abs(round(sqrt(dwn_mod(5)^2+dwn_mod(9)^2),3))) ' urad/d'];...
    [num2str(abs(round(sqrt(dwn_mod(6)^2+dwn_mod(10)^2),3))) ' urad/d'];...
    [num2str(abs(round(sqrt(dwn_mod(7)^2+dwn_mod(11)^2),3))) ' urad/d'];...
    [num2str(abs(round(sqrt(dwn_mod(8)^2+dwn_mod(12)^2),3))) ' urad/d']};

text(stalocs(1,:)-1.5,stalocs(2,:)+0.6,lbl1,'color','b','fontsize',18)
text(stalocs(1,:)-1.5,stalocs(2,:)+0.3,lbl2,'color','k','fontsize',18)

title(['"Post" (Modeled) ' load_dir])
set(gca,'fontsize',18)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print(['../tiltcompare/inflation_reversal/in_progress_matfiles/' load_dir '/post_signal'],'-dtiff')

% residuals
figure(14)
clf
hold on
plot(stalocs(1,:),stalocs(2,:),'^k','markerfacecolor','k','markersize',10)
plot(calxy(1,:),calxy(2,:),'k')

quiver(stalocs(1,:),stalocs(2,:),non_sym2(5:8)'*2,non_sym2(9:end)'*2,'k','linewidth',2,'autoscale','off')
quiver(stalocs(1,:),stalocs(2,:),zeros(1,4)*20^1,non_sym2(1:4)'*20^1,'b','linewidth',2,'autoscale','off')

lbl1={[num2str(abs(round(non_sym2(1),3))) ' cm/d'];[num2str(abs(round(non_sym2(2),3))) ' cm/d'];...
    [num2str(abs(round(non_sym2(3),3))) ' cm/d'];[num2str(abs(round(non_sym2(4),3))) ' cm/d']};
lbl2={[num2str(abs(round(sqrt(non_sym2(5)^2+non_sym2(9)^2),3))) ' urad/d'];...
    [num2str(abs(round(sqrt(non_sym2(6)^2+non_sym2(10)^2),3))) ' urad/d'];...
    [num2str(abs(round(sqrt(non_sym2(7)^2+non_sym2(11)^2),3))) ' urad/d'];...
    [num2str(abs(round(sqrt(non_sym2(8)^2+non_sym2(12)^2),3))) ' urad/d']};

text(stalocs(1,:)-1.5,stalocs(2,:)+0.6,lbl1,'color','b','fontsize',18)
text(stalocs(1,:)-1.5,stalocs(2,:)+0.3,lbl2,'color','k','fontsize',18)

title(['Residual (Observed - Modeled) ' load_dir])
set(gca,'fontsize',18)

fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print(['../tiltcompare/inflation_reversal/in_progress_matfiles/' load_dir '/non_sym_signal'],'-dtiff')

save(['../tiltcompare/inflation_reversal/in_progress_matfiles/' load_dir '/non_sym_signal'],'up','dwn','dwn_mod','non_sym2','U','a','Cu','c','Cd')