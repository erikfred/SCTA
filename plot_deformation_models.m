% plot_deformation_models.m
%
% Takes output from USGS volcanic deformation code and generates plots of
% uplift/subsidence and tilt for Axial Seamount
%

clear all; close all

%%%%%CONFIG%%%%%
model_name='def_pre_1.xyzuvw'; %text file of model output
%%%END CONFIG%%%

origin=[-130.1,45.9];
mod_loc=['../axial_deformation_models/dmodels/a_matlab_functions/Axial_Deformation/' model_name];
save_loc=['../axial_deformation_models/' model_name(1:end-7)];
load('../calderacoordinates.mat')

if ~exist(save_loc,'dir')
    mkdir(save_loc)
end

%instrument location info
ASHES.lat=45.93363;ASHES.lon=-130.01368;ASHES.depth=1542;
AXCC1.lat=45.954681;AXCC1.lon=-130.008896;AXCC1.depth=-1528;
AXEC2.lat=45.939671;AXEC2.lon=-129.973801;AXEC2.depth=-1519;
AXID1.lat=45.925732;AXID1.lon=-129.977997;AXID1.depth=-1527.5;

caldera_loc=llh2local(caldera_rim',origin)';
crds=llh2local([AXCC1.lon,AXCC1.lat]',origin);
AXCC1.x=crds(1); AXCC1.y=crds(2);
crds=llh2local([AXEC2.lon,AXEC2.lat]',origin);
AXEC2.x=crds(1); AXEC2.y=crds(2);
crds=llh2local([AXID1.lon,AXID1.lat]',origin);
AXID1.x=crds(1); AXID1.y=crds(2);
crds=llh2local([ASHES.lon,ASHES.lat]',origin);
ASHES.x=crds(1); ASHES.y=crds(2);

%import deformation data
if exist([save_loc '/def.mat'],'file')
    %load existing structure
    load([save_loc '/def.mat'])
else
    %load model data
    fid=fopen(mod_loc);
    u=textscan(fid,'%f %f %f %f %f %f','headerlines',2);
    
    U.x=u{1};
    U.y=u{2};
    U.z=u{3};
    U.ux=u{4};
    U.uy=u{5};
    U.uz=u{6};
    
    %grid data
    [U.xg,U.yg]=meshgrid(floor(min(U.x)):0.5:max(U.x),floor(min(U.y)):0.5:max(U.y));
    U.ug=griddata(U.x,U.y,U.uz,U.xg,U.yg);
    
    %calculate tilts
    [U.xgrad,U.ygrad]=gradient(U.ug);
    
    save([save_loc '/def'],'U')
end

%plot vertical deformation field
figure(7)
clf
contourf(U.xg,U.yg,U.ug);
hold on
plot(caldera_loc(:,1),caldera_loc(:,2),'k','linewidth',1)
plot(AXCC1.x,AXCC1.y,'k^','markersize',15,'linewidth',2)
plot(AXID1.x,AXID1.y,'k^','markersize',15,'linewidth',2)
plot(AXEC2.x,AXEC2.y,'k^','markersize',15,'linewidth',2)
plot(ASHES.x,ASHES.y,'k^','markersize',15,'linewidth',2)
axis equal
ch=colorbar;
ylabel(ch,'vertical deformation (m)','fontsize',15)
ylabel('north (km)')
xlabel('east (km)')
set(gca,'fontsize',15)

orient portrait
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print([save_loc '/vertical_deformation'],'-dtiff')

%plot tilt directions
figure(8)
clf
hold on
plot(caldera_loc(:,1),caldera_loc(:,2),'k','linewidth',1)
plot(AXCC1.x,AXCC1.y,'k^','markersize',15,'linewidth',2)
plot(AXID1.x,AXID1.y,'k^','markersize',15,'linewidth',2)
plot(AXEC2.x,AXEC2.y,'k^','markersize',15,'linewidth',2)
plot(ASHES.x,ASHES.y,'k^','markersize',15,'linewidth',2)
quiver(U.xg,U.yg,-U.xgrad,-U.ygrad,'k')
axis equal
ylabel('north (km)')
xlabel('east (km)')
set(gca,'fontsize',15)
box on

orient portrait
fh=gcf;
fh.PaperUnits='inches';
fh.PaperPosition=[0 0 8.5 11];
print([save_loc '/tilt'],'-dtiff')