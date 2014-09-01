%%Make plots comparing isotropic and anisotropic results

LX = 150000;
LY = 75000;
t_profile_depth = 35000;

%% Read texture from binary markers file
% texdir = '/work/dogmatix/max/convection_runs/ani_from_gk_runs/1mm/75km/';
% % atexturefilename = '/work/dogmatix/max/convection_runs/ani_from_gk_runs/1mm/75km/Texture.0.103500';
% amarkerfilename = '/work/dogmatix/max/convection_runs/ani_from_gk_runs/1mm/75km/Markers.0.114000';
% loadgrid
% 
% isodir = '/work/dogmatix/max/convection_runs/gk_spinup_runs/1mm/75km_2/'
% isodir='/work/dogmatix/max/convection_runs/gk_spinup_runs/1mm/75km_2_lowres/'

%% for 75 km case
anidir = '/Users/max/Documents/Berkeley/Research with Michael/Convection_Anisotropic_Ice/Data/1mm_75km/'
amarkerfilename = [anidir 'Markers.0.135500'];
isodir = '/Users/max/Documents/Berkeley/Research with Michael/Convection_Anisotropic_Ice/Data/gk_spinup_runs/1.0_75km/';
imarkerfilename = [isodir 'Markers.0.88000'];

%% for 15km case 
% anidir = '/Users/max/Documents/Berkeley/Research with Michael/Convection_Anisotropic_Ice/Data/0.1mm_15km/'
% isodir = '/Users/max/Documents/Berkeley/Research with Michael/Convection_Anisotropic_Ice/Data/gk_spinup_runs/0.1_15km/';
% amarkerfilename = [anidir 'Markers.0.117500'];
% imarkerfilename = [isodir 'Markers.0.94000'];

%% get texture and markers
% texture=readMarkerTextureBinary(atexturefilename,24);
ma = getBinaryMarkers(amarkerfilename,1);
mi = getBinaryMarkers(imarkerfilename,0);


%% grid temperature data for plotting
lx = LX;
ly = LY;
nx = 400; ny = 100;
X0 = linspace(0,lx,nx);
Y0 = linspace(0,ly,ny);
[X Y] = meshgrid(X0,Y0);
mask = ma.x > 0 & ma.y > 0;
Tfa = TriScatteredInterp(ma.x(mask),ma.y(mask),ma.T(mask),'linear');
syyfa = TriScatteredInterp(ma.x(mask),ma.y(mask),ma.syy(mask));
vyfa = TriScatteredInterp(ma.x(mask),ma.y(mask),ma.vy(mask));
mask = mi.x > 0 & mi.y > 0;
Tfi = TriScatteredInterp(mi.x(mask),mi.y(mask),mi.T(mask),'linear');
syyfi = TriScatteredInterp(mi.x(mask),mi.y(mask),mi.syy(mask));
vyfi = TriScatteredInterp(mi.x(mask),mi.y(mask),mi.vy(mask));

ga.T=Tfa(X,Y);
ga.syy=syyfa(X,Y);
ga.vy=vyfa(X,Y);
ga.X = X0;
ga.Y = Y0;
gi.T = Tfi(X,Y);
gi.syy = syyfi(X,Y);
gi.vy=vyfi(X,Y);
gi.X = X0;
gi.Y = Y0;
clear mask X Y 

%% Plots
depth_slice = 65;
secondsInYear = 365*24*3600;
figure, plot(X0/1000,ga.T(depth_slice,:)), hold on, plot(X0/1000,gi.T(depth_slice,:),'r');
ylabel('T(K)','FontSize',16); xlabel('km','FontSize',16);
title(sprintf('Temperature Profile, %.0f km depth',Y0(depth_slice)/1000),'FontSize',16);
legend(gca,'Anisotropic','Isotropic');

secondsInYear = 365*24*3600;
figure, plot(X0/1000,-ga.syy(2,:)/1.3/930), hold on, plot(X0/1000,-gi.syy(2,:)/1.3/930,'r');
ylabel('m','FontSize',16); xlabel('km','FontSize',16);
title(sprintf('Dynamic Topography'),'FontSize',16);
legend(gca,'Anisotropic','Isotropic','Location','SouthEast');

figure, plot(X0/1000,-100*ga.vy(depth_slice,:)*secondsInYear), hold on, plot(X0/1000,-100*gi.vy(depth_slice,:)*secondsInYear,'r')
title(sprintf('Velocity Profile, %.0f km depth',Y0(depth_slice)/1000),'FontSize',16);
ylabel('Velocity (cm/yr)','FontSize',16);
xlabel('km','FontSize',16);
legend(gca,'Anisotropic','Isotropic','Location','SouthEast');

%% Plot two temperature fields
figure, subplot(2,1,1);
imagesc(gi.X,gi.Y,gi.T); caxis([100 260]); colorbar;
axis equal tight
set(gca,'Visible','Off');
subplot(2,1,2);
imagesc(ga.X,ga.Y,ga.T); caxis([100 260]); colorbar;
set(gca,'Visible','Off');

axis equal tight
