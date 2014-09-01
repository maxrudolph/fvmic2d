%make a movie of convection simulation output
clear, close all
%% Set Output Options
avi_file_name = 'yielding_gk.avi';
avi_fps = 4; % frame rate
aviobj = avifile(avi_file_name,'fps',avi_fps);

%% Generate file list
fileNamePrefix='Markers.0';
files=dir([fileNamePrefix '.*']);
[filelist ltezero] = sort_files(files);
filelist = filelist(~ltezero);
LX=300000;
LY=30000;
NX=600;
NY=300;
nfiles = size(filelist,1);
secondsInYear=365*24*3600;

%% Open a figure window for plotting
fig=figure;

%% Begin iteration over marker files
for ifile=1:nfiles;
markerfilename = filelist(ifile);
m = getBinaryMarkers(markerfilename{1},0);


%% grid temperature data for plotting
lx = LX;
ly = LY;
nx = 700; ny = 350;
X0 = linspace(0,lx,nx);
Y0 = linspace(0,ly,ny);
[X Y] = meshgrid(X0,Y0);
mask = m.x > 0 & m.y > 0;
Tf = TriScatteredInterp(m.x(mask),m.y(mask),m.T(mask));
g.T=Tf(X,Y);
g.X = X0;
g.Y = Y0;
clear mask X Y X0 Y0

%% calculate Nusselt number
%grid is uniformly spaced, so nusselt number = total heat flow/conducted
%heat flow = dtdz_top/dtdz total
dtdz_top = (g.T(4,:)-g.T(3,:))/(g.Y(4)-g.Y(3));
dtdz_top = mean(dtdz_top(~isnan(dtdz_top)));
dtdz_total = (260-100)/LY;
Nu = dtdz_top/dtdz_total;

%% plot temperature image
imagesc(g.X,g.Y,g.T), axis equal tight;
fontname = 'Helvetica';
fontsize=16;
title(sprintf('Elapsed Time %3.f kyr, Nu = %.2f',m.elapsedTime/secondsInYear/1000,Nu),'FontName',fontname,'FontSize',fontsize);
xlabel([num2str(LX/1000) ' km'],'FontName',fontname,'FontSize',16);
ylabel([num2str(LY/1000) ' km'],'FontName',fontname,'FontSize',16);
set(gca,'XTick',[])
set(gca,'YTick',[])
colorbar
caxis([100 260])
drawnow;
F = getframe(fig);
aviobj = addframe(aviobj,F);

end

%% close avi file
aviobj = close(aviobj);