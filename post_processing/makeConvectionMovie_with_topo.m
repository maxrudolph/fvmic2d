%make a movie of convection simulation output
clear, close all, fclose all
%% Set Output Options
avi_file_name = 'Titan_gk_30km_0.1mm_aspect2.avi';
avi_fps = 4; % frame rate
aviobj = avifile(avi_file_name,'fps',avi_fps);

%% Generate file list
fileNamePrefix='Markers.0';
files=dir([fileNamePrefix '.*']);
[filelist ltezero] = sort_files(files);
filelist = filelist(~ltezero);
LX=60000;
LY=30000;
NX=600;
NY=300;
nfiles = size(filelist,1);
secondsInYear=365*24*3600;

%% Open a figure window for plotting
fig=figure;

%% Begin iteration over marker files
for ifile=1:2:nfiles;
disp([num2str( (ifile-1)/nfiles*100 ) ' % complete']);
markerfilename = filelist(ifile);
m = getBinaryMarkers(markerfilename{1},0);
nf_filename = ['loadNodalFields_0_' sscanf(markerfilename{1},[fileNamePrefix '.%s']) '.petscbin'];
nf=loadNodalFieldsPetscBin2( nf_filename );

%% grid temperature data for plotting
lx = LX;
ly = LY;
nx = 400; ny = 200;
X0 = linspace(0,lx,nx);
Y0 = linspace(0,ly,ny);
[X Y] = meshgrid(X0,Y0);
mask = m.x > 0 & m.y > 0;
Tf = TriScatteredInterp(m.x(mask),m.y(mask),m.T(mask));
Topof = TriScatteredInterp(m.x(mask),m.y(mask),m.syy(mask));
g.T=Tf(X,Y);
g.topo=Topof(X(:,:),Y(:,:));
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
fontname = 'Helvetica';
fontsize=16;

subplot(2,1,1)
plot(linspace(0,LX,nf.NX)/1000,-nf.syy(2,:)/1.3/930)
axis([0 LX/1000 -30 30])
ylabel(['Uplift (m)'],'FontName',fontname,'FontSize',16);
title(sprintf('Elapsed Time %3.f kyr, Nu = %.2f',m.elapsedTime/secondsInYear/1000,Nu),'FontName',fontname,'FontSize',fontsize);
set(gca,'pos',[0.1 0.66 0.75 0.33-0.1]); set(gca,'XTickLabel',[]);
subplot(2,1,2)
imagesc(g.X,g.Y,g.T), axis equal tight;

xlabel([num2str(LX/1000) ' km'],'FontName',fontname,'FontSize',16);
ylabel([num2str(LY/1000) ' km'],'FontName',fontname,'FontSize',16);
set(gca,'XTick',[])
set(gca,'YTick',[])
colorbar
caxis([100 260])
set(gca,'pos',[0.1 0.1 0.75 0.66-0.1])
drawnow;
F = getframe(fig);
aviobj = addframe(aviobj,F);

end

%% close avi file
aviobj = close(aviobj);