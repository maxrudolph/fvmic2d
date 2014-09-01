%function markers = getGriddedMarkers(fn)
%reads markers from binary file
% close all
%clear markers;
iMonte=0;
iTime=200;
fn=sprintf('/mnt/cluster/markercode/vep/ridges/output/GriddedMarkers.%d.%d',iMonte,iTime);
%fn = sprintf('/Users/max/projects/markercode/vep/ridges/output/GriddedMarkers.%d.%d',iMonte,iTime);
% fn = sprintf('/Users/max/projects/markercode/vep/ridges/output.promising/GriddedMarkers.%d.%d',iMonte,iTime);
% fn = sprintf('/Users/max/projects/markercode/vep/ridges/runs/GriddedMarkers.%d.%d',iMonte,iTime);
% /Users/max/projects/markercode/vep/ridges/runs
% fn = '/tmp/runs/GriddedMarkers.1.200'
fh = fopen(fn,'r');

%get number of markers
nx = fread(fh,1, 'int');
ny = fread(fh,1, 'int');
% test= fread(fh,1, 'int');
%get X

% markers.X = fread(fh,nMark,'double');
% markers.Y = 
% fread(fh,nMark,'double');
markers.rho = zeros(ny,nx);
markers.rho(:) = fread(fh,nx*ny,'double');
markers.vx = zeros(ny,nx);
markers.vx(:) = fread(fh,nx*ny,'double');
markers.vy = zeros(ny,nx);
markers.vy(:) = fread(fh,nx*ny,'double');
markers.vz = zeros(ny,nx);
markers.vz(:) = fread(fh,nx*ny,'double');
markers.T = zeros(ny,nx);
markers.T(:) = fread(fh,nx*ny,'double');
markers.D = zeros(ny,nx);
markers.D(:) = fread(fh,nx*ny,'double');
markers.Ddot = zeros(ny,nx);
markers.Ddot(:) = fread(fh,nx*ny,'double');
markers.p = fread(fh,[ny nx],'double');
markers.eta = fread(fh,[ny nx],'double');

%figure(2);
%plotProfile(markers);

doplot=0;
if(doplot)
    figure(1);
    npx=4;
    npy=2;
    subplot(npy,npx,1)
    fieldToPlot = markers.rho;
    clims = [min(fieldToPlot(fieldToPlot >=100)) max(fieldToPlot(fieldToPlot >= 100))];
    imagesc(fieldToPlot,clims); colorbar;
    title('rho')
    subplot(npy,npx,2)
    fieldToPlot = markers.D;
    clims = [min(fieldToPlot(fieldToPlot ~=-1)) max(fieldToPlot(fieldToPlot ~= -1))];
    if(clims(1) == clims(2)) clims(2) = clims(1)+0.01; end
    imagesc(fieldToPlot,clims); colorbar;
    title('D')
    subplot(npy,npx,3)
    fieldToPlot = markers.Ddot;
    clims = [min(fieldToPlot(fieldToPlot ~=-1)) max(fieldToPlot(fieldToPlot ~= -1))];
    if(clims(1) == clims(2)) clims(2) = clims(1)+.01; end
    imagesc(fieldToPlot,clims); colorbar;
    title('Ddot')
    subplot(npy,npx,4)
    fieldToPlot = markers.vx;
    clims = [min(fieldToPlot(fieldToPlot ~=-1)) max(fieldToPlot(fieldToPlot ~= -1))];
    if(clims(1) == clims(2)) clims(2) = clims(1)+.01; end
    imagesc(fieldToPlot,clims); colorbar;
    title('vx')
    subplot(npy,npx,5)
    fieldToPlot = markers.vy;
    clims = [min(fieldToPlot(fieldToPlot ~=-1)) max(fieldToPlot(fieldToPlot ~= -1))];
    if(clims(1) == clims(2)) clims(2) = clims(1)+.01; end
    imagesc(fieldToPlot,clims); colorbar;
    title('vy')
    subplot(npy,npx,6)
    fieldToPlot = markers.vz;
    clims = [min(fieldToPlot(fieldToPlot ~=-1)) max(fieldToPlot(fieldToPlot ~= -1))];
    if(clims(1) == clims(2)) clims(2) = clims(1)+.01; end
    imagesc(fieldToPlot,clims); colorbar;
    title('vz')
    subplot(npy,npx,7)
    fieldToPlot = markers.T;
    clims = [min(fieldToPlot(fieldToPlot ~=-1)) max(fieldToPlot(fieldToPlot ~= -1))];
    if(clims(1) == clims(2)) clims(2) = clims(1)+.01; end
    imagesc(fieldToPlot,clims); colorbar;
    title('T')
    subplot(npy,npx,8)
    fieldToPlot = markers.eta;
    clims = [min(fieldToPlot(fieldToPlot ~=-1)) max(fieldToPlot(fieldToPlot ~= -1))];
    if(clims(1) == clims(2)) clims(2) = clims(1)+.01; end
    imagesc(fieldToPlot,clims); colorbar;
    title('eta')
end


% markers.VY = fread(fh,nMark,'double');
% markers.VZ = fread(fh,nMark,'double');
% markers.Mat = fread(fh,nMark,'integer*1');
% markers.T = fread(fh,nMark,'double');
% markers.eta = fread(fh,nMark,'double');
% markers.D = fread(fh,nMark,'double');
% markers.exx = fread(fh,nMark,'double');
% markers.exy = fread(fh,nMark,'double');
% markers.sxx = fread(fh,nMark,'double');
% markers.sxy = fread(fh,nMark,'double');
% markers.p = fread(fh,nMark,'double');
% markers.wxy = fread(fh,nMark,'double');
% markers.wxz = fread(fh,nMark,'double');
% markers.wyz = fread(fh,nMark,'double');
% markers.rho = fread(fh,nMark,'double');
fclose(fh);



% [XI,YI] = meshgrid(X,Y);
%X	Y	VX	VY	T	Eta	exx	exy	sxx	sxy	MaterialId
%X	Y	VX	VY	T	Eta	exx	exy	sxx	sxy	p	w MaterialId
% mmat = griddata(markers.X,markers.Y,markers.Mat,XI,YI);
% mexx = griddata(markers.X,markers.Y,markers.exx,XI,YI);
% mexy = griddata(markers.X,markers.Y,markers.exy,XI,YI);
% mD = griddata(markers.X,markers.Y,markers.D,XI,YI);
% mvx = griddata(markers.X,markers.Y,markers.VX,XI,YI);
% mvy = griddata(markers.X,markers.Y,markers.VY,XI,YI);
% mvz = griddata(markers.X,markers.Y,markers.VZ,XI,YI);
% mrho = griddata(markers.X,markers.Y,markers.rho,XI,YI);
% mpr = griddata(markers.X,markers.Y,markers.p,XI,YI);
% mt = griddata(markers.X,markers.Y,markers.T,XI,YI);
% figure;
% nplot=4;
% subplot(1,nplot,1);
% imagesc(mrho)
% axis image ij
% %
% subplot(1,nplot,2);
% imagesc(mvx)
% axis image ij
% colorbar
% %
% subplot(1,nplot,3);
% imagesc(mvy)
% axis image ij
% colorbar
% 
% subplot(1,nplot,4);
% imagesc(mt)
% axis image ij
% colorbar
% 
% 
% title(sprintf('material, iTime=%d',iTime));
% drawnow;