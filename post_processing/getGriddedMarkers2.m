function markers = getGriddedMarkers2( fn , texture )
%reads markers from binary file


% clear markers;
%  fn =
% fn=sprintf('/Users/max/projects/markercode/vep/ridges/output/GriddedMarkers.%d.%d',iMonte,iTime);
%  fn =
%  sprintf('/Users/max/projects/markercode/vep/ridges/output/GriddedMarkers.%d.%d',iMonte,iTime);
% fn='/Users/max/projects/markercode/vep/ridges/output/GriddedMarkers.10.4000';
% fn= 'otest';
% fn = sprintf('/Users/max/projects/markercode/vep/ridges/output.promising/GriddedMarkers.%d.%d',iMonte,iTime);
% fn = sprintf('/Users/max/projects/markercode/vep/ridges/runs/GriddedMarkers.%d.%d',iMonte,iTime);
% /Users/max/projects/markercode/vep/ridges/runs
% fn = '/tmp/runs/GriddedMarkers.1.200'
fh = fopen(fn,'r');

%get number of markers
nx = fread(fh,1, 'int');
ny = fread(fh,1, 'int');
markers.nx = nx;
markers.ny = ny;

% test= fread(fh,1, 'int');
%get X
markers.elapsedTime = fread(fh,1,'double');
% markers.X = fread(fh,nMark,'double');
% markers.Y = 
% fread(fh,nMark,'double');
% markers.rho = zeros(ny,nx);
% markers.rho(:) = fread(fh,nx*ny,'double');
markers.Z = fread(fh,[nx ny],'double');
markers.vx = zeros(nx,ny);
markers.vx(:) = fread(fh,nx*ny,'double');
markers.vy = zeros(nx,ny);
markers.vy(:) = fread(fh,nx*ny,'double');
markers.vz = zeros(nx,ny);
markers.vz(:) = fread(fh,nx*ny,'double');
markers.mat = zeros(nx,ny);
markers.mat(:) = fread(fh,nx*ny,'double');

markers.T = zeros(nx,ny);
markers.T(:) = fread(fh,nx*ny,'double');
markers.Tdot = zeros(nx,ny);
markers.Tdot(:) = fread(fh,nx*ny,'double');
markers.eta = fread(fh,[nx ny],'double');
markers.mu  = fread(fh,[nx ny],'double');
%texture goes here
if(texture)
markers.N11 = fread(fh,[nx ny],'double');
markers.N12 = fread(fh,[nx ny],'double');
markers.N21 = fread(fh,[nx ny],'double');
markers.N22 = fread(fh,[nx ny],'double');
end
markers.D = zeros(nx,ny);
markers.D(:) = fread(fh,nx*ny,'double');
markers.Ddot = zeros(nx,ny);
markers.Ddot(:) = fread(fh,nx*ny,'double');
markers.exx = fread(fh,[nx ny],'double');
markers.exy = fread(fh,[nx ny],'double');
markers.Eii = fread(fh,[nx ny],'double');
markers.sxx = fread(fh,[nx ny],'double');
markers.syy = fread(fh,[nx ny],'double');
markers.szz = fread(fh,[nx ny],'double');
markers.sxz = fread(fh,[nx ny],'double');
markers.syz = fread(fh,[nx ny],'double');
markers.sxy = fread(fh,[nx ny],'double');

markers.sii = sqrt(0.5*(markers.sxx.^2+markers.syy.^2+markers.szz.^2)+markers.sxy.^2+markers.syz.^2 + markers.sxz.^2);

markers.p = fread(fh,[nx ny],'double');
markers.rho = fread(fh,[nx ny],'double');
markers.rhodot = fread(fh,[nx ny],'double');
markers.cpu = fread(fh,[nx ny],'double');
% markers.eta = fread(fh,[ny nx],'double');
% figure(1);
% plotProfile(markers);
% 
% LX = 10000;
% LY = 3000;
% X = linspace(0,LX,nx);
% Y = linspace(0,LY,ny);
% figure(10);
% npx=1;
% npy=8;
% 
% mask = (markers.rho > 10);
% 
% %special colormap that has white at very end
% 
% 
% % subplot(npy,npx,1)
% figure(1)
% fieldToPlot = markers.rho;
% clims = [min(fieldToPlot(fieldToPlot >=100)) max(fieldToPlot(fieldToPlot >= 100))];
% % imagesc(X,Y,fieldToPlot,'AlphaData',mask,clims); colorbar; axis equal image;
% plotWithMask(X,Y,fieldToPlot,mask);  axis equal image;
% title('rho')
% % subplot(npy,npx,2)
% figure(2)
% fieldToPlot = markers.D;
% clims = [min(fieldToPlot(fieldToPlot ~=-1)) max(fieldToPlot(fieldToPlot ~= -1))];
% if(clims(1) == clims(2)) clims(2) = clims(1)+0.01; end
% plotWithMask(X,Y,fieldToPlot,mask);  axis equal image;
% 
% title('D')
% % subplot(npy,npx,3)
% figure(3)
% fieldToPlot = markers.Ddot;
% clims = [min(fieldToPlot(fieldToPlot ~=-1)) max(fieldToPlot(fieldToPlot ~= -1))];
% if(clims(1) == clims(2)) clims(2) = clims(1)+.01; end
% plotWithMask(X,Y,fieldToPlot,mask);  axis equal image;
% title('Ddot')
% % subplot(npy,npx,4)
% figure(4)
% fieldToPlot = markers.vx;
% clims = [min(fieldToPlot(fieldToPlot ~=-1)) max(fieldToPlot(fieldToPlot ~= -1))];
% if(clims(1) == clims(2)) clims(2) = clims(1)+.01; end
% plotWithMask(X,Y,fieldToPlot,mask);  axis equal image;
% title('vx')
% % subplot(npy,npx,5)
% figure(5)
% fieldToPlot = markers.vy;
% clims = [min(fieldToPlot(fieldToPlot ~=-1)) max(fieldToPlot(fieldToPlot ~= -1))];
% if(clims(1) == clims(2)) clims(2) = clims(1)+.01; end
% plotWithMask(X,Y,fieldToPlot,mask);  axis equal image;
% title('vy')
% % subplot(npy,npx,6)
% figure(6)
% fieldToPlot = markers.vz;
% clims = [min(fieldToPlot(fieldToPlot ~=-1)) max(fieldToPlot(fieldToPlot ~= -1))];
% if(clims(1) == clims(2)) clims(2) = clims(1)+.01; end
% plotWithMask(X,Y,fieldToPlot,mask);  axis equal image;
% title('vz')
% % subplot(npy,npx,7)
% figure(7)
% fieldToPlot = markers.T;
% clims = [min(fieldToPlot(fieldToPlot ~=-1)) max(fieldToPlot(fieldToPlot ~= -1))];
% if(clims(1) == clims(2)) clims(2) = clims(1)+.01; end
% plotWithMask(X,Y,fieldToPlot,mask);  axis equal image;
% title('T')
% % subplot(npy,npx,8)
% figure(8)
% fieldToPlot = markers.eta;
% clims = [min(fieldToPlot(fieldToPlot ~=-1)) max(fieldToPlot(fieldToPlot ~= -1))];
% if(clims(1) == clims(2)) clims(2) = clims(1)+.01; end
% plotWithMask(X,Y,fieldToPlot,mask);  axis equal image;
% title('eta')
% drawnow;
% 
% ofn = sprintf('eta.%d.%d.eps',iMonte,iTime);
% print(gcf,'-depsc2','-painters',ofn); 
% 
% % markers.VY = fread(fh,nMark,'double');
% % markers.VZ = fread(fh,nMark,'double');
% % markers.Mat = fread(fh,nMark,'integer*1');
% % markers.T = fread(fh,nMark,'double');
% % markers.eta = fread(fh,nMark,'double');
% % markers.D = fread(fh,nMark,'double');
% % markers.exx = fread(fh,nMark,'double');
% % markers.exy = fread(fh,nMark,'double');
% % markers.sxx = fread(fh,nMark,'double');
% % markers.sxy = fread(fh,nMark,'double');
% % markers.p = fread(fh,nMark,'double');
% % markers.wxy = fread(fh,nMark,'double');
% % markers.wxz = fread(fh,nMark,'double');
% % markers.wyz = fread(fh,nMark,'double');
% % markers.rho = fread(fh,nMark,'double');
fclose(fh);
% 
% 
% 
% % [XI,YI] = meshgrid(X,Y);
% %X	Y	VX	VY	T	Eta	exx	exy	sxx	sxy	MaterialId
% %X	Y	VX	VY	T	Eta	exx	exy	sxx	sxy	p	w MaterialId
% % mmat = griddata(markers.X,markers.Y,markers.Mat,XI,YI);
% % mexx = griddata(markers.X,markers.Y,markers.exx,XI,YI);
% % mexy = griddata(markers.X,markers.Y,markers.exy,XI,YI);
% % mD = griddata(markers.X,markers.Y,markers.D,XI,YI);
% % mvx = griddata(markers.X,markers.Y,markers.VX,XI,YI);
% % mvy = griddata(markers.X,markers.Y,markers.VY,XI,YI);
% % mvz = griddata(markers.X,markers.Y,markers.VZ,XI,YI);
% % mrho = griddata(markers.X,markers.Y,markers.rho,XI,YI);
% % mpr = griddata(markers.X,markers.Y,markers.p,XI,YI);
% % mt = griddata(markers.X,markers.Y,markers.T,XI,YI);
% % figure;
% % nplot=4;
% % subplot(1,nplot,1);
% % imagesc(mrho)
% % axis image ij
% % %
% % subplot(1,nplot,2);
% % imagesc(mvx)
% % axis image ij
% % colorbar
% % %
% % subplot(1,nplot,3);
% % imagesc(mvy)
% % axis image ij
% % colorbar
% % 
% % subplot(1,nplot,4);
% % imagesc(mt)
% % axis image ij
% % colorbar
% % 
% % 
% % title(sprintf('material, iTime=%d',iTime));
% drawnow;
