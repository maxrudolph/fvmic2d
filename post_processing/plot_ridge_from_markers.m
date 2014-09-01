% plot some ridge results from the markers

filename = 'Markers.0.0.gridded';

g = getGriddedMarkers2(filename,0);

NY = g.ny;
NX = g.nx;

loadgrid;
LX = max(grid.x);
LY = max(grid.y);

X = linspace(0,LX,NX);
Y = linspace(0,LY,NY);
%%
% figure, imagesc(X,Y,g.vx'), axis image, title('vx'), colorbar
% figure, imagesc(X,Y,g.vy'), axis image, title('vy'), colorbar
% figure, imagesc(X,Y,g.vz'), axis image, title('vz'), colorbar
% figure, imagesc(X,Y,g.D'), axis image, title('D'), colorbar
% figure, imagesc(X,Y,g.Ddot'), axis image, title('Ddot'), colorbar
% 
% figure, imagesc(X,Y,g.T'), axis image, title('T'), colorbar
figure, imagesc(X,Y,g.p'), axis image, title('p'), colorbar
figure, imagesc(X,Y,log10(g.eta')), axis image, title('eta'), colorbar
% figure, imagesc(X,Y,g.rho'), axis image, title('rho'), colorbar
% figure, imagesc(X,Y,g.rhodot'), axis image, title('rhodot'), colorbar
% figure, imagesc(X,Y,g.sii), axis image, title('sii'), colorbar

%% extract ridge topography
tmp = g.rho>10;
[A I] = max(diff(tmp'));
dy = diff(Y(1:2));
topo = -I*dy;
slope = 180*atan(diff(smooth(topo(:),10))./diff(X(:)))/pi;
figure,subplot(2,1,1), plot(X, topo), axis equal
subplot(2,1,2), plot(X(1:end-1), slope)