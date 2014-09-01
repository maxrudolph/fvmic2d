%plot results from anisotropic convection experiments

filename = 'Markers.0.1600';
clear m g;

m = getBinaryMarkers(filename,1);

%%
mask = m.x>0 & m.y>0;
NX=200;
NY=200;
xmin =0;
xmax = max(m.x);
ymin=0;
ymax = max(m.y);
X = linspace(xmin,xmax,NX);
Y = linspace(ymin,ymax,NY);
[X1 Y1]=meshgrid(X,Y);

%%
g.vx = griddata(m.x(mask),m.y(mask),m.vx(mask),X1,Y1);
g.vy = griddata(m.x(mask),m.y(mask),m.vy(mask),X1,Y1);
g.T = griddata(m.x(mask),m.y(mask),m.T(mask),X1,Y1);

%%
figure, imagesc(X,Y,g.T), axis image, title('T'), colorbar;
figure, imagesc(X,Y,g.vx), axis image, title('vx'), colorbar;
figure, imagesc(X,Y,g.vy), axis image, title('vy'), colorbar;
%%
g.N11 = griddata(m.x(mask),m.y(mask),m.N11(mask),X1,Y1);
g.N22 = griddata(m.x(mask),m.y(mask),m.N22(mask),X1,Y1);
g.N12 = griddata(m.x(mask),m.y(mask),m.N12(mask),X1,Y1);
g.N21 = griddata(m.x(mask),m.y(mask),m.N21(mask),X1,Y1);
figure, imagesc(log10(abs(g.N11))), colorbar, title('N11'), axis image

%%
g.srr = griddata(m.x(mask),m.y(mask),m.srr(mask),X1,Y1);
figure, imagesc(log10(g.srr)), title('srr'), colorbar;