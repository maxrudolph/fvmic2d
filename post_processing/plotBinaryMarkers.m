clear
filename = './Markers.-5.-5';
texture =0;

%open file for reading
fh = fopen(filename,'r');

%first get number of markers
clear m;
nMark = fread(fh,1,'int32');
elapsedTime = fread(fh,1,'double')
%read marker fields
m.cpu = fread(fh,nMark,'int32');
m.cellX = fread(fh,nMark,'int32');
m.cellY = fread(fh,nMark,'int32');

m.x = fread(fh,nMark,'double');
m.y = fread(fh,nMark,'double');
m.z = fread(fh,nMark,'double');
m.vx = fread(fh,nMark,'double');
m.vy = fread(fh,nMark,'double');
m.vz = fread(fh,nMark,'double');
m.Mat = fread(fh,nMark,'int32');
m.T = fread(fh,nMark,'double');
m.Tdot = fread(fh,nMark,'double');
m.eta = fread(fh,nMark,'double');
if(texture)
    m.N11 = fread(fh,nMark,'double');
    m.N12 = fread(fh,nMark,'double');
    m.N21 = fread(fh,nMark,'double');
    m.N22 = fread(fh,nMark,'double');
end

m.D = fread(fh,nMark,'double');
m.Ddot = fread(fh,nMark,'double');

m.exx = fread(fh,nMark,'double');
m.exy = fread(fh,nMark,'double');
m.Eii = fread(fh,nMark,'double');

m.sxx = fread(fh,nMark,'double');
m.syy = fread(fh,nMark,'double');
m.szz = fread(fh,nMark,'double');
m.sxz = fread(fh,nMark,'double');
m.syz = fread(fh,nMark,'double');
m.sxy = fread(fh,nMark,'double');
m.p = fread(fh,nMark,'double');
m.rho = fread(fh,nMark,'double');
m.rhodot = fread(fh,nMark,'double');


if(texture)
m.res = fread(fh,nMark,'double');
m.srr = fread(fh,nMark,'double');
end
fclose(fh);
mask = m.x > 0 & m.y > 0;

%%
N=hist3([m.x(mask) m.y(mask)],[70 70]);
figure, imagesc(N), colorbar, title('Nspace')
nx1=max(m.cellX);
ny1=max(m.cellY);
Ncell=zeros(nx1+1,ny1+1);
for ix=0:nx1
    for jy=0:ny1
        Ncell(ix+1,jy+1)=sum( m.cellX == ix & m.cellY == jy);
    end
end
figure, imagesc(Ncell), colorbar, title('Ncell');

%% plot the temperature field
nx = 200;
ny = 100;

xmax = max(m.x);
ymax = max(m.y);
X0 = linspace(0,xmax,nx);
Y0 = linspace(0,ymax,ny);
[X Y] = meshgrid(X0,Y0);
 matf = TriScatteredInterp(m.x(mask),m.y(mask),m.Mat(mask));
 g.mat=matf(X,Y);
 rhof = TriScatteredInterp(m.x,m.y,m.rho); 
 g.rho = rhof(X,Y);
% g.T = griddata(m.x(mask),m.y(mask),m.T(mask),X,Y);
%g.rho = griddata(m.x(mask),m.y(mask),m.rho(mask),X,Y);
g.vx = griddata(m.x(mask),m.y(mask),m.vx(mask),X,Y);
g.vy = griddata(m.x(mask),m.y(mask),m.vy(mask),X,Y);
g.v = (g.vx.^2 + g.vy.^2).^(1/2);
 Df = TriScatteredInterp(m.x,m.y,m.D);
 g.D = Df(X,Y);


%%
g.sxx = griddata(m.x(mask),m.y(mask),m.sxx(mask),X,Y);
g.sxy = griddata(m.x(mask),m.y(mask),m.sxy(mask),X,Y);
% g.sxz = griddata(m.x(mask),m.y(mask),m.sxz(mask),X,Y);
g.exx = griddata(m.x(mask),m.y(mask),m.exx(mask),X,Y);
g.exy = griddata(m.x(mask),m.y(mask),m.exy(mask),X,Y);
g.p = griddata(m.x(mask),m.y(mask),m.p(mask),X,Y);

if(texture)
    g.res = griddata(m.x(mask),m.y(mask),m.res(mask),X,Y);  
    g.srr = griddata(m.x(mask),m.y(mask),m.srr(mask),X,Y);
g.N11 = griddata(m.x(mask),m.y(mask),m.N11(mask),X,Y);
g.N12 = griddata(m.x(mask),m.y(mask),m.N12(mask),X,Y);
g.N21 = griddata(m.x(mask),m.y(mask),m.N21(mask),X,Y);
g.N22 = griddata(m.x(mask),m.y(mask),m.N22(mask),X,Y);
end
[X1 Y1] = meshgrid(linspace(0,xmax,15),linspace(0,ymax,15));
%% plot streamlines
figure, streamline(X,Y,g.vx,g.vy,X1,Y1)