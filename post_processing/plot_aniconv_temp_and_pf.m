%% Read texture from binary markers file
% texturefilename = 'Texture.0.110500';
% markerfilename = 'Markers.0.110500';

 texturefilename = 'Texture.0.121000';
 markerfilename = 'Markers.0.121000';
loadgrid

%% get texture and markers
texture=readMarkerTextureBinary(texturefilename,24);
m = getBinaryMarkers(markerfilename,1);

if(length(m.T) ~= length(texture))
    disp('WARNING: length of texture not equal to number of markers. There may be a problem')
end

%% grid temperature data for plotting
lx = max(grid.x);
ly = max(grid.y);
nx = 600; ny = 300;
X0 = linspace(0,lx,nx);
Y0 = linspace(0,ly,ny);
[X Y] = meshgrid(X0,Y0);
mask = m.x > 0 & m.y > 0;
Tf = TriScatteredInterp(m.x(mask),m.y(mask),m.T(mask));
g.T=Tf(X,Y);
g.X = X0;
g.Y = Y0;
clear mask X Y X0 Y0

%% define points near where we would like to plot temperature, strain rate
npx = 6;
npy = 4;
% left side
xt = linspace(0.0,1.0,npx+1);
xt1 = xt(1:end-1) + diff(xt)/2;
yt = linspace(0.4,1,npy+1);
yt1 = yt(1:end-1) + diff(yt)/2;
X = ones(npy,1) * xt1;
Y = yt1'*ones(1,npx);
mask = false(npy,npx);
mask(1,:) = true;
mask(end,:) = true;
mask(:,1) = true;
mask(:,end) = true;
xp = X*lx;
yp = Y*ly;

%% build list of marker indices near xp and yp
clear midx
for ix = 1:npx;
    for jy = 1:npy
        if( mask(jy,ix) )
            %loop over xp
            %mask out markers with x,y -xp,yp > lx/dx,ly.dy
            dxmax = lx/npx/2;
            dymax = ly/npy/2;
            mm = (abs(m.x-xp(jy,ix)) < dxmax & abs(m.y-yp(jy,ix)) < dymax);
            %now loop over mm, find minimum distance
            minr = Inf;
            minidx = 0;
            
            r = sqrt( (m.x(mm)-xp(jy,ix)).^2 + (m.y(mm)-yp(jy,ix)).^2);
            minr = min(r);
            minidx = find(r  == minr ); % position within mm
            mmidx = find(mm);
            minidx = mmidx(minidx); %global marker number
            midx(jy,ix) = minidx;
        end
    end
end


%% make figure with temperature in middle and pole figures around edges
h=figure,
% subplot(npy,npx,1);
% calculate plot numbers asssociated with center plot
idx = 1:npx*npy;
idx = reshape(idx,[npx npy])';
idxnm = idx(~mask);

%% make figure
%geometry
sps = 0.05; % spacing between subplots
spw = (1.0-sps*(npx-1))/npx;
sph = (1.0-sps*(npy-1))/npy;
tpw = spw*(npx-2)+sps*(npx-2-1);
tph = sph*(npy-2)+sps*(npy-2-1);
flag = 0;

% select outer plots in sequence
% idxm = idx(mask);

for jy=1:npy
    for ix=1:npx
        if( mask(jy,ix) )
            figure(h)
            pos=[(ix-1)*(spw+sps) (npy-jy)*(sph+sps) spw sph];
            hsp = subplot('Position',pos); hold all;            
            %     sp = subplot(npy,npx,idxm(i));
            plotTexturePoleFigureMTEX(texture,midx(jy,ix)); set(gca,'Visible','off','DataAspectRatio',[1 1 1]);            
            p=get(gca,'pos')
        end
    end
end

subplot('Position',[spw+sps sph+sps tpw tph]);
imagesc(g.X,g.Y,g.T), hold on, set(gca,'Visible','Off'); hold all;
scatter(xp(mask),yp(mask),'ks','MarkerFaceColor','k')
set(gca,'DataAspectRatio',[1 1 1])
colorbar
%%
saveas(gcf,'test.eps','psc2')


