%% Read texture from binary markers file
anidir = '/Users/max/Documents/Berkeley/Research with Michael/Convection_Anisotropic_Ice/Data/1mm_75km/'
cd(anidir)
texturefilename = [anidir 'Texture.0.135500'];
markerfilename = [anidir 'Markers.0.135500'];
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
figure, imagesc(g.X,g.Y,g.T);
npx = 6;
npy = 4;
% left side
[xp yp p] = impixel;

np = length(xp);
if( np> 2*npx + 2*(npy-2) )
    np = 2*npx + 2*(npy-2);
end
close
%%
mask = zeros(npy,npx);
mask(1,:) = 1;
mask(:,1) = 1;
mask(end,:) = 1;
mask(:,end) = 1;

%% build list of marker indices near xp and yp
clear midx
for ix = 1:np;
    %loop over xp
    %mask out markers with x,y -xp,yp > lx/dx,ly.dy
    dxmax = lx/npx/2;
    dymax = ly/npy/2;
    mm = (abs(m.x-xp(ix)) < dxmax & abs(m.y-yp(ix)) < dymax);
    %now loop over mm, find minimum distance
    minr = Inf;
    minidx = 0;
    
    r = sqrt( (m.x(mm)-xp(ix)).^2 + (m.y(mm)-yp(ix)).^2);
    minr = min(r);
    minidx = find(r  == minr ); % position within mm
    mmidx = find(mm);
    minidx = mmidx(minidx); %global marker number
    midx(ix) = minidx;
end


%% make figure with temperature in middle and pole figures around edges
h=figure,
% subplot(npy,npx,1);
% calculate plot numbers asssociated with center plot
idx = 1:npx*npy;
idx = reshape(idx,[npx npy])';


% make figure
%geometry
sps = 0.05; % spacing between subplots
spw = (1.0-sps*(npx-1))/npx;
sph = (1.0-sps*(npy-1))/npy;
tpw = spw*(npx-2)+sps*(npx-2-1);
tph = sph*(npy-2)+sps*(npy-2-1);
flag = 0;

% select outer plots in sequence
% idxm = idx(mask);

ix1 = 1;
for ix=1:npx
    for jy=1:npy
        if(mask(jy,ix))
            figure(h)
            pos=[(ix-1)*(spw+sps) (npy-jy)*(sph+sps) spw sph];
            hsp = subplot('Position',pos); hold all;
            %     sp = subplot(npy,npx,idxm(i));
            if(ix1 <= np )
                plotTexturePoleFigureMTEX(texture,midx(ix1)); set(gca,'Visible','off','DataAspectRatio',[1 1 1]);
                text(1.4,1.4,num2str(ix1),'FontName','Times','FontSize',16)
            end           
            p=get(gca,'pos')
            ix1 = ix1+1;
        end
    end
end

subplot('Position',[spw+sps sph+sps tpw tph]);
imagesc(g.X,g.Y,g.T), hold on, set(gca,'Visible','Off'); hold all;
scatter(xp(1:np),yp(1:np),'ks','MarkerFaceColor','k')
for ix=1:np    
    text(xp(ix)+0.02*lx,yp(ix),num2str(ix))
end
set(gca,'DataAspectRatio',[1 1 1])
colorbar





%%
saveas(gcf,'test.eps','psc2')


