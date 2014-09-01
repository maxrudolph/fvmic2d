%% Read texture from binary markers file
anidir = '/Users/max/Documents/Berkeley/Research with Michael/Convection_Anisotropic_Ice/Data/1mm_75km/'
cd(anidir)
texturefilename = [anidir 'Texture.0.121000'];
markerfilename = [anidir 'Markers.0.121000'];
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


% left side
%% 
idxplot_start = 3;
idxplot = idxplot_start;
figure(2);
clf, imagesc(g.X,g.Y,g.T); hold on;
[xp yp p] = impixel;

while(true)
    if(idxplot>idxplot_start)
    [xp yp p] = impixel;
    end
    
    
    
    % build list of marker indices near xp and yp
    clear midx
    ix = 1;
    %loop over xp
    %mask out markers with x,y -xp,yp > lx/dx,ly.dy
    dxmax = lx/20;
    dymax = ly/20;
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
    
    
    
    %figure for pole figure
    

    % select outer plots in sequence
    % idxm = idx(mask);
    
    figure(idxplot)
    
                
    plotTexturePoleFigureMTEX(texture,midx); set(gca,'Visible','off','DataAspectRatio',[1 1 1]);
    text(1.4,1.4,num2str(idxplot),'FontName','Times','FontSize',16)
                
    p=get(gca,'pos')

    figure(2); hold on;
%     subplot('Position',[spw+sps sph+sps tpw tph]);
%     if(idxplot == idxplot_start)
%         imagesc(g.X,g.Y,g.T), hold on, set(gca,'Visible','Off'); hold all;
%     end
    scatter(xp,yp,'ks','MarkerFaceColor','k')
    %     for ix=1:np
    text(xp+0.02*lx,yp,num2str(idxplot))
    %     end
    %     set(gca,'DataAspectRatio',[1 1 1])
    %     colorbar
    idxplot = idxplot+1;
end



%%
% saveas(gcf,'test.eps','psc2')


