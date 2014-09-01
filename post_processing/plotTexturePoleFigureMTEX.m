function plotTexturePoleFigureMTEX( varargin )
texture = varargin{1};
if(nargin == 1 )
    imark = 1;
else
    imark = varargin{2};
end
if(nargin < 3)
    cmax = 800;
else
    cmax = varargin{3};
end

sp = gca; % get current subplot

% column 1 is theta 2 is phi
tex = [texture(imark).ctheta texture(imark).cphi];

x = cos(tex(:,1)) .* sin(tex(:,2));
y = sin(tex(:,1)) .* sin(tex(:,2));
z = cos(tex(:,2));

v = vector3d(y,-x,z);

sx = 0;
sy = 0;
sz = 1;
vs = vector3d(sx,sy,sz);
ss = symmetry('triclinic');
cs = symmetry('hexagonal');
% odf = fibreODF( Miller(v), cs);
odf = fibreODF( Miller(vs),Miller(v),cs,ss,'halfwidth',20*degree);
ftmp = figure; % make a temporary figure
plotpdf(odf,Miller(0,0,0,1),'antipodal')
cmap = colormap;
hold on
plot(v,'MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',1.0)

props = get(gca);

tmp = get(gca,'children');
%copy contents of temporary figure and axes to sp;
copyobj(allchild(gca),sp);
axes(sp);
colormap(cmap);
caxis([0 cmax])

% set(gca,props);
close(ftmp);

% plot(v)
% t=linspace(0,2*pi,360);
% xc = cos(t); yc = sin(t);
% scatter(x(idxout),y(idxout),'r.'), hold on, scatter(x(idxin),y(idxin),'k.'),plot(xc,yc)
