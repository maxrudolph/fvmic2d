% read texture info from ascii file and plot a pole figure
%  close all
% filename = '../output/texture.49.dat';

% tex=load(filename);
%% Read texture from binary markers file
filename = '../output/texture.10000.dat';
texture=readMarkerTextureASCII(filename)

%%
% column 1 is theta 2 is phi
imark = 1;
tex = [texture(imark).ctheta texture(imark).cphi];

x = cos(tex(:,1)) .* sin(tex(:,2));
y = sin(tex(:,1)) .* sin(tex(:,2));
z = cos(tex(:,2));

idxout = z>0;
idxin = z<0;

t=linspace(0,2*pi,360);
xc = cos(t); yc = sin(t);
figure, scatter(x(idxout),y(idxout),'rx'), hold on, scatter(x(idxin),y(idxin),'kx'),plot(xc,yc),axis equal
% figure, plotTexturePoleFigureMTEX_simple( texture )
texture_index( texture.cphi, texture.ctheta )
%second method - make r,theta bins
nr = 10;
nt = 50;

rbins = linspace(0,1,nr+1);
tbins = linspace(-pi,pi,50+1);

Rvals = histc(sqrt(x.*x+y.*y),rbins);
Tvals = histc( atan2(y,x),tbins ); % comes out between -Pi and Pi with -x axis =0
Rvals = Rvals(1:end-1);
Tvals = Tvals(1:end-1);
%figure out patch geometery for each r, theta bin

