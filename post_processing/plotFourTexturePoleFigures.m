% read texture info from ascii file and plot a pole figure
% close all
% filename = '../output/texture.49.dat';

% tex=load(filename);
%% Read texture from binary markers file

step = [0 1000 2000 3000 4000 4980];
nplot = 5;
figure
for i=1:nplot
    filename = ['texture.' num2str(step(i)) '.dat'];
% filename = 'texture.4000.dat';
texture=readMarkerTextureASCII(filename,8)

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
subplot(1,nplot,i);
scatter(x(idxout),y(idxout),'r.'), hold on, scatter(x(idxin),y(idxin),'k.'),plot(xc,yc)
set(gca,'Visible','off','DataAspectRatio',[1 1 1]);
h=text(0,-1.5,sprintf('E_{tot}=%1.1e',texture.Eii),'HorizontalAlignment','Center')
set(h,'FontName','Times')
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
end
