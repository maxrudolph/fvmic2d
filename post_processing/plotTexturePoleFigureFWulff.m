function plotTexturePoleFigureFWulff( texture, imark )

% column 1 is theta 2 is phi
tex = [texture(imark).ctheta texture(imark).cphi];

x = cos(tex(:,1)) .* sin(tex(:,2));
y = sin(tex(:,1)) .* sin(tex(:,2));
z = cos(tex(:,2));

idxout = z>0;
idxin = z<0;

t=linspace(0,2*pi,360);
xc = cos(t); yc = sin(t); 

XO = x(idxout)./(1-z(idxout));
YO = y(idxout)./(1-z(idxout));

XI = x(idxin)./(1+z(idxin));
YI = y(idxin)./(1+z(idxin));
scatter(XO,YO,'r.'), hold on, scatter(XI,YI,'k.'),plot(xc,yc)
