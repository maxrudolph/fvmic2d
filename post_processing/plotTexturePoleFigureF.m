function plotTexturePoleFigureF( texture, imark )

% column 1 is theta 2 is phi
tex = [texture(imark).ctheta texture(imark).cphi];

x = cos(tex(:,1)) .* sin(tex(:,2));
y = sin(tex(:,1)) .* sin(tex(:,2));
z = cos(tex(:,2));

idxout = z>0;
idxin = z<0;

t=linspace(0,2*pi,360);
xc = cos(t); yc = sin(t); 
scatter(x(idxout),y(idxout),'r.'), hold on, scatter(x(idxin),y(idxin),'k.'),plot(xc,yc)
