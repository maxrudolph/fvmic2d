function result = plotProfileAndDensity(mark)

%get ridge topo from density

%size of domain, for scale
lx = 10000;
ly = 3000;%meters
dy = ly/size(mark.rho,1);
dx = lx/size(mark.rho,2);

t1 = diff(mark.rho);

nx1 = size(mark.rho,2);

maxyidx = 120;

profile = zeros(nx1,1);
for i=1:nx1
%     maxyidx = max(find(mark.rho(:,i) < 100));
    profile(i) = -dy*find(abs(t1(1:maxyidx,i))==max(abs(t1(1:maxyidx,i))),1,'first');
end
xvals = (0+dx/2):dx:(lx-dx/2);

filterwidth = 5;
%smooth the profile
s1 = zeros(1,nx1);
for i=1:nx1
    s1(i) = mean(profile( max(1,(i-filterwidth)) : min(i+filterwidth,nx1) ) );
end
s1 = s1-min(s1);
subplot(2,1,1);
plot(xvals,s1);
subplot(2,1,2);
imagesc(mark.rho);

slope = diff(s1)./diff(xvals);
slope = atan(slope)/pi*180;
% figure; plot(slope);

%extract relief
%define mean boundary elevation
bw=5;
bel=mean(s1([1:bw nx1-bw:nx1]))
relief=max( abs(s1-bel) )

result.s1 =s1;
result.relief = relief;
result.xvals=xvals;