clear all; close all;
xH2O = 0.1;
numP = 201;
numT = 101;
F = zeros(numP, numT);
P = zeros(numP, numT);
T = zeros(numP, numT);
minP = 0.09; %GPa
maxP = 7;
minT = 1200;  %K
maxT = 1600;

Plist = linspace(minP,maxP,numP-1);
Tlist = linspace(minT,maxT,numT-1);

dp = Plist(2)-Plist(1);
dT = Tlist(2)-Tlist(1);
Plist = [Plist Plist(end)+dp];
Tlist = [Tlist Tlist(end)+dT];

for i = 1:numP;
    for j = 1:numT;
        F(i,j) = KatzMelt(Plist(i), Tlist(j), xH2O);
        P(i,j)=Plist(i);
        T(i,j)=Tlist(j);
    end
end

melt_table.P = P(1:end-1,1:end-1);
melt_table.T = T(1:end-1,1:end-1);
melt_table.F = F(1:end-1,1:end-1);
dFdp = zeros(numP-1,numT-1);
dFdT = zeros(numP-1,numT-1);

for i = 1:numP-1
    for j = 1:numT-1
        dFdp(i,j) = (F(i+1,j)-F(i,j))/dp;
        dFdT(i,j) = (F(i,j+1)-F(i,j))/dT;
    end
end
melt_table.dFdP = dFdp;
melt_table.dFdT = dFdT;
save(['my_melt_table_' num2str(xH2O) '.mat'],'melt_table');