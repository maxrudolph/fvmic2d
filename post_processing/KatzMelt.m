function meltF = KatzMelt(P,T,xBulkH2O)


% P = 1;
% T = linspace(900,1400,501);
% xBulkH2O = [0, 0.02, 0.05, 0.1, 0.3];



A1 = 1085.7;  % Tsolid parameters
A2 = 132.9;
A3 = -5.1;
B1 = 1475;  % TliqidLherz parameters
B2 = 80;
B3 = -3.2;
C1 = 1780;  %Tliquid parameters
C2 = 45;
C3 = -2;
k = 43;  % other parameters
gamma = 0.75;
DH2O = 0.01;
beta1 = 1.50;
beta2 = 1.50;
M = .17;  % weight fraction of CPX in solid periodotite
r0 = 0.5;
r1 = 0.08;
x1 = 12;
x2 = 1;
lambda = 0.6;


Fcpxout = M/(r0 + r1*P);
% xSatH2O = x1*P^lambda + x2*P;
% for h = 1:length(xBulkH2O);


j = 1;
for i = 1:length(T);
    %             if xBulkH2O(h) >= 0.3;
    %                 meltF(j) = bisection3(0,1,1e-4,20,@(F) ((T(i) - (A1 + A2*P + A3*P^2 - k*((xBulkH2O(h))/(DH2O + F*(1 - DH2O)))^gamma))...
    %                     /(B1 + B2*P + B3*P^2 - (A1 + A2*P + A3*P^2)))^beta1 - F);
    %                 if meltF(j) > Fcpxout
    %                     deltaT = k*xSatH2O^gamma;  % temperature decrease in solidus due to xH2O in the melt
    %                     Tsolid = A1 + A2*P + A3*P^2 - deltaT;  % hydrous solidus
    %                     TliquidLherz = B1 + B2*P + B3*P^2 - deltaT;  % hydrous lherzolite liquidus
    %                     Tcpxout = Fcpxout^(1/beta1)*(TliquidLherz - Tsolid) + Tsolid;
    %                     Tliquid = C1 + C2*P +C3*P^2 - deltaT;
    %                     Fopx = Fcpxout + (1 - Fcpxout)*((T(i)-Tcpxout)/(Tliquid - Tcpxout))^beta2;
    %                     meltF(j) = Fopx;
    %                 end
    %             else
    meltF(j) = bisection3(0,1,1e-4,20,@(F) ((T(i) - (A1 + A2*P + A3*P^2 - k*((xBulkH2O)/(DH2O + F*(1 - DH2O)))^gamma))...
        /(B1 + B2*P + B3*P^2 - (A1 + A2*P + A3*P^2)))^beta1 - F);
    if meltF(j) > Fcpxout
        deltaT = k*xBulkH2O^gamma;  % temperature decrease in solidus due to xH2O in the melt
        Tsolid = A1 + A2*P + A3*P^2 - deltaT;  % hydrous solidus
        TliquidLherz = B1 + B2*P + B3*P^2 - deltaT;  % hydrous lherzolite liquidus
        Tcpxout = Fcpxout^(1/beta1)*(TliquidLherz - Tsolid) + Tsolid;
        Tliquid = C1 + C2*P +C3*P^2 - deltaT;
        Fopx = Fcpxout + (1 - Fcpxout)*((T(i)-Tcpxout)/(Tliquid - Tcpxout))^beta2;
        meltF(j) = Fopx;
    end
    j = j + 1;
end
end
%     plot(T,meltF,'LineWidth',3); hold on;

