function y1=f1(x,P,T,XH2Obulk) %cpx function

A1 = 1085.7;
A2 = 132.9;
A3 = -5.1;
B1 = 1475;
B2 = 80;
B3 = -3.2;
C1 = 1780;
C2 = 45;
C3 = -2;
r0 = 0.5;
r1 = .08;
beta1 = 1.5;
beta2 = 1.5;
k = 43;
gamma = .75;
x1=12;
x2=1;
delta=0.6;

DH2O = 0.01;

M = .15;
%XH2Obulk=0.1; %units in wt %

Tsolidus = A1 + A2*P + A3*P^2;
Tlhrzliq = B1 + B2*P + B3*P^2;
Fcpxout = M/(r0 + r1*P);
Tliquidus = C1 + C2*P + C3*P^2;
Tcpxout = Fcpxout^(1/beta1)*(Tlhrzliq - Tsolidus) + Tsolidus;

Fcpx = ((T - Tsolidus)/(Tlhrzliq - Tsolidus))^beta1;
Fopx = Fcpxout + (1 - Fcpxout ) *((T - Tcpxout)/(Tliquidus - Tcpxout))^beta2;

XH2Ocpx = XH2Obulk/(DH2O + x*(1 - DH2O));
XH2Osat = x1*P^delta + x2*P;

%XH2Oopx = XH2Obulk/(DH2O + Fhopx*(1 - DH2O));

DTcpx = k*XH2Ocpx^gamma;
%DTopx = k*XH2Oopx^gamma;
if(XH2Ocpx>XH2Osat)
    DTcpx=k*XH2Osat^gamma;
end

y1= x - ((T - (Tsolidus - DTcpx))/(Tlhrzliq - Tsolidus))^beta1;
end