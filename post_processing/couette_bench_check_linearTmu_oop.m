% compare benchmark poiseuille flow solution with analytic solution

filename = 'loadNodalFields_0_0.petscbin';

loadgrid;

nf = loadNodalFieldsPetscBin(filename);
rho = nf.rho(2,2);
mu = nf.etaS(2,2);
g=10;
LX=grid.x(end);

xv = grid.xc(2:end);

% close all

figure, pcolor(grid.x,grid.y,nf.T), colorbar, title('Temperature')
figure, plot(grid.x,nf.T(5,:),'.'), title('Temperature profile')

%This is a solution for couette flow with mu=mu0*T

u1 = 1e-2; %velocity of top plate wrt bottom
QonR = 5000;
T1=1000;
T0=1;
mu0=1.0;

%%
% analytic_soln = (yv*L-yv.^2)*rho*g/mu;
h=LX;
tau=-(((mu0*T0 - mu0*T1)*u1)/(h*(log(h*mu0*T0) - log(h*mu0*T1))));

soln=(h*tau*log(h*mu0*T1))/(mu0*(T0 - T1)) -  (h*tau*log(mu0*(h*T0 - T0*xv + T1*xv)))/(mu0*T0 - mu0*T1)

figure, hold on
plot( grid.xc(2:end), nf.vz(5,2:end), 'g.')
hold on
plot( xv,soln,'rx')


%% calculate shear stress sxz based on actual solution

NX = length(grid.x);
tauNum = zeros(NX-2,1);
for i=2:NX-1
   tauNum(i-1) = nf.etaS(3,i-1) * ( (nf.vz(3,i)-nf.vz(3,i-1))/(grid.xc(i)-grid.xc(i-1)) );
end
figure, plot(tauNum)

