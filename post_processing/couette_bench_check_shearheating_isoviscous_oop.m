% compare numerical solution to Turcotte and Schubert example 6-23: Couette
% flow with visous dissipation. Here, flow is taken to be out-of-plane.

filename = 'loadNodalFields_0_100.petscbin';

loadgrid;

nf = loadNodalFieldsPetscBin(filename);
rho = nf.rho(2,2);
mu = nf.etaS(2,2);

LX=grid.x(end);

xv = grid.x;

% close all

figure, pcolor(grid.x,grid.y,nf.T), colorbar, title('Temperature')
figure, plot(grid.x,nf.T(5,:),'.'), title('Temperature profile')

%This is a solution for couette flow with mu=mu0*T

u1 = 5; %velocity of top plate wrt bottom
T1=10;
T0=1;
k= 5.0;

%%
% analytic_soln = (yv*L-yv.^2)*rho*g/mu;
PrE = mu*u1^2/k/(T1-T0);
h=LX;

tnd = xv/h*(1+PrE/2)-xv.^2/h^2*PrE/2;
soln = (T1-T0)*tnd+T0;

figure, hold on
plot( grid.x(1:end), nf.T(5,:), 'g.')
hold on
plot( xv,soln,'rx')


%% calculate shear stress sxy based on actual solution

NX = length(grid.x);
tauNum = zeros(NX-2,1);
for i=2:NX-1
   tauNum(i-1) = nf.etaS(3,i-1) * ( (nf.vz(3,i)-nf.vz(3,i-1))/(grid.xc(i)-grid.xc(i-1)) );
end
% figure, plot(tauNum)

