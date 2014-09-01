% compare benchmark poiseuille flow solution with analytic solution

filename = 'loadNodalFields_0_5.petscbin';

loadgrid;

nf = loadNodalFieldsPetscBin(filename);
rho = nf.rho(2,2);
mu = nf.etaS(2,2);
g=10;
LY=grid.y(end);

yv = grid.yc(2:end);

% close all

figure, pcolor(grid.x,grid.y,nf.T), colorbar, title('Temperature')
figure, plot(grid.y,nf.T(:,1),'.'), title('Temperature profile')

%This is a solution for couette flow with mu=mu0*T

u1 = 1e-2; %velocity of top plate wrt bottom
QonR = 5000;
T1=100;
T0=1;
mu0=1.0;

%%
% analytic_soln = (yv*L-yv.^2)*rho*g/mu;
h=LY;
tau=-(((mu0*T0 - mu0*T1)*u1)/(h*(log(h*mu0*T0) - log(h*mu0*T1))));

soln=(h*tau*log(h*mu0*T1))/(mu0*(T0 - T1)) -  (h*tau*log(mu0*(h*T0 - T0*yv + T1*yv)))/(mu0*T0 - mu0*T1)

figure, hold on
plot( grid.yc(2:end), nf.vx(1:end-1,10), 'g.')
hold on
plot( yv,soln,'rx')


%% calculate shear stress sxy based on actual solution

NY = length(grid.y);
tauNum = zeros(NY-2,1);
for i=2:NY-1
   tauNum(i-1) = nf.etaS(i,2) * ( (nf.vx(i,2)-nf.vx(i-1,2))/(grid.yc(i+1)-grid.yc(i)) );
end
% figure, plot(tauNum)

