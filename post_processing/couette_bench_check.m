% compare benchmark poiseuille flow solution with analytic solution

filename = 'loadNodalFields_0_5.petscbin';

loadgrid;

nf = loadNodalFieldsPetscBin(filename);
rho = nf.rho(2,2);
mu = nf.etaS(2,2);
g=10;
LY=grid.y(end);

yv = grid.yc(2:end);

close all

figure, pcolor(grid.x,grid.y,nf.T), colorbar, title('Temperature')
figure, plot(grid.y,nf.T(:,1)), title('Temperature profile')



u1 = 1e-9; %velocity of top plate wrt bottom
QonR = 5000;
T1=1300;
T0=1000;

% analytic_soln = (yv*L-yv.^2)*rho*g/mu;
analytic_soln=u1* (exp(-QonR*(T1-T0)/T0^2*(1-yv/LY))-1)./(exp(-QonR*(T1-T0)/T0^2)-1);
figure
plot( grid.yc(2:end), nf.vx(1:end-1,1), 'b')
hold on
plot( yv,analytic_soln,'rx')