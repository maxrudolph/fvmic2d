% compare benchmark poiseuille flow solution with analytic solution

ifile = 'loadNodalFields_0_10.petscbin';

loadgrid;

nf = loadNodalFieldsPetscBin('loadNodalFields_0_10.petscbin');
rho = nf.rho(2,2);
mu = nf.etaS(2,2);
g=10;
L=grid.y(end);

yv = grid.yc(2:end);
analytic_soln = (yv*L-yv.^2)*rho*g/mu;

figure
plot( grid.yc(2:end), nf.vx(1:end-1,1), 'b')
hold on
plot( grid.yc(2:end), analytic_soln/2, 'rx')