%calculate nusselt number

clear

loadgrid
loadNodalFields_0_700

%%


NX=length(grid.x);
NY=length(grid.y);
LX=max(grid.x);
LY=max(grid.y);
T=reshape(lastT,[NX NY])';

for i=1:NX-1
    %compute cell-centered dT/dz
    dtdz(i) = 1/2*(T(2,i)-T(1,i) + T(2,i+1)-T(1,i+1));
    dtdz(i) = dtdz(i)/(grid.y(2)-grid.y(1));
    %weight by cell width
%     dtdz(i) = dtdz(i)*
    dx(i)=(grid.x(i+1)-grid.x(i));
end
%gerya's definition, table 16.1
Nu = LY/(LX*(T(end,10)-T(1,10)))*sum(dtdz.*dx)

%calculate normalized rms velocity
vx=reshape(vx,[NX NY])';
vy=reshape(vy,[NX NY])';
integrand=0;
for i=1:NX-1
    for j=1:NY-1
        %compute cell-center vx and vy
        vx1 = (vx(j,i)+vx(j,i+1))/2;
        vy1 = (vy(j,i)+vy(j+1,i))/2;
        
        integrand=integrand + (vx1^2+vy1^2)*(grid.x(i+1)-grid.x(i))*(grid.y(j+1)-grid.y(j));
    end
end
rho0=4000;
Cp=1250;
k=5;
vrms = sqrt(integrand/LX/LY)*LY/(k/Cp/rho0)