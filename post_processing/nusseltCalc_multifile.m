%calculate nusselt number

clear
loadgrid
filelist=dir('loadNodalFields_0*.m');
nskip=5;
%%
NX=length(grid.x);
NY=length(grid.y);
LX=max(grid.x);
LY=max(grid.y);

%sort filelist
for i=1:length(filelist)
    stepnum(i) = str2double( filelist(i).name(find(filelist(i).name == '_',1,'last')+1:end-2));
end
[num idx] = sort(stepnum);
filelist = filelist(idx);

clear stepnum

files=1:nskip:length(filelist);
if(files(end) ~= length(filelist))
files=[files length(filelist)];
end

for iFile=1:length(files)
    %strip off .m at end of filename, load results
    disp(filelist(files(iFile)).name(1:end-2))
    run(filelist(files(iFile)).name(1:end-2));
    %get temperature field
    T=reshape(lastT,[NX NY])';
    %compute Nu
    for i=1:NX-1
        %compute cell-centered dT/dz
        dtdz(i) = 1/2*(T(2,i)-T(1,i) + T(2,i+1)-T(1,i+1));
        dtdz(i) = dtdz(i)/(grid.y(2)-grid.y(1));
        %weight by cell width
        %     dtdz(i) = dtdz(i)*
        dx(i)=(grid.x(i+1)-grid.x(i));
    end
    %gerya's definition, table 16.1
    Nu(iFile) = LY/(LX*(T(end,10)-T(1,10)))*sum(dtdz.*dx)
    
    %blankenbach 1989 definition
    num=0;
    denom=0;
    for i=1:NX-1
        num=num+ 1/2*(T(2,i)-T(1,i) + T(2,i+1)-T(1,i+1))/(grid.y(2)-grid.y(1))*(grid.x(i+1)-grid.x(i));
        denom=denom+ 1/2*(T(end,i)+T(end,i+1))*(grid.x(i+1)-grid.x(i));
    end
    
    Nu1(iFile) = LY*num/denom;
    
    
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
    vrms(iFile) = sqrt(integrand/LX/LY)*LY/(k/Cp/rho0)
    time(iFile) = elapsedTime;
    stepnum(iFile) = str2double( filelist(files(iFile)).name(find(filelist(files(iFile)).name == '_',1,'last')+1:end-2)); 
end
%%
siny = 60*60*24*7*365.25;
timey=time/siny;
nu_bench = 4.8844;

figure, plot(timey,Nu,'rx');
hold on
plot([0 max(timey)],[1 1]*nu_bench,'k');
title(sprintf('res %dx%d, Nu=%e\n',NX,NY,mean(Nu(end-3:end))));
print(gcf,'-depsc','Nu_time_plot.eps')



