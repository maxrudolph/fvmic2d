%calculate nusselt number

clear
fclose all
[status pd] = unix('echo $PETSC_DIR');
% pd = '/da/'
PETSC_DIR='/usr/local/petsc';
setenv('PETSC_DIR',PETSC_DIR);

addpath([PETSC_DIR '/share/petsc/matlab']);

loadgrid
filelist=dir('loadNodalFields_0*.petscbin');
nskip=1;
%%
NX=length(grid.x);
NY=length(grid.y);
LX=max(grid.x);
LY=max(grid.y);

%sort filelist
for i=1:length(filelist)
    stepnum(i) = str2double( filelist(i).name(find(filelist(i).name == '_',1,'last')+1:end-9));
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
    disp(filelist(files(iFile)).name)
 
    nf=loadNodalFieldsPetscBin2(filelist(files(iFile)).name);
    if( isfield(nf,'gridx') )
        clear grid;
        grid.x = nf.gridx(1,:);
        grid.y = nf.gridy(:,1);
        NX = length(nf.gridx(1,:));
        NY = length(nf.gridy(:,1));
    else
        loadgrid
    end
    nxs(iFile) = NX;
    nys(iFile) = NY;
    T=nf.T;
    vx=nf.vx;
    vy=nf.vy;
    elapsedTime = nf.elapsedTime;
    
    %compute Nu
    for i=1:nf.NX-1
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

    
    
       
    %calculate temperature profile through center of model
    T1 = nf.T(:,floor((NX+1)/2));
    %nondimensionalize this
    Ttop = nf.T(1,floor((NX+1)/2));
    Tbtm= nf.T(end,floor((NX+1)/2));
    T1 = (T1-Ttop)/(Tbtm-Ttop);
    
        
    %calculate non-dimensional temperature gradients in upper left and
    %right corners
    q1(iFile) = LY/(Tbtm-Ttop)*(nf.T(2,1)-nf.T(1,1))/(grid.y(2)-grid.y(1))
    q2(iFile) = LY/(Tbtm-Ttop)*(nf.T(2,end)-nf.T(1,end))/(grid.y(2)-grid.y(1))
    
    
    %calculate normalized rms velocity
    vx=reshape(vx,[NX NY])';
    vy=reshape(vy,[NX NY])';
    integrand=0;
    for i=1:NX-1
        for j=1:NY-1
            %compute cell-center vx and vy
            vx1 = (vx(j,i)^2+vx(j,i+1)^2)/2;
            vy1 = (vy(j,i)^2+vy(j+1,i)^2)/2;
            
            integrand=integrand + (vx1+vy1)*(grid.x(i+1)-grid.x(i))*(grid.y(j+1)-grid.y(j));
        end
    end
    rho0=4000;
    Cp=1250;
    k=5;
    vrms(iFile) = LY*rho0*Cp/k *sqrt(1/(LX*LY)*integrand)
    time(iFile) = elapsedTime;
    stepnum(iFile) = str2double( filelist(files(iFile)).name(find(filelist(files(iFile)).name == '_',1,'last')+1:end-9)); 
%      figure, pcolor(grid.x,grid.y,nf.T), axis ij
end

%% set benchmark values
%case 1a
% nu_bench = 4.8844;
% vrms_bench = 42.865;
% q1_bench = 8.0593;
% q2_bench = 0.5888;


%case 1c
% nu_bench = 21.972;
% vrms_bench = 833.99;
% q1_bench = 45.964;
% q2_bench = 0.8772;

%case 2a
 nu_bench =  10.066;
 vrms_bench = 480.43;
 q1_bench = 17.531;
 q2_bench = 1.0085;



%%
siny = 60*60*24*7*365.25;
timey=time/siny;

figure, subplot(1,4,1);
hold all
plot(timey,Nu,'kx');
hold on
plot([0 max(timey)],[1 1]*nu_bench,'k');
subplot(1,4,2);
plot(timey,vrms,'k.'), hold on;
plot([0 max(timey)],[1 1]*vrms_bench,'k');
title(sprintf('res %dx%d, Nu=%e\n',NX,NY,mean(Nu(end-3:end))));
subplot(1,4,3);
plot(stepnum,timey,'.')
print(gcf,'-depsc','Nu_time_plot.eps')
subplot(1,4,4);
plot(timey,q1,'g.'), hold on, plot(timey,q2,'r.'), 
plot([0 max(timey)],[1 1]*q1_bench,'k');
plot([0 max(timey)],[1 1]*q2_bench,'k');
title(sprintf('q1=%fm q2=%f',mean(q1(end-3:end)),mean(q2(end-3:end))))

%% Plot dimensionless temperature profile and benchmark values
H=max(grid.y);
figure, plot(T1,(H-grid.y)/H);
hold on
%Case 1a:
    scatter(.4222,.2249);
    scatter(.5778,.7751);
%Case 1b:
% scatter(.4284,.1118);
% scatter(.5716,.8882);
%Case 1c:
%    scatter(0.4322,0.0577);
%    scatter(0.5678,0.9423);
%case 2a:
%       scatter(.7405,.0623);
%       scatter(.8323,.8243);

