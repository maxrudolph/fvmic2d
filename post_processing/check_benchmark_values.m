%calculate nusselt number

clear
close all
fclose all
[status pd] = unix('echo $PETSC_DIR');
% pd = '/da/'
PETSC_DIR='/opt/petsc';
setenv('PETSC_DIR',PETSC_DIR);

addpath([PETSC_DIR '/share/petsc/matlab']);

[bvx,bvy,bp,bT] = load_vankeken_benchmark_results();

% loadgrid

filelist=dir('../output/loadNodalFields_0_*.petscbin');
snums = zeros(size(filelist));
for i=1:length(snums)
    undpos = find(filelist(i).name == '_',1,'last');
    dotpos = find(filelist(i).name == '.',1,'first');
    snums(i) = str2num(filelist(i).name(undpos+1:dotpos-1));
end
[a,i] = sort(snums);
filelist = filelist(i);
snums = snums(i);

markfile = '../output/Markers.0.0';
% mark = getBinaryMarkers(markfile,0);
nskip=1;
iFile=1;
nfiles = length(filelist);
for iFile = 1:nfiles
    nf=loadNodalFieldsPetscBin2(['../output/' filelist(iFile).name]);
    
    % nf1=loadNodalFieldsPetscBin2(['../output/' filelist(end-1).name]);
    % nf2=loadNodalFieldsPetscBin2(['../output/' filelist(end).name]);
    grid.x = nf.gridx(1,:);
    grid.y = nf.gridy(:,1);
    NX = length(nf.gridx(1,:));
    NY = length(nf.gridy(:,1));
    
    % get cell-centered vx, vy
    slabv = 1.58e-9;
    
    % vxc = (nf2.vx(1:end-1,1:end-1) + nf2.vx(1:end-1,2:end))/2;
    % vyc = (nf2.vy(1:end-1,1:end-1) + nf2.vy(2:end,1:end-1))/2;
    
    % figure, imagesc(nf2.vx); title('vx'); colorbar; caxis([-slabv slabv]);
    % figure, imagesc(nf2.vy); title('vy'); colorbar; caxis([-slabv slabv]);
    % figure, imagesc(nf2.p); title('p'); colorbar;
    
    LX=max(grid.x);
    LY=max(grid.y);
    % xc = nf1.gridx(1,1:end-1) + diff(nf1.gridx(1,:))/2;
    % yc = nf1.gridy(1:end-1,1) + diff(nf1.gridy(:,1))/2;
    % [X,Y] = meshgrid(xc,yc);
    % figure, pcolor(xc,yc,sqrt(vxc.^2+vyc.^2)); shading flat;colorbar; caxis([-slabv slabv]);
    % set(gca,'YDir','reverse');
    nsl = 50;
    
    slx = rand(nsl,1)*LX;
    sly = rand(nsl,1)*LY;
    
    % hold on
    % streamline(xc,yc,vxc,vyc,slx,sly);
    
    % figure, pcolor(nf2.gridx/1e3,nf2.gridy/1e3,nf2.T); title('T'); hold on;
    % set(gca,'YDir','reverse');
    % colorbar;
    % h = streamline(xc/1e3,yc/1e3,vxc,vyc,xc(end-1)*ones(nsl,1)/1e3,linspace(51000,599000,nsl)'/1e3);
    % shading interp
    % set(h,'Color','k')
    % xlabel('Distance km)')
    % ylabel('Depth (km)');
    % axis equal tight
    
    % figure, pcolor(nf2.gridx/1e3,nf2.gridy/1e3,nf2.T-nf1.T); title('T'); hold on; caxis([-25 25]);
    % set(gca,'YDir','reverse');
    % colorbar;
    % h = streamline(xc/1e3,yc/1e3,vxc,vyc,xc(end-1)*ones(nsl,1)/1e3,linspace(51000,599000,nsl)'/1e3);
    % shading interp
    % set(h,'Color','k')
    % xlabel('Distance km)')
    % ylabel('Depth (km)');
    % axis equal tight
    
    % re-sample T onto a 6x6 km grid
    newx = 0:6000:660000;
    newy = 0:6000:600000;
    [X,Y] = meshgrid(newx,newy);
    newT = interp2(nf.gridx,nf.gridy,nf.T,X,Y,'linear')-273;
   
    newvx = interp2(nf.gridx,nf.gridy,nf.vx,X,Y,'linear');
    newvy = interp2(nf.gridx,nf.gridy,nf.vy,X,Y,'linear');
    newp = interp2(nf.gridx,nf.gridy,nf.p,X,Y,'linear');
    
    Tslab = 0;
    for i=1:36
        Tslab = Tslab + newT(i,i)^2;
    end
    Tslab = sqrt(Tslab / 36);
    
    Twedge = 0;
%     mask = zeros(size(newT));
    for i=10:21
        for j=10:i
            Twedge = Twedge + newT(j,i)^2;
%             mask(j,i) = 1;
        end
    end
    Twedge = sqrt(Twedge/78);
%     figure, imagesc(mask);
    
    Tslabs(iFile) = Tslab;
    Twedges(iFile) = Twedge;
    T1111(iFile) = newT(11,11);
    times(iFile) = nf.elapsedTime/3.15e7/1e6;
end
% disp(sprintf('T_11,11 = %e\n||Tslab || = %e\n||Twedge|| = %e\n',newT(11,11),Tslab,Twedge));

figure;
vscale = 3.156e9;
subplot(2,2,1);
imagesc(newvx*vscale);
title('vx - my calc');
subplot(2,2,2);
imagesc(-newvy*vscale);
subplot(2,2,3);
imagesc(newp*vscale);
subplot(2,2,4);
imagesc(newT*vscale);
title('Max temperature')

figure, subplot(1,2,1);
imagesc(bT); colorbar;
subplot(1,2,2);
imagesc(newT); colorbar;

mask = times>10;
[fit1,gof1] = extrapolate_values(times(mask),T1111(mask));
[fit2,gof2] = extrapolate_values(times(mask),Tslabs(mask));
[fit3,gof3] = extrapolate_values(times(mask),Twedges(mask));




figure, subplot(3,1,1);
plot(times,T1111);
hold on
plot([times(1) times(end)],[1 1]*fit1.c,'r');
ylabel('T_{11,11}');
subplot(3,1,2);
plot(times,Tslabs);
hold on
plot([times(1) times(end)],[1 1]*fit2.c,'r');
ylabel('||Tslab||');
subplot(3,1,3);
plot(times,Twedges);
hold on
plot([times(1) times(end)],[1 1]*fit3.c,'r');
ylabel('||Twedge||');
xlabel('Time (Myr)');



%% check divergence field
divv = zeros(NY-1,NX-1);
vx = nf.vx;
vy = nf.vy;
for i=1:NX-1
    for j=1:NY-1
        divv(j,i) = (vx(j,i+1)-vx(j,i))/(grid.x(i+1)-grid.x(i)) + (vy(j+1,i)-vy(j,i))/(grid.y(j+1)-grid.y(j));
    end
end
figure, imagesc(divv)



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

