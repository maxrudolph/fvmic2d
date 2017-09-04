clear
close all
fclose all;
[status pd] = unix('echo $PETSC_DIR');
% pd = '/da/'
PETSC_DIR='/opt/petsc-3.5';
%PETSC_DIR='/usr/local/petsc';
setenv('PETSC_DIR',PETSC_DIR);

addpath([PETSC_DIR '/share/petsc/matlab']);

% [bvx,bvy,bp,bT] = load_vankeken_benchmark_results();


s_in_yr = 3.156e7;
ro = 3300;
g = 10.0;


% loadgrid
output_dir = '~/subduction_test66/output';
plate_thickness=50000;          %PLATE THICKNESS!!!!
%slab_angle=40;                  % SLAB ANGLE !!!!
slab_angle=20;                  % SLAB ANGLE !!!!
xH2O=0.1;
%output_dir ='~/fvmic2d/output';
% output_dir = '../case2';


filelist=dir([output_dir '/loadNodalFields_0_*.petscbin']);
snums = zeros(size(filelist));
for i=1:length(snums)
    undpos = find(filelist(i).name == '_',1,'last');
    dotpos = find(filelist(i).name == '.',1,'first');
    snums(i) = str2num(filelist(i).name(undpos+1:dotpos-1));
end
[a,i] = sort(snums);
filelist = filelist(i);
snums = snums(i);

% markfile = '../output/Markers.0.0';
% mark = getBinaryMarkers(markfile,0);
nskip=1;
iFile=1;
nfiles = length(filelist);
for iFile = 1:6
    nf=loadNodalFieldsPetscBin2([output_dir '/' filelist(iFile).name]);
    
    % nf1=loadNodalFieldsPetscBin2(['../output/' filelist(end-1).name]);
    % nf2=loadNodalFieldsPetscBin2(['../output/' filelist(end).name]);
    grid.x = nf.gridx(1,:);
    grid.y = nf.gridy(:,1);
    NX = length(nf.gridx(1,:));
    NY = length(nf.gridy(:,1));
    
    nf.gridxc = zeros(size(nf.gridx));
    nf.gridyc = zeros(size(nf.gridy));
    for i=1:NY
        for j=1:NX
            if( j == NX)
                nf.gridxc(i,j) = nf.gridx(i,j) + 0.5*(nf.gridx(i,j)-nf.gridx(i,j-1));
            else
                nf.gridxc(i,j) = (nf.gridx(i,j+1)+nf.gridx(i,j))/2;
            end
            if( i == NY )
                nf.gridyc(i,j) = nf.gridy(i,j) + 0.5*(nf.gridy(i,j) - nf.gridy(i-1,j));
            else
                nf.gridyc(i,j) = (nf.gridy(i+1,j)+nf.gridy(i,j))/2;
            end
        end
    end
    
    
    
    %%
%        reduced dynamic pressure when eta = 3e19 instead of 1e21
    nf.p = 3e19/1e21*nf.p; 
    %%
    npres=300;
    
    LY = max(nf.gridy(:,1));
    slabx=linspace(plate_thickness/tand(slab_angle),LY/tand(slab_angle),npres);
    slaby=tand(slab_angle)*slabx;
    slabx = slabx+1860; 
    pintern = nf.p(2:end,2:end);    % throw out ghost cells

    
    xc = nf.gridxc(1:end-1,1:end-1);
    yc = nf.gridyc(1:end-1,1:end-1);
    
    %     figure, pcolor(xc(1,:),yc(:,1),pintern);
    p = interp2(xc,yc,pintern,slabx,slaby);
%     hold on;
%     plot(slabx,slaby);
    %     figure, plot(slabx,p);
    r = sqrt( (slabx-slabx(1)).^2 + (slaby-slaby(1)).^2 );
    ptot = trapz(r,p)
    
    Ttot = trapz(r,p.*r)
end

%%
% Calculate P and T at Vx and Vy nodes
for i=2:NY-1
    for j=2:NX-1
        Tvx(i,j) = (nf.T(i,j) + nf.T(i+1,j))/2;
        Pvx(i,j) = (nf.p(i+1,j+1) + nf.p(i+1,j))/2 + nf.gridyc(i,j)*ro*g;
    end
end

for i=2:NY-1
    for j=2:NX-1
        Tvy(i,j) = (nf.T(i,j) + nf.T(i,j+1))/2;
        Pvy(i,j) = (nf.p(i+1,j+1) + nf.p(i,j+1))/2 + nf.gridy(i,j)*ro*g;
    end
end

% truncate pressure
Pvx(Pvx<0) = 0;
Pvy(Pvy<0) = 0;

% calculate F at vx and vy nodes
for i=2:NY-1
    for j=2:NX-1
        meltfvx(i,j) = KatzMelt(Pvx(i,j)/1e9, Tvx(i,j) - 273, xH2O);
    end
end

for i=2:NY-1
    for j=2:NX-1
        meltfvy(i,j) = KatzMelt(Pvy(i,j)/1e9, Tvy(i,j) - 273, xH2O);
    end
end
% calculate T at cell centers
for i=1:NY-2
    for j=1:NX-2
        Tc(i,j) = (Tvx(i,j) + Tvx(i,j+1) + Tvy(i,j) + Tvy(i+1,j))/4;
    end
end
meltfc = zeros(249,262);
% calculate P at cell centers
for i=1:NY-2
    for j=1:NX-2
        Pc(i,j) = nf.p(i+1,j+1) + nf.gridyc(i,j)*ro*g;
    end
end
Pc(Pc<0) = 0;
% calculate F at cell centers
for i=2:NY-2
    for j=2:NX-2
        meltfc(i,j) = KatzMelt(Pc(i,j)/1e9, Tc(i,j) - 273, xH2O);
    end
end

meltfvx(imag(meltfvx) ~= 0) = NaN;
meltfvy(imag(meltfvy) ~= 0) = NaN;
meltfc(imag(meltfc) ~= 0) = NaN;

% calculate df/dx and df/dy at cell centers
for i=1:NY-2
    for j=1:NX-2
        dfdx(i+1,j+1) = (meltfvx(i,j+1) - meltfvx(i,j))/(nf.gridx(i,j+1) - nf.gridx(i,j));
        dfdy(i+1,j+1) = (meltfvy(i+1,j) - meltfvy(i,j))/(nf.gridy(i+1,j) - nf.gridy(i,j));
    end
end

% calculate veolocity at cell centers
for i=1:NY-2;
    for j=1:NX-2;
        vxc(i+1,j+1) = (nf.vx(i,j) + nf.vx(i,j+1))/2;
        vyc(i+1,j+1) = (nf.vy(i,j) + nf.vy(i+1,j))/2;
    end
end
% calculate cell areas
for i=1:NY-2
    for j=1:NX-2
        dA(i+1,j+1) = (nf.gridx(i,j+1)-nf.gridx(i,j))*(nf.gridy(i+1,j)-nf.gridy(i,j));
        dX(j+1) = nf.gridx(i,j+1)-nf.gridx(i,j);
    end
end

% calculate melt production df/dt
melt_production = (vxc.*dfdx + vyc.*dfdy).*dA;
melt_production(melt_production<0) = 0;
melt_production(isnan(melt_production)) = 0;
melt_production(:,1:2) = 0;
melt_production = melt_production(2:end,2:end); % this accounts for the fact that the first row and column are ghost cells

melt_rate = (vxc.*dfdx + vyc.*dfdy).*s_in_yr;
melt_rate(melt_rate<0) = 0;
melt_rate(isnan(melt_rate)) = 0;
melt_rate(:,1:2) = 0;
melt_rate = melt_rate(2:end,2:end);

mxc = nf.gridxc(1:end-2,1:end-2);
myc = nf.gridyc(1:end-2,1:end-2);
figure, pcolor(mxc,myc, meltfc); shading flat; colorbar

%% Make figure showing temperature field and streamlines
figure;
set(gcf,'Position',[560   558  522*2   390*2]);
set(gca,'FontSize',11,'FontName','Helvetica');
%pcolor(nf.gridx/1e3,nf.gridy/1e3,nf.T); shading flat
contourf(nf.gridx/1e3,nf.gridy/1e3,nf.T, 100, 'color', 'none');
hold on;
hcb=colorbar;
colormap hot
axis equal tight
hcb.Label.String='Temperature (K)';
hcb.Label.FontSize=12;
set(gca,'YDir','reverse');
nsl = 50;
LX = max(nf.gridx(1,:))/1e3;
LY = max(nf.gridy(:,1))/1e3;
%     slx = 0.99*LX*ones(nsl,1);
%     slx = 440*ones(nsl,1);
slx = 440*ones(nsl,1);
sly = linspace(0,120,nsl)';
%     sly = linspace(0,120,nsl)';
% vxc = (nf.vx(1:end-1,1:end-1) + nf.vx(1:end-1,2:end))/2;
% vyc = (nf.vy(1:end-1,1:end-1) + nf.vy(2:end,1:end-1))/2;
plate_mask = sqrt(vxc.^2+vyc.^2)<1e-14;
% plot overriding plate polygon
bwb=bwboundaries(plate_mask);
bwb=bwb{1}; %assume bwb structure only contains one polygon
hold on;
hp = plot( xc(1,bwb(:,2))/1e3, yc(bwb(:,1),1)/1e3);
hp.LineWidth=3;
hp.Color=0.7*[1 1 1];
% alpha(double(plate_mask));

hold on
mask = xc(1,:) <= 450*1000;
hsl=streamline(xc(:,mask)/1e3,yc(:,mask)/1e3,vxc(:,mask),vyc(:,mask),slx,sly);
set(hsl,'Color','k')

%     set(gca,'XLim',[plate_thickness/tand(slab_angle)/1e3 440]);
%     set(gca,'YLim',[0 160]);
set(gca,'XLim',[0 440]);
set(gca,'YLim',[0 160]);
%     xlabel('Distance from wedge corner (km)');
xlabel('Distance from trench (km)');
ylabel('Depth (km)');

% figure;
%     set(gca,'XTick',[plate_thickness/tand(slab_angle)/1e3, plate_thickness/tand(slab_angle)/1e3+50, plate_thickness/tand(slab_angle)/1e3+100, plate_thickness/tand(slab_angle)/1e3+150,plate_thickness/tand(slab_angle)/1e3+200,plate_thickness/tand(slab_angle)/1e3+250, plate_thickness/tand(slab_angle)/1e3+300]);
%     set(gca,'XTickLabel',[0,50,100,150,200,250,300]);
set(gca,'YTick',[0,30,60,90,120,150]);
%comment out above if needed for distance from trench!!!
set(gca,'YTick',[0,50,100,150,200]);

a1=gca;
a2=axes();
a2.Position = a1.Position;
[cd,hc]=contour(mxc/1e3,myc/1e3,melt_rate/1e-11,9); shading interp;
colormap(a2,'Jet');
hcb2=colorbar('South');
hcb2.Label.String = 'Melt Production (year^{-1}) x10^{-11}';
hcb2.Label.FontSize=11;
hcb2.FontSize=11;
a2.YTick=[];
a2.XTick=[];
set(gca,'Color','none');
set(gcf, 'Color', 'w');
a2.PlotBoxAspectRatio = a1.PlotBoxAspectRatio;
a2.PlotBoxAspectRatioMode = a1.PlotBoxAspectRatioMode;
a2.XLim = a1.XLim;
a2.YLim = a1.YLim;
a2.YDir = a1.YDir;
set(gca,'FontSize',11,'FontName','Arial');
% caxis([7 50])
% set(hcb2,'YTick',[1,2,3,4,5])
%% sum melt production vertically
meltsum = sum(melt_production,1)*s_in_yr;
meltsum(:,1:2) = 0;

% figure, plot(nf.gridxc(1,2:end-1),meltsum(2:end)./dX(2:end))

melttot = sum(meltsum(1:230))
% export_fig 

% %%
% [m,n]=size(nf.p);
% % for i = 1:m
% %     for j = 1:n
% %         if nf.p(i,j) < -1e9;
% %             nf.p(i,j) = 0.72*1e9;
% %         end
% %         if nf.p(i,j) > 0;
% %             nf.p(i,j) = 0;
% %         end
% %     end
% % end
% 
% figure;
% set(gcf,'Position',[560   558  522*2   390*2]);
% set(gca,'FontSize',11,'FontName','Helvetica');
% %pcolor(nf.gridx/1e3,nf.gridy/1e3,nf.T); shading flat
% contourf(nf.gridx/1e3,nf.gridy/1e3,nf.p/1e9, 100, 'color', 'none');
% hold on;
% hcb=colorbar;
% colormap hot
% axis equal tight
% hcb.Label.String='Dynamic Pressure (GPa)';
% hcb.Label.FontSize=12;
% set(gca,'YDir','reverse');
% nsl = 50;
% LX = max(nf.gridx(1,:))/1e3;
% LY = max(nf.gridy(:,1))/1e3;
% %     slx = 0.99*LX*ones(nsl,1);
% %     slx = 440*ones(nsl,1);
% slx = 440*ones(nsl,1);
% sly = linspace(0,120,nsl)';
% %     sly = linspace(0,120,nsl)';
% % vxc = (nf.vx(1:end-1,1:end-1) + nf.vx(1:end-1,2:end))/2;
% % vyc = (nf.vy(1:end-1,1:end-1) + nf.vy(2:end,1:end-1))/2;
% plate_mask = sqrt(vxc.^2+vyc.^2)<1e-14;
% % plot overriding plate polygon
% bwb=bwboundaries(plate_mask);
% bwb=bwb{1}; %assume bwb structure only contains one polygon
% hold on;
% hp = plot( xc(1,bwb(:,2))/1e3, yc(bwb(:,1),1)/1e3);
% hp.LineWidth=3;
% hp.Color=0.7*[1 1 1];
% % alpha(double(plate_mask));
% 
% hold on
% mask = xc(1,:) <= 450*1000;
% hsl=streamline(xc(:,mask)/1e3,yc(:,mask)/1e3,vxc(:,mask),vyc(:,mask),slx,sly);
% set(hsl,'Color','k')
% 
% %     set(gca,'XLim',[plate_thickness/tand(slab_angle)/1e3 440]);
% %     set(gca,'YLim',[0 160]);
% set(gca,'XLim',[0 440]);
% set(gca,'YLim',[0 160]);
% %     xlabel('Distance from wedge corner (km)');
% xlabel('Distance from trench (km)');
% ylabel('Depth (km)');

% export_fig test.png -m2.5
