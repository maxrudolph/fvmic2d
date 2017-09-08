
%calculate nusselt number

%clear
close all
fclose all;
% load('norootnodeTPV.mat');
[status pd] = unix('echo $PETSC_DIR');
% pd = '/da/'
 PETSC_DIR='~/sw/petsc-3.5.4';
%PETSC_DIR='/usr/local/petsc';
setenv('PETSC_DIR',PETSC_DIR);

addpath([PETSC_DIR '/share/petsc/matlab']);

% [bvx,bvy,bp,bT] = load_vankeken_benchmark_results();

load('melt_table_0.1.mat');

vscale = 3.156e9;
s_in_yr = 3.156e7;


% loadgrid
output_dir = '../case66';
plate_thickness=50000;          %PLATE THICKNESS!!!!
%slab_angle=40;                  % SLAB ANGLE !!!!
slab_angle=20;                  % SLAB ANGLE !!!!

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
for iFile = nfiles:nfiles
    nf=loadNodalFieldsPetscBin2([output_dir '/' filelist(iFile).name]);
    [output_dir '/' filelist(iFile).name]
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
       % reduced dynamic pressure when eta = 3e19 instead of 1e21
    nf.p = 3e19/1e21*nf.p; 
    %%
    npres=300;
    
    LY = max(nf.gridy(:,1));
    slabx=linspace(plate_thickness/tand(slab_angle),LY/tand(slab_angle),npres);
    slaby=tand(slab_angle)*slabx;
    slabx = slabx+1860;
    pintern = nf.p(2:end,2:end);    % throw out ghost cells
    xc = nf.gridxc(2:end,2:end);
    yc = nf.gridyc(2:end,2:end);
    
    figure, pcolor(xc(1,:),yc(:,1),pintern);
    p = interp2(xc,yc,pintern,slabx,slaby);
    hold on;
    plot(slabx,slaby);
    figure, plot(slabx,p);
    r = sqrt( (slabx-slabx(1)).^2 + (slaby-slaby(1)).^2 );
    ptot = trapz(r,p)
    
    Ttot = trapz(r,p.*r)
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
    vxc = (nf.vx(1:end-1,1:end-1) + nf.vx(1:end-1,2:end))/2;
    vyc = (nf.vy(1:end-1,1:end-1) + nf.vy(2:end,1:end-1))/2;
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
    
    % re-sample T onto a 6x6 km grid
    newx = linspace(min(nf.gridx(:)),max(nf.gridx(:)),440);
    newy = linspace(min(nf.gridy(:)),max(nf.gridy(:)),300);
    [X,Y] = meshgrid(newx,newy);
    newT = interp2(nf.gridx,nf.gridy,nf.T,X,Y,'linear')-273.0;
    
    newvx = interp2(nf.gridx,nf.gridyc,nf.vx,X,Y,'linear');
    newvy = interp2(nf.gridxc,nf.gridy,nf.vy,X,Y,'linear');
    newp = interp2(nf.gridxc,nf.gridyc,nf.p,X,Y,'linear');
    %
    newtotp = 3300*10*Y+newp;
    % calculate melt fraction
    meltf= interp2( melt_table.T, melt_table.P*1e9, melt_table.F, newT, newtotp ,'linear',0.0);
    % dfdp is in units of 1/GPa. Convert this into 1/Pa by multiplying by
    % (1 GPa/1e9 Pa)
    dfdp = interp2( melt_table.T, melt_table.P*1e9, melt_table.dFdP/1e9, newT, newtotp,'linear',0.0);
    dfdT = interp2( melt_table.T, melt_table.P*1e9, melt_table.dFdT, newT, newtotp,'linear',0.0);
    %% calculate gradient of pressure and temperature
    dx = newx(2)-newx(1);
    dy = newy(2)-newy(1);
    [dpdx,dpdy] = gradient(newtotp,dx,dy);
    [dTdx,dTdy] = gradient(newT,dx,dy);

    dPdt = (dpdx.*newvx + dpdy.*newvy)*s_in_yr;
    dTdt = (dTdx.*newvx + dTdy.*newvy)*s_in_yr;
    
    melt_production = dPdt.*dfdp + dTdt.*dfdT;
    melt_production(melt_production<0) = 0;
    melt_production(isnan(melt_production)) = 0;
    
    % figure;
%     set(gca,'XTick',[plate_thickness/tand(slab_angle)/1e3, plate_thickness/tand(slab_angle)/1e3+50, plate_thickness/tand(slab_angle)/1e3+100, plate_thickness/tand(slab_angle)/1e3+150,plate_thickness/tand(slab_angle)/1e3+200,plate_thickness/tand(slab_angle)/1e3+250, plate_thickness/tand(slab_angle)/1e3+300]);
%     set(gca,'XTickLabel',[0,50,100,150,200,250,300]);
    set(gca,'YTick',[0,30,60,90,120,150]);
    %comment out above if needed for distance from trench!!!
    set(gca,'YTick',[0,50,100,150,200]);
    set(gca,'FontSize',11,'FontName','Helvetica');
    a1=gca;
    a2=axes();
    a2.Position = a1.Position;
    [cd,hc]=contour(X/1e3,Y/1e3,(meltf),10); shading interp;
    colormap(a2,'Jet');
    hcb2=colorbar('South');
    hcb2.Label.String = 'Melt Production (year^{-1}) x10^{-11}';
    hcb2.Label.FontSize=12;
    hcb2.FontSize=12;
    a2.YTick=[];
    a2.XTick=[];
    set(gca,'Color','none');
    a2.PlotBoxAspectRatio = a1.PlotBoxAspectRatio;
    a2.PlotBoxAspectRatioMode = a1.PlotBoxAspectRatioMode;
    a2.XLim = a1.XLim;
    a2.YLim = a1.YLim;
    a2.YDir = a1.YDir;
    set(gca,'Color','none');
    set(gcf, 'Color', 'w');
%     caxis([3 40])
end

% %% make plot showing melt volume vs distance 
% ij = size(melt_production);
% melt_vol=zeros(ij(2),1);
% 
% for j=1:ij(2)
%     for i=1:ij(1)
%         melt_vol(j) = melt_vol(j) + melt_production(i,j)*max(newx)/(ij(2)-1)*max(newy)/(ij(1)-1);
%     end
% end
% 
% figure(301);
% plot(newx/1e3,melt_vol)
% set(gca,'XLim',[0, 560]);    set(hsl,'Color','k')
% 
% xlabel('Distance from trench (km)');
% ylabel('Melt production (kg/m/yr)');

%% make figure showing dynamic pressure
% nf.p(nf.p < -3e7 ) = -3e7;
% nf.p(isnan(p)) = 0;
% figure, pcolor(nf.gridx,nf.gridy, nf.p); axis ij; shading flat; colorbar

%% make figure showing upward velocity vy or Temperature difference
% figure;
%     set(gcf,'Position',[560   558  522*2   390*2]);
%     set(gca,'FontSize',11,'FontName','Helvetica');
%     %pcolor(nf.gridx/1e3,nf.gridy/1e3,nf.T); shading flat
% %     contourf(nf.gridx/1e3,nf.gridy/1e3, nf.vy, 100, 'color', 'none');
%     tempdiff = newT-rootnodenewT;
%     contourf(X/1e3,Y/1e3,tempdiff,200,'color','none');
% 
%     hold on;
%     hcb=colorbar;
%     colormap(flipud(brewermap([],'RdBu')));
%     axis equal tight
%     hcb.Label.String='\Delta T (K)';
%     hcb.Label.FontSize=12;
%     hcb.Direction='normal';
%     set(gca,'YDir','reverse');
%     nsl = 50;
%     LX = max(nf.gridx(1,:))/1e3;
%     LY = max(nf.gridy(:,1))/1e3;
% %     slx = 0.99*LX*ones(nsl,1);
% %     slx = 440*ones(nsl,1);
%     slx = 440*ones(nsl,1);
%     sly = linspace(0,120,nsl)';
% %     sly = linspace(0,120,nsl)';
%     vxc = (nf.vx(1:end-1,1:end-1) + nf.vx(1:end-1,2:end))/2;
%     vyc = (nf.vy(1:end-1,1:end-1) + nf.vy(2:end,1:end-1))/2;
%     plate_mask = sqrt(vxc.^2+vyc.^2)<1e-14;
%     % plot overriding plate polygon
%     bwb=bwboundaries(plate_mask);
%     bwb=bwb{1}; %assume bwb structure only contains one polygon
%     hold on;
%     hp = plot( xc(1,bwb(:,2))/1e3, yc(bwb(:,1),1)/1e3);
%     hp.LineWidth=3;
%     hp.Color=0.7*[1 1 1];
%     % alpha(double(plate_mask));
%     
%     hold on
%     mask = xc(1,:) <= 450*1000;
%     hsl=streamline(xc(:,mask)/1e3,yc(:,mask)/1e3,vxc(:,mask),vyc(:,mask),slx,sly);
%     set(hsl,'Color','k')
% %     quiver(X(1:8:end,1:8:end)/1e3,Y(1:8:end,1:8:end)/1e3,newvx(1:8:end,1:8:end),newvy(1:8:end,1:8:end),'color','k');
% 
%     
% %     set(gca,'XLim',[plate_thickness/tand(slab_angle)/1e3 440]);
% %     set(gca,'YLim',[0 160]);    
%     set(gca,'XLim',[0 440]);
%     set(gca,'YLim',[0 160]);
%     caxis ([-130 130]);
% %     caxis ([-2e-10 3.3789e-10]);
% %     caxis ([-max(max(nf.vx)) max(max(nf.vx))]);
% %     xlabel('Distance from wedge corner (km)');
%     xlabel('Distance from trench (km)');   
%     ylabel('Depth (km)');
%     
%     % re-sample T onto a 6x6 km grid
%     newx = 0:3000:660000;
%     newy = 0:3000:600000;
%     [X,Y] = meshgrid(newx,newy);
%     newT = interp2(nf.gridx,nf.gridy,nf.T,X,Y,'linear')-273;
%     
%     newvx = interp2(nf.gridx,nf.gridyc,nf.vx,X,Y,'linear');
%     newvy = interp2(nf.gridxc,nf.gridy,nf.vy,X,Y,'linear');
%     newp = interp2(nf.gridxc,nf.gridyc,nf.p,X,Y,'linear');
% 
%     %
%     newtotp = 3300*10*Y+newp;
%     % calculate melt fraction
%     meltf = interp2( melt_table.T, melt_table.P*1e9, melt_table.F, newT, newtotp ,'linear',0.0);
%     dfdp = interp2( melt_table.T, melt_table.P*1e9, melt_table.dFdP, newT, newtotp,'linear',0.0);
%     dfdT = interp2( melt_table.T, melt_table.P*1e9, melt_table.dFdT, newT, newtotp,'linear',0.0);
%     %% calculate gradient of pressure and temperature
%     dx = newx(2)-newx(1);
%     dy = newy(2)-newy(1);
%     [dpdx,dpdy] = gradient(newtotp,dx,dy);
%     [dTdx,dTdy] = gradient(newT,dx,dy);
% 
%     dPdt = (dpdx.*newvx + dpdy.*newvy)*s_in_yr;
%     dTdt = (dTdx.*newvx + dTdy.*newvy)*s_in_yr;
%     
%     melt_production = dPdt.*dfdp/1e9 + dTdt.*dfdT;
%     melt_production(melt_production<0) = 0;
%     melt_production(isnan(melt_production)) = 0;
%     
%     % figure;
% %     set(gca,'XTick',[plate_thickness/tand(slab_angle)/1e3, plate_thickness/tand(slab_angle)/1e3+50, plate_thickness/tand(slab_angle)/1e3+100, plate_thickness/tand(slab_angle)/1e3+150,plate_thickness/tand(slab_angle)/1e3+200,plate_thickness/tand(slab_angle)/1e3+250, plate_thickness/tand(slab_angle)/1e3+300]);
% %     set(gca,'XTickLabel',[0,50,100,150,200,250,300]);
%     set(gca,'YTick',[0,30,60,90,120,150]);
%     %comment out above if needed for distance from trench!!!
%     set(gca,'YTick',[0,50,100,150,200]);
%     set(gca,'FontSize',11,'FontName','Helvetica');
%     a1=gca;
% %     a2=axes();
% %     a2.Position = a1.Position;
% %     pcolor(X,Y,newvy); shading interp;
% %     colormap(a2,'Jet');
% %     hcb2=colorbar('South');
% %     hcb2.Label.String = 'Melt Production (year^{-1}) x10^{-11}';
% %     hcb2.Label.FontSize=12;
% %     hcb2.FontSize=12;
% %     a2.YTick=[];
% %     a2.XTick=[];
% %     set(gca,'Color','none');
% %     a2.PlotBoxAspectRatio = a1.PlotBoxAspectRatio;
% %     a2.PlotBoxAspectRatioMode = a1.PlotBoxAspectRatioMode;
% %     a2.XLim = a1.XLim;
% %     a2.YLim = a1.YLim;
%     a2.YDir = a1.YDir;
%     set(gca,'Color','none');
%     set(gcf, 'Color', 'w');
%     
%     
% %% make figure showing dynamic pressure
% nf.p(isnan(nf.p)) = 0;
% % newp(newp < newp(18,55)) = newp(18,55);
% figure;
%     set(gcf,'Position',[560   558  522*2   390*2]);
%     set(gca,'FontSize',11,'FontName','Helvetica');
%     %pcolor(nf.gridx/1e3,nf.gridy/1e3,nf.T); shading flat
%     pcolor(nf.gridx/1e3, nf.gridy/1e3, (nf.p)/1e6);shading flat
%     hold on;
%     hcb=colorbar;
%     colormap jet
%     axis equal tight
%     hcb.Label.String='Dynamic pressure (MPa)';
%     hcb.Label.FontSize=12;
% %     hcb.Direction='reverse';
%     set(gca,'YDir','reverse');
%     nsl = 50;
%     LX = max(nf.gridx(1,:))/1e3;
%     LY = max(nf.gridy(:,1))/1e3;
% %     slx = 0.99*LX*ones(nsl,1);
% %     slx = 440*ones(nsl,1);
%     slx = 440*ones(nsl,1);
%     sly = linspace(0,120,nsl)';
% %     sly = linspace(0,120,nsl)';
%     vxc = (nf.vx(1:end-1,1:end-1) + nf.vx(1:end-1,2:end))/2;
%     vyc = (nf.vy(1:end-1,1:end-1) + nf.vy(2:end,1:end-1))/2;
%     plate_mask = sqrt(vxc.^2+vyc.^2)<1e-14;
%     % plot overriding plate polygon
%     bwb=bwboundaries(plate_mask);
%     bwb=bwb{1}; %assume bwb structure only contains one polygon
%     hold on;
%     hp = plot( xc(1,bwb(:,2))/1e3, yc(bwb(:,1),1)/1e3);
%     hp.LineWidth=3;
%     hp.Color=0.7*[1 1 1];
%     % alpha(double(plate_mask));
%     
%     hold on
%     mask = xc(1,:) <= 450*1000;
%     hsl=streamline(xc(:,mask)/1e3,yc(:,mask)/1e3,vxc(:,mask),vyc(:,mask),slx,sly);
%     set(hsl,'Color','k')
% %     quiver(X(1:10:end,1:10:end)/1e3,Y(1:10:end,1:10:end)/1e3,newvx(1:10:end,1:10:end),newvy(1:10:end,1:10:end),'color','k');
% 
%     
% %     set(gca,'XLim',[plate_thickness/tand(slab_angle)/1e3 440]);
% %     set(gca,'YLim',[0 160]);    
%     set(gca,'XLim',[0 440]);
%     set(gca,'YLim',[0 160]);
%     caxis([-22 0])
% %     xlabel('Distance from wedge corner (km)');
%     xlabel('Distance from trench (km)');   
%     ylabel('Depth (km)');
%     
% %     % re-sample T onto a 6x6 km grid
% %     newx = 0:3000:660000;
% %     newy = 0:3000:600000;
% %     [X,Y] = meshgrid(newx,newy);
% %     newT = interp2(nf.gridx,nf.gridy,nf.T,X,Y,'linear')-273;
% %     
% %     newvx = interp2(nf.gridx,nf.gridyc,nf.vx,X,Y,'linear');
% %     newvy = interp2(nf.gridxc,nf.gridy,nf.vy,X,Y,'linear');
% %     newp = interp2(nf.gridxc,nf.gridyc,nf.p,X,Y,'linear');
% % 
% %     %
% %     newtotp = 3300*10*Y+newp;
% %     % calculate melt fraction
% %     meltf = interp2( melt_table.T, melt_table.P*1e9, melt_table.F, newT, newtotp ,'linear',0.0);
% %     dfdp = interp2( melt_table.T, melt_table.P*1e9, melt_table.dFdP, newT, newtotp,'linear',0.0);
% %     dfdT = interp2( melt_table.T, melt_table.P*1e9, melt_table.dFdT, newT, newtotp,'linear',0.0);
% %     %% calculate gradient of pressure and temperature
% %     dx = newx(2)-newx(1);
% %     dy = newy(2)-newy(1);
% %     [dpdx,dpdy] = gradient(newtotp,dx,dy);
% %     [dTdx,dTdy] = gradient(newT,dx,dy);
% % 
% %     dPdt = (dpdx.*newvx + dpdy.*newvy)*s_in_yr;
% %     dTdt = (dTdx.*newvx + dTdy.*newvy)*s_in_yr;
% %     
% %     melt_production = dPdt.*dfdp/1e9 + dTdt.*dfdT;
% %     melt_production(melt_production<0) = 0;
% %     melt_production(isnan(melt_production)) = 0;
% %     
%     % figure;
% %     set(gca,'XTick',[plate_thickness/tand(slab_angle)/1e3, plate_thickness/tand(slab_angle)/1e3+50, plate_thickness/tand(slab_angle)/1e3+100, plate_thickness/tand(slab_angle)/1e3+150,plate_thickness/tand(slab_angle)/1e3+200,plate_thickness/tand(slab_angle)/1e3+250, plate_thickness/tand(slab_angle)/1e3+300]);
% %     set(gca,'XTickLabel',[0,50,100,150,200,250,300]);
%     set(gca,'YTick',[0,30,60,90,120,150]);
%     %comment out above if needed for distance from trench!!!
%     set(gca,'YTick',[0,50,100,150,200]);
%     set(gca,'FontSize',11,'FontName','Helvetica');
%     a1=gca;
% %     a2=axes();
% %     a2.Position = a1.Position;
% %     pcolor(X,Y,newvy); shading interp;
% %     colormap(a2,'Jet');
% %     hcb2=colorbar('South');
% %     hcb2.Label.String = 'Melt Production (year^{-1}) x10^{-11}';
% %     hcb2.Label.FontSize=12;
% %     hcb2.FontSize=12;
% %     a2.YTick=[];
% %     a2.XTick=[];
% %     set(gca,'Color','none');
% %     a2.PlotBoxAspectRatio = a1.PlotBoxAspectRatio;
% %     a2.PlotBoxAspectRatioMode = a1.PlotBoxAspectRatioMode;
% %     a2.XLim = a1.XLim;
% %     a2.YLim = a1.YLim;
% %     a2.YDir = a1.YDir;
%     set(gca,'Color','none');
%     set(gcf, 'Color', 'w');
% 
% % tipx = (newx - plate_thickness/tand(slab_angle))/1e3;
% % figure(301);
% % plot(tipx,melt_vol)
% % set(gca,'XLim',[0, (560000 - plate_thickness/tand(slab_angle))/1e3]);
% % xlabel('Distance from wedge corner (km)');
% % ylabel('Melt production (kg/m/yr)');
%%
%
% %     [grid.xc grid.yc] = meshgrid(grid.xc,grid.yc);
%     % get cell-centered vx, vy
%     slabv = 1.58e-9 ;
%
%     % vxc = (nf2.vx(1:end-1,1:end-1) + nf2.vx(1:end-1,2:end))/2;
%     % vyc = (nf2.vy(1:end-1,1:end-1) + nf2.vy(2:end,1:end-1))/2;
%
%     % figure, imagesc(nf2.vx); title('vx'); colorbar; caxis([-slabv slabv]);
%     % figure, imagesc(nf2.vy); title('vy'); colorbar; caxis([-slabv slabv]);
%     % figure, imagesc(nf2.p); title('p'); colorbar;
%
%     LX=max(grid.x);
%     LY=max(grid.y);
%     % xc = nf1.gridx(1,1:end-1) + diff(nf1.gridx(1,:))/2;
%     % yc = nf1.gridy(1:end-1,1) + diff(nf1.gridy(:,1))/2;
%     % [X,Y] = meshgrid(xc,yc);
%     % figure, pcolor(xc,yc,sqrt(vxc.^2+vyc.^2)); shading flat;colorbar; caxis([-slabv slabv]);
%     % set(gca,'YDir','reverse');
%     nsl = 50;
%
%     slx = rand(nsl,1)*LX;
%     sly = rand(nsl,1)*LY;
%
%     % hold on
%     % streamline(xc,yc,vxc,vyc,slx,sly);
%
%     % figure, pcolor(nf2.gridx/1e3,nf2.gridy/1e3,nf2.T); title('T'); hold on;
%     % set(gca,'YDir','reverse');
%     % colorbar;
%     % h = streamline(xc/1e3,yc/1e3,vxc,vyc,xc(end-1)*ones(nsl,1)/1e3,linspace(51000,599000,nsl)'/1e3);
%     % shading interp
%     % set(h,'Color','k')
%     % xlabel('Distance km)')
%     % ylabel('Depth (km)');
%     % axis equal tight
%
%     % figure, pcolor(nf2.gridx/1e3,nf2.gridy/1e3,nf2.T-nf1.T); title('T'); hold on; caxis([-25 25]);
%     % set(gca,'YDir','reverse');
%     % colorbar;
%     % h = streamline(xc/1e3,yc/1e3,vxc,vyc,xc(end-1)*ones(nsl,1)/1e3,linspace(51000,599000,nsl)'/1e3);
%     % shading interp
%     % set(h,'Color','k')
%     % xlabel('Distance km)')
%     % ylabel('Depth (km)');
%     % axis equal tight
%
%     % re-sample T onto a 6x6 km grid
%     newx = 0:3000:660000;
%     newy = 0:3000:600000;
%     [X,Y] = meshgrid(newx,newy);
%     newT = interp2(nf.gridx,nf.gridy,nf.T,X,Y,'linear')-273;
%
%     newvx = interp2(nf.gridx,nf.gridyc,nf.vx,X,Y,'linear');
%     newvy = interp2(nf.gridxc,nf.gridy,nf.vy,X,Y,'linear');
%     newp = interp2(nf.gridxc,nf.gridyc,nf.p,X,Y,'linear');
%
%     newtotp = 3300*10*Y+newp;
%     % calculate melt fraction
%     meltf = interp2( melt_table.T, melt_table.P*1e9, melt_table.F, newT, newtotp ,'linear',0.0);
%     dfdp = interp2( melt_table.T, melt_table.P*1e9, melt_table.dFdP, newT, newtotp,'linear',0.0);
%     dfdT = interp2( melt_table.T, melt_table.P*1e9, melt_table.dFdT, newT, newtotp,'linear',0.0);
%
%
%     Tslab = 0;
%     for i=1:36
%         Tslab = Tslab + newT(i,i)^2;
%     end
%     Tslab = sqrt(Tslab / 36);
%
%     Twedge = 0;
% %     mask = zeros(size(newT));
%     for i=10:21
%         for j=10:i
%             Twedge = Twedge + newT(j,i)^2;
% %             mask(j,i) = 1;
%         end
%     end
%     Twedge = sqrt(Twedge/78);
% %     figure, imagesc(mask);
%
%     Tslabs(iFile) = Tslab;
%     Twedges(iFile) = Twedge;
%     T1111(iFile) = newT(11,11);
%     times(iFile) = nf.elapsedTime/3.15e7/1e6;
% end
% % disp(sprintf('T_11,11 = %e\n||Tslab || = %e\n||Twedge|| = %e\n',newT(11,11),Tslab,Twedge));
% %% calculate relative magnitudes of gradients in total pressure and dynamic pressure
% dx = newx(2)-newx(1);
% dy = newy(2)-newy(1);
% [dpdx,dpdy] = gradient(newtotp,dx,dy);
% [dpdyndx,dpdyndy] = gradient(bp*1e6,dx,dy);
%
% pratio = sqrt( dpdyndx.^2 + dpdyndy.^2 )./(4000*10);
%
% %% calculate gradient of pressure and temperature
% dx = newx(2)-newx(1);
% dy = newy(2)-newy(1);
% [dpdx,dpdy] = gradient(newtotp,dx,dy);
% [dTdx,dTdy] = gradient(newT,dx,dy);
% dPdt = (dpdx.*newvx + dpdy.*newvy)*s_in_yr;
% dTdt = (dTdx.*newvx + dTdy.*newvy)*s_in_yr;
%
% melt_production = dPdt.*dfdp/1e9 + dTdt.*dfdT;
% melt_production(melt_production<0) = 0;
% melt_production(isnan(melt_production)) = 0;
% %%
% figure;
% % subplot(1,2,1);
% fs=16;
% fn = 'Helvetica'
% pcolor(X/1e3,Y/1e3,newT);hold on;
% % add streamlines
% nsl = 30;
% sx = ones(nsl,1)*(LX-1000)/1e3;
% sy = linspace(60,600,nsl+2);
% sy = sy(2:end-1);
% h=streamline(X/1e3,Y/1e3,newvx,newvy,sx,sy);
% set(gca,'FontName',fn);
% set(h,'Color','k');
% set(gca,'YDir','reverse');
% shading interp;
% % title('Temperature and Streamlines');
% axis equal tight
% colorbar
% % xlabel('km');
% % ylabel('km');
% colormap(redblue());
% set(gca,'FontSize',fs,'FontName',fn);
% set(gcf,'Renderer','Painters')
% saveas(gcf,'temperature.eps','psc2');
%
% figure;
% % subplot(1,2,2);
% % pcolor(X/1e3,Y/1e3,log10(meltf));
% pcolor(X/1e3,Y/1e3,(melt_production/1e-11));
% caxis([0 10]);
% set(gca,'YDir','reverse');
% shading interp
% axis equal tight
% % set(gca,'XLim',[0 400]);
% % set(gca,'YLim',[0 400]);
% hold on;
% plot([0 600],[0 600],'w');
% plot([50 660],[50 50],'w');
% colorbar
% % xlabel('km');
% % ylabel('km');
% set(gca,'FontSize',fs,'FontName',fn);
% set(gca,'FontName',fn);
% colormap(hot);
% set(gcf,'Renderer','Painters')
% saveas(gcf,'melt_production.eps','psc2');
%
%
% %%
% figure;
% pcolor(X/1e3,Y/1e3,bT); shading interp;
% hold on;
% h=streamline(X/1e3,Y/1e3,bvx,-bvy,sx,sy);
% set(gca,'YDir','reverse');
%
% %%
%
% figure;
% vscale = 3.156e9;
% subplot(2,2,1);
% imagesc(newvx*vscale);
% subplot(2,2,2);
% imagesc(-newvy*vscale);
% subplot(2,2,3);
% imagesc(newp*vscale);
% subplot(2,2,4);
% imagesc(newT*vscale);
% title('Max temperature')
%
% figure, subplot(1,2,1);
% imagesc(bT); colorbar;
% subplot(1,2,2);
% imagesc(newT); colorbar;
%
% mask = times>10;
% [fit1,gof1] = extrapolate_values(times(mask),T1111(mask));
% [fit2,gof2] = extrapolate_values(times(mask),Tslabs(mask));
% [fit3,gof3] = extrapolate_values(times(mask),Twedges(mask));
%
%
%
%
% figure, subplot(3,1,1);
% plot(times,T1111);
% hold on
% plot([times(1) times(end)],[1 1]*fit1.c,'r');
% ylabel('T_{11,11}');
% subplot(3,1,2);
% plot(times,Tslabs);
% hold on
% plot([times(1) times(end)],[1 1]*fit2.c,'r');
% ylabel('||Tslab||');
% subplot(3,1,3);
% plot(times,Twedges);
% hold on
% plot([times(1) times(end)],[1 1]*fit3.c,'r');
% ylabel('||Twedge||');
% xlabel('Time (Myr)');
%
%
%
% %% check divergence field
% divv = zeros(NY-1,NX-1);
% vx = nf.vx;
% vy = nf.vy;
% for i=1:NX-1
%     for j=1:NY-1
%         divv(j,i) = (vx(j,i+1)-vx(j,i))/(grid.x(i+1)-grid.x(i)) + (vy(j+1,i)-vy(j,i))/(grid.y(j+1)-grid.y(j));
%     end
% end
% figure, imagesc(divv)
%
%
%
% %%
% NX=length(grid.x);
% NY=length(grid.y);
% LX=max(grid.x);
% LY=max(grid.y);
%
% %sort filelist
% for i=1:length(filelist)
%     stepnum(i) = str2double( filelist(i).name(find(filelist(i).name == '_',1,'last')+1:end-9));
% end
% [num idx] = sort(stepnum);
% filelist = filelist(idx);
%
% clear stepnum
%
% files=1:nskip:length(filelist);
%
% if(files(end) ~= length(filelist))
%     files=[files length(filelist)];
% end
%
% for iFile=1:length(files)
%     %strip off .m at end of filename, load results
%     disp(filelist(files(iFile)).name)
%
%     nf=loadNodalFieldsPetscBin2(filelist(files(iFile)).name);
%     if( isfield(nf,'gridx') )
%         clear grid;
%         grid.x = nf.gridx(1,:);
%         grid.y = nf.gridy(:,1);
%         NX = length(nf.gridx(1,:));
%         NY = length(nf.gridy(:,1));
%     else
%         loadgrid
%     end
%     nxs(iFile) = NX;
%     nys(iFile) = NY;
%     T=nf.T;
%     vx=nf.vx;
%     vy=nf.vy;
%     elapsedTime = nf.elapsedTime;
%
%     %compute Nu
%     for i=1:nf.NX-1
%         %compute cell-centered dT/dz
%         dtdz(i) = 1/2*(T(2,i)-T(1,i) + T(2,i+1)-T(1,i+1));
%         dtdz(i) = dtdz(i)/(grid.y(2)-grid.y(1));
%         %weight by cell width
%         %     dtdz(i) = dtdz(i)*
%         dx(i)=(grid.x(i+1)-grid.x(i));
%     end
%     %gerya's definition, table 16.1
%     Nu(iFile) = LY/(LX*(T(end,10)-T(1,10)))*sum(dtdz.*dx)
%
%     %blankenbach 1989 definition
%     num=0;
%     denom=0;
%     for i=1:NX-1
%         num=num+ 1/2*(T(2,i)-T(1,i) + T(2,i+1)-T(1,i+1))/(grid.y(2)-grid.y(1))*(grid.x(i+1)-grid.x(i));
%         denom=denom+ 1/2*(T(end,i)+T(end,i+1))*(grid.x(i+1)-grid.x(i));
%     end
%
%     Nu1(iFile) = LY*num/denom;
%
%
%
%
%     %calculate temperature profile through center of model
%     T1 = nf.T(:,floor((NX+1)/2));
%     %nondimensionalize this
%     Ttop = nf.T(1,floor((NX+1)/2));
%     Tbtm= nf.T(end,floor((NX+1)/2));
%     T1 = (T1-Ttop)/(Tbtm-Ttop);
%
%
%     %calculate non-dimensional temperature gradients in upper left and
%     %right corners
%     q1(iFile) = LY/(Tbtm-Ttop)*(nf.T(2,1)-nf.T(1,1))/(grid.y(2)-grid.y(1))
%     q2(iFile) = LY/(Tbtm-Ttop)*(nf.T(2,end)-nf.T(1,end))/(grid.y(2)-grid.y(1))
%
%
%     %calculate normalized rms velocity
%     vx=reshape(vx,[NX NY])';
%     vy=reshape(vy,[NX NY])';
%     integrand=0;
%     for i=1:NX-1
%         for j=1:NY-1
%             %compute cell-center vx and vy
%             vx1 = (vx(j,i)^2+vx(j,i+1)^2)/2;
%             vy1 = (vy(j,i)^2+vy(j+1,i)^2)/2;
%
%             integrand=integrand + (vx1+vy1)*(grid.x(i+1)-grid.x(i))*(grid.y(j+1)-grid.y(j));
%         end
%     end
%     rho0=4000;
%     Cp=1250;
%     k=5;
%     vrms(iFile) = LY*rho0*Cp/k *sqrt(1/(LX*LY)*integrand)
%     time(iFile) = elapsedTime;
%     stepnum(iFile) = str2double( filelist(files(iFile)).name(find(filelist(files(iFile)).name == '_',1,'last')+1:end-9));
%     %      figure, pcolor(grid.x,grid.y,nf.T), axis ij
% end
%
% %% set benchmark values
% %case 1a
% % nu_bench = 4.8844;
% % vrms_bench = 42.865;
% % q1_bench = 8.0593;
% % q2_bench = 0.5888;
%
%
% %case 1c
% % nu_bench = 21.972;
% % vrms_bench = 833.99;
% % q1_bench = 45.964;
% % q2_bench = 0.8772;
%
% %case 2a
% nu_bench =  10.066;
% vrms_bench = 480.43;
% q1_bench = 17.531;
% q2_bench = 1.0085;
%
%
%
% %%
% siny = 60*60*24*7*365.25;
% timey=time/siny;
%
% figure, subplot(1,4,1);
% hold all
% plot(timey,Nu,'kx');
% hold on
% plot([0 max(timey)],[1 1]*nu_bench,'k');
% subplot(1,4,2);
% plot(timey,vrms,'k.'), hold on;
% plot([0 max(timey)],[1 1]*vrms_bench,'k');
% title(sprintf('res %dx%d, Nu=%e\n',NX,NY,mean(Nu(end-3:end))));
% subplot(1,4,3);
% plot(stepnum,timey,'.')
% print(gcf,'-depsc','Nu_time_plot.eps')
% subplot(1,4,4);
% plot(timey,q1,'g.'), hold on, plot(timey,q2,'r.'),
% plot([0 max(timey)],[1 1]*q1_bench,'k');
% plot([0 max(timey)],[1 1]*q2_bench,'k');
% title(sprintf('q1=%fm q2=%f',mean(q1(end-3:end)),mean(q2(end-3:end))))
%
% %% Plot dimensionless temperature profile and benchmark values
% H=max(grid.y);
% figure, plot(T1,(H-grid.y)/H);
% hold on
% %Case 1a:
% scatter(.4222,.2249);
% scatter(.5778,.7751);
% %Case 1b:
% % scatter(.4284,.1118);
% % scatter(.5716,.8882);
% %Case 1c:
% %    scatter(0.4322,0.0577);
% %    scatter(0.5678,0.9423);
% %case 2a:
% %       scatter(.7405,.0623);
% %       scatter(.8323,.8243);
%
