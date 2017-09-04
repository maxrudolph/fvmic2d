%calculate nusselt number

clear
% close all
% fclose all;
[status pd] = unix('echo $PETSC_DIR');
% pd = '/da/'
 PETSC_DIR='/opt/petsc-3.5';
%PETSC_DIR='/usr/local/petsc';
setenv('PETSC_DIR',PETSC_DIR);

addpath([PETSC_DIR '/share/petsc/matlab']);

% [bvx,bvy,bp,bT] = load_vankeken_benchmark_results();

load('melt_table_0.1.mat');

vscale = 3.156e9;
s_in_yr = 3.156e7;


% loadgrid
output_dir = '~/subduction_test47/output';
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
    npres=300;
    
    LY = max(nf.gridy(:,1));
    slabx=linspace(plate_thickness/tand(slab_angle),LY/tand(slab_angle),npres);
    slaby=tand(slab_angle)*slabx;
    slabx = slabx+1860;
    pintern = nf.p(2:end,2:end);    % throw out ghost cells

    
    xc = nf.gridxc(2:end,2:end);
    yc = nf.gridyc(2:end,2:end);
    
%     figure, pcolor(xc(1,:),yc(:,1),pintern);
    p = interp2(xc,yc,pintern,slabx,slaby);
    hold on;
%     plot(slabx,slaby);
%     figure, plot(slabx,p);
    r = sqrt( (slabx-slabx(1)).^2 + (slaby-slaby(1)).^2 );
    ptot = trapz(r,p)
    
    Ttot = trapz(r,p.*r)
    %% Make figure showing temperature field and streamlines
%     figure;
%     set(gcf,'Position',[560   558  522*2   390*2]);
%     set(gca,'FontSize',11,'FontName','Helvetica');
%     %pcolor(nf.gridx/1e3,nf.gridy/1e3,nf.T); shading flat
%     contourf(nf.gridx/1e3,nf.gridy/1e3,nf.T, 100, 'color', 'none');
%     hold on;
%     hcb=colorbar;
%     colormap hot
%     axis equal tight
%     hcb.Label.String='Temperature (K)';
%     hcb.Label.FontSize=12;
%     set(gca,'YDir','reverse');
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
%     hold on;
%     hp = plot( xc(1,bwb(:,2))/1e3, yc(bwb(:,1),1)/1e3);
%     hp.LineWidth=3;
%     hp.Color=0.7*[1 1 1];
    % alpha(double(plate_mask));
    
%     hold on
    mask = xc(1,:) <= 450*1000;
%     hsl=streamline(xc(:,mask)/1e3,yc(:,mask)/1e3,vxc(:,mask),vyc(:,mask),slx,sly);
%     set(hsl,'Color','k')
%     
%     set(gca,'XLim',[plate_thickness/tand(slab_angle)/1e3 440]);
%     set(gca,'YLim',[0 160]);    
%     set(gca,'XLim',[0 440]);
%     set(gca,'YLim',[0 160]);
% %     xlabel('Distance from wedge corner (km)');
%     xlabel('Distance from trench (km)');   
%     ylabel('Depth (km)');
%     
    % re-sample T onto a 6x6 km grid
    newx = 0:3000:660000;
    newy = 0:3000:600000;
%     newx = linspace(0,560000,264);
%     newy = linspace(0,200000,251);
    [X,Y] = meshgrid(newx,newy);
    newT = interp2(nf.gridx,nf.gridy,nf.T,X,Y,'linear')-273;
    
    newvx = interp2(nf.gridx,nf.gridyc,nf.vx,X,Y,'linear');
    newvy = interp2(nf.gridxc,nf.gridy,nf.vy,X,Y,'linear');
    newp = interp2(nf.gridxc,nf.gridyc,nf.p,X,Y,'linear');
        % reduced dynamic pressure when eta = 3e19 instead of 1e21
    newp = 3e19/1e21*newp;
    %
    newtotp = 3300*10*Y+newp;
    % calculate melt fraction
    meltf = interp2( melt_table.T, melt_table.P*1e9, melt_table.F, newT, newtotp ,'linear',0.0);
    dfdp = interp2( melt_table.T, melt_table.P*1e9, melt_table.dFdP, newT, newtotp,'linear',0.0);
    dfdT = interp2( melt_table.T, melt_table.P*1e9, melt_table.dFdT, newT, newtotp,'linear',0.0);
    %% calculate gradient of pressure and temperature
    dx = newx(2)-newx(1);
    dy = newy(2)-newy(1);
    [dpdx,dpdy] = gradient(newtotp,dx,dy);
    [dTdx,dTdy] = gradient(newT,dx,dy);

    dPdt = (dpdx.*newvx + dpdy.*newvy)*s_in_yr;
    dTdt = (dTdx.*newvx + dTdy.*newvy)*s_in_yr;
    
    melt_production = dPdt.*dfdp/1e9 + dTdt.*dfdT;
    melt_production(melt_production<0) = 0;
    melt_production(isnan(melt_production)) = 0;
    
    % figure;
%     set(gca,'XTick',[plate_thickness/tand(slab_angle)/1e3, plate_thickness/tand(slab_angle)/1e3+50, plate_thickness/tand(slab_angle)/1e3+100, plate_thickness/tand(slab_angle)/1e3+150,plate_thickness/tand(slab_angle)/1e3+200,plate_thickness/tand(slab_angle)/1e3+250, plate_thickness/tand(slab_angle)/1e3+300]);
%     set(gca,'XTickLabel',[0,50,100,150,200,250,300]);
%     set(gca,'YTick',[0,30,60,90,120,150]);
%     %comment out above if needed for distance from trench!!!
%     set(gca,'YTick',[0,50,100,150,200]);
%     set(gca,'FontSize',11,'FontName','Helvetica');
%     a1=gca;
%     a2=axes();
%     a2.Position = a1.Position;
%     [cd,hc]=contour(X/1e3,Y/1e3,(melt_production/1e-11)); shading interp;
%     colormap(a2,'Jet');
%     hcb2=colorbar('South');
%     hcb2.Label.String = 'Melt Production (year^{-1}) x10^{-11}';
%     hcb2.Label.FontSize=12;
%     hcb2.FontSize=12;
%     a2.YTick=[];
%     a2.XTick=[];
%     set(gca,'Color','none');
%     a2.PlotBoxAspectRatio = a1.PlotBoxAspectRatio;
%     a2.PlotBoxAspectRatioMode = a1.PlotBoxAspectRatioMode;
%     a2.XLim = a1.XLim;
%     a2.YLim = a1.YLim;
%     a2.YDir = a1.YDir;
end

%% make plot showing melt volume vs distance 
ij = size(melt_production);
melt_vol=zeros(ij(2),1);

for j=1:ij(2)
    for i=1:ij(1)
        melt_vol(j) = melt_vol(j) + melt_production(i,j)*max(newx)/(ij(2)-1)*max(newy)/(ij(1)-1);
    end
end
tipx = (newx - plate_thickness/tand(slab_angle))/1e3;
% figure;
plot(newx/1e3,melt_vol,'LineWidth',2)
hold on;
% set(gca,'XLim',[0, (440000 - plate_thickness/tand(slab_angle))/1e3]);
set(gca,'XLim',[0,480]);
xlabel('Distance from wedge corner (km)');
ylabel('Melt production (kg/m/yr)');
% hleg = legend('Reference case','No thickened root','Thickened root','Double dip angle','Location','NE');
hleg = legend('0 km','5 km','10 km','15 km','20 km','Location','NE');
set(hleg,'FontSize',12)
% v = get(hleg,'title');
% set(v,'string','Position of crustal root');

%% upward velocity vs melt volume
min_vy=min(nf.vy);
dist=linspace(0,660,264);
figure,plot(dist,melt_vol);
figure,plot(dist,min_vy);


