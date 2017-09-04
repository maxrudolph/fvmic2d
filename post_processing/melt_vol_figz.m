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


vscale = 3.156e9;
s_in_yr = 3.156e7;
ro = 3300;
g = 10.0;


% loadgrid
output_dir = '~/subduction_test50/output';
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
meltfvx(imag(meltfvx) ~= 0) = NaN;
meltfvy(imag(meltfvy) ~= 0) = NaN;

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

melt_rate = (vxc.*dfdx + vyc.*dfdy).*s_in_yr.*dA;
melt_rate(melt_rate<0) = 0;
melt_rate(isnan(melt_rate)) = 0;
melt_rate(:,1:2) = 0;

% sum melt production vertically
meltsum = sum(melt_production,1);
meltsum2 = sum(melt_rate,1);
%%
newx = nf.gridxc(1,2:end-1);
tipx = (newx - plate_thickness/tand(slab_angle))/1e3;
melt_vol = meltsum(2:end)./dX(2:end)*s_in_yr;
melt_vol2 = meltsum2(2:end)./dX(2:end);
melt_vol(:,1:10) = 0;
melt_vol2(:,1:10) = 0;
% figure;
plot(newx/1e3,melt_vol,'LineWidth',2)
% plot(newx/1e3,melt_vol,'LineWidth',2,'Color',[0    0.5000         0])
hold on;
% set(gca,'XLim',[0, (440000 - plate_thickness/tand(slab_angle))/1e3]);
set(gca,'XLim',[0,480]);
% set(gca,'YLim',[0,7e-5]);
xlabel('Distance from trench (km)','FontSize',12);
ylabel('Melt production (m^3/m^2/yr)','FontSize',12);
% hleg = legend('No thickened root','Thickened root (15 km)','Thickened root (narrow)','Thickened root (50 K hotter)','Thickened root (deep decouple)','Location','NE');
hleg = legend('Reference case','Decoupling no root','Thickened root','Thickened root w/ decoupling','Small root','Thickened root deep decoupling', 'Location','NE');
% hleg = legend('No thickened root','5 km','10 km','15 km','20 km','Location','NE');
% hleg = legend('290 km','240 km','190 km','Location','NE');
set(hleg,'FontSize',12)
set(gca,'FontSize',12)
set(gcf, 'Color', 'w');

% export_fig test.png -m2

% v = get(hleg,'title');
% set(v,'string','Position of crustal root');
% %% make plot showing melt volume vs distance
% ij = size(melt_production);
% melt_vol=zeros(ij(2),1);
% newx = nf.gridx(1,:);
% newy = nf.gridy(:,1);
% newx(:,end) = [];
% newy(end,:) = [];
% [X,Y] = meshgrid(newx,newy);
% 
% for j=1:ij(2)
%     for i=1:ij(1)
%         melt_vol(j) = melt_vol(j) + melt_production(i,j)*max(newx)/(ij(2)-1)*max(newy)/(ij(1)-1);
%     end
% end
% 
% figure(301);
% plot(newx/1e3,melt_vol)
% set(gca,'XLim',[0, 560]);
% xlabel('Distance from trench (km)');
% ylabel('Melt production (kg/m/yr)');

