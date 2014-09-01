%make a movie of convection simulation output
clear, close all
%% Set Output Options
avi_file_name = 'yielding_gk.avi';
avi_fps = 4; % frame rate
aviobj = avifile(avi_file_name,'fps',avi_fps);

%% Generate file list
fileNamePrefix='Markers.0';
files=dir([fileNamePrefix '.*']);
[filelist ltezero] = sort_files(files);
filelist = filelist(~ltezero);
LX=150000;
LY=75000;

nfiles = size(filelist,1);
secondsInYear=365*24*3600;

%% Open a figure window for plotting
fig=figure;

%% Begin iteration over marker files
for ifile=1:nfiles;
    markerfilename = filelist(ifile);
    g = getGriddedMarkers2(markerfilename{1},0);
    g.X = linspace(0,LX,g.nx);
    g.Y = linspace(0,LY,g.ny);
    
    %% calculate Nusselt number
    %grid is uniformly spaced, so nusselt number = total heat flow/conducted
    %heat flow = dtdz_top/dtdz total
    dtdz_top = (g.T(4,:)-g.T(3,:))/(g.Y(4)-g.Y(3));
    dtdz_top = mean(dtdz_top(~isnan(dtdz_top)));
    dtdz_total = (260-100)/LY;
    Nu = dtdz_top/dtdz_total;
    
    %% plot temperature image
    imagesc(g.X,g.Y,g.p'), axis equal tight;
%     caxis([100 260])
% caxis([0 5])
    fontname = 'Helvetica';
    fontsize=16;
    title(sprintf('Elapsed Time %3.f kyr, Nu = %.2f',g.elapsedTime/secondsInYear/1000,Nu),'FontName',fontname,'FontSize',fontsize);
    xlabel([num2str(LX/1000) ' km'],'FontName',fontname,'FontSize',16);
    ylabel([num2str(LY/1000) ' km'],'FontName',fontname,'FontSize',16);
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    colorbar
    
    drawnow;
    F = getframe(fig);
    aviobj = addframe(aviobj,F);
    
end

%% close avi file
aviobj = close(aviobj);