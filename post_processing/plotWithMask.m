function dataToPlot =plotWithMask(X,Y,im1,mask)

%idea is to do something like imagesc
minval=min(im1(mask));
maxval=max(im1(mask));
if(minval==maxval)
    minval=minval-sqrt(eps(minval));
    maxval=maxval+sqrt(eps(maxval));
end
nx=size(im1,1);
ny=size(im1,2);
dataToPlot=zeros(nx,ny,3);

cm1 = jet;
cidx = linspace(minval,maxval,size(cm1,1));

%map image to cdata

dataToPlot(:,:,1) = interp1(cidx,cm1(:,1),im1);
dataToPlot(:,:,2) = interp1(cidx,cm1(:,2),im1);
dataToPlot(:,:,3) = interp1(cidx,cm1(:,3),im1);
% dataToPlot(mask,:) = [1 1 1];
for i=1:nx
    for j=1:ny
        if(~mask(i,j))
            dataToPlot(i,j,1) = 1;
            dataToPlot(i,j,2) = 1;
            dataToPlot(i,j,3) = 1;
        end
    end
end

image(X/1000,Y/1000,dataToPlot);
set(gca,'FontName','Times New Roman','FontSize',12)
hcb=colorbar;

set(hcb,'Ytick',[1 16 32 48 64],'YTickLabel',{linspace(minval,maxval,5)});
xlabel('X (km)');
ylabel('Y (km)');
