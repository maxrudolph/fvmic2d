% plot some ridge results from the markers

% close all

filename = 'Markers.0.1000.gridded';
%  filename = 'test.gridded';

lx = 4000;
ly = 2000;
g = getGriddedMarkers2(filename,0);
secondsInYear = 365.25*24*60*60;
disp(['Elapsed time ' num2str(g.elapsedTime/secondsInYear) ' yr'])

NY = g.ny;
NX = g.nx;
X = linspace(0,lx,NX);
Y = linspace(0,ly,NY);
%%
figure, imagesc(X,Y,g.vx'), axis image, title('vx'), colorbar
figure, imagesc(X,Y,g.vy'), axis image, title('vy'), colorbar
figure, imagesc(X,Y,g.vz'), axis image, title('vz'), colorbar
figure, imagesc(X,Y,g.D'), axis image, title('D'), colorbar
figure, imagesc(X,Y,g.Ddot'), axis image, title('Ddot'), colorbar

figure, imagesc(X,Y,g.T'), axis image, title('T'), colorbar
figure, imagesc(X,Y,g.p'), axis image, title('p'), colorbar
figure, imagesc(X,Y,log10(g.eta')), axis image, title('eta'), colorbar
figure, imagesc(X,Y,log10(g.mu')), axis image, title('mu'), colorbar
figure, imagesc(X,Y,g.rho'), axis image, title('rho'), colorbar
figure, imagesc(X,Y,g.rhodot'), axis image, title('rhodot'), colorbar
%%
figure, imagesc(X,Y,g.sii'), axis image, title('sii'), colorbar
figure, imagesc(X,Y,g.sxx'), axis image, title('sxx'), colorbar
figure, imagesc(X,Y,g.sxy'), axis image, title('sxy'), colorbar
%%
% figure, imagesc(X,Y,g.eii), axis image, title('sii'), colorbar
figure, imagesc(X,Y,g.exx'), axis image, title('exx'), colorbar
figure, imagesc(X,Y,g.exy'), axis image, title('exy'), colorbar

%% plot topography
grho = diff(g.rho,1,2);
tprof = zeros(NX,1);
for i=1:NX
    i1=find( grho(i,:) > sqrt(1000),1,'first' );
    tprof(i) = -Y(i1);
end
figure, plot(X,tprof);