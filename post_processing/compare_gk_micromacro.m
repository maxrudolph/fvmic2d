%% read file
filename = '../log1';

fh = fopen(filename,'r');
% T        sii             eii_mm          eii_gk          eta_mm          eta_gk
fgetl(fh);
fgetl(fh)
fgetl(fh);

data=fscanf(fh,'%e %e %e %e %e %e\n');

fclose(fh);

%% reshape to get useful quantities
data1=reshape(data,[6 length(data)/6]);
data1 = data1';

%% make a meshgrid of temperature and stress
tmin = min(data1(:,1));
tmax = max(data1(:,1));
smin = min(data1(:,2));
smax = max(data1(:,2));

T = linspace(tmin,tmax,20);
S = linspace(log10(smin),log10(smax),20);

[T1 S1] = meshgrid(T,S);

%% grid the data
etamm = griddata(data1(:,1),log10(data1(:,2)),log10(data1(:,5)),T1,S1);
etagk = griddata(data1(:,1),log10(data1(:,2)),log10(data1(:,6)),T1,S1);

%%
figure, hold on
[C h]=contour(T1,S1,etamm,'k');
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
[C h]=contour(T1,S1,etagk,'r--');
set(gca,'FontName','Times')
xlabel('Temperature (K)')
ylabel('log_{10} s')