function [texture] = readMarkerTextureASCII(fn)

fh = fopen(fn,'r');

NT = fscanf(fh,'NT=%e\n',1);
texture.eii=fscanf(fh,'eii=%e\n',1);
texture.Eii=fscanf(fh,'Eii=%e\n',1);
texture.eta=fscanf(fh,'eta=%e\n',1);

data=fscanf(fh,'%e\t%e\n',[2 NT^3]);
data = data';
texture.ctheta = data(:,1);
texture.cphi = data(:,2);
fclose(fh);