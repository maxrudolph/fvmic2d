function [texture] = readConsolidatedMarkerTextureBinary(fnb)
idx = 0;


fn = fnb;
fh = fopen(fn,'rb');

nCPU = fread(fh,1,'int')
%     for iCPU=1:nCPU
nMark = fread(fh,1,'int')
NT = fread(fh,1,'int')

for i=idx+1:idx+nMark
    texture(i).ctheta = fread(fh,NT^3,'double');
    texture(i).cphi = fread(fh,NT^3,'double');
end
idx = idx+nMark;
%     end
fclose(fh);

