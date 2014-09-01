function [m] = getBinaryMarkers(filename,texture)
%open file for reading
fh = fopen(filename,'r');

%first get number of markers
clear m;
nMark = fread(fh,1,'int32');
m.elapsedTime = fread(fh,1,'double')
%read marker fields
m.cpu = fread(fh,nMark,'int32');
m.cellX = fread(fh,nMark,'int32');
m.cellY = fread(fh,nMark,'int32');

m.x = fread(fh,nMark,'double');
m.y = fread(fh,nMark,'double');
m.z = fread(fh,nMark,'double');
m.vx = fread(fh,nMark,'double');
m.vy = fread(fh,nMark,'double');
m.vz = fread(fh,nMark,'double');
m.Mat = fread(fh,nMark,'int32');
m.T = fread(fh,nMark,'double');
m.Tdot = fread(fh,nMark,'double');
m.eta = fread(fh,nMark,'double');
m.mu = fread(fh,nMark,'double');
if(texture)
    m.N11 = fread(fh,nMark,'double');
    m.N12 = fread(fh,nMark,'double');
    m.N21 = fread(fh,nMark,'double');
    m.N22 = fread(fh,nMark,'double');
end

m.D = fread(fh,nMark,'double');
m.Ddot = fread(fh,nMark,'double');

m.exx = fread(fh,nMark,'double');
m.exy = fread(fh,nMark,'double');
m.Eii = fread(fh,nMark,'double');

m.sxx = fread(fh,nMark,'double');
m.syy = fread(fh,nMark,'double');
m.szz = fread(fh,nMark,'double');
m.sxz = fread(fh,nMark,'double');
m.syz = fread(fh,nMark,'double');
m.sxy = fread(fh,nMark,'double');
m.p = fread(fh,nMark,'double');
m.rho = fread(fh,nMark,'double');
m.rhodot = fread(fh,nMark,'double');


if(texture)
m.res = fread(fh,nMark,'double');
m.srr = fread(fh,nMark,'double');
end
fclose(fh);
mask = m.x > 0 & m.y > 0;

names = fieldnames(m);
