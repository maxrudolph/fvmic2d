function [vx vy p T] = load_vankeken_benchmark_results()

nx = 111;
ny = 101;

bmfilepath = './';

fh = fopen([bmfilepath 'vx.dat'],'r');
data = fscanf(fh,'%le');
fclose(fh);

vx = reshape(data,[nx ny])';

fh = fopen([bmfilepath 'vy.dat'],'r');
data = fscanf(fh,'%le');
fclose(fh);

vy = reshape(data,[nx ny])';

fh = fopen([bmfilepath 'p.dat'],'r');
data = fscanf(fh,'%le');
fclose(fh);

p = reshape(data,[nx ny])';

fh = fopen([bmfilepath 'T.dat'],'r');
data = fscanf(fh,'%le');
fclose(fh);

T = reshape(data,[nx ny])';

figure;
subplot(2,2,1);
imagesc(vx);
title('vankeken - vx');
subplot(2,2,2);
imagesc(vy);
title('vy');
subplot(2,2,3);
imagesc(p);
title('P');
subplot(2,2,4);
imagesc(T);
title('T');