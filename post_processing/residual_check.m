% examine the residual and solution from an iterative solution attempt
close all
stokes = PetscBinaryRead('stokes_linearsystem.petscbin','cell',10000);

if( length(stokes) == 5 )
    lhs = stokes{1};
    prec = stokes{2};
    x = stokes{3};
    rhs = stokes{4};
    r = stokes{5};
elseif(length(stokes) == 4)
    lhs = stokes{1};    
    prec = lhs;
    x = stokes{2};
    rhs = stokes{3};
    r = stokes{4};
end
n = size(lhs,1);

pdof = 3:3:n;
vxdof = 1:3:n-2;
vydof = 2:3:n-1;

vdof = union(vxdof,vydof);

a00=lhs(vdof,vdof);
a01=lhs(vdof,pdof);
a10=lhs(pdof,vdof);
a11=lhs(pdof,pdof);
p00=prec(vdof,vdof);
p01=prec(vdof,pdof);
p10=prec(pdof,vdof);
p11=prec(pdof,pdof);

nx = 21;
ny = 21;

pres = reshape(r(pdof),[nx ny])';
vxres = reshape(r(vxdof),[nx ny])';
vyres = reshape(r(vydof),[nx ny])';
figure, imagesc(pres), title('p residual')
figure, imagesc(vxres), title('vx residual')
figure, imagesc(vyres), title('vy residual')

p =  reshape(x(pdof),[nx ny])';
vx =  reshape(x(vxdof),[nx ny])';
vy =  reshape(x(vydof),[nx ny])';
figure, imagesc(p), title('p field'), colorbar

%% try scaling lhs matrix
d1 = max(abs(lhs'));
d = sparse(diag(d1.^-1));
lhs_s = d*lhs;
condest(lhs)
condest(lhs_s)