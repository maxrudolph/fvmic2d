% read linear system for debugging
clear;
close all;
filename = 'stokes_linearsystem.petscbin';
ls = PetscBinaryRead(filename,'cell',10000);
LHS = ls{1};

lhsf = full(LHS);
% look for rows with only zeros
for i=1:size(lhsf,1)
   if( max(abs(lhsf(i,:))) < 1e-10 )
       error(['small value at row ' num2str(i)]);
   end
end