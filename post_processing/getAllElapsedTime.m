% print elapsed time for all files in directory

flist = dir('*.petscbin');

for i=1:length(flist)
   nf=loadNodalFieldsPetscBin2(flist(i).name);
   disp( [flist(i).name ' ' sprintf('%e',nf.elapsedTime)])
end