%process csv files
clear;
dirnames=['data/vep-ridge1';
    'data/vep-ridge2';
    'data/vep-ridge3';
    'data/vep-ridge4'];

i1=1;
for idir =1:size(dirnames,1)
    dn = dirnames(idir,:);
    fns = dir([dn '/runInfo.*.csv']);

    for ifile = 1:size(fns,1)
        
        fn = [dn '/' fns(ifile).name];

        %begin reading the file
        fid = fopen(fn,'r');
        clear mydata
        iRow  =1;
        while(~feof(fid))
            l=fgetl(fid);
            if(l(1) == '#')
                %ignore this line
            else
                cidx=find(l == ',');
                label1 =l(1:cidx-1);
                num = sscanf(l(cidx+1:length(l)),'%e\n',1);
                mydata{iRow,1}=label1;
                mydata{iRow,2} =num;
                iRow=iRow+1;
            end
        end
        fclose(fid);
        %end reading the file

        %pull out any desired fields
        idx1=find(strcmp(mydata(:,1),'vbz'));
        vbzs(i1) = mydata{idx1,2};

        idx1=find(strcmp(mydata(:,1),'damageB0'));
        damageB0s(i1) = mydata{idx1,2};
        
       
        i1 = i1+1;
    end
end

figure;
loglog(vbzs,damageB0s,'LineStyle','none','Marker','x');