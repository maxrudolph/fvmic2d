function result=getParamsFromCSV(fn)
%process csv files



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
vbz = mydata{idx1,2};

idx1=find(strcmp(mydata(:,1),'vbx'));
vbx = mydata{idx1,2};

idx1=find(strcmp(mydata(:,1),'damageB0'));
damageB0 = mydata{idx1,2};

idx1=find(strcmp(mydata(:,1),'damageAlpha3_0'));
damageAlpha3 = mydata{idx1,2};

idx1=find(strcmp(mydata(:,1),'damagem0'));
damagem = mydata{idx1,2};

result.vbz=vbz;
result.vbx=vbx;
result.damageB0=damageB0;
result.damageAlpha3=damageAlpha3;
result.damagem=damagem;