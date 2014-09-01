clear
iMonte=0;
% iTime=0;
dTime = 100;%interval between write-outs

stopt=0; 
stopm=0;

while(~stopm)
    iMonte = iMonte+1;
    iTime = 0;
    stopt=0;
    
    itmax=0;
    while(~stopt)
        fn = sprintf('/mnt/cluster/runs/vep-ridge3/output/GriddedMarkers.%d.%d',iMonte,iTime);
        if(isempty(dir(fn)) && iTime > 0)
            iTime = iTime-dTime;
            disp(sprintf('last output iMonte=%d, iTime=%d\n',iMonte,iTime))
            itmax = iTime;
            stopt=1;
            %get the markers
            %iTime = iTime - dTime;
            fn = sprintf('/mnt/cluster/runs/vep-ridge3/output/GriddedMarkers.%d.%d',iMonte,iTime);
            %markers=getGriddedMarkers(fn);
            %pause;
            
        elseif(isempty(dir(fn)) && iTime == 0)
            stopm=1;
            stopt=1;
            disp(sprintf('last montecarlo iteration present is %d\n',iMonte-1))
        end
        iTime=iTime+dTime;
    end
    %delete unwanted files
    stopt=0;
    if(itmax>0)
        for i=1:itmax-dTime
           if( mod(i,50) && ~mod(i,dTime) )
               fn = sprintf('/mnt/cluster/runs/vep-ridge3/output/GriddedMarkers.%d.%d',iMonte,i);              
               unix(['rm ' fn]);
               fn = sprintf('/mnt/cluster/runs/vep-ridge3/output/Markers.%d.%d',iMonte,i);
               unix(['rm ' fn]);
                fn = sprintf('/mnt/cluster/runs/vep-ridge3/output/loadNodalFields.%d.%d.m',iMonte,i);
               unix(['rm ' fn]);
               
           end
        end

    end
end
fn = sprintf('/mnt/cluster/runs/vep-ridge3/output/GriddedMarkers.%d.%d',iMonte,iTime);



% markers=getGriddedMarkers(fn);