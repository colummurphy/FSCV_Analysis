function sitenames=getSiteNames(nlxFileNames,csc_map)
%provide csc_map with csc ids in odd cells, and site names in even cells
%provide nlxFileNames
%output sitenames for the filenames provided
sitenames={};
mapcsc=csc_map(1:2:end);
mapnames=csc_map(2:2:end);
label='\d';
for ii=1:length(nlxFileNames)
    numids=regexp(nlxFileNames{ii},label);
    cscid=nlxFileNames{ii}(numids);
    aa=startsWith(mapcsc,cscid);
    bb=endsWith(mapcsc,cscid);      %has to start and end with csc num so not 2 channels for instance
    aa1=find(aa);
    bb1=find(bb);
    singleChID=intersect(aa1,bb1);
    %check also same # digits, so not 3 and 33 both
    if length(singleChID)>1
        singleChIDs=singleChID;
        singleChID=[];
        for xx=1:length(singleChIDs)
            if isequal(length(mapcsc{singleChIDs(xx)}),length(cscid))
                singleChID=singleChIDs(xx);
            end
        end
    end
    sitenames{ii}=mapnames{singleChID};
end

end