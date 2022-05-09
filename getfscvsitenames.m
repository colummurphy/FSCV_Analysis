function sites=getfscvsitenames(pathname)
sites=repmat('',1,4);
origpath=pathname(1:strfind(pathname,'cvtotxt')-1);
filesorig=dir(origpath);
filenamesorig={filesorig.name};
targfiles=strfind(filenamesorig,'1dr_');
labeledfiles=find(~cellfun(@isempty,targfiles));
%contains function does not work for 2013 for cells 
filesorig=filenamesorig(labeledfiles);
if isempty(filesorig)
    %try '1dr2_'
    targfiles=strfind(filenamesorig,'1dr2_');
    labeledfiles=find(~cellfun(@isempty,targfiles));
    if isempty(labeledfiles)
        %still cannot find, try 2dr
        targfiles=strfind(filenamesorig,'2dr_');
        labeledfiles=find(~cellfun(@isempty,targfiles));
    end
    %contains function does not work for 2013 for cells 
    filesorig=filenamesorig(labeledfiles);
end
if ~isempty(filesorig)
sitenamessp=strsplit(filesorig{1},'_');
if length(sitenamessp)>5
    %should be 6 for 4 channels
sites=sitenamessp(2:5);
else
    sites=sitenamessp(2:end-1);
    while length(sites)<4
        sites{end+1}='xx';
    end
end
end
        
end
