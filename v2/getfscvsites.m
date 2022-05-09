function fscvsites=getfscvsites(pathname)
%05/2022
%Get FSCV site names (e.g. pl1, pl2, based on how the file was named during
%recording in CVtoTXT folder), assuming Tarheel format
filesorig=dir(pathname);
filenamesorig={filesorig.name};
targfiles=strfind(filenamesorig,'1dr_');
labeledfiles=find(~cellfun(@isempty,targfiles));
if isempty(labeledfiles) 
targfiles=strfind(filenamesorig,'1dr2_');
labeledfiles=find(~cellfun(@isempty,targfiles));
end
if isempty(labeledfiles)    
targfiles=strfind(filenamesorig,'1dr3_');
labeledfiles=find(~cellfun(@isempty,targfiles));
end
filesorig=filenamesorig(labeledfiles);
sitenamessp=strsplit(filesorig{1},'_');
fscvsites=sitenamessp(2:5);
