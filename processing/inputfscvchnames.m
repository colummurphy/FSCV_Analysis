function inputfscvchnames(pathdir,fscvchnames)
%provide directory with fscv_multi trial by trial fscv data
%provide fscv ch name labels and resave all files with this var
%eg) fscvchnames={'p3','pl2','pl3','cl5'}
d=dir(pathdir);
cd(pathdir);
filenames={d.name};
targfiles=strfind(filenames,'fscv_multi');
processfiles=find(~cellfun(@isempty,targfiles));
%contains function does not work for 2013 for cells 
files=filenames(processfiles);

for ifile=1:size(files,2)
    disp(['processing file # ' ...
        num2str(ifile) '/' num2str(size(files,2))]);
    load(files{ifile});
    filename=files{ifile};

    %resave all loaded data in same paths with fscvchnames var
    save([pathdir, filename],'fscv','samplesLFP','tsLFP','NlxEventTTL','NlxEventTS','nlxFileNames','fscvchnames')
    %disp(['saved to: ' pathSave filename]);
end



end