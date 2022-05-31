function logbadtrials(subj,sess,trtype,trial,chs)
%write bad trials based on ch# etc. to excel sheet (TXT file)
%eg. logbadtrials(hgui.subject,hgui.sessionnum,hgui.trialtype,hgui.trialnum,1:4)
%trialnum is the label of trial num eg. trial num 100 = trial 1
%xlswrite(FILE,ARRAY,SHEET,RANGE) writes to the specified SHEET and RANGE.
pathsave=fullfile('Y:', 'data_MIT', 'analysis', subj, 'triallogs', filesep);
if ~isdir(pathsave)
    mkdir(pathsave);
end
trialtypeidx=strfind(trtype,'_');
trialtype=trtype(1:trialtypeidx-1);
sessid='';
if isnumeric(sess)
    sessid=num2str(sess);
else
    sessid=sess;
end
filename=['chronic' sessid '_' trialtype '_trialselection.txt'];

manlist=[pathsave filename];

A = repmat(trial,1,4);      %make row of 4 for all 4 channels
AB=zeros(1,4);
AB(chs)=1;
rowwrite=A.*AB;
fileid=0;
if exist(manlist)>0
    fileID = fopen(manlist,'a');        %append if txt file already exists
else
    fileID = fopen(manlist,'w');
    fprintf(fileID,'%s \t %s \t%s \t%s\t\r\n','ch 1','ch 2','ch 3','ch 4');
end
fprintf(fileID,'%d\t %d \t%d \t%d\t\r\n',rowwrite);
disp(['writing : ' num2str(rowwrite)]);
fclose(fileID);

end