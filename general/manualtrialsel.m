function goodtrials=manualtrialsel(manlist,ich,alltrials,type)
%05/2022 updated so that it looks at actual column #'s nto just where
%columns are populated, which could be inaccurate esepcially if no ch1
%numbers
%open bad trials if exists in homedir from excel sheet
%chronicXX_trialselection
%use targeted ich for selecting column
%must be saved in original data location
chbadtrials=[];
if exist(manlist)>0
    [~,sheetsav] = xlsfinfo(manlist); %find # sheets in worksheet
    sheet=find(contains(sheetsav,type));
    if ~isempty(sheet)
       % xlsdata = xlsread(manlist,sheet);
       xlsdata=readtable(manlist,'Sheet',sheet);    %read xls as table (table format)
       %NEED TO CONVERT TO NUMBERS FOR EACH CHANNEL NOW
       % chbadtrials=xlsdata-100;        % 101 first trial
        chbadtrials=xlsdata-99;        %default, 100 first trial
        %CONDITIONAL IF A CHANNEL IS EMPTY MAKE chbad trials=alltrials
    end
end
goodtrials=alltrials;         %default, all trials goodtrials
goodtrials={};
for ich=1:4
    %add 5/2019 all channels, not just #2
if ~isempty(chbadtrials)        
    %use manual selected trials
    if size(chbadtrials,2)>1
        goodtrials{ich}=alltrials(find(~ismember(alltrials,chbadtrials(:,ich))));
    else
        %only single ch in excel sheet
        goodtrials{ich}=alltrials(find(~ismember(alltrials,chbadtrials)));
    end        
end
end
end