function goodtrials=manualtrialsel(manlist,ich,alltrials,type)
%open bad trials if exists in homedir from excel sheet
%chronicXX_trialselection
%use targeted ich for selecting column
%must be saved in original data location
chbadtrials=[];
if exist(manlist)>0
    [~,sheetsav] = xlsfinfo(manlist); %find # sheets in worksheet
    sheet=find(contains(sheetsav,type));
    if ~isempty(sheet)
        xlsdata = xlsread(manlist,sheet);
       % chbadtrials=xlsdata-100;        % 101 first trial
        chbadtrials=xlsdata-99;        %default, 100 first trial
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