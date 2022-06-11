function goodtrials=manualtrialsel(manlist,ich,alltrials,type)
%05/2022 updated so that it looks at actual column #'s nto just where
%columns are populated, which could be inaccurate esepcially if no ch1
%numbers
%open bad trials if exists in homedir from excel sheet
%chronicXX_trialselection
%use targeted ich for selecting column
%must be saved in original data location
chbadtrials=[];
final=[];
if exist(manlist)>0
    [~,sheetsav] = xlsfinfo(manlist); %find # sheets in worksheet
    sheet=find(contains(sheetsav,type));
    if ~isempty(sheet)
       xlsdata=readcell(manlist);    %read xls as table (table format)
        channels={}
        for i= 1:size(xlsdata,2)
            if ischar(xlsdata{1,i})==0
                break
            end
            if contains(xlsdata{1,i}, 'Ch')
                channels=[channels, xlsdata(:,i)]
            elseif contains(xlsdata{1,i}, 'ch')
                channels=[channels, xlsdata(:,i)]
            end
        end
    
        for i= 1:size(channels,2)
            if contains(channels{1,i}, '1')
                channels{1,i}=1
            elseif contains(channels{1,i}, '2')
                channels{1,i}=2
            elseif contains(channels{1,i}, '3')
                channels{1,i}=3
            elseif contains(channels{1,i}, '4')
                channels{1,i}=4
            end
        end
        final=str2double(string(channels(2:end,:)))-99;   
    end
end

names=sheetnames(manlist); %added this
if sum(contains(names,type))==0; %added this
    final=0; %added this
end; %added this
for ch=1:4 
    if final==0; %added this
        chbadtrials(:,ch)=alltrials; %added this
        goodtrials{ch}=[]; %added this
    elseif sum(~isnan(final(:,ch)))==0 % adding comment to allow for commit
        chbadtrials(:,ch)=alltrials;
        goodtrials{ch}=[]
    else 
        goodtrials{ch}=alltrials(~ismember(alltrials,final(:,ch)));
    end
end
%change