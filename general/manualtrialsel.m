function goodtrials=manualtrialsel(manlist,ich,alltrials,type)
%JWC / USA made changes 06/01/2022 to look at headers for channnel #s
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


       %NEED TO CONVERT TO NUMBERS FOR EACH CHANNEL NOW
       % chbadtrials=xlsdata-100;        % 101 first trial
%         chbadtrials=[final(1,:); final(2:size(final,1),:)-99];        %default, 100 first trial
        %CONDITIONAL IF A CHANNEL IS EMPTY MAKE chbad trials=alltrials
    end
end

for ch=1:4 %added these loops but i dont think this will work unless we remove/change lines 28, and 29 since that sets goodtrials=alltrials
    if sum(~isnan(final(:,ch)))==0
        chbadtrials(:,ch)=alltrials;
        goodtrials{ch}=[];
    else 
        goodtrials{ch}=alltrials(~ismember(alltrials,final(:,ch)));
    end
end

% badtrials=zeros(size(final,1),size(final,2));
% for i = 1:size(final,1)
%     for j = 1:size(final,2)
%         if isnan(final(i,j))
%             badtrials(i,j)=alltrials(i,j);
%         end
%     end
% end

% goodtrials=alltrials;         %default, all trials goodtrials
% goodtrials={};
% for ich=1:4
%     %add 5/2019 all channels, not just #2
%     if ~isempty(chbadtrials)        
%         %use manual selected trials
%         if size(chbadtrials,2)>1
%             goodtrials{ich}=alltrials(find(~ismember(alltrials,chbadtrials(:,ich))));
%         else
%             %only single ch in excel sheet
%             goodtrials{ich}=alltrials(find(~ismember(alltrials,chbadtrials)));
%         end        
%     end
% end