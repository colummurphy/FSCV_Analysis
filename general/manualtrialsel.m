function goodtrials=manualtrialsel(manlist,ich,alltrials,type)
%line by line comments need to be added, what is final?
%JWC / USA made changes 06/01/2022 to look at headers for channnel #s
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

       xlsdata=readcell(manlist);    %read xls in cell format
        channels={};
        for i= 1:size(xlsdata,2)
            if ischar(xlsdata{1,i})==0 % Once the for loop gets to a column header with no words the loop exits, since notes is always one column away it catches the four channels
                break
            end
            if contains(xlsdata{1,i}, 'Ch') || contains(xlsdata{1,i}, 'ch') %only captures columns with channel in the header
                channels=[channels, xlsdata(:,i)];
            end
        end
    
        for i= 1:size(channels,2) % converts the trial header from words to numeric
            if contains(channels{1,i}, '1')
                channels{1,i}=1;
            elseif contains(channels{1,i}, '2')
                channels{1,i}=2;
            elseif contains(channels{1,i}, '3')
                channels{1,i}=3;
            elseif contains(channels{1,i}, '4')
                channels{1,i}=4;
            end
        end
        final=str2double(string(channels(2:end,:)))-99;   %final is the channel output in matrix output to itterate in the following loops. 
    end
end

names=sheetnames(manlist); %saves a variable with the names of all the sheets in the xls doc
if sum(contains(names,type))==0 %if the type isnt in the sheet names set the final matrix output to 0
    final=0; 
end 
for ch=1:4 
    if final==0 %makes all trials bad for sheets not in the xls doc
        chbadtrials(:,ch)=alltrials; 
        goodtrials{ch}=[]; 
    elseif sum(~isnan(final(:,ch)))==0  %if all the columns are empty then make all the channels bad channels
        chbadtrials(:,ch)=alltrials;
        goodtrials{ch}=[];
    else  %if the matrix has filled columns make them good trials for a given channel. 
        goodtrials{ch}=alltrials(~ismember(alltrials,final(:,ch)));
    end
end