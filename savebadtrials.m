function savebadtrials(subj,sess,varargin)
%save bad trials from text document in triallogs folder to excelsheet
%trialnum is the label of trial num eg. trial num 100 = trial 1
%xlswrite(FILE,ARRAY,SHEET,RANGE) writes to the specified SHEET and RANGE.
pathsave=fullfile('Z:', 'inj-monkey2', 'analysis', subj, 'triallogs', filesep);
if ~isdir(pathsave)
    mkdir(pathsave);
end
argnum=1;
ephys=0;    %ephys log, l
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'ephys'
            ephys=1;
    end
    argnum=argnum+1;
end
chcols={'B','C','D','E'};
xlsexist=0;
dirpath=dir(pathsave);
xlsname=['chronic' num2str(sess) '_trialselection.xlsx'];
if ephys
    xlsname=['chronic' num2str(sess) 'l' '_trialselection.xlsx'];
end
manlist=[pathsave xlsname];
trtypes={'bigreward','smallreward','targetbreak','fixbreak'};
for ix=1:length(trtypes)
    if ix==1
        %check if xls exists first, since otherwise will check for newly
        %written file during this loop & want to completely write anew new
        %sheets for different types if generated from scratch
        if exist(manlist)==2
            xlsexist=1;
        end
    end
    trialtype=trtypes{ix};
    filename=['chronic' num2str(sess) '_' trialtype '_trialselection.txt'];
    if ephys
        filename=['chronic' num2str(sess) 'l' '_' trialtype '_trialselection.txt'];
    end
    readfile={};
    if any(contains({dirpath.name},filename,'IgnoreCase',true))
        %txt file exists
        readfile=importdata([pathsave filename],'\t');
    end
    %open pre-existing xls file
    xlsdata=[];
    %spreadsheet already exists
    %read bad trials from excel spreadsheet
    sheet = 1;      %sheet 1 default big
    if ~isempty(strfind(trialtype,'small'))
    sheet=2;    %small reward sheet is #2
    elseif ~isempty(strfind(trialtype,'target'))
    sheet=3;    %target error sheet
    elseif ~isempty(strfind(trialtype,'fix'))
    sheet=4;    %target error sheet
    end
    if ~isempty(readfile)
        for ich=1:4
        trialnums=sort(unique(readfile.data(:,ich)));
        if xlsexist 
            [~,sheetsav] = xlsfinfo(manlist); %find # sheets in worksheet
            if sheet<=length(sheetsav)
                xlsdata = xlsread(manlist,sheet);            
                [ndata,txtdata,xlsdata] = xlsread(manlist,sheet);  
                ch1col=find(contains(txtdata(1,:),'ch 1','ignorecase',true));
                lastrow=1;
                finishfind=0;
                while finishfind<=0
                    lastrow=lastrow+1;
                    numsexist=find(~isnan([xlsdata{lastrow,:}]));
                    if isempty(numsexist)
                        finishfind=1;
                    end
                end
                xlswrite(manlist,trialnums,sheet,[chcols{ich} num2str(lastrow) ':' chcols{ich} num2str(lastrow+size(trialnums,1)-1)]);
            end
        else
            %generate new spreadsheet
            if ich==1
                xlswrite(manlist,readfile.textdata,trtypes{ix},['B1:E1']);
            end
            xlswrite(manlist,trialnums,trtypes{ix},[chcols{ich} '2:' chcols{ich} num2str(2+size(trialnums,1)-1)]);
        end
        end
    end
end

end