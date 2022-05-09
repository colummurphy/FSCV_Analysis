function sites=getlfpsites(sessnums,chids, varargin)
%Some modifications 2022
%Unlike getsites, chids here is {'cl4-cl5','cl1-cl5',...}
%for each sess num, get corresponding site id for specified chids
%open site map from excel
sitelist='C:\Users\putamen\Dropbox (MIT)\monkey_fscv\patra\patra_map_04132019.xlsx';
argnum=1;
cleoflag=0;
sheetname='patra_map_spikes2021.xlsx';
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'cleo'
            sitelist='Z:\inj-monkey2\analysis\cleo\cleo_map_05102019.xlsx';
            cleoflag=1;
        case 'path'
            %supplied path name
            argnum=argnum+1;
            sitelist=[varargin{argnum} sheetname];
    end
    argnum=argnum+1;
end
%read bad trials from excel spreadsheet
sheet = 1;      %sheet 1 default site map
[~,sheetsav] = xlsfinfo(sitelist); %find # sheets in worksheet
sheetnum=find(strcmp(sheetsav,'sitemap')==1);
[ndata,txtdata,xlsdata] = xlsread(sitelist,sheetnum);
sesscol=find(strcmp(xlsdata(2,:),'day')==1);
chcols=3:6;
sitecols=7:10;
chcol1=find(strcmp(xlsdata(2,:),'ch1')==1);
if cleoflag
    chcols=chcol1:chcol1+3;
    sitecols=chcols(end)+1:chcols(end)+4;
    sesscol=sesscol-1;
end
%validrows=find(~isnan(ndata(:,2)));
%validrows=find(isnumeric(regexp(txtdata(:,sitecols),'.\d')));
[validrows,~]=find(~cellfun(@isempty,regexp(txtdata(:,sitecols),'.\d')));       %get rows with defined site names
sites={};
count=1;
for ises=1:length(sessnums)
    sessrow=find(ndata(:,sesscol)==sessnums(ises));
    if cleoflag
        sessrow=sessrow+1;
    end
    for ich=1:length(chids)
        %find most recent row where lfp ch not recorded by da system
        %then find previous to that the siteID for the given ch name
        chnames{1}=chids{ich};
        num=regexp(chnames{1},'\d+','match');%check if ch name contains digit, otherwise a phys channel
        if isempty(num)
            %If phys channel (e.g. eyex), then just assign and continue
            sitenamefull=chnames{1};
            sites(count).sessnum=sessnums(ises);
            sites(count).probeid=chids{ich};
            sites(count).site=sitenamefull;
            sites(count).ch=sitecol;        %chnum original
            count=count+1;
            continue
        end
        if contains(chids{ich},'-')
            %2 part name
            mid=strfind(chids{ich},'-');
            chnames{1}=chids{ich}(1:mid-1);
            chnames{2}=chids{ich}(mid+1:end);
        end
        if ~any(contains(txtdata(sessrow,chcols),chnames{1})) && ...
                ~any(contains(txtdata(sessrow,chcols),chnames{2}))
            %if not a da recording site for current session        
        sitename={};
        for iich=1:length(chnames)
           % [siterows,sitecols2]=find(strcmp(txtdata,chnames{iich}));
            [siterows sitecols2]=find(contains(txtdata(:,1:6),chnames{iich}));     %switch to contains to include pl1x also as contians depth info 05/22
            siterows=intersect(siterows,validrows);
            if ~isempty(siterows)
                rowdiffs=sessrow-siterows;
                [~,nearestrowid]=min(rowdiffs(rowdiffs>0));       %find nearest past session to get id
                nearestrow=siterows(nearestrowid);
                %sitecol=find(strcmp(txtdata(nearestrow,chcols),chnames{iich})==1);
                sitecol=find(contains(txtdata(nearestrow,chcols),chnames{iich})==1);    %Also modified 2022
                sitename{iich}=txtdata(nearestrow,sitecols(sitecol));
            end
        end
        sitenamefull=sitename{1}{:};
        if length(chnames)>1
            sitenamefull=[sitename{1}{:} '-' sitename{2}{:}];
        end
        sites(count).sessnum=sessnums(ises);
        sites(count).probeid=chids{ich};
        sites(count).site=sitenamefull;
        sites(count).ch=sitecol;        %chnum original
        count=count+1;
        end
    end

end

end
