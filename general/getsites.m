function sites=getsites(sessnums,chids,varargin)
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
[~,sheetsav] = xlsfinfo(sitelist); %find # sheets in worksheet
sheetnum=find(strcmp(sheetsav,'sitemap')==1);
[ndata,txtdata,xlsdata] = xlsread(sitelist,sheetnum);
chcols=3:6;
sitecols=7:10;
datecol=1;
sesscol=2;
initialdayrow=4;
if cleoflag
    chcol1=find(strcmp(xlsdata(2,:),'ch1')==1);
    chcols=chcol1:chcol1+3;
    sitecols=chcols(end)+1:chcols(end)+4;
    initialdayrow=5;
end
sites={};
count=1;
cleotoggle=0;
for ises=1:length(sessnums)
    if sessnums(ises)<30 && ises>1 && cleotoggle==0
        %switch to cleo, assume rest of #'s cleo
        sitelist='Z:\inj-monkey2\analysis\cleo\cleo_map_11152019.xlsx';
        [~,sheetsav] = xlsfinfo(sitelist); %find # sheets in worksheet
        sheetnum=find(strcmp(sheetsav,'sitemap')==1);
        [ndata,txtdata,xlsdata] = xlsread(sitelist,sheetnum);
        chcol1=find(strcmp(xlsdata(2,:),'ch1')==1);
        chcols=chcol1:chcol1+3;
        sitecols=chcols(end)+1:chcols(end)+4;
        initialdayrow=5;
        cleotoggle=1;
        cleoflag=1;
    end
    
    sessrow=find(ndata(:,sesscol)==sessnums(ises));
    if cleoflag
        sessrow=sessrow+1;
    end
    
    for ich=1:length(chids)
        sitecol=find(strcmp(txtdata(sessrow,chcols),chids{ich})==1);
        sensortype='stp';
        if ~isempty(sitecol)
        sitename=txtdata(sessrow,sitecols(sitecol));
        origdate='';
        if cleoflag
        if contains(sitename,'p10') 
            origdate='04/10/2017';
            sensortype='uip';
        elseif contains(sitename,'p11')
            origdate='04/12/2017';
            sensortype='uip';
        elseif contains(sitename,'p12')
            origdate='04/14/2017';
            sensortype='uip';
        end
        end
        curdate=txtdata(sessrow,datecol);
        daynum=datenum(curdate,'mm/dd/yyyy');
        if isempty(origdate)
            origdate=txtdata(initialdayrow,datecol);
        end
        if ~cleoflag && (strcmp(chids{ich},'p3') || strcmp(chids{ich},'p5') || strcmp(chids{ich},'cl3'))
            sensortype='uip';
        end
        relday=daynum-datenum(origdate,'mm/dd/yyyy');
        if ~contains(sitename,'x')
        sites(count).sessnum=sessnums(ises);
        sites(count).probeid=chids{ich};
        sites(count).site=sitename{:};
        sites(count).ch=sitecol;        %chnum original
        sites(count).date=curdate;
        sites(count).day=relday;
        sites(count).type=sensortype;
        count=count+1;
        end
        end
    end
end

end
