function CM_getfileinfo(varargin)
global hgui
hgui.sessionnum='';

% CM - separator for mac OR pc
namy = strfind(hgui.PathName, ("/"|"\") );

hgui.ephysid='';
hgui.trialtype=hgui.PathName(namy(end-1)+1:namy(end)-1);
if strcmpi(hgui.trialtype,'cvtotxt')
    %not behavior file, maybe raw data recorded so no name
    hgui.trialtype=hgui.FileName;
end
s1=strfind(hgui.PathName,'chronic');
ss1=strfind(hgui.PathName,'_');
if ~isempty(s1)
    s2=ss1(ss1>s1);
    hgui.session=hgui.PathName(s1:s2(1)-1);
    hgui.sessionnum=str2num(hgui.session(regexpi(hgui.session,'\d')));
    if isempty(str2num(hgui.session(end)))
        %character at end of sessid meaning ephys
        hgui.ephysid=hgui.session(end);
    end
    s3=namy(namy>s2(1));
    hgui.date=hgui.PathName(s2(1)+1:s3(1)-1);
end
hgui.dirIndex=[]; isFSCV=contains(hgui.FileName,'fscv','ignorecase',1);
isTXT=contains(hgui.FileName,'TXT','ignorecase',1);
hgui.trialnum=0;
%check if converted behavior format
if isTXT==0
   % hgui.dirIndex=dir([hgui.PathName '*.mat']);
    hgui.dirIndex=dir([hgui.PathName 'fscv_multi_*']);
    if isempty(hgui.dirIndex)
        %look for other mat files besides fscv_multi
        hgui.dirIndex=dir([hgui.PathName '*.mat']);
    end
    iscsc=contains(hgui.FileName,'csc','ignorecase',1);
    %if is csc split file in directory search for fscv-multi w/ same trial
    if iscsc
        %replace csc name with matching fscv_multi name/file in directory
        trialid=strfind(hgui.FileName,'_');
        trialnum=hgui.FileName(trialid+1:end);
        matchfileid=regexpi({hgui.dirIndex.name},['fscv_multi_' trialnum]);
        matchfileid=~cellfun('isempty',matchfileid);       %make array of zeros/1 for matching
        hgui.FileName=hgui.dirIndex(find(matchfileid==1)).name; %replace name
    end
        namx=strfind(hgui.FileName,'.mat');
        hgui.trialnum=hgui.FileName(namx-3:namx-1); 
else
%check if raw cvtotxt data recorded so no name
    hgui.dirIndex=dir([hgui.PathName '*_CVtoTXT']);
    hgui.trialnum=[];
    hgui.trialtype=hgui.FileName;
end

%check if original tarheel format file, need to modify pathname/dirindex
D=dir(hgui.PathName);
%determine if original tarheel file format selected
%if there is a .txt file in same folder with same name
istar=find(ismember({D.name},[hgui.FileName '.txt'])==1);
if ~isempty(istar)
    %go to \cvtotxt folder
    hgui.PathName=[hgui.PathName 'cvtotxt\'];
    D=dir(hgui.PathName);
    dirid=find(ismember({D.name},['0_' hgui.FileName '_cvtotxt'])==1);
    if isempty(dirid)
        %try with caps
        dirid=find(ismember({D.name},['0_' hgui.FileName '_CVtoTXT'])==1);
        hgui.FileName=['0_' hgui.FileName '_CVtoTXT'];
    else
        hgui.FileName=['0_' hgui.FileName '_cvtotxt'];
    end
    hgui.dirIndex=dir([hgui.PathName '*_CVtoTXT']);
    hgui.trialnum=[];
    hgui.trialtype=hgui.FileName;
end

hgui.currentIndex=1;
for ii=1:length(hgui.dirIndex)
    aa=strcmp(hgui.dirIndex(ii).name,hgui.FileName);
    if aa==1
        hgui.currentIndex=ii;
    end
end

%get subject name
iscleo=~isempty(strfind(hgui.PathName,'cleo'));
ispatra=~isempty(strfind(hgui.PathName,'patra'));
iscfmea=~isempty(strfind(hgui.PathName,'cfmea'));
if iscleo
    hgui.subject='cleo';
elseif ispatra
    hgui.subject='patra';
elseif iscfmea
    hgui.subject='cfmea';
    endnameid=strfind(hgui.FileName,'_');
    hgui.session=hgui.FileName(1:endnameid-1);
    hgui.sessionnum=hgui.session;
else
    hgui.subject='';
    hgui.subject='patra';   %07/02/2021
end
chmap=[];
if strcmp(hgui.subject,'cleo') && ~isempty(hgui.sessionnum)    
    configname=['cleo_chronic' num2str(hgui.sessionnum) '.m'];
    if exist(configname,'file')>0
        run(configname);    %get ncschannels for this session
        chmap=ncschannels;
    end
elseif strcmp(hgui.subject,'patra') && ~isempty(hgui.sessionnum)  
    configname=['chronic' num2str(hgui.sessionnum) 'chconfig.m'];
    if exist(configname,'file')>0
        run(configname);    %get ncschannels for this session
        chmap=ncschannels;
    end
end
hgui.chmap=chmap;
end