function synclfp(subject,sessnum,align,varargin)
%just lfp signals, make artificial fscv folder with letter after num
%not made for whole file recordings or fixed intervals--option still needs
%to be made/fixed noticed 05/2021

log=[];
subjectname='patra';
if strcmp(subject,'cleo')
    subjectname='cleo';
end
sessid='';
if ~isnumeric(sessnum)
    sessid=sessnum;
else
    sessid=num2str(sessnum);
end
%save path in separate directory from data
pctype=computer;
ispc=strcmpi(pctype,'pcwin64');
homedir='Z:';
if ~ispc
    homedir=[filesep];
end
graserver='inj-monkey2';
fscvdir='patra_fscv2';
analysispath=fullfile(homedir,graserver,'analysis',subjectname,['chronic' sessid],filesep);
if ~isdir(analysispath)
    mkdir(analysispath);
end
%default on putamen pc in lab, get data directory from config file
configdir='A:\mit\injectrode\experiments\fscv\matlab\analysis\analysis\config\';
if ~ispc
    %chunky dir
    configdir=fullfile(filesep,'home','schwerdt','matlab','analysis','analysis','config',filesep);
end
ncschannels={};
csc_map={};
paths={};
event_codes={};
switch subjectname
    case 'patra'
        run([configdir ['chronic' sessid 'chconfig.m']]);        %get 'paths'
        run([configdir 'patra_map_bipolar.m']);      %get csc_map
    case 'cleo'
        run([configdir ['cleo_chronic' sessid '.m']]);        %get 'paths'
        run([configdir ['cleo_map_bipolar' sessid '.m']]);      %get csc_map & event_codes
end

%selectcsc=cscchs;
alignname='';
intervalsplit=30;       %30 minutes default split for fixedintervals
splitfiles=0;           %default merge files
mergeflag=0;
aligns={};      %multiple alignment events
switch align
    case 'smallrew'
        alignname='smallreward';
    case 'bigrew'
        alignname='bigreward';
    case 'targetbreak'
        alignname='targetbreak';
    case 'fixbreak'
        alignname='fixbreak';
    case 'fixedint'
        alignname='fixedintervals';
    case 'all'
        alignname='bigreward';       %do all
        aligns={'bigreward','smallreward','targetbreak','fixbreak'};
    otherwise
        alignname='bigreward';       %do all
        aligns={'bigreward','smallreward','targetbreak','fixbreak'};
        %default
end
if strcmpi(alignname,'fixedintervals')
    %get interval minutes that want to split into
    %next argument after should be given
    if ~isnumeric(varargin{1})
        error('interval not given in argument');
    end
    intervalsplit=varargin{1};
end
if isempty(intervalsplit)
    %if empty argument, then means entire file
    mergeflag=1;
end
argnum=1;
cleolickflag=0;
addlfps=0;
cscadd={};
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'merge'
            %merge entire recording into 1 file
            mergeflag=1;
        case 'cleolick37'
            cleolickflag=1;     %for chronic 25 where lick rec
        case 'addlfps'
        %just add csc channels to existing folder already sync sigs da
        addlfps=1;
        argnum=argnum+1;
        cscadd=varargin{argnum};
    end
    argnum=argnum+1;
end
aligns{1}=alignname;
[fscvSyncTTL, nlxSyncID, alignmentIDs]=getAlign(alignname);
%fscvSyncTTL/nlsSyncID is universal sync event trigger between nlx & fscv independent
%of behavioral alignment event

if cleolickflag==0
    [lfpchs,bipolarchs,otherchs,sumchs]=getcscids(ncschannels,csc_map);
else
    [lfpchs,bipolarchs,otherchs,sumchs]=getcscids(ncschannels,csc_map,'cleolick37');
end
selectcsc=[lfpchs otherchs];
%selectcsc=[33:35];          %nlx channels that want to store, only eye
preOffset=30;           %seconds from event ID to align to
durationTrial=60;       %seconds duration of file past start trial index, if zero, then up to next trial event
fscvRate=10;            %fscv sampling rate

%get fscv file names from defined path
if ~isdir(paths{1})
    %fscv folder may not exist since only lfps
    mkdir(paths{1});
end
dirfscv=dir([paths{1} '*_cvtotxt']);       %get cvtotxt files in path
filesfscv={dirfscv.name};
pathfscv=paths{1};
sep=findstr(filesep,pathfscv);      %indexes of separator '\' or '/' in pathname

%get events file from path{2}
direvents=dir([paths{2} '*.nev']);
fileevents=[paths{2} direvents.name];       %get events file name
%load neuralynx events file (mat converted from dg_nlx2mat.m
%load dg_Nlx2Mat_EventStrings, dg_Nlx2Mat_Timestamps, dg_Nlx2Mat_TTL
dg_Nlx2Mat(fileevents); 
load([[fileevents(1:end-4)] '.mat']);
nlx_ts=dg_Nlx2Mat_Timestamps.*1e-6;      %convert to sec (timestamps for events)

%load lfp samples & eye data
pathnlx=paths{2};
sep2=findstr(filesep,pathnlx);
[samples, nlxFileNames, TS]=ephys_getCSC(selectcsc,pathnlx);
%TS here is vector nlx timestamps in seconds
ratelfp=getSampleRate(TS);

%% get filenames
savefolder=fullfile(pathfscv,'matlab',filesep);
if ~isdir(savefolder)
    status = mkdir(savefolder);
end


%find Nlx IDs when evt signals alignment ID for trigger
idsTrigNlx=find(dg_Nlx2Mat_TTL==nlxSyncID);   
trigNlxTTLTS=nlx_ts(idsTrigNlx);        %convert NLx data point to Nlx TS

%%
%repeat this section of script for different behavior alignment events
%skipping previous sections so that do not need to reprocess csc/fscv raw
for ialign=1:length(aligns)
    alignname=aligns{ialign};
disp('splitting & saving data')
disp(alignname)
[fscvSyncTTL, nlxSyncID, alignmentIDs]=getAlign(alignname);

fscvTargetTrialsIDs=[];
nlxTargetTSs=[];
lfp_tID=[];
if ~strcmp(alignname,'fixedintervals')
    %split event data and recorded data based on specific trial event id
    %find fscv TS IDs of targeted event IDs
    [nlxTargetIDs,nlxTargetIDs2]=find(ismember(dg_Nlx2Mat_TTL, alignmentIDs )~=0);
    nlxTargetIDs2=sort(nlxTargetIDs2);  nlxTargetIDs=nlxTargetIDs2(1);
    %remove repeats
    for ii=2:length(nlxTargetIDs2)
        if nlxTargetIDs2(ii)~=nlxTargetIDs2(ii-1)+1
            nlxTargetIDs=[nlxTargetIDs nlxTargetIDs2(ii)];
        end
    end
    nlxTargetTSs2=nlx_ts(nlxTargetIDs);   
    nlxTargetTSs=nlxTargetTSs2;
    %Get corresponding TS ID for Nlx ephys/eye recording samples that will be
    %extracted as needed from samples variable loaded previously
    for ii=1:length(nlxTargetTSs)
        [mins,mintsdif]=min(abs(TS-nlxTargetTSs(ii)));      %based on minimum distance from targeted TS from all TS's to get actual TSid
        if ~isempty(mintsdif)
        %tTemp=find(TS<=nlxTargetTSs(ii)+.0001 & TS>=nlxTargetTSs(ii)-0.0001);
        lfp_tID(ii)=mintsdif;
        else
            lfp_tID(ii)=nan;
        end
    end

end

%%
%Get samples offset from trial targeted Sync event (ie. how much to store)
nlxSamplesOffset=ceil(preOffset*ratelfp);
numLFPch=length(samples);
%%%%%%%%%%%%%%%%%%%%%%%%%
%Split fscv & ephys data based on trial sequences
syncData={};    
syncLFP={};
tsLFP=[];
NlxEventTTL=[];
NlxEventTS=[];
fscvData=[];
fscvEvents=[];
%newfolder=[pathfscv 'matlab\' alignname '\' ];
if mergeflag==1
    alignname='merged';
end
newfolder=fullfile(pathfscv,'matlab',alignname,filesep);
if ~isdir(newfolder)
    status = mkdir(newfolder);
end
PathName4=newfolder;
initialNum=100;
previousID=0;
log=[];
totaltrials=length(nlxTargetTSs);
totaltrials=length(lfp_tID);

if durationTrial==0
    %fixed interval condition, just 1 less
    totaltrials=totaltrials-1;
end
for ii=1:totaltrials
    disp(['split trial # ' ...
        num2str(ii) '/' num2str(totaltrials)]);
    if lfp_tID(ii)==previousID
        continue
    end
    %get start IDs for fscv & nlx samples to capture
    currentLFPIndex=lfp_tID(ii)-nlxSamplesOffset;
    %check if currentLFPIndex>0 07/31/2018
    if currentLFPIndex<1 
        warning('first trial not included because index <1');
        continue
    end
    if durationTrial==0
        %fixed interval condition
        if ii==length(lfp_tID)
            if ii>length(lfp_tID)
                nextLFPIndex=length(TS)-1;
            else
                nextLFPIndex=lfp_tID(ii+1)-1;
            end
        else
        nextLFPIndex=lfp_tID(ii+1)-1;
        end
    else
        nextLFPIndex=currentLFPIndex+ceil(durationTrial*ratelfp);
    end
    if nextLFPIndex<length(TS)        
        for LFPch=1:numLFPch
            syncLFP{LFPch}=samples{LFPch}(currentLFPIndex:nextLFPIndex);
        end
        tsLFP=TS(currentLFPIndex:nextLFPIndex);
        nlxTargetIDsWithinTrial=find(nlx_ts>=tsLFP(1) & nlx_ts<=tsLFP(end));
        NlxEventTTL=dg_Nlx2Mat_TTL(nlxTargetIDsWithinTrial);
        NlxEventTS=nlx_ts(nlxTargetIDsWithinTrial);        
    else
        %get to end of recording       
            nextidxlfp=nextLFPIndex;
            if nextidxlfp>length(samples{1})
                nextidxlfp=length(samples{1});
            end
            for LFPch=1:numLFPch
                syncLFP{LFPch}=samples{LFPch}(currentLFPIndex:nextidxlfp);
            end
            tsLFP=TS(currentLFPIndex:nextidxlfp);
            nexttargetid=nlxTargetIDs(ii)+100;
            if nexttargetid>length(dg_Nlx2Mat_TTL)
                nexttargetid=length(dg_Nlx2Mat_TTL);
            end
            NlxEventTTL=dg_Nlx2Mat_TTL(nlxTargetIDs(ii):nexttargetid);
            NlxEventTS=nlx_ts(nlxTargetIDs(ii):nexttargetid);
    end
    if lfp_tID(ii)<=length(TS)
        saveName=['fscv_multi_' num2str(initialNum)];        
        samplesLFP=syncLFP;
        if splitfiles==0
            %merge all signals into one file
            save([PathName4, saveName],'samplesLFP','tsLFP','NlxEventTTL','NlxEventTS','nlxFileNames','-v7.3');
            prompt = ['saved as: ' saveName];
        end
        initialNum=initialNum+1;
    end
   % fprintf('%d \r', floor(ii/length(fscvTargetTrialsIDs)*100)); 
    previousID=lfp_tID(ii);
    if mergeflag==1
        %only do one run for entire recording
        break
    end
end
savelogname=['log_getFSCV'];
save([PathName4, savelogname],'log');
end
%save log of latencies
end