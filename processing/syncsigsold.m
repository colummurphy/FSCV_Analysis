function syncsigs(paths,ncschannels,csc_map,align,varargin)
%getFSCV 09/14/2018
%combine CSC/FSCV recordings and organize by trials
%load LFP, fscv, eye data and split according to behavior nlx events
%according to targeted IDs defined in getAlign
%Get targeted alignment ID's, ie event codes that each exported trial will
%be synchronized to - preOffset with entire trial length = durationTrial
%paths{1}= directory with cvtotxt files to merge and convert
%paths{2}= directory with nlx events and csc files to merge & convert
%paths might be stored in chronicXXchconfig
%align = name ' ' of align event ie. 'smallrew' 'bigrew' 'fixedint'
%'targetbreak' 
%ncs channels has the names of sites that want to retrieve for a given
%session, get by loading eg. chronic67chconfig
%csc_map (from patra_map_bipolar_chunky) has the site defined by csc ids
%need csc ids to load appropriate ncs files

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
        aligns={'bigreward','smallreward','targetbreak','fixbreak','fixedintervals'};

    otherwise
        alignname='fixedintervals';
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
        case 'split'
            %split files
            splitfiles=1;       %DO NOT USE, reconvert will split
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
%alignname='smallReward';
%alignname='bigReward';
%alignname='targetBreak';
%alignname='fixedintervals';  
aligns{1}=alignname;
[fscvSyncTTL, nlxSyncID, alignmentIDs]=getAlign(alignname);
%fscvSyncTTL/nlsSyncID is universal sync event trigger between nlx & fscv independent
%of behavioral alignment event

%selectcsc=[1:16, 33:35];    %nlx channels that want to store cleo
%selectcsc=[1:24, 33:37];    %nlx channels that want to store patra
%selectcsc=[1:24, 33:41];    %nlx channels that want to store patra > day28, after new lick fix dc input
%selectcsc=[33:41];    %just eye/lick data
%selectcsc=[1:16, 26:32, 33:35, 39:41, 42];    %Corrected 03/14 nlx ephys, eye, lick, pulse
%lfpCh=[2:15, 26:32];    %corrected 3/14
%lfpCh=[3, 4, 5, 11, 14, 15, 26:32];    %chronic30
%lfpCh=[3, 4, 5, 8,9, 11, 14, 15, 27:32];    %chronic38
%lfpCh=[3, 5, 9, 10, 11, 14, 15, 26:32];    %chronic54
%lfpCh=[2:5, 7, 9, 10, 11, 14, 15, 27:32];    %chronic58
%lfpCh=[2:5, 7, 8,  11, 12, 14, 26:32];    %chronic60-67, 83
%lfpCh=[2,3,5, 8,9,  11, 14, 15, 27:32];    %chronic31
%lfpCh=[2:4, 9, 10, 11, 14, 15, 27:32];    %chronic71
%lfpCh=[2,3, 7, 8, 10, 11, 14, 15, 26:32];    %chronic46-45
if cleolickflag==0
    [lfpchs,bipolarchs,otherchs,sumchs]=getcscids(ncschannels,csc_map);
else
    [lfpchs,bipolarchs,otherchs,sumchs]=getcscids(ncschannels,csc_map,'cleolick37');
end
%otherCh=[33:35,39:41,42];   %eye, lick, pulse
selectcsc=[lfpchs otherchs];
%selectcsc=[33:35];          %nlx channels that want to store, only eye
preOffset=30;           %seconds from event ID to align to
durationTrial=60;       %seconds duration of file past start trial index, if zero, then up to next trial event
fscvRate=10;            %fscv sampling rate

%get fscv file names from defined path
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

Iread=[];
fscvStartTrialTS=[];        %fscv domain TS of cheetah start TTLs
%find first file & only include text files for each channel
%find first recording in sequence (ie. 100)
%ispresent=contains(filesfscv,'100_cvtotxt','IgnoreCase',true);
%contains does not work in 2013 version matlab
ispresent=strfind(filesfscv,'100_cvtotxt');     %case-sensitive
firstFileIDs=find(~cellfun(@isempty,ispresent));
%contains function does not work for 2013 for cells 
%files=filenames(processfiles);

%firstFileIDs=find(ispresent>0);         %find first rec files
if isempty(firstFileIDs)
    warning('file with 100 not found, finding 200 for first idx file');
    %maybe first file named was 200, check
   % ispresent=contains(files,'200_cvtotxt','IgnoreCase',true);
    ispresent=strfind(filesfscv,'200_cvtotxt');
    firstFileIDs=find(~cellfun(@isempty,ispresent));
    %firstFileIDs=find(ispresent>0);
end
firstFileName=filesfscv{1};
%check to make sure this first file in list is also the first channel
%default naming convention is '0_1dr_xx_100_cvtotxt'
correspondsCh0=regexpi(firstFileName,'0_');
firstch=0;      %default 1st channel for 4 channel fscv recording
if correspondsCh0(1)~=1
    warning('first channel not included');
    if find(regexpi(firstFileName,'\d')==1)
        %check that first character of name is a number corresponding to ch
        firstch=str2num(firstFileName(1));  %get first ch #
        display(['first fscv ch is #' num2str(firstch)]);
    else
        error('incorrect naming or format of fscv cvtotxt files')
    end
        
end

%% get filenames

%firstChannelFiles=contains(filesfscv,[num2str(firstch) '_1dr'],'IgnoreCase',true);
    firstChannelFiles=strfind(filesfscv,[num2str(firstch) '_1dr']);
    fileNamesids=find(~cellfun(@isempty,firstChannelFiles));
fileNames=filesfscv(fileNamesids);
%fileNames=filesfscv(firstChannelFiles>0);       %cvtotxt files of 1st channel only
arrayedFiles =vertcat(fileNames{:});        %only works when length of strings same
arrayedNames=arrayedFiles;
endNameIdx=regexpi(arrayedNames(1,:),'_cvtotxt');
targetNames=arrayedNames(:,3:endNameIdx-1);     %name without channel #
endName=arrayedNames(1,endNameIdx:end);
fileNums=targetNames(:,end-2:end);
%file #'s of continuous fscv recordings, should start from 100 typically
fileNums2=str2num(fileNums);            

%sort filenames
[sorted,sortID]=sort(fileNums2,'ascend');
sortedNames=targetNames(sortID,:);

%check if already sorted (not really needed..)
xx=sorted-fileNums2;
if sum(xx)==0
    display('names already sorted as loaded');
end

%get .txt files up one directory with TTls from neuralynx
pathNameUp=pathfscv(1:sep(end-1));      %directory with fscv ttl files

fscvEvents=[];
%load and concatenate .txt files containing TTLs from cheetah/vcortex
disp('concatenating fscv TTL files');
lastTS=0;
for fileidx=1:size(sortedNames,1)
    %load txt file with fscv ttls
    txtLoad=load([pathNameUp sortedNames(fileidx,:) '.txt']);
    
    %store in fscvEvents
    if isempty(fscvEvents)
        fscvEvents=txtLoad;
    else
        %shift TS's upward depending on last TS recorded from last file
        % and fscv time interval (eg. 1/fscvRate)
        txtLoad2=txtLoad;
        txtLoad2(:,1)=txtLoad2(:,1)+lastTS+1/fscvRate;
        fscvEvents=[fscvEvents; txtLoad2];
    end
    
    %get last recorded TS for this file
    lastTS=fscvEvents(end,1);
end
disp(['recorded ' num2str(lastTS/60/60) ' hours']);

%Find fscv data point when TTL goes to 241 in FSCV, corresponding to Nlx first trigger
fscvTTLIdx=find(fscvEvents(:,3)==fscvSyncTTL);  % get TTL on Ids
fscvTTLTS=fscvEvents(fscvTTLIdx,1);         %get fscv TS of TTL on periods
%nlx_ts is the timestamps of all TTL's inputted/outputted in Neuralynx
%actual TTL values in dg_Nlx2Mat_TTL. Neuralynx has much faster sampling 
%rate but only generates # once so the first encoded TS in Neuralynx for 
%outputted TTL represents several successive TTL TS's in fscvEvents FSCV 
%TSs (since TTL is several hundred ms long)

%find Nlx IDs when evt signals alignment ID for trigger
idsTrigNlx=find(dg_Nlx2Mat_TTL==nlxSyncID);   
trigNlxTTLTS=nlx_ts(idsTrigNlx);        %convert NLx data point to Nlx TS

%calculate time difference between Nlx first TTL TS (ie. cheetah time base)
% and fscv first TTL TS (ie. fscv time base)
fscvToNlxOffset=trigNlxTTLTS(1)-fscvTTLTS(1);
fscvNlxEvents=fscvEvents;       
%fscvNlxEvents has timestamps for fscv signals and cheetah events
%create same events file with fscv TS shifted to synchronize with Neuralynx TS
%fscvNlxEvents(:,1)=fscvEvents(:,1)+fscvToNlxOffset;
%remove successive TTL's (ie. hold period) in fscv TTLs corresponding 
%to same sync ttl pulse from nlx (~ 200 - 300 ms?, so repeated 2 - 3 times)
%to make sure same amount of TTL's in cheetah, can check delays
fscvTTLIdxCompressed=[];
IDcount=1;
previdx=-1;
for idx=1:length(fscvTTLIdx)
    if fscvTTLIdx(idx)-1~=previdx
        fscvTTLIdxCompressed(IDcount)=fscvTTLIdx(idx);
        IDcount=IDcount+1;
    end
    previdx=fscvTTLIdx(idx);
end
%get fscv TS of TTL compressed IDs  matching length of idsTrigNlx
syncfscvTTLTS=fscvNlxEvents(fscvTTLIdxCompressed,1);    
if length(syncfscvTTLTS)~=length(idsTrigNlx)
    warning('length of TTLs activated in fscv do not match TTLs in Nlx')
    %check if period between first two ttl's match for fscv/nlx systems
    interinterval_fscv=syncfscvTTLTS(2)-syncfscvTTLTS(1);
    interinterval_nlx=trigNlxTTLTS(2)-trigNlxTTLTS(1);
    %3/1/2019, add checks for making sure signal events aligned, based on
    %errors in cleo fscv 17, 20 lines below added
    difffscv=diff(syncfscvTTLTS');
    diffnlx=diff(trigNlxTTLTS);
    difffscv2=difffscv(1:min(length(difffscv),length(diffnlx)));
    diffnlx2=diffnlx(1:min(length(difffscv),length(diffnlx)));
    difboth=difffscv2-diffnlx2;
    [r,p]=corr(difffscv2',diffnlx2');           %correlatiosn betweetn intervals between both
    obstaclefound=find(abs(difboth)>5);     %discrepancy between interinterval of fscv ttl's and nlx ttl's > 5 s (1 s maybe understandable for prolonged breaks)
    if length(obstaclefound)>10 && r<0.95
        %will create multiple insychronies afterwards
        %def marker of something wrong
        if length(idsTrigNlx)>length(fscvTTLIdxCompressed)  
            disp('nlx trig events more than fscv events captured')
            diffnlx3=diffnlx;
            diffnlx3(obstaclefound(1))=[];      %remove proble nlx ts trigger that was not captured in fscv & make sure signal events now aligned based on corr below
            [r,p]=corr(difffscv2',diffnlx3');
            if r>0.95
                disp(['removing obstacle at ' num2str(obstaclefound(1)) ' trigger #'])
                %if removal of obstalce restores signal alignments                
                idsTrigNlx(obstaclefound(1))=[];        %remove problem nlx trigger that was not captured in fscv
                trigNlxTTLTS=nlx_ts(idsTrigNlx);
                interinterval_nlx=trigNlxTTLTS(2)-trigNlxTTLTS(1);
                if length(syncfscvTTLTS)==length(idsTrigNlx)
                    disp(['length of triggers equal, fscv ttls = ' num2str(length(syncfscvTTLTS)) ...
                        ', nlx ttls = ' num2str(length(idsTrigNlx))]);
                end
            end
        else
            warning('fscv trig events more than nlx events triggered, nothing to do, need to manually fix');
        end
    end
    if abs(interinterval_fscv-interinterval_nlx)>0.2
        warning(['discrepancy between period from trial 1 to trial 2 for 2 systems is > 0.2 s = ' num2str(abs(interinterval_fscv-interinterval_nlx))]);
    end    
end
%Because of slow fscv sampling rate & long pulse can be large discrepancy
%between TTL TS's since fscv system has inherent delay that gets progressively
%larger as recording continues.....................longer.............
%Therefore, Nlx referred TS's created in fscvNlxEvents are useless
%need to refer each TTL in fscv itself to corresponding TTL in nlx to
%synchronize trial task events accurately
%% 
%now incorporate NLx events into fscvEvents around each sync trigger TTL
%as well as nonlinear merged time Nlx TS's
fscvNlxEvents(10,12)=0;     %expand to have 10 columns (ie 5 columns for events)
eventsininterval=[];
eventspreinterval=[];
eventsinterval=[];
%scan at each triggered FSCV TTL
for ii=1:length(fscvTTLIdxCompressed)
    %get # of events encoded in Nlx between trials
    %events (idxs) from current trigger to preOffset period prior 
    eventspreinterval=find(nlx_ts...
        <=nlx_ts(idsTrigNlx(ii))-preOffset);
    if ~isempty(eventspreinterval)
        eventspreinterval=eventspreinterval(end):(idsTrigNlx(ii)-1);
    else
        eventspreinterval=1:(idsTrigNlx(ii)-1);
    end
    %events (idxs) from current trigger until next trigger
    if ii~=length(idsTrigNlx)
        eventspostinterval=idsTrigNlx(ii):(idsTrigNlx(ii+1));  
    else
        eventspostinterval=idsTrigNlx(ii):length(dg_Nlx2Mat_TTL);
    end
    %eventspreinterval redundant unless > preOffset period since last
    %trigger (duration that could produce mismatch) or first trigger
    timelasttrigger=60;
    if ii~=1
        timelasttrigger=nlx_ts(idsTrigNlx(ii))-...
            nlx_ts(idsTrigNlx(ii-1));
    end
    if timelasttrigger>preOffset && ii~=1
        %only use preinterval events if long delay between triggers
        eventsinterval=[eventspreinterval eventspostinterval];
    else
        %otherwise only use post events
        eventsinterval=eventspostinterval;
    end
    if ii==1
        %if first trigge,r get everything before
        eventsinterval=[eventspreinterval eventspostinterval];
    end
    %get range of nlx TS's between trigger intervals
    firstTS=nlx_ts(eventsinterval(1)); 
    lastTS=nlx_ts(eventsinterval(end));
    %get reference idx of start TTL trigger in fscv
    reffscvID=fscvTTLIdxCompressed(ii);     %first idx of trigger in fscv
    refNlxTS=nlx_ts(idsTrigNlx(ii));     %Nlx TS of trigger
    %get firstTS id relative to fscv trigger
    td=firstTS-refNlxTS;        
    tdinsamplesFSCV=round(td*fscvRate);        %convert to FSCV samples
    firstTSfscvIDx=reffscvID+tdinsamplesFSCV;
    TSrangesamplesFSCV=round(firstTS*fscvRate):round(lastTS*fscvRate);
    TSrangesamplesFSCV=TSrangesamplesFSCV./fscvRate;
    %store nlx TS range in appropriate row and column 4 in fscv data
    if firstTSfscvIDx<1
        %event ts's in nlx before trigger may occur prior to start of fscv
        %recording
        TSrangesamplesFSCV=TSrangesamplesFSCV(-(firstTSfscvIDx-2):end);
        firstTSfscvIDx=1;
       % eventsinterval=
    end
    %store predicted timestamps for current interval
    fscvNlxEvents(firstTSfscvIDx:firstTSfscvIDx+length(TSrangesamplesFSCV)-1,4)=TSrangesamplesFSCV';
    if ii==length(fscvTTLIdxCompressed)
        %store predicted ts's to end of fscv file
        finalstoredidx=firstTSfscvIDx+length(TSrangesamplesFSCV)-1;
        if finalstoredidx<size(fscvNlxEvents,1)
            prelasttss=fscvNlxEvents(finalstoredidx,4)+1/fscvRate;
            prelastfidx=round(prelasttss*fscvRate);
            remaininglength=size(fscvNlxEvents,1)-finalstoredidx-1;
            idstoend=prelastfidx:prelastfidx+remaininglength;
            tsstoend=idstoend./fscvRate;
            fscvNlxEvents(finalstoredidx+1:end,4)=tsstoend';
        end
    end
    fscvNlxEvents(reffscvID,4)=refNlxTS;       %store Nlx TS for trigger in col 4
    refCode=dg_Nlx2Mat_TTL(idsTrigNlx(ii));         %refCode = nlxSyncID, redundant
    fscvNlxEvents(reffscvID,5)=refCode;          %store event id for trigger Should be 3 in col 5
    %scan all events in trial interval encoded from nlx until next 'trial'
    %write event codes and timestamps into fscvNlxEvents
    %added pre trial events to preOffset period 05/01/2018
    for jj=1:length(eventsinterval)-1
        %get next NLX event code in interval
        eventCodetostore=dg_Nlx2Mat_TTL(eventsinterval(jj));   
        %get its NLX TS
        eventTStostore=nlx_ts(eventsinterval(jj)); 
        %calc time difference between event to store & trigger ref TS
        td=eventTStostore-refNlxTS;        
        tdinsamplesFSCV=round(td*fscvRate);        %convert to FSCV samples
        %store TS of event in appropriate row and column 4
        if reffscvID+tdinsamplesFSCV<1
            continue
        end
        fscvNlxEvents(reffscvID+tdinsamplesFSCV,4)=eventTStostore;
        %find unfilled column to store event code for given TS row
        firstCheck=5;
        if fscvNlxEvents(reffscvID+tdinsamplesFSCV,firstCheck)==0
            fscvNlxEvents(reffscvID+tdinsamplesFSCV,firstCheck)=eventCodetostore;
        else 
        %store event code for same TS in fscvNlxEvents in another column 
        %if targeted column already taken
            while fscvNlxEvents(reffscvID+tdinsamplesFSCV,firstCheck)~=0 && firstCheck<10
                firstCheck=firstCheck+1;
            end
            fscvNlxEvents(reffscvID+tdinsamplesFSCV,firstCheck)=eventCodetostore; 
            %store event in empty space
        end
    end
   
    
end        

%Get FSCV echem data now for each chanel
%break down files again based on start trial, load tarheel text files,
%concatanete data, and break apart based on fscvStartTrialsIDs & include
%all event IDs in exported files
%find # of channels based on selected text files
display('finding # channels')
numFilesPerCh=size(sortedNames,1);
foundch0=strfind(filesfscv,'0_1dr');     %case-sensitive
foundch0=find(~cellfun(@isempty,foundch0));
foundch1=strfind(filesfscv,'1_1dr');     %case-sensitive
foundch1=find(~cellfun(@isempty,foundch1));
foundch2=strfind(filesfscv,'2_1dr');     %case-sensitive
foundch2=find(~cellfun(@isempty,foundch2));
foundch3=strfind(filesfscv,'3_1dr');     %case-sensitive
foundch3=find(~cellfun(@isempty,foundch3));

%foundch0=contains(filesfscv,'0_1dr');  foundch0=find(foundch0>0);
%foundch1=contains(filesfscv,'1_1dr');  foundch1=find(foundch1>0);
%foundch2=contains(filesfscv,'2_1dr');  foundch2=find(foundch2>0);
%foundch3=contains(filesfscv,'3_1dr');  foundch3=find(foundch3>0);
numch=0;
chnums=[];
if length(foundch0)>=numFilesPerCh
    numch=numch+1;
    chnums=[chnums 0];
end
if length(foundch1)>=numFilesPerCh
    numch=numch+1;
    chnums=[chnums 1];
end
if length(foundch2)>=numFilesPerCh
    numch=numch+1;
    chnums=[chnums 2];
end
if length(foundch3)>=numFilesPerCh
    numch=numch+1;
    chnums=[chnums 3];
end
%fixed above 1/18/2018, previosu code meant even if ch0/1 nonexistent but 
%ch3/4 exist counts as 4 channels.. now look at individual chs

%merge FSCV data for each channel and store in concatenated array
%fixed to load based on actual existing ch nums in chnums 01/18/2018
%channels when loaded below normalized so if ch 2 & 3 only recorded,
%becomes ch0 & 1
txtDataCh0=[];
txtDataCh1=[];
txtDataCh2=[];
txtDataCh3=[];
for jj=1:numch
    disp(['merging ch ' num2str(chnums(jj)) ' fscv txt data'])
    for ii=1:numFilesPerCh
        %load cvtotxt file based on chnums
        tempLoad=load([pathfscv num2str(chnums(jj)) '_' sortedNames(ii,:) endName]);
        switch jj
            case 1
                if ii==1
                    txtDataCh0=tempLoad;
                else
                    txtDataCh0=[txtDataCh0; tempLoad];
                end
            case 2
                if ii==1
                    txtDataCh1=tempLoad;
                else
                    txtDataCh1=[txtDataCh1; tempLoad];
                end
           case 3
                if ii==1
                    txtDataCh2=tempLoad;
                else
                    txtDataCh2=[txtDataCh2; tempLoad];
                end
            case 4
                if ii==1
                    txtDataCh3=tempLoad;
                else
                    txtDataCh3=[txtDataCh3; tempLoad];
                end
        end
    end
end
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
    indicesTXT=fscvNlxEvents(:,5:end);      %just to look at event codes
    %find fscv TS IDs of targeted event IDs
    [fscvTargetTrialsIDs,fscvStartTrialsIDs2]=find(ismember(indicesTXT, alignmentIDs )~=0);
    fscvTargetTrialsIDs=sort(fscvTargetTrialsIDs); 
    fscvTargetTrialsIDs2=fscvTargetTrialsIDs(1);
    %remove repeats
    for ii=2:length(fscvTargetTrialsIDs)
        if fscvTargetTrialsIDs(ii)>fscvTargetTrialsIDs(ii-1)+4
            %only include if next id is > 4 samples of previous, otherwise it is repeat
            %03/08/18 changed from 2 samples barrier to 4 samples barrier (ie 400 ms),
            %because event 45 is programmed to appear 200ms after 18/19 ?
            fscvTargetTrialsIDs2=[fscvTargetTrialsIDs2 fscvTargetTrialsIDs(ii)];
        end
    end
    fscvTargetTrialsIDs=fscvTargetTrialsIDs2;
    %find Nlx TS IDs of targeted event IDs
    %nlxTargetIDs will be bigger than fscvTarget.. because dg_Nlx2Mat is
    %searched serially whereas fscvNlxEvents combines several different
    %target IDs in a single fscv "TS" row IDx, as with redundancy above
    [nlxTargetIDs,nlxTargetIDs2]=find(ismember(dg_Nlx2Mat_TTL, alignmentIDs )~=0);
    nlxTargetIDs2=sort(nlxTargetIDs2);  nlxTargetIDs=nlxTargetIDs2(1);
    %remove repeats
    for ii=2:length(nlxTargetIDs2)
        if nlxTargetIDs2(ii)~=nlxTargetIDs2(ii-1)+1
            nlxTargetIDs=[nlxTargetIDs nlxTargetIDs2(ii)];
        end
    end
    nlxTargetTSs2=nlx_ts(nlxTargetIDs);
    %get Neuralynx TS for same Target IDs directly as stored in fscvNlxEvents
    nlxTargetTSs=fscvNlxEvents(fscvTargetTrialsIDs,4);
    if length(nlxTargetIDs)~=length(fscvTargetTrialsIDs)
        warning('# of nlx target ids does not match with fscv target ids');
    end
    discrepancy=0;
    %check if # of target sync id's in nlx is same as fscv system, otherwise 
    %break in recording or something wrong with communication between systems
    if length(nlxTargetTSs)==length(nlxTargetTSs2)
        discrepancy=nlxTargetTSs-nlxTargetTSs2'; %check if Ts's match up, maybe some error < 1ms
    else
        lengthFSCVIDs=length(nlxTargetTSs);
        lengthNlxIDs=length(nlxTargetTSs2);
        discrepancy=nlxTargetTSs(1:min(lengthFSCVIDs,lengthNlxIDs),:)-nlxTargetTSs2(:,1:min(lengthFSCVIDs,lengthNlxIDs))';
        warning(['# of nlx target ids is different from fscv target ids by ' num2str(length(nlxTargetTSs2)-length(nlxTargetTSs))]);
    end
    if any(abs(discrepancy)>0.1)
        trialsLagging=find(abs(discrepancy)>0.1);   %get IDs for trials with significant sync lag
        %check if any ts's stored are more than 0.1 s different between systems
        warning('time discrepancies between saved nlx ID ts and synchronized nlx ids in fscvNlxEvents > 0.1 s');
        display(['lagging trials: ' num2str(trialsLagging','%10.0f\n')]);
        if length(trialsLagging)>15
            disp(['attempt to realign'])
            intervalstsFSCV=diff(nlxTargetTSs(1:min(length(nlxTargetTSs),length(nlxTargetTSs2))));
            intervalstsNlx=diff(nlxTargetTSs2(1:min(length(nlxTargetTSs),length(nlxTargetTSs2))));
            [r,p]=corr(intervalstsFSCV,intervalstsNlx');
            if length(nlxTargetTSs)>length(nlxTargetTSs2) && r<0.95
                %more fscv targ ieds than nlx
                %remove one of fscv target ids causing desychcronization
                if length(nlxTargetTSs)-1==length(nlxTargetTSs2)
                    disp(['more fscv events than nlx, removing trial # ' num2str(trialsLagging(1)) ]);
                    fscvTargetTrialsIDs(trialsLagging(1))=[];
                    nlxTargetTSs=fscvNlxEvents(fscvTargetTrialsIDs,4);
                    intervalstsFSCV=diff(nlxTargetTSs(1:min(length(nlxTargetTSs),length(nlxTargetTSs2))));
                    intervalstsNlx=diff(nlxTargetTSs2(1:min(length(nlxTargetTSs),length(nlxTargetTSs2))));
                    [r,p]=corr(intervalstsFSCV,intervalstsNlx');
                    if r<0.9
                        error(['fscv and nlx still descyrhonized after removing potential fscv obstacle event, need to manually fix']);
                    end
                else
                    error(['more than one fscv event causing desychronization, need to manually fix']);
                end        
            else
                error(['more nlx target events than fscv target events, & more than 15 trials desychronized, need to manually fix']);
            end
        end
    end

    %Get corresponding TS ID for Nlx ephys/eye recording samples that will be
    %extracted as needed from samples variable loaded previously
    %MAY NOT BE EFFICIENT, CHECK
    for ii=1:length(nlxTargetTSs)
        tTemp=find(TS<=nlxTargetTSs(ii)+.0001 & TS>=nlxTargetTSs(ii)-0.0001);
        lfp_tID(ii)=tTemp(1);
    end
else
    preOffset=0;
    durationTrial=0;
    if mergeflag==1
        %whole file
        fscvTargetTrialsIDs=[1 size(fscvNlxEvents,1)];
    else
    %split based on fixed time intervals irregardless of task events    
    intervalinsamples=fscvRate.*intervalsplit.*60;
    totalrecordedsamples=size(fscvNlxEvents,1);
    numsections=ceil(totalrecordedsamples/intervalinsamples);
    firstid=find(fscvNlxEvents(:,4)>=TS(1));
    firstfscvid=firstid(1);     %if fscv started before nlx then will be offset
    firsteffective=ceil(firstfscvid/10)*10;
    if firsteffective <=0
        firsteffective=1;
    end
    intervalids=firsteffective:intervalinsamples-1:totalrecordedsamples;    %offset 1s
    if intervalids(end)+5*fscvRate*60<totalrecordedsamples
        %if more than 5 minutes remaining not accounted for need to keep
        intervalids(end+1)=totalrecordedsamples;
    end
    fscvTargetTrialsIDs=intervalids;
    end
    
    %looking at all data, not just around task events, so 
    %need to fill out all "zero" nlx ts's that were not interpolated
    %around task events
    nlxtss=fscvNlxEvents(:,4);
    disp('check nlx ts integrated w/ fscv');
    
    %aaaa=diff(nlxtss);          %find sig excursions in previously interpolated TS's
    %xxx=find(aaaa>5);
    %nlxtss(xxx)=nan;
    nlxtss(nlxtss==0)=nan;      %any ts's equal to zero are those that were not filled, make nan
    verm=version('-release');
    if str2num(verm(regexp(verm,'\d')))>2016
       [nlxtssfill,tf]=fillmissing(nlxtss,'linear');   %linear interpolate those points
        fscvNlxEvents(:,4)=nlxtssfill;      %replace fscvNlxEvents TSss
        nlxTargetTSs=fscvNlxEvents(fscvTargetTrialsIDs,4);
        %figure; plot(nlxtss,'*')        %PLOT   if not chunky (2013 version)
        %hold on; plot(nlxtssfill)
        %legend('before re-interpolation','after filling');
    else
        error('must use fillmissing from new version of matlab otherwise may cause delays with interp1');
        %DO NOT USE< MAY CAUSE DELAYS IN SYNC 12/2019
        %fillmissing function not available pre 2016b
        nlxtssii = 1:length(nlxtss);
        m = isnan(nlxtss);        
        %s = spline(nlxtssii(~m),nlxtss(~m),nlxtssii(m));
        s = interp1(nlxtssii(~m),nlxtss(~m),nlxtssii(m),'linear','extrap');
        % replace NaN values with interpolated values; plot to see results
        nlxtssfill=nlxtss;
        nlxtssfill(m) = s; 
        fscvNlxEvents(:,4)=nlxtssfill;      %replace fscvNlxEvents TSss
        nlxTargetTSs=fscvNlxEvents(fscvTargetTrialsIDs,4);
    end
    disp(['check re-interpolation to fill gaps not accounted for between'...
        'usually long delays in task']);
    targetTSoverrecord=find(nlxTargetTSs>max(TS));  %if find that target TS over max recorded TS
    %shorten this target TS, assume at end:
    if targetTSoverrecord==length(nlxTargetTSs)
        %nlxTargetTSs(end)=[];
        if nlxTargetTSs(end-1)<TS(end)-ratelfp*120 && fscvNlxEvents(end,4)<TS(end)-ratelfp*120
        nlxTargetTSs(end)=TS(end)-ratelfp*10;
        else
            nlxTargetTSs(end)=[];
        end
        fscvTargetTrialsIDs(end)=[];
    end
    %check if nlxTargetTSs(1) < first timestamp of nlx
    %ie started fscv recording before nlx
    %not needed anymore since ceil above for firsteffective
    %add extra index at end to account for remaining time not covered by
    %intervals split above
    if nlxTargetTSs(end)~=fscvNlxEvents(end,4) && fscvNlxEvents(end,4)<TS(end)
        nlxTargetTSs(end+1)=fscvNlxEvents(end,4);
        fscvTargetTrialsIDs(end+1)=size(fscvNlxEvents,1);
    elseif nlxTargetTSs(end)<TS(end)-ratelfp*120 && fscvNlxEvents(end,4)<TS(end)-ratelfp*120
        %at least 120 s remaining
        nlxTargetTSs(end+1)=TS(end)-ratelfp*10; %final target 10s before end recording
    end
    %Get corresponding TS ID for Nlx samples
    %MAY NOT BE EFFICIENT, CHECK
    for ii=1:length(nlxTargetTSs)
        tTemp=find(TS<=nlxTargetTSs(ii)+.0001 & TS>=nlxTargetTSs(ii)-0.0001);       
        lfp_tID(ii)=tTemp(1);
    end
    
end
%%
%Get samples offset from trial targeted Sync event (ie. how much to store)
fscvSamplesOffset=preOffset*fscvRate;
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
totaltrials=length(fscvTargetTrialsIDs);
if durationTrial==0
    %fixed interval condition, just 1 less
    totaltrials=totaltrials-1;
end
for ii=1:totaltrials
    disp(['split trial # ' ...
        num2str(ii) '/' num2str(totaltrials)]);
    if fscvTargetTrialsIDs(ii)==previousID
        continue
    end
    %get start IDs for fscv & nlx samples to capture
    currentFSCVIndex=fscvTargetTrialsIDs(ii)-fscvSamplesOffset;
    currentLFPIndex=lfp_tID(ii)-nlxSamplesOffset;
    %check if currentLFPIndex>0 07/31/2018
    if currentLFPIndex<1 || currentFSCVIndex<1
        warning('first trial not included because index <1');
        continue
    end
    if durationTrial==0
        %fixed interval condition
        if ii==length(fscvTargetTrialsIDs)
            nextFSCVIndex=length(fscvNlxEvents)-1;
            if ii>length(lfp_tID)
                nextLFPIndex=length(TS)-1;
            else
                nextLFPIndex=lfp_tID(ii+1)-1;
            end
        else
        nextFSCVIndex=(fscvTargetTrialsIDs(ii+1)-1);
        nextLFPIndex=lfp_tID(ii+1)-1;
        end
    else
        nextFSCVIndex=currentFSCVIndex+durationTrial*fscvRate;
        nextLFPIndex=currentLFPIndex+ceil(durationTrial*ratelfp);
    end
    if currentFSCVIndex<1
        currentFSCVIndex=1;
    end
    %Extract fscv data for appropriate trial indeces from merged txt data
    if ii<length(fscvTargetTrialsIDs) && nextLFPIndex<length(TS)
        %if not last target ID, get to next trial end
        if nextFSCVIndex>size(txtDataCh0,1)
            nextFSCVIndex=size(txtDataCh0,1);
        end
        if currentFSCVIndex>size(txtDataCh0,1)
            display(['start sync trigger ID = ' num2str(currentFSCVIndex)]);
            display(['length of samples = ' num2str(size(txtDataCh0,1))]);
            error(['start sync trigger > length of FSCV recording' +sprintf('\n') ...
                'check if FSCV recording experienced break or data drop'])
        end
        syncData(ii).ch0=txtDataCh0(currentFSCVIndex:nextFSCVIndex,:);
        if numch>1
            syncData(ii).ch1=txtDataCh1(currentFSCVIndex:nextFSCVIndex,:);
        end
        if numch>2
            syncData(ii).ch2=txtDataCh2(currentFSCVIndex:nextFSCVIndex,:);
        end
        if numch>3
            syncData(ii).ch3=txtDataCh3(currentFSCVIndex:nextFSCVIndex,:);
        end
        syncData(ii).events=fscvNlxEvents(currentFSCVIndex:nextFSCVIndex,:);
        for LFPch=1:numLFPch
            syncLFP{LFPch}=samples{LFPch}(currentLFPIndex:nextLFPIndex);
        end
        tsLFP=TS(currentLFPIndex:nextLFPIndex);
        nlxTargetIDsWithinTrial=find(nlx_ts>=tsLFP(1) & nlx_ts<=tsLFP(end));
        NlxEventTTL=dg_Nlx2Mat_TTL(nlxTargetIDsWithinTrial);
        NlxEventTS=nlx_ts(nlxTargetIDsWithinTrial);
        %check how much apparent deviation in time from fscv trial & nlx
        %trial start ID, if over 500 ms, then too much
        if tsLFP(1)>syncData(ii).events(1,4)+0.5 || tsLFP(1)<syncData(ii).events(1,4)-0.5
            warning(['trial ' num2str(ii) ' has > 0.5 s latency'])
            display(['latency first idx in saved fscv trial data and ' ...
                'nlx is: ' num2str(tsLFP(1)-syncData(ii).events(1,4))]);
            writelog=['trial ' num2str(ii) ' (file ' num2str(initialNum) ...
                ') sample 1 latency: ' ...
                    num2str(tsLFP(1)-syncData(ii).events(1,4)) ' s'];
            %record in log file
            if isempty(log)
                log=writelog;
            else
                log=[log sprintf('\n') writelog];
            end
                
        end
    else
        %get to end of recording
        if fscvTargetTrialsIDs(ii)<=length(fscvNlxEvents)
            nextfidx=length(fscvNlxEvents);
            if nextfidx>size(txtDataCh0,1)
                nextfidx=size(txtDataCh0,1);
            end
            syncData(ii).ch0=txtDataCh0(currentFSCVIndex:nextfidx,:);
            if numch>1
            syncData(ii).ch1=txtDataCh1(currentFSCVIndex:nextfidx,:);
            syncData(ii).ch2=txtDataCh2(currentFSCVIndex:nextfidx,:);
            syncData(ii).ch3=txtDataCh3(currentFSCVIndex:nextfidx,:);
            end
            syncData(ii).events=fscvNlxEvents(currentFSCVIndex:length(fscvNlxEvents),:);
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
    end
    if fscvTargetTrialsIDs(ii)<=length(fscvNlxEvents)
        saveName=['fscv_multi_' num2str(initialNum)];
        fscv.ch0=syncData(ii).ch0;
        if numch>1
        fscv.ch1=syncData(ii).ch1;
        end
        if numch>2
        fscv.ch2=syncData(ii).ch2;
        end
        if numch>3
        fscv.ch3=syncData(ii).ch3;
        end
        fscv.events=syncData(ii).events;
        samplesLFP=syncLFP;
        if splitfiles==0
            %merge all signals into one file
            save([PathName4, saveName],'fscv','samplesLFP','tsLFP','NlxEventTTL','NlxEventTS','nlxFileNames','-v7.3');
            prompt = ['saved as: ' saveName];
        end
        if splitfiles==1
            %save nlx signals as separate files and fscv as separate file
            %save in nlx format compatible with dg / lfp lib
            dg_Nlx2Mat_Timestamps=tsLFP.*1e6;     %convert to Microseconds
            dg_Nlx2Mat_Samples={};
            dg_Nlx2Mat_SamplesUnits='V';
            for nlxch=1:numLFPch
                %save each nlx channel as separate file
                dg_Nlx2Mat_Samples=samplesLFP{nlxch};   %already in volts
                savecsc=['csc' num2str(selectcsc(nlxch)) '_' num2str(initialNum)];
                save([PathName4 savecsc],'dg_Nlx2Mat_Samples','dg_Nlx2Mat_SamplesUnits','dg_Nlx2Mat_Timestamps','-v7.3');
            end
            %save fscv file without nlx signals, but with events
            save([PathName4, saveName],'fscv','tsLFP','NlxEventTTL','NlxEventTS','nlxFileNames','-v7.3');
            prompt = ['saved as, split ' num2str(numLFPch) ' chs: ' saveName];
        end        
        initialNum=initialNum+1;
    end
   % fprintf('%d \r', floor(ii/length(fscvTargetTrialsIDs)*100)); 
    previousID=fscvTargetTrialsIDs(ii);
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