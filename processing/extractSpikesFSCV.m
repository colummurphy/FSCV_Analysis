function extractSpikesFSCV(subject,sessnum,varargin)
%extract spikes using dg_fitFSCVart dg_rmFSCVart
%made for entire csc file recording, also output interpolated data as
%temporary storage unless checked otherwise varargin
%use chronicXXchconfig.m to include
%ncsnoartifacts={'cl1', 's4','s3'};  %to denote what channels don't have
%artifacts

log=[];
flagNLX=0;
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
homedir='Z:';                           %NAS MAPPED NETWORK DRIVE CHANGE HERE
if ~ispc
    homedir=[filesep];
end
graserver='inj-monkey2';        %change here and fscvdir after move to Pitt 05/2021
graserver='data_MIT';
fscvdir='patra_fscv2';          %now used in "chronicXXchconfig" to load paths          
fscvdir='patra_fscv';    
patra_ephys='patra_ephys';       %add for use in "chronicXXchconfig" to load paths         

analysispath=fullfile(homedir,graserver,'analysis',subjectname,['chronic' sessid],filesep);
if ~isdir(analysispath)
    mkdir(analysispath);
end
%default on putamen pc in lab, get data directory from config file
configdir='C:\Users\putamen\Documents\MATLAB\fscv\analysis\config\';
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
intervalsplit=[];       %default full file for fixed interval default
splitfiles=1;           %default split channels into separate files as original files
mergeflag=1;            %default merge flag (i.e. not task separated)
aligns={};      %multiple alignment events
align='fixedint';       %default, no longer in argument of function
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

argnum=1;
cleolickflag=0;
cscadd={};
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'nomerge'
            %don't merge entire recording into 1 file
            mergeflag=0;
        case 'cleolick37'
            cleolickflag=1;     %for chronic 25 where lick rec
        case 'interval'
            argnum=argnum+1;
            intervalsplit=varargin{argnum};
        case 'nlx'
            %convert to nlx
            flagNLX=1;
           % argnum=argnum+1;
           % cscadd=varargin{argnum};
        case 'sites'
            %explicitly indicate sites to convert, not from
            %chronicXXchconfig file
            %should be in format {'sitename1','sitename2',...}
            argnum=argnum+1;
            ncschannels=varargin{argnum};
    end
    argnum=argnum+1;
end
if isempty(intervalsplit)
    %if empty argument, then means entire file
    mergeflag=1;
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
noartcsc=[];
if exist('ncsnoartifacts')
    %channels listing no artifacts from chronicXXchconfig.m
    noartcsc=getcscids(ncsnoartifacts,csc_map);
end

%selectcsc=[33:35];          %nlx channels that want to store, only eye
preOffset=30;           %seconds from event ID to align to
durationTrial=60;       %seconds duration of file past start trial index, if zero, then up to next trial event

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
[samples, nlxFileNames, TS,header,nlxjunk]=ephys_getCSCwheader(selectcsc,pathnlx);
%TS here is vector nlx timestamps in seconds
ratelfp=getSampleRate(TS);

ratelfp=30e3;           %DEFAULT

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
nlxTargetIDs=[];
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
else
     preOffset=0;
    durationTrial=0;
    if mergeflag==1
        %whole file
       % fscvTargetTrialsIDs=[1 size(fscvNlxEvents,1)];
        lfp_tID(1)=1;
        nlxTargetIDs=1:length(dg_Nlx2Mat_TTL);
    else
      warning('not merging all trials, not single whole file interval, program not programmed for this option');
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
if mergeflag==0
    alignname='tasksplit';
end
subfolder=alignname;
if mergeflag==1
    subfolder='';
end
savepath=fullfile(pathfscv,'matlab','spikes',subfolder,filesep);
if ~isdir(savepath)
    status = mkdir(savepath);
end
temppath=fullfile(pathfscv,'matlab','spikes','interp',subfolder,filesep);
nospikepath=fullfile(pathfscv,'matlab','spikes',subfolder,filesep);
%storage for interpolated raw data
if ~isdir(temppath)
    status = mkdir(temppath);
end
initialNum=100;
previousID=0;
log=[];
totaltrials=length(nlxTargetTSs);
totaltrials=length(lfp_tID);

if durationTrial==0
    %fixed interval condition, just 1 less?????
   % totaltrials=totaltrials-1;
end

%get header information for writing CSC
ADBitVoltstr={};
ADBitVolts=[];
 for k = 1:length(header)
    if regexp(header{k}, '^\s*-ADBitVolts\s+')
        ADBitVoltstr = regexprep(header{k}, ...
            '^\s*-ADBitVolts\s+', '');
        ADBitVolts = str2num(ADBitVoltstr);                           
    end
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
           nextLFPIndex=length(TS);  %to end of recording          
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
            %nexttargetid=nlxTargetIDs(ii)+100;
            %if nexttargetid>length(dg_Nlx2Mat_TTL)
                nexttargetid=length(dg_Nlx2Mat_TTL);
           % end
            NlxEventTTL=dg_Nlx2Mat_TTL(nlxTargetIDs(ii):nexttargetid);
            NlxEventTS=nlx_ts(nlxTargetIDs(ii):nexttargetid);
    end
    if lfp_tID(ii)<=length(TS)
        saveName=['fscv_multi_' num2str(initialNum)];        
        samplesLFP=syncLFP;
        
        %Spike Extraction by interpolating FSCV artifacts 08/2021
        %get threshold by averaging all samples all chs with artifacts
        %{
        samplesall=cell2mat(samplesLFP); %out of memory error
        if size(samplesLFP,2)> size(samplesLFP,1)
            samplesall=cell2mat(samplesLFP');
        end
        %}
        %ids of only artifact containing and neural signals
        artchs=find(~ismember(selectcsc,otherchs) & ~ismember(selectcsc,noartcsc));
        %samplesmean=nanmean(samplesall(artchs,:),1);
        sumsamples=samplesLFP{artchs(1)};
        for iachs=2:length(artchs)
            sumsamples=sumsamples+samplesLFP{artchs(iachs)};
        end
        samplesmean=sumsamples./length(artchs);
        
        %some samples are nan if frames have excessively variable duration
        %as found in ephys_getCSCwheader
        thres=nanmean(rms2(samplesmean,'omitnan'))*2;
       % xingidx1 = dg_fitFSCVart(tsLFP, samplesmean, thres); 
        xingidx = dg_fitFSCVart2(tsLFP, samplesmean, thres); 
        disp(['found ' num2str(length(xingidx)) ' artifacts']);   
        for iachs=1:length(artchs)
            %interpolate artifacts away with fixed 11 ms window for all chs, 0.7
            %ms backward and 10.3 ms forward
            samplesinterp=dg_rmFSCVart(samplesLFP{artchs(iachs)},xingidx,1e-3*30e3,10e-3*30e3); 
            samplesLFP{artchs(iachs)}=samplesinterp;
        end
        
        %find channels that have spikes physiological 09/10/2021 using
        %interpolated data, or non-interpolated (if didn't have fscv
        %artifacts), scanning all neural channels now
        sampleperiod=1/ratelfp;
        filtsamples={};
        for iachs=1:length(lfpchs)
            filtsamples{iachs}=filterLFP(samplesLFP{iachs},ratelfp,[100 inf]); %filter for spike detection below
        end
        samples_spikes={};      %normal spikes detected
        samples_bispikes={};    %biphasic spikes detected
        fileID = fopen([nospikepath 'spikes.txt'],'w');       %log # spikes per channel detected in separate file
        fprintf(fileID,'%s\t','csc');
        fprintf(fileID,'%s\t','# spikes');
        fprintf(fileID,'%s\n','# biphasic spikes');
        for iachs=1:length(lfpchs)
            %filtsamples=filterLFP(samplesLFP{iachs},ratelfp,[100 inf]); %filter for spike detection below
            thress=4*std(filtsamples{iachs});
            [sS, sS2] = lfp_findSpikesFSCV(iachs, thress, filtsamples,sampleperiod,'invert');   %invert so spikes are positive going
           tsS=TS(sS);    %convert to timestmsaamps for detected spikes
           tsS2=TS(sS2);    %convert to timestmsaamps for detected spikes
           samples_spikes{iachs}=sS;
           samples_bispikes{iachs}=sS2;
           sname=['csc' num2str(lfpchs(iachs))];
           fprintf(fileID,'%s\t',sname);
            fprintf(fileID,'%d\t',length(sS));
            fprintf(fileID,'%d\n',length(sS2));
        end       
        fclose(fileID);     %close text file writing

        if splitfiles==0
            %merge all signals into one file
            save([temppath, saveName],'samplesLFP','tsLFP','NlxEventTTL','NlxEventTS','nlxFileNames','-v7.3');
            prompt = ['saved as: ' saveName];
        end
        if splitfiles==1
            %save nlx signals as separate files and fscv as separate file
            %save in nlx format compatible with dg / lfp lib
            dg_Nlx2Mat_Timestamps=round(tsLFP.*1e6);     %convert to Microseconds * MAKE SURE TO ROUND OTHERWISE GET ROUNDOFF ERROR
            dg_Nlx2Mat_Samples={};
            dg_Nlx2Mat_SamplesUnits='microVolts'; %convert to microvolts for Plexon offline sorter
            for nlxch=1:numLFPch
                %save each nlx channel as separate file
                dg_Nlx2Mat_Samples=samplesLFP{nlxch}.*1e6;   %convert to microvolts for Plexon offline sorter
                savecsc=['csc' num2str(selectcsc(nlxch)) '_' num2str(initialNum)];
                if ~flagNLX
                    save([temppath savecsc],'dg_Nlx2Mat_Samples','dg_Nlx2Mat_SamplesUnits','dg_Nlx2Mat_Timestamps','-v7.3');
                else
                    %save spike detection data, whether spikes detected or
                    %not
                    spikes=samples_spikes{nlxch};
                    spikesbi=samples_bispikes{nlxch};
                    save([nospikepath 'spikeids_' savecsc],'spikes','spikesbi');
                    %convert to NLX using dg_writeCSC
                    %first convert to frames of 512 samples per frame,
                    %timestamp for last sample
                    %save in nlx format compatible with dg / lfp lib
                    dg_Nlx2Mat_Timestamps=round(tsLFP.*1e6);     %in microseconds for nlx
                    %convert to ad units for nlx
                    samplesLFP2={};                    
                    samplesLFP2{nlxch} = samplesLFP{nlxch}./ADBitVolts; %in ad units for nlx, see header
                    dg_Nlx2Mat_Samples=round(samplesLFP2{nlxch});   
                    dg_Nlx2Mat_SamplesUnits='AD'; 
                    nlx_TS=dg_Nlx2Mat_Timestamps(1:512:end);    %get TS every 512 timestamps to put into frames
                    samplesnlx=reshape(dg_Nlx2Mat_Samples,512,length(dg_Nlx2Mat_Samples)/512);   %reshape samples into 512 length frames
                    dg_writeCSC([temppath savecsc '.ncs'], nlx_TS, samplesnlx, header,'nlxjunk',nlxjunk);   %need to write junk that came with origianl nlx file in order to read into plexon
                end                
            end
            prompt = ['saved as, split ' num2str(numLFPch) ' chs: ' saveName];
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
save([savepath, savelogname],'log');
end
%save log of latencies
end