function compiletrialsdir(pathdir,varargin)
%1/2019, fix to output trialbytrial(1).ttarg, tfix, tfixeye, etc with exact
%timestamps relative to ts(1) for nlx signals, rather than just samples for
%fscv
%06/28/2018 make function call for running in chunky
%modified 04/08/2018, 
%saves trialbytrial(ich) data, with samplesfix, samplestarg times,
%& selected after fail trials, after switch side trials
%saves max da conc calculated from
%maximum value detected until next trial 
%maximum change (ie. negative/positive peak) immediately following target
%until align, 
%& maximum change until 2 s post align time (ie. reward period)
%04/12/2018 open directory
%07/27/2018 change prefix baseline to 1 s before 
%average before actual fix ttl (4)
%instead of left/right condition
argnum=1;
configmode='patrabipolar';     %for pc, default patra, just for event-codes
%configmode='patrabichunky';     %for chunky, needs to specify dir for pcr temp
baseline='fix';     %subtract to fixation period
baseline='prefix';  %subtract to prior to fix cue - 1 second
%
nofscv=0;
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'fscvchs'
            %user provided fscv selected channels / NOT USED detected
            %instead
            argnum=argnum+1;
            selch=varargin{argnum};
        case 'config'
            %user provided config file
            argnum=argnum+1;
            configmode=varargin{argnum};        
        case 'nofscv'
            nofscv=1;
    end
    argnum=argnum+1;
end

plotparam.baseline=baseline;
avgperiod=1;        %period for averaging baseline if prefix
eventIDsFix=[29 30];        %
eventIDsFixCue=4;           %fix apperance
eventIDsFixOn=5;            %start fixation on fix cue
eventIDsTarget=10;          %target apperance, 11 is when start fix on targ
eventIDsTargetOn=11;            %eye on target, start target fixation 1/10/19 add

%[dParam,param,csc_map,event_codes]=getparams('patra');
%[dParam,param,csc_map,event_codes]=getparams('patrabipolar');
[parameters,csc_map,event_codes]=getparams(configmode);

plotparam.cminshow=0;

time=[20 40];

Vox_ID=125;         %when plotting raw Iox
interval=2.5;
stdthres=3;
%cscale2=[-0.2 0.4];

d=dir(pathdir);
savepath=pathdir;
namn=strfind(savepath,filesep);
namy=savepath(namn(end-2)+1:namn(end-1)-1);

pathdirs=ismember({d.name},['ch1'; 'ch2'; 'ch3'; 'ch4']);   %get dirs with ch folders
pathdirid=find(pathdirs==1);
folder=fullfile(pathdir, d(pathdirid(1)).name, filesep); %get first folder
chnum=str2num(d(pathdirid(1)).name(3:end)); %get chnum from folder name
numchs=length(pathdirid);       %number of ch folders to analyze
if nofscv
    folder=pathdir;
    numchs=1;
end
ds=dir(folder);
ds2=find([ds.isdir]==0);            %get actual file id's
nFiles=size(ds2,2);
dsfilenames={ds(ds2).name};
PathName=folder;
FileName=dsfilenames{1};
dataLoaded=load([PathName FileName]);
samplesperscan=dataLoaded.parameters.samplesperscan;
samplespersec=dataLoaded.parameters.samplerate;
chnum=[];
Vrange=[];
Vrange_cathodal=[];
Vrange_all=[];
if ~nofscv
chnum=dataLoaded.chnum;
Vrange=dataLoaded.Vrange;
Vrange_cathodal=dataLoaded.Vrange_cathodal;
Vrange_all=[Vrange Vrange_cathodal];
end
trialSortTS=zeros(1,nFiles);
samples=zeros(samplesperscan,nFiles);       %best cv's each file
win=time(1)*samplespersec:1/samplespersec:time(2)*samplespersec;
plotparam.trialtype=namy;
plotparam.samplespersec=samplespersec;
plotparam.samplesperscan=samplesperscan;
plotparam.win=win;
plotparam.interval=2.5;
alignTime=median(win(1):win(end));
%initialize structure arrays
trialbytrial={};
switchtrials={};
afterfailtrials={};
maxDAs={};
samplesfix={};
samplesfixeye={};
samplestarg={};
samplestargeye={};  %added 1/10/19
%save also ts's for nlx signals 
tfix={};
tfixeye={};
ttarg={};
ttargeye={};
samplesfixnexttemp=[];
fscvchs=[];
logfile=[];         %log errors/warnings
for ich=1:numchs
    %scroll through channels
    if ~nofscv
        chnum=str2num(d(pathdirid(ich)).name(3:end)); %get chnum from folder name
        folder=fullfile(pathdir,d(pathdirid(ich)).name,filesep); %get first folder
        fscvchs(ich)=chnum;
    end
    ds=dir(folder);
    ds2=find([ds.isdir]==0);            %get actual file id's
    dsfilenames={ds(ds2).name};
    PathName=folder;
    disp(['processing fscv ch' num2str(chnum) ]);
    status=0;
    sumCount=1;         %for averaging, flag
    sumcolors=[];
for ifile=1:nFiles
    %scroll files for given ch
    dataLoaded=load([PathName dsfilenames{ifile}]);

   % dataLoaded=load(files{ifile});
    %only analyze data within time window exclusive to specific trial
    %load fscv data
    Itemp=[];
    IM=[];
    IBG=[];
    IPH=[];
    Isub=[];
    Iox=[];
    if ~nofscv
        Itemp=dataLoaded.Ipcr.DAiso;
      %  Itemp(310)=dataLoaded.Ipcr.DAproj(310);      %***glitch in saving, all 310 samples are nan for no reason
        IM=dataLoaded.Ipcr.Mplot;
        IBG=dataLoaded.Ipcr.BGplot;
        IPH=dataLoaded.Ipcr.pHplot;
        if ifile==nFiles && length(Itemp)~=length(trialbytrial(ich).da(ifile-1,:))
            warning(['last file # ' num2str(ifile) ' skipped'])
            %last file, time not total duration, skip
            continue
        end
        Isub=[];
        if isfield(dataLoaded.Isub,'data')
            Isub=dataLoaded.Isub.data;
        else
            Isub=dataLoaded.Isub;
        end
        Iox=Isub(Vox_ID,:);
        samples1=dataLoaded.Idata.anodal; samples2=dataLoaded.Idata.cathodal;   
        %load cv's, if any correlated cv's found
        if ~isempty(samples1)
            if size(samples1,1)>size(samples1,2)
                %default
            samples(:,ifile)=flipud([samples1; samples2]');        %store best CV's in each trial
            else
                samples(:,ifile)=flipud([samples1'; samples2']');
            end
        else
            samples(:,ifile)=nan(1,samplesperscan);    %if no cv exists in file, save as nan's
        end
    end
    %get relevant event ts's
    ts1=dataLoaded.ncsread.LFPts(1);        %first time point in recording
    events=dataLoaded.ncsread.LFPeventTTL;      %get event TTLs from NCS record
    idxfixTTLs=find(ismember(events,eventIDsFixCue)==1);   
    %get targeted TTL idx %07/272018 changed eventIDsFix to eventIDsFixCue for actual fix cue
    if ich==1
        %find evets closest to alignment index alignTime
        tssfix=(dataLoaded.ncsread.LFPeventTS(idxfixTTLs)-ts1);    %in seconds
        tssfixtoaln=tssfix-alignTime./samplespersec;
        [~,closest]=min(abs(tssfixtoaln));
        samplesfixnexttemp(ifile)=win(end);        %default end window
        if ~isempty(closest)
            trialbytrial(1).tfix(ifile)=tssfix(closest);    %get ts for closest to alignTime
            tsfix=round(trialbytrial(1).tfix(ifile).*10)./10;   %rounded to 100 ms
            [samplesfixtemp, idxfixtemp, bb]=intersect(tsfix*samplespersec,win);         %only within window
            if ~isempty(samplesfixtemp)
                samplesfix{ich}(ifile)=samplesfixtemp(1);
                if closest<length(tssfix)
                    samplesfixnexttemp(ifile)=round(tssfix(closest+1)*10)./10.*samplespersec;
                end
            else
                samplesfix{ich}(ifile)=nan;
                trialbytrial(1).tfix(ifile)=nan;
            end
        else
           samplesfix{ich}(ifile)=nan;       %default if not found
           trialbytrial(1).tfix(ifile)=nan;
            %write log file if empty
            warning(['fix cue idx not recorded, check file ' num2str(ifile)])
                writelog=['trial ' num2str(ifile) ' / no appropriate tsfix '];
                %record in log file
                if isempty(logfile)
                    logfile=writelog;
                else
                    logfile=[logfile char(10) writelog];
                end                
        end        
        
        idxeyefixTTLs=find(ismember(events,eventIDsFixOn)==1);  %eye on fix cue
        tssfixeye=(dataLoaded.ncsread.LFPeventTS(idxeyefixTTLs)-ts1);    %in seconds
        tssfixeyetoaln=tssfixeye-alignTime./samplespersec;
        [~,closest]=min(abs(tssfixeyetoaln));
        if ~isempty(closest)
            trialbytrial(1).tfixeye(ifile)=tssfixeye(closest);    %get ts for closest to alignTime
            tsfixeye=round(trialbytrial(1).tfixeye(ifile).*10)./10;   %rounded to 100 ms
            [samplesfixeyetemp, idxfixeyetemp, bb]=intersect(tsfixeye*samplespersec,win);         %only within window
            if ~isempty(samplesfixeyetemp)                
                samplesfixeye{ich}(ifile)=samplesfixeyetemp(1);
            else
                samplesfixeye{ich}(ifile)=nan;
                trialbytrial(1).tfixeye(ifile)=nan;
            end
        else           
           samplesfixeye{ich}(ifile)=nan;         %default if not found
           trialbytrial(1).tfixeye(ifile)=nan; 
            warning(['fix cue eye idx not recorded, check file ' num2str(ifile)])
            writelog=['trial ' num2str(ifile) ' / no appropriate tsfixeye '];
            if isempty(logfile)
                logfile=writelog;
            else
                logfile=[logfile char(10) writelog];
            end
        end
       
    idxtargTTLs=find(ismember(events,eventIDsTarget)==1);
    tsstarg=(dataLoaded.ncsread.LFPeventTS(idxtargTTLs)-ts1);   %in seconds
    tsstargtoaln=tsstarg-alignTime./samplespersec;
    [~,closest]=min(abs(tsstargtoaln));
    if ~isempty(closest)
        trialbytrial(1).ttarg(ifile)=tsstarg(closest);
        tstarg=round(trialbytrial(1).ttarg(ifile).*10)./10;   %rounded to 100 ms
        samplestargtemp=intersect(tstarg*samplespersec,win);       %only within window
        if ~isempty(samplestargtemp)
            samplestarg{ich}(ifile)=samplestargtemp(1);
        else
            samplestarg{ich}(ifile)=nan;
            trialbytrial(1).ttarg(ifile)=nan;
        end
    else
        samplestarg{ich}(ifile)=nan;
        trialbytrial(1).ttarg(ifile)=nan;
        warning(['target cue idx not recorded, check file ' num2str(ifile)])
        writelog=['trial ' num2str(ifile) ' / no appropriate tstarg ' ];
        %record in log file
        if isempty(logfile)
            logfile=writelog;
        else
            logfile=[logfile char(10) writelog];
        end
    end

    idxtargonTTLs=find(ismember(events,eventIDsTargetOn)==1);
    tsstargeye=(dataLoaded.ncsread.LFPeventTS(idxtargonTTLs)-ts1);   %in seconds
    tsstargeyetoaln=tsstargeye-alignTime./samplespersec;
    [~,closest]=min(abs(tsstargeyetoaln));
    if ~isempty(closest)
        trialbytrial(1).ttargeye(ifile)=tsstargeye(closest);
        tstargeye=round(trialbytrial(1).ttargeye(ifile).*10)./10;   %rounded to 100 ms
        samplestargontemp=intersect(tstargeye*samplespersec,win);       %only within window
        if ~isempty(samplestargontemp)            
            samplestargeye{ich}(ifile)=samplestargontemp(1);
        else
            samplestargeye{ich}(ifile)=nan;  
            trialbytrial(1).ttargeye(ifile)=nan;
        end
    else
        samplestargeye{ich}(ifile)=nan;  
        trialbytrial(1).ttargeye(ifile)=nan;
        warning(['target cue on idx not recorded, check file ' num2str(ifile)])
        writelog=['trial ' num2str(ifile) ' / no appropriate tstargon ' ];
        if isempty(logfile)
            logfile=writelog;
        else
            logfile=[logfile char(10) writelog];
        end
    end
    end
   
    %load lfp data if exists
    if isfield(dataLoaded.nlx,'resampled')
        %DON"T USE ANYMORE
        if ifile==1 && ich==1
            %only load these for first file & first ch since same across files
            trialbytrial(ich).cscNames=dataLoaded.nlx.cscNames;
            %trialbytrial(ich).lfpch=dataLoaded.nlx.lfpCh;    %csc original ch names
           % trialbytrial(ich).eyech=dataLoaded.nlx.eyeCh;
           % trialbytrial(ich).lickch=dataLoaded.nlx.lickCh;
            %trialbytrial(ich).pulsech=dataLoaded.nlx.pulseCh;
            trialbytrial(ich).eyeid=dataLoaded.nlx.eyeid;    %id's for csc names
            trialbytrial(ich).lfpid=dataLoaded.nlx.lfpid;
            trialbytrial(ich).lickid=dataLoaded.nlx.lickid;
            trialbytrial(ich).pulseid=dataLoaded.nlx.pulseid;

        end  
        if ich==1
            %only load for first ch since LFP same across channels
            %slightly varies from trial to trial but pretty much same
            %for sake of aligning below update for eacht rial
            

            %load & store csc data
            cscdata=dataLoaded.nlx.resampled;

            %reorganize csc signals by id --> trials
            %each {ich} cell contains all data recorded for channel for all
            %trials concatenated in the cell
            for icsc=1:size(cscdata,1)
               trialbytrial(ich).lfp{icsc}(ifile,:)=cscdata(icsc,:);               
            end
        end
    end
    if ich==1
       trialbytrial(ich).relts=dataLoaded.ncsread.LFPts-dataLoaded.ncsread.LFPts(1);
      trialbytrial(ich).tscsc=dataLoaded.ncsread.LFPts;
            trialbytrial(ich).ratelfp=dataLoaded.ncsread.LFPsamplingfreq;%time stamp relative to 0 for each trial
        %get absolute time stamps for each trial
        aidx=round(alignTime*trialbytrial(ich).ratelfp/samplespersec);
        if ~isfield(dataLoaded.ncsread,'nolfp')
            trialbytrial(ich).alignts(ifile)=dataLoaded.ncsread.LFPts(aidx);    
        else
            aidx=round(alignTime);
            %if just dopamine recording no lfps
            trialbytrial(ich).alignts(ifile)=dataLoaded.ncsread.LFPts(aidx);      
        end
    end
    
    %continue finding relevent event times
    %relts=trialbytrial(1).relts;
    infob=computevents(dataLoaded.ncsread,[],event_codes);
    preveventidx=round(infob.preveventts*10)/10.*samplespersec; %fscv idx units
    trialbytrial(ich).prevoutcomeidx(ifile)=preveventidx; %fscv idx units
    nexteventidx=round(infob.nexttrialts*10)/10.*samplespersec; %fscv idx units
    if infob.currside~=-1
        %if was able to retrieve events properly for trials around current
        trialbytrial(ich).prevtargcueidx(ifile)=round(infob.prevtargcuets*10)/10.*samplespersec; %fscv idx units
        trialbytrial(ich).prevfixcueidx(ifile)=round(infob.prevfixcuets*10)/10.*samplespersec; %fscv idx units
        trialbytrial(ich).prevoutcome{ifile}=infob.prevevent;      %store string of previous trial outcome
        trialbytrial(ich).nextfixcueidx(ifile)=round(infob.nextfixcuets*10)/10.*samplespersec; %fscv idx units
    else
        trialbytrial(ich).prevtargcueidx(ifile)=nan; %fscv idx units
        trialbytrial(ich).prevfixcueidx(ifile)=nan; %fscv idx units
        trialbytrial(ich).prevoutcome{ifile}=nan;      %store string of previous trial outcome
        trialbytrial(ich).nextfixcueidx(ifile)=nan; %fscv idx units
    end

    if ~nofscv
    %use trial min as baseline
    baselineDA=nanmin(Itemp);
    idbaseline=[samplesfix{1}(ifile) samplestarg{1}(ifile)];
    idbaseline=idbaseline(idbaseline>1 & idbaseline<=length(Itemp));
    if length(idbaseline)~=2
        idbaseline=[1 length(Itemp)];
    end
    switch baseline
        case 'prefix'
            idbaseline=[samplesfix{1}(ifile)-1-avgperiod.*samplespersec samplesfix{1}(ifile)-1];
            idbaseline=idbaseline(idbaseline>1 & idbaseline<=length(Itemp));
            if length(idbaseline)~=2
                idbaseline=[1 length(Itemp)];
            end
    end
    %for calculating average changes large/small reward
    baselineDA=nanmean(Itemp(idbaseline(1):idbaseline(2)));     
    baselineM=nanmean(IM(idbaseline(1):idbaseline(2)));     
    baselinePH=nanmean(IPH(idbaseline(1):idbaseline(2)));     
    baselineBG=nanmean(IBG(idbaseline(1):idbaseline(2)));     
    baselineOx=nanmean(Iox(idbaseline(1):idbaseline(2)));     
    
    %subtract baseline from signals
    Itemp=Itemp-baselineDA;
    IM=IM-baselineM;
    IBG=IBG-baselineBG;
    IPH=IPH-baselinePH;
    Iox=Iox-baselineOx;
    
    %da change metrics
    %max amplitude (abs) after target until next trial cue
    %{
    [maxwithintrial, maxid]=nanmax(Itemp(samplestarg{1}(ifile):samplesfixnexttemp(ifile)));
    maxDAs{ich}(ifile).positiveposttarget=Itemp(samplestarg{1}(ifile)+maxid-1);    
    %max change (ie decrease/increase) from target to align (rew) time
    [maxwithintrial, maxid]=nanmax(abs(Itemp(samplestarg{1}(ifile):alignTime)));
    maxDAs{ich}(ifile).posttargetprealigndelta=Itemp(samplestarg{1}(ifile)+maxid-1);  
    %max change (ie decrease/increase) from target to 2s after reward/align time
    [maxwithintrial, maxid]=nanmax(abs(Itemp(samplestarg{1}(ifile):alignTime+20)));
    maxDAs{ich}(ifile).posttarget=Itemp(samplestarg{1}(ifile)+maxid-1);  
    %average immediate amplitude (ie decrease/increase) from target to reward/align time 
    meantarget=nanmean(Itemp(samplestarg{1}(ifile):alignTime));
    maxDAs{ich}(ifile).posttargetmean=meantarget; 
    %max change from fix appearance to prior last trial
    if isnan(preveventidx)
        preveventidx=1;
    end
    [maxwithintrial, maxid]=nanmax(abs(Itemp(preveventidx:samplesfix{1}(ifile))));
    if ~any(isnan(maxid)) && ~isempty(maxid)
        maxDAs{ich}(ifile).pretrial=Itemp(preveventidx+maxid(1)-1); 
    else
        maxDAs{ich}(ifile).pretrial=nan; 
    end
    
    %max change from curr reward/align to next trial
    [maxwithintrial, maxid]=nanmax(abs(Itemp(alignTime:samplesfixnexttemp(ifile))));
    maxDAs{ich}(ifile).posttrial=Itemp(alignTime+maxid-1);   
  %}
    %if good trial, average color plot
    BGvectors=Isub(:,idbaseline(1):idbaseline(2));   
    bgdata=mean(BGvectors,2);
    %tempdata=Isub-bgdata;
    tempdata=Isub;
    if ifile>1
        %not first file so have stored sumcolors
    if size(tempdata,2)~=size(sumcolors,2)
        %not equal lengths
        tempdata=tempdata(:,1:size(tempdata,2));
    end
    end
   % if trialbytrial(ich).badtrials(ifile)==0
        if sumCount==1
            sumcolors=tempdata;        %first data stored
            sumCount=sumCount+1;        %increment count        
        elseif ifile<nFiles
            if size(tempdata,2)>size(sumcolors,2)
                %ignore first data set if different length
                sumcolors=tempdata;
                sumCount=1;     %reset counter
            else
                sumcolors=tempdata+sumcolors;      %add to existing data
            end
            sumCount=sumCount+1;        %increment count        
        end
    %end
    %average all if last file
    if ifile==nFiles
        %if trialbytrial(ich).badtrials(ifile)==0
            %sum with existing data if good trial
            if isequal(size(tempdata,2),size(sumcolors,2))
            sumcolors=tempdata+sumcolors;
            sumCount=sumCount+1;
            trialbytrial(ich).avgcolor=sumcolors./sumCount;
            end
       % else
            %if bad trial sum use existing data wihtout adding current
         %   trialbytrial(ich).avgcolor=sumcolors./sumCount;
       % end
    end
    if ifile>1
        %not first file so have stored data
        if size(Itemp,2)~=size(trialbytrial(ich).da,2)
            %not equal lengths
            if size(Itemp,2)>=size(trialbytrial(ich).da,2)
            Itemp=Itemp(1:size(trialbytrial(ich).da,2));
            IM=IM(1:size(trialbytrial(ich).da,2));
            IPH=IPH(1:size(trialbytrial(ich).da,2));
            IBG=IBG(1:size(trialbytrial(ich).da,2));
            Iox=Iox(1:size(trialbytrial(ich).da,2));
            else
            Itemp=repmat(nan,1,size(trialbytrial(ich).da,2));
            IM=repmat(nan,1,size(trialbytrial(ich).da,2));
            IPH=repmat(nan,1,size(trialbytrial(ich).da,2));
            IBG=repmat(nan,1,size(trialbytrial(ich).da,2));
            Iox=repmat(nan,1,size(trialbytrial(ich).da,2));
            end
        end
    end
    trialbytrial(ich).da(ifile,:)=Itemp;
    trialbytrial(ich).m(ifile,:)=IM;
    trialbytrial(ich).ph(ifile,:)=IPH;
    trialbytrial(ich).bg(ifile,:)=IBG;
    trialbytrial(ich).iox(ifile,:)=Iox;
    
    %sorting trials by threshold crossing
    passesThresholdTime=find(Itemp>stdthres*std(Itemp(~isnan(Itemp)))); 
    % passesThresholdTime=find(Itemp>stdthres*std(ITdata(~isnan(Itemp)))); %all data;
    if ~isempty(passesThresholdTime)

        passesThresholdTime=passesThresholdTime(1);
    else
       passesThresholdTime=max(win(1):win(end));
    end
    trialSortTS(ifile)=passesThresholdTime;
    
    %da processing finished if fscv exists
    end
    if isfield(dataLoaded,'info')
   if ~isfield(dataLoaded.info,'behav')
       %behavior not calculated in previously
       dataLoaded.info.behav=infob;              
   end
    else
       dataLoaded.info.behav=infob;              
    end
    %compute behavioral metrics
    if ~isempty(dataLoaded.info.behav.target_rt)
        trialbytrial(ich).target_rt(ifile)=dataLoaded.info.behav.target_rt;
    else
        trialbytrial(ich).target_rt(ifile)=nan;
    end
    if ~isempty(dataLoaded.info.behav.fix_rt)
        trialbytrial(ich).fix_rt(ifile)=dataLoaded.info.behav.fix_rt;
    else
        trialbytrial(ich).fix_rt(ifile)=nan;
    end    
    %current reward side
    trialbytrial(ich).rewardside(ifile)=0;          %0 is left side big reward
    if isfield(dataLoaded.info.behav,'bigside')
        if strncmp(dataLoaded.info.behav.bigside,'right',4)
            trialbytrial(ich).rewardside(ifile)=1;      %1 is right side big reward
        end
    else
        if strncmp(dataLoaded.info.behav.currside,'right',4)
            trialbytrial(ich).rewardside(ifile)=1;      %1 is right side big reward
        end
    end
    %big reward side (ie. block)
    trialbytrial(ich).bigside(ifile)=nan;
    if isfield(dataLoaded.info.behav,'currbigside')
        if strncmp(dataLoaded.info.behav.currbigside,'right',4)
            trialbytrial(ich).bigside(ifile)=1;
        else
            trialbytrial(ich).bigside(ifile)=0;
        end
    end
    %previous trial successful or failure
    trialbytrial(ich).prevsuccess(ifile)=0;         %0 previous trial not successful unrewarded
    if strncmp(dataLoaded.info.behav.prevevent,'reward',6)
        trialbytrial(ich).prevsuccess(ifile)=1;     %1 previous trial successful rewarded
    end
    %previous trial reward side
    trialbytrial(ich).prevside(ifile)=0;            %0 previous trial was left, rewarded
    if strncmp(dataLoaded.info.behav.prevevent(8:end),'big',3)
        trialbytrial(ich).prevside(ifile)=1; %1 previous trial was right, rewarded
    end
    %next trial side
    trialbytrial(ich).nextside(ifile)=1;
    if strncmp(infob.nexttrial,'left',4)
        trialbytrial(ich).nextside(ifile)=0;
    end
    %next trial outcome sucess/failure
    trialbytrial(ich).nextsuccess(ifile)=0;
    if strncmp(infob.nextevent,'reward',6)
        trialbytrial(ich).nextsuccess(ifile)=1;
    end
    %next trial reward size
    trialbytrial(ich).nextreward(ifile)=nan;
    if strncmp(infob.nextevent(8:end),'big',6)
        trialbytrial(ich).nextreward(ifile)=1;
    end
    if strncmp(infob.nextevent(8:end),'small',5)
        trialbytrial(ich).nextreward(ifile)=0;
    end
    %previous trial condition of big rewarded side ;
    trialbytrial(ich).prevrewside(ifile)=nan;
    if isfield(infob,'prevrewside')
        trialbytrial(ich).prevrewside(ifile)=1;
    if strncmp(infob.prevrewside,'left',4)
        trialbytrial(ich).prevrewside(ifile)=0;
    end
    end
    %previous trial reward size
    trialbytrial(ich).prevrewtype(ifile)=nan;
    if isfield(infob,'prevrewtype')
        
    if strncmp(infob.prevrewtype(8:end),'big',6)
        trialbytrial(ich).prevrewtype(ifile)=1;
    end
    if strncmp(infob.prevrewtype(8:end),'small',5)
        trialbytrial(ich).prevrewtype(ifile)=0;
    end
    end
    
    
    
    
    if ifile==1
        switchtrials{ich}=[];
        afterfailtrials{ich}=[];
    end
    %get list of trial id's where previous trial has different reward side
    if trialbytrial(ich).prevside(ifile)~=trialbytrial(ich).rewardside(ifile)
        %if different direction than current trial
        if ~isempty(switchtrials{ich})
            switchtrials{ich}=[switchtrials{ich} ifile];
        else
            switchtrials{ich}=ifile;            
        end
    end
    
    %get list of trial id's where previous trial was error
    if trialbytrial(ich).prevsuccess(ifile)==0
        %if different direction than current trial
        if ~isempty(afterfailtrials{ich})
            afterfailtrials{ich}=[afterfailtrials{ich} ifile];
        else 
            afterfailtrials{ich}=ifile;
        end
    end
    
    %display status
    percentcomplete=round(ifile/nFiles*1000)/10;
    if percentcomplete>=25 && status==0
        disp([num2str(percentcomplete) '%']);
        status=1;
    elseif percentcomplete>=50 && status==1
        disp([num2str(percentcomplete) '%']);
        status=2;
    elseif percentcomplete>=75 && status==2
        disp([num2str(percentcomplete) '%']);
        status=3;
    end

end

end

%%%%%%%
%save data
subjectlim=strfind(savepath,'_chronic');
pathlim=strfind(savepath(1:subjectlim),filesep);
subjectname=savepath(pathlim(end)+1:subjectlim-1);
pathlim2=strfind(savepath(subjectlim+1:end),'_');
sessid=savepath(subjectlim+8:subjectlim+pathlim2-1);
sesstype='bigreward';
if contains(savepath,'reward')
    isbigrew=contains(savepath,'bigreward');
    if ~isbigrew
        sesstype='smallreward';
    end
else
    %target or fix break type
    istargbreak=contains(savepath,'targetbreak');
    sesstype='targetbreak';
    if ~istargbreak
        isfixbreak=contains(savepath,'fixbreak');
        if isfixbreak
            sesstype='fixbreak';
        end
    end
end

pctype=computer;
ispc=strcmpi(pctype,'pcwin64');
homedir='Z:';
if ~ispc
    homedir=[filesep 'smbshare'];
end
graserver='inj-monkey2';
fscvdir='patra_fscv2';
savepath=fullfile(homedir,graserver,'analysis',subjectname,['chronic' sessid],sesstype,filesep);
if ~isdir(savepath)
    mkdir(savepath);
end
saveName=['trials_chronic' sessid '_' sesstype ];
disp(['saving ' saveName]);
save([savepath saveName], 'fscvchs','trialbytrial','samplesfix','samplesfixeye','samplestarg','samplestargeye','switchtrials','afterfailtrials','maxDAs','-v7.3');
savelogname=['log_compile'];
save([savepath savelogname],'logfile');
end