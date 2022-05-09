function compiletrialsdir(pathdir,varargin)
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
%configmode='patrabipolar';     %for pc, default patra
configmode='patrabichunky';     %for chunky, needs to specify dir for pcr temp
baseline='fix';     %subtract to fixation period
baseline='prefix';  %subtract to prior to fix cue - 1 second
%

while argnum<=length(varargin)
    switch varargin{argnum}
        case 'fscvchs'
            %user provided fscv selected channels
            argnum=argnum+1;
            selch=varargin{argnum};
        case 'config'
            %user provided config file
            argnum=argnum+1;
            configmode=varargin{argnum};            
    end
    argnum=argnum+1;
end

plotparam.baseline=baseline;
avgperiod=1;        %period for averaging baseline if prefix
eventIDsFix=[29 30];        %
eventIDsFixCue=4;           %fix apperance
eventIDsFixOn=5;            %start fixation on fix cue
eventIDsTarget=10;          %target apperance, 11 is when start fix on targ

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
ds=dir(folder);
ds2=find([ds.isdir]==0);            %get actual file id's
nFiles=size(ds2,2);
dsfilenames={ds(ds2).name};
PathName=folder;
FileName=dsfilenames{1};
dataLoaded=load([PathName FileName]);
samplesperscan=dataLoaded.parameters.samplesperscan;
samplespersec=dataLoaded.parameters.samplerate;
chnum=dataLoaded.chnum;

trialSortTS=zeros(1,nFiles);
samples=zeros(samplesperscan,nFiles);       %best cv's each file
Vrange=dataLoaded.Vrange;
Vrange_cathodal=dataLoaded.Vrange_cathodal;
Vrange_all=[Vrange Vrange_cathodal];

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
samplestarg={};
fscvchs=[];
logfile=[];         %log errors/warnings
for ich=1:numchs
    %scroll through channels
    chnum=str2num(d(pathdirid(ich)).name(3:end)); %get chnum from folder name
    folder=fullfile(pathdir,d(pathdirid(ich)).name,filesep); %get first folder
    ds=dir(folder);
    ds2=find([ds.isdir]==0);            %get actual file id's
    dsfilenames={ds(ds2).name};
    PathName=folder;
    disp(['processing fscv ch' num2str(chnum) ]);
    status=0;
    fscvchs(ich)=chnum;
    sumCount=1;         %for averaging, flag
    sumcolors=[];
for ifile=1:nFiles
    %scroll files for given ch
    dataLoaded=load([PathName dsfilenames{ifile}]);

   % dataLoaded=load(files{ifile});
    %only analyze data within time window exclusive to specific trial
    %load fscv data
    Itemp=dataLoaded.Ipcr.DAiso;
  %  Itemp(310)=dataLoaded.Ipcr.DAproj(310);      %***glitch in saving, all 310 samples are nan for no reason
    IM=dataLoaded.Ipcr.Mplot;
    IBG=dataLoaded.Ipcr.BGplot;
    IPH=dataLoaded.Ipcr.pHplot;
    Iox=dataLoaded.Isub.data(Vox_ID,:);
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
    
    %get relevant event ts's
    ts1=dataLoaded.ncsread.LFPts(1);        %first time point in recording
    events=dataLoaded.ncsread.LFPeventTTL;      %get event TTLs from NCS record
    idxfixTTLs=find(ismember(events,eventIDsFixCue)==1);   
    %get targeted TTL idx %07/272018 changed eventIDsFix to eventIDsFixCue for actual fix cue
    tsfix=round((dataLoaded.ncsread.LFPeventTS(idxfixTTLs)-ts1).*10)./10;   %rounded to 100 ms
    [samplesfixtemp, idxfixtemp, bb]=intersect(tsfix*samplespersec,win);         %only within window
    samplesfix{ich}(ifile)=1;       %default if not found
    if ~isempty(samplesfixtemp)        
        samplesfix{ich}(ifile)=samplesfixtemp(1);
    else
        %write log file if empty
        warning(['fix cue idx not recorded, check file ' num2str(ifile)])
            writelog=['trial ' num2str(ifile) ' / tsfix: ' num2str(tsfix)];
            %record in log file
            if isempty(logfile)
                logfile=writelog;
            else
                logfile=[logfile char(10) writelog];
            end
    end
    idxeyefixTTLs=find(ismember(events,eventIDsFixOn)==1);  %eye on fix cue
    tsfixeye=round((dataLoaded.ncsread.LFPeventTS(idxeyefixTTLs)-ts1).*10)./10;   %rounded to 100 ms
    [samplesfixeyetemp, idxfixeyetemp, bb]=intersect(tsfixeye*samplespersec,win);         %only within window
    if isempty(samplesfixeyetemp)           
        %outside window, search before 7 s
        [samplesfixeyetemp, idxfixeyetemp, bb]=intersect(tsfixeye*samplespersec,win(1)-70:win(1)+30); 
    end
   % assignin('base','tsfixeye',tsfixeye);
   samplesfixeye{ich}(ifile)=1;         %default if not found
   if ~isempty(samplesfixeyetemp)
        samplesfixeye{ich}(ifile)=samplesfixeyetemp(1);
   end
       %for some reason idx of fixation cue not recorded
       if isempty(samplesfixeyetemp)
            warning(['fix cue eye idx not recorded, check file ' num2str(ifile)])
            writelog=['trial ' num2str(ifile) ' / tsfixeye: ' num2str(tsfixeye)];
    %record in log file
            if isempty(logfile)
                logfile=writelog;
            else
                logfile=[logfile char(10) writelog];
            end
       end
    idxtargTTLs=find(ismember(events,eventIDsTarget)==1);
    tstarg=round((dataLoaded.ncsread.LFPeventTS(idxtargTTLs)-ts1).*10)./10;   %rounded to 100 ms
    samplestargtemp=intersect(tstarg*samplespersec,win);       %only within window
           %for some reason idx of target cue not recorded
    if isempty(samplestargtemp)
        warning(['target cue idx not recorded, check file ' num2str(ifile)])
        writelog=['trial ' num2str(ifile) ' / tstarg: ' num2str(tstarg)];
        %record in log file
        if isempty(logfile)
            logfile=writelog;
        else
            logfile=[logfile char(10) writelog];
        end
    end
   samplestarg{ich}(ifile)=1;         %default if not found
   if ~isempty(samplestargtemp)
        samplestarg{ich}(ifile)=samplestargtemp(1);
   end
    %samplestarg{ich}(ifile)=samplestargtemp(1);
    

    %load lfp data if exists
    if isfield(dataLoaded.nlx,'resampled')
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
            trialbytrial(ich).envwinlfp=dataLoaded.nlx.envwin;   %smooth window used for lfp
            trialbytrial(ich).filtbetah=dataLoaded.nlx.filtbetah;    %filter used for lfp
            if isfield(dataLoaded.nlx,'filtbetal')
                trialbytrial(ich).filtbetal=dataLoaded.nlx.filtbetal;    %filter used for lfp
            end
            if isfield(dataLoaded.nlx,'filtgammal')
                trialbytrial(ich).filtgammal=dataLoaded.nlx.filtgammal;    %filter used for lfp
            end
            if isfield(dataLoaded.nlx,'filttheta')
                trialbytrial(ich).filttheta=dataLoaded.nlx.filttheta;    %filter used for lfp
            end
            if isfield(dataLoaded.nlx,'filtdelta')
                trialbytrial(ich).filttheta=dataLoaded.nlx.filttheta;    %filter used for lfp
            end
            trialbytrial(ich).tscsc=dataLoaded.ncsread.LFPts;
            trialbytrial(ich).ratelfp=dataLoaded.ncsread.LFPsamplingfreq;
        end  
        if ich==1
            %only load for first ch since LFP same across channels
            %slightly varies from trial to trial but pretty much same
            %for sake of aligning below update for eacht rial
            trialbytrial(ich).relts=dataLoaded.ncsread.LFPts-dataLoaded.ncsread.LFPts(1);    %time stamp relative to 0 for each trial

            %load & store csc data
            cscdata=dataLoaded.nlx.resampled;

            %reorganize csc signals by id --> trials
            %each {ich} cell contains all data recorded for channel for all
            %trials concatenated in the cell
            for icsc=1:size(cscdata,1)
               trialbytrial(ich).lfp{icsc}(ifile,:)=cscdata(icsc,:);
                if isfield(dataLoaded.nlx,'resampledbl')
                    %other band of filtered lfp beta low
                    if icsc<=size(dataLoaded.nlx.resampledbl,1)
                    trialbytrial(ich).lfp2{icsc}(ifile,:)=dataLoaded.nlx.resampledbl(icsc,:);
                    end
                end
                if isfield(dataLoaded.nlx,'resampledgl')
                    %other band of filtered lfp beta low
                    if icsc<=size(dataLoaded.nlx.resampledgl,1)
                    trialbytrial(ich).lfp3{icsc}(ifile,:)=dataLoaded.nlx.resampledgl(icsc,:);
                    end
                end
                if isfield(dataLoaded.nlx,'resampledt')
                    %other band of filtered lfp theta
                    if icsc<=size(dataLoaded.nlx.resampledt,1)
                    trialbytrial(ich).lfp4{icsc}(ifile,:)=dataLoaded.nlx.resampledt(icsc,:);
                    end
                end
                if isfield(dataLoaded.nlx,'resampledd')
                    %other band of filtered lfp delta
                    if icsc<=size(dataLoaded.nlx.resampledd,1)
                    trialbytrial(ich).lfp5{icsc}(ifile,:)=dataLoaded.nlx.resampledd(icsc,:);
                    end
                end
            end
        end
    end
    
    %continue finding relevent event times
    relts=trialbytrial(1).relts;
    infob=computevents(dataLoaded.ncsread,[],event_codes);
    trialbytrial(ich).prevoutcomeidx=round(infob.preveventts*10)/10.*samplespersec; %fscv idx units
    nexteventidx=round(infob.nexttrialts*10)/10.*samplespersec; %fscv idx units
    trialbytrial(ich).prevtargcueidx=round(infob.prevtargcuets*10)/10.*samplespersec; %fscv idx units
    trialbytrial(ich).prevfixcueidx=round(infob.prevfixcuets*10)/10.*samplespersec; %fscv idx units
    trialbytrial(ich).prevoutcome=infob.prevevent;      %store string of previous trial outcome
    trialbytrial(ich).nextfixcueidx=round(infob.nextfixcuets*10)/10.*samplespersec; %fscv idx units
    
    %event time for fix cue appearance next window
    samplesfixnexttemp=win(end);        %default end window
    if length(tsfix)>=idxfixtemp(1)+1
        samplesfixnexttemp=tsfix(idxfixtemp(1)+1).*samplespersec;
    end
    
    %use trial min as baseline
    baselineDA=nanmin(Itemp);

    idbaseline=[samplesfix{ich}(ifile) samplestarg{ich}(ifile)];
    switch baseline
        case 'prefix'
            idbaseline=[samplesfix{ich}(ifile)-1-avgperiod.*samplespersec samplesfix{ich}(ifile)-1];
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
    [maxwithintrial, maxid]=nanmax(Itemp(samplestarg{ich}(ifile):samplesfixnexttemp));
    maxDAs{ich}(ifile).positiveposttarget=Itemp(samplestarg{ich}(ifile)+maxid-1);    
    %max change (ie decrease/increase) from target to align (rew) time
    [maxwithintrial, maxid]=nanmax(abs(Itemp(samplestarg{ich}(ifile):alignTime)));
    maxDAs{ich}(ifile).posttargetprealigndelta=Itemp(samplestarg{ich}(ifile)+maxid-1);  
    %max change (ie decrease/increase) from target to 2s after reward/align time
    [maxwithintrial, maxid]=nanmax(abs(Itemp(samplestarg{ich}(ifile):alignTime+20)));
    maxDAs{ich}(ifile).posttarget=Itemp(samplestarg{ich}(ifile)+maxid-1);  
    %average immediate amplitude (ie decrease/increase) from target to reward/align time 
    meantarget=nanmean(Itemp(samplestarg{ich}(ifile):alignTime));
    maxDAs{ich}(ifile).posttargetmean=meantarget; 
    %max change from fix appearance to prior last trial
    [maxwithintrial, maxid]=nanmax(abs(Itemp(preveventidx:samplesfix{ich}(ifile))));
    if ~any(isnan(maxid)) && ~isempty(maxid)
        maxDAs{ich}(ifile).pretrial=Itemp(preveventidx+maxid(1)-1); 
    else
        maxDAs{ich}(ifile).pretrial=nan; 
    end
    
    %max change from curr reward/align to next trial
    [maxwithintrial, maxid]=nanmax(abs(Itemp(alignTime:samplesfixnexttemp)));
    maxDAs{ich}(ifile).posttrial=Itemp(alignTime+maxid-1);     
    %if good trial, average color plot
    BGvectors=dataLoaded.Isub.data(:,idbaseline(1):idbaseline(2));   
    bgdata=mean(BGvectors,2);
    %tempdata=dataLoaded.Isub.data-bgdata;
    tempdata=dataLoaded.Isub.data;
    if ifile>1
        %not first file so have stored sumcolors
    if size(tempdata,2)~=size(sumcolors,2)
        %not equal lengths
        tempdata=tempdata(:,1:size(sumcolors,2));
    end
    end
   % if trialbytrial(ich).badtrials(ifile)==0
        if sumCount==1
            sumcolors=tempdata;        %first data stored
            sumCount=sumCount+1;        %increment count        
        elseif ifile<nFiles
            sumcolors=tempdata+sumcolors;      %add to existing data
            sumCount=sumCount+1;        %increment count        
        end
    %end
    %average all if last file
    if ifile==nFiles
        if trialbytrial(ich).badtrials(ifile)==0
            %sum with existing data if good trial
            sumcolors=tempdata+sumcolors;
            sumCount=sumCount+1;
            trialbytrial(ich).avgcolor=sumcolors./sumCount;
        else
            %if bad trial sum use existing data wihtout adding current
            trialbytrial(ich).avgcolor=sumcolors./sumCount;
        end
    end
    if ifile>1
        %not first file so have stored data
        if size(Itemp,2)~=size(trialbytrial(ich).da,2)
            %not equal lengths
            Itemp=Itemp(1:size(trialbytrial(ich).da,2));
            IM=IM(1:size(trialbytrial(ich).da,2));
            IPH=IPH(1:size(trialbytrial(ich).da,2));
            IBG=IBG(1:size(trialbytrial(ich).da,2));
            Iox=Iox(1:size(trialbytrial(ich).da,2));
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
namee=strfind(savepath,'chronic');
namm=savepath(namee:namee+8);
saveName=['trials_' namm '_' namy ];
disp(['saving ' saveName]);
save([savepath filesep saveName], 'fscvchs','trialbytrial','samplesfix','samplesfixeye','samplestarg','switchtrials','afterfailtrials','maxDAs','-v7.3');
savelogname=['log_compile'];
save([savepath filesep savelogname],'logfile');
end