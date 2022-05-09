function autotrialdir(pathdir,varargin)
global plotParam processed parameters
%10/14/2018 no longer re-integrate lfp data, will get separately in split
%files during plotTrials4
%use getparams to call settings_default for "default" recordsettings ie no
%argument from user
%so far (<07/02/2018) we have been using previous events because in beginning 
%of loop we opened events for Iread not processed.Iread, 
%NEED TO REDO ALL TRIALS
%06/30 added infobehav.prevrewside & infobehav.prevrewtype in calcBehav so
%that we can determine whether when switched big sides, if preceding trial
%informed that next trial would already be big or not (ie complete or
%partial surprise).
%need patra_map_xx config file in config folder
%path dir should be pointed to dir of bigreward_pro or other dir with fscv
%split by trials aligned to mid of file (30s default)
%automatically load saved trial multifscv files 
%create ipca data automatically and find correlated cv's within time window
%supply glitch width (default +/-4 maybe too wide, especially if hf glitch
%short)
glitchwidth=3;      %default pre 06/29/2018 was 4, for HF glitch only
nanwidth=5;      %default pre 06/29/2018 was 4, for movement signals above rthresout
timewin=[];
argnum=1;
selch=[ 1 2 3 4];            %default fscv channels to plot
configmode='patrabipolar';     %for pc, default patra
%configmode='patrabichunky';     %for chunky, needs to specify dir for pcr temp
recordsettings='default';       %default parameters
nlxchs={};
numpcs=3;
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
        case 'glitchwidth'
            %user provided glitch width (for hf electircal glitch)
            argnum=argnum+1;
            glitchwidth=varargin{argnum};
        case 'nanwidth'
            %user provide padding out for movement artifacts
            argnum=argnum+1;
            nanwidth=varargin{argnum};
        case 'timewin'
            %user provides window of time in each data to store
            argnum=argnum+1;
            timewin=varargin{argnum};
        case 'settings'
            argnum=argnum+1;
            recordsettings=varargin{argnum};
        case 'nlxchs'
            %user provides ncs channel "names" mapped by patra_map_bipolar
            argnum=argnum+1;
            nlxchs=varargin{argnum};
        case 'numpcs'
            %how many pcs to use, da, ph, m, bg
            argnum=argnum+1;
            numpcs=varargin{argnum};
            
    end
    argnum=argnum+1;
end
d=dir(pathdir);
cd(pathdir);
filenames={d.name};
targfiles=strfind(filenames,'fscv_multi');
processfiles=find(~cellfun(@isempty,targfiles));
%contains function does not work for 2013 for cells 
files=filenames(processfiles);
nfiles=length(files);

settingsfile=[];                    %load settings file (plotParam), if empty default
filtLFP=[17 34];  
csc_map={};

tParam.timeWin=[];
tParam.timeWin=[20 40];       %window to search for cv's if empty, will search all file

tParam.samplingRate=10;        %10 samples per sec
t_start=0;  t_end=60;           %in seconds
BGavg=3;                        %average BG +/- samples
event_ref_code=4;               %fix on appear cue (4) for bg averaging.
offsetbg=10;                    %use 5 s offset before event ref for bg averaging
%initialize variables that will be saved
BG=[];
Idata={};
Ipcr={};
Isub={};
Vrange=[];
Vrange_cathodal=[];
chnum=selch;
events=[];
info={};
ncsread={};
parameters={};
plotParam={};


%[parameters,csc_map,event_codes]=getparams('patra');
%[parameters,csc_map,event_codes]=getparams('patrabipolar');
[parameters,csc_map,event_codes]=getparams(configmode,recordsettings,nlxchs);

Vrange=parameters.Vrange;
Vrange_cathodal=parameters.Vrange_cathodal;

parameters.NCSchannels=[];           %load all csc data
parameters.badChannels=[];

sep=findstr(filesep,files{1});
filename=files{1};
%load single file to get default parameters
%pathdir=pathdir;
[Iread,LFPread,samplesNCS]=loadall(pathdir, filename,parameters,selch);
loadedCSCchs=[];
containslfps=~isempty(LFPread.LFPchNames);
containsda=~isempty(Iread(1).data);
removeIDs=[];
if containslfps
for ii=1:length(LFPread.LFPchNames)
    getnum=strfind(LFPread.LFPchNames{ii},'.ncs');
    loadedCSCchs(ii)=str2num(LFPread.LFPchNames{ii}(4:getnum-1));
end
nlx.cscNames=getcscnames(loadedCSCchs,csc_map);     
%remove licky & lickz, since lickx has magnitude (generated in
%reconvertFSCV), and remove eyeY (do not use)
removeIDs=find(ismember(nlx.cscNames,{'licky','lickz','eyey'})==1);
nlx.cscNames(removeIDs)=[];
%Need to remove above everytime file loaded in directory
else
    %make dummy variable in parameters, so we can load LFPread data for
    %ncsread for event ts's
    %find csc file in path and its channel id just as a reference
    matchcsc=regexpi({d.name},['csc']);
    idscsc=~cellfun('isempty',matchcsc);       %make array of zeros/1 for matching   
    fileid=find(idscsc==1);     %get file id that matches both target trial/cscnum
    if ~isempty(fileid)
        cscnumid=strfind(d(fileid(1)).name,'_');
        cscnum=str2num(d(fileid(1)).name(4:cscnumid-1));
        parameters.NCSchannels=cscnum;
        [Iread,LFPread,samplesNCS]=loadall(pathdir, filename,parameters,selch);
    end
    nlx.cscNames=[];
end

pathdirNew={};
for ich=1:length(selch)
pathdirNew{ich}=[pathdir 'analyzed' filesep 'ch' num2str(selch(ich)) filesep ];
if ~isdir(pathdirNew{ich})
    status = mkdir(pathdirNew{ich});
end
end

if ~containsda
pathdirNew=[pathdir 'analyzed' filesep ];
if ~isdir(pathdirNew)
    status = mkdir(pathdirNew);
end
end

sampleratencs=LFPread.LFPsamplingfreq;


%load ch names from patra_map based on loaded csc chs
%populate plotParam
getplotsettings(filtLFP,nlx.cscNames,event_codes,sampleratencs,settingsfile);       
nlx=plotParam;      %nlx.lfpid nlx.eyeid, nlx.. have the id's of categorized chs
nlx.selectCSC=[nlx.lfpid nlx.eyeid nlx.lickid nlx.pulseid];
parameters.BGavg=nlx.BGavg;

parameters.sampleratencs=sampleratencs;

rateLFP=LFPread.LFPsamplingfreq;
tsLFP=LFPread.LFPts;
nlx.filtlick=[0 10];        %freq bands of lick filter
    %nlx.filtbetal=[12 18];        %low beta%    changed from 13 to 12 07/07/18
    %nlx.filtbetah=[18 33];        %high beta
   % nlx.filtbeta=[13 33];           %broad beta%ake out 10/14/2018
    %nlx.filtgammal=[40 58];           %low gamma
  %  nlx.filttheta=[6 10];           %theta
   % nlx.filtdelta=[0.5 4];           %delta
status=0;
%param.NCSchannels=[];
for ifile=1:nfiles
    filename=files{ifile};
    load([pathdir, filename]);       %for non-lfp recorded signals
    [processed.Iread,LFPread,samplesNCS]=loadall(pathdir, filename,parameters,selch);
    samplesNCS(removeIDs,:)=[];     %remove unused CSC channels, defined above
    processed.LFPread=LFPread;
    ncsread=processed.LFPread;
    ncsread.LFPchNames(removeIDs)=[];
    if isempty(ncsread.LFPeventTS)
        %manually assign event ts ttl since no lfps recorded
        ncsread.LFPeventTS=NlxEventTS;
        ncsread.LFPeventTTL=NlxEventTTL;
        ncsread.LFPsamplingfreq=30e3;
        ncsread.LFPts=fscv.events(:,4);
        ncsread.nolfp=1;
        processed.LFPread=ncsread;
    end
    %get duration info on fscv data loaded
    sizeData=size(processed.Iread(selch(1)).data); 
    samplesperscan=sizeData(1);
    lengthData=sizeData(2);
    t_end=lengthData;
    nlx_events=[];
    events=[];
    Idata.anodal=[];
    Idata.cathodal=[];
    %for each channel of data (in selch) get detected da 
    %process nlx data
    %get nlx events into fscv domain
    if size(processed.Iread(selch(1)).events,1)>0
        [nlx_events,events]=readEvents(processed.Iread(selch(1)).events,event_ref_code);
    end
    %process CSC signals with input filter, if exists
    if containslfps
    for ich=1:length(nlx.selectCSC)
        %lfp signal --> filter, square, envelope
        if ismember(ich,nlx.lfpid)
            %filter at "beta-band" as defined in nlx.filtlfp
            %default filter is filtbetah now
            nlx.resampled(ich,:)=samplesNCS(ich,:);   %store unfiltered data            
            if isfield(nlx,'filtbeta')
                %if another filter band specified 07/07/18
                nlx.resampledb(ich,:)=filterLFP(samplesNCS(ich,:),rateLFP,nlx.filtbeta);
                %square & envelope signal
                winlength=round(rateLFP*.5/mean(nlx.filtbeta));
                nlx.resampledb(ich,:)=nlx.resampledb(ich,:).^2;   %get power V^2
                nlx.resampledb(ich,:)=smoothwin(nlx.resampledb(ich,:),winlength);   %smoothing 
            end
        %lick signal --> filter 
        elseif ismember(ich,nlx.lickid)
            nlx.resampled(ich,:)=filterLFP(samplesNCS(ich,:),rateLFP,nlx.filtlick);
        %eye/pulse signal --> dont do anything already filt in reconvertFSCV.m            
        else
            nlx.resampled(ich,:)=samplesNCS(ich,:);
        end
    end
    end
    
    %for each channel of fscv (in selch) get detected da 
    samples={};         %temporary storage of samples detected
    if containsda
    for ich=1:length(selch)
        info={};
        chnum=selch(ich);
        %filter all fscv data
        processed.Iread(chnum).data = filterDisplay([],[],processed.Iread(chnum).rawdata);
        %prepare process fscv data in time window selected
        [aa, closeevent]=min(abs(events-mean(tParam.timeWin).*tParam.samplingRate));
        idxbgref=events(closeevent)-offsetbg;     %1 s before first appropriate event (display fix)
        [processed.Isub(chnum).data,processed.BG(chnum).data,BGrefmatrix]=...
            refreshPlots(processed.Iread(chnum).data,tParam.samplingRate,idxbgref,nlx.BGavg,t_end); 
        [processed.Isub(chnum).rawdata,processed.BG(chnum).rawdata,bb]=...
            refreshPlots(processed.Iread(chnum).rawdata,tParam.samplingRate,idxbgref,nlx.BGavg,t_end); 
      %compute pca components of data
    %Ipcr = getITPCR(Isub,parameters,1);
 %   Ipcr = getpct(processed, parameters,chnum);
        info{chnum}.behav=calcBehav(processed,tParam.timeWin,event_codes);    
        
    %    badids=findsigartifacts(allchdata,parameters);  %timestamps
      %  badids=round(badids.*parameters.samplerate)+1;  %convert to samples

        %find ids where HF glitches, skip these data points
        glitchids=findglitches(processed.Iread(chnum).rawdata,...
            parameters.glitchThres,'glitchwidth',glitchwidth);
    %combine hf glitches w/ movement artifact periods & store
%    processed.glitchids{selch(ii)}=sort(glitchids')';
        %find CV maximally correlated to DA within time window of trial
   % samplesTemp=detectDA(dataFSCV,parameters,parameters,tParam);       %input dataFSCV, output samplesTemp
        %get pct dopamine signal to remove M/Q thres & glitches
        %{
        Ipcr = getpctdirect(processed.Isub(chnum).data,...
            'ats',parameters.CV,'cts',parameters.Cts,...
            'qthres',parameters.QaDAPH,'rthresm',parameters.RThresOut,...
            'imaxbg',max(abs(processed.Iread(chnum).data(:,idxbgref))),...
            'glitchids',glitchids,'nanwidth',nanwidth,'numpcs',2);  
        %}
        %01/2019, get 3 component PCA signal with DA/PH/M (NO BG)
        %previously just 2 above
        Ipcr = getpctdirect(processed.Isub(chnum).data,...
            'ats',parameters.CV,'cts',parameters.Cts,...
            'qthres',parameters.QaDAPH,'rthresm',parameters.RThresOut,...
            'imaxbg',max(abs(processed.Iread(chnum).data(:,idxbgref))),...
            'glitchids',glitchids,'nanwidth',nanwidth,'numpcs',numpcs);  
        %09/06/2018, added 'numpcs' to getpctdirect argument to reduce #
        %pcs from 3 to just da & ph to increase SNR
        detectedda=detectdatransientschunks(...
            processed.Iread(chnum).data,parameters,...
            glitchids,plotParam,'timestamp',...
            processed.LFPread.LFPts(1)); 
        if ~isempty(detectedda.maxTS)
        if ~isempty(detectedda.cv)
            rr=abs(corr(detectedda.cv,parameters.K(:,1)));         %correlate retrieved CV's to da template & keep best
            maxCV=detectedda.cv(:,find(rr==max(rr)));  
            Idata.anodal=maxCV(1:round(samplesperscan/2));
            Idata.cathodal=maxCV(round(samplesperscan/2+1):round(samplesperscan));          
        end
        end
        Iread=processed.Iread(chnum);           %fixed 06/30/2018, 
        %so far we have been using previous events because in beginning 
        %of loop we opened events for Iread not processed.Iread, 
        %NEED TO REDO ALL TRIALS
        Isub=processed.Isub(chnum);
        info=info{chnum};
        BG=processed.BG(chnum).data;
                saveName=['CVITdata_ch' num2str(chnum) '_' filename];
        save([pathdirNew{ich} saveName],'Ipcr','nlx','ncsread','parameters','Idata','BG',...
            'Vrange','Vrange_cathodal','Isub','info','events','chnum');
    end    
    else
        %no da recorded, just lfps
                        saveName=['CVITdata_ch' '_' filename];
                save([pathdirNew saveName],'nlx','ncsread','parameters','Idata',...
            'events');
        
    end
    %display status
    percentcomplete=round(ifile/nfiles*1000)/10;
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
