function getTRdata(sessnum,varargin)
%Get trial by trial data and save into format easily loadable for analysis
%Saved in \analysis\patra\chronicXX\tr\
%04/30/2022 - attempt to simplify analysis code from extracted trial by
%trial FSCV and LFP
%used code from plotsession 
%assume Patra single subject for now
% Major changes:
%   Make indiscriminant of trial type (i.e. big and small together, ideally with error trials all consolidated into single data file together, but most essions error trials have not been extracted)
%   Separate FSCV analysis and variable properties from LFP, too mixed up
%   in original code
%   Unified time stamps in (s) for FSCV and LFP rather than going to FSCV
%   sample periods (100 ms) and then converting from LFP to this range
% Minor changes:
%   Output analysis data into single folder in local directory (rather than
%   server that takes more time to transfer/upload) and input analysis from
%   local directory if possible unless otherwise specified

%%NEED CHECKING
%While making this new function realized PL1-p1 not created in chronic 83
%in chunky copied folders now in our folder, need to check if in external?
%%DEV NOTES
%With 10 ephys channels and 4 DA channels embedded in trlists =5.27 GB trlists file
%Save each channel separately in the trlist format?
%Saving just DA channels embedded in list format = 3.24 GB (???) The total
%"" for 10 ephys channels = 2 GB (~200 MB when separated)
%of the original separate files in big reward = 370 MB ~1.5 GB for all 4 types(Struct really makes it big just because of extra array for good/bad?)
%Save just all DA channel in list format = 10 MB, save ephys channels

%%TESTING FUNCTIONS / VARS
%sessnum=83; varargin={'fscvchs',1:4};


%default analysis folder in matlab user path
analysisdir=fullfile(userpath,'fscv','analysis',filesep);    %default analysis dir
configdir=fullfile(analysisdir,'config',filesep);   %path with config files for data directory and channels to be analyzed
%Configdir also stores excel sheet for patra
%REMOVE THESE WHEN IMPROVE to indicate on spreadsheet
homedirsplit=split(userpath,filesep);%Loaded in chronicXXchconfig file
homedir=fullfile(homedirsplit{1:end-1},filesep);    %C:\Users....\Documents\"

%default fscv channels to analyze (up to 4)
maxfscvchs=4;
fscvchs=1:4;
subjectname='patra';
ratelfp=1e3;

%Analysis window
trialalnTS=30;    %all fscv data made to align at sample 300 (30 s) which is outcome TS, replaces plotparam.alignidx=300
win_an=[-10 5];                   %Time window for analysis, instead of time=[25 35] around rew outcome, now relative to actual targeted alignment event
alnevt='targeye';       %aln to when eye on target, base sub before



sessid=num2str(sessnum);        %session number converted to string
argnum=1;
%select type of analysis based on input argument
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'fscvchs'
            %user provided fscv selected channels
            argnum=argnum+1;
            fscvchs=varargin{argnum};
            disp(['fscv chs ' num2str(fscvchs)]);
        case 'analysisdir'
            %provide directory for fscv analysis scripts
            analysisdir=varargin{argnum};
        
    end
    argnum=argnum+1;
end

%load config file for session to get directory and nlx channels to load
%get paths{1} & paths{2} and ncschannels
run([configdir ['chronic' sessid 'chconfigsimple.m']]);            %get ncschannels
run([configdir 'patra_map_bipolar']);   %run patra_map_bipolar for ch settings

%%%%IMPROVEMENT NEEDED HERE
%get paths{1} just use spreadsheet that has list of sessions recorded and
    %their directories, right now comes from chronicXXchconfig loaded above                                                     

%run config file to get ncschannels & paths e.g. chronic83chconfigsimple.m
homedirf=strfind(paths{1},'1dr');
subjectdir=paths{1}(1:homedirf-1);
manlist=[subjectdir 'chronic' sessid '_trialselection.xlsx'];      %list for manually selected bad trials
%assignin('base','fcnStatus',targpathtemp)   %store targpath in workspace
manlistlfp=[subjectdir 'chronic' sessid 'l_trialselection.xlsx'];      %list for manually selected bad trials

%Data PATHS for trial by trial fscv data/nlx parameters for each trial type
fscvpathtemp.bigreward=fullfile(paths{1}, 'matlab','bigreward_pro',filesep);
fscvpathtemp.smallreward=fullfile(paths{1}, 'matlab','smallreward_pro',filesep);
fscvpathtemp.targetbreak=fullfile(paths{1}, 'matlab','targetbreak_pro',filesep);
fscvpathtemp.fixbreak=fullfile(paths{1}, 'matlab','fixbreak_pro',filesep);

%Original analysis paths (with data analyzed from compiletrialsdir and
%extractsession), mainly used to just get good trials, not for actual
%analysis data, except trlist.mat
rootOG=strfind(subjectdir,'patra_fscv');    %Get root index to re-index to analysis folder
OGdir=subjectdir(1:rootOG-1);
analysispathOG.home=fullfile(OGdir,'analysis',subjectname,['chronic' sessid],filesep);
if ~isfolder(analysispathOG.home)
   error([analysispathOG ' does not exist']);
end
analysispathOG.bigreward=fullfile(OGdir,'analysis',subjectname,['chronic' sessid],'bigreward',filesep);
analysispathOG.smallreward=fullfile(OGdir,'analysis',subjectname,['chronic' sessid],'smallreward',filesep);
analysispathOG.targetbreak=fullfile(OGdir,'analysis',subjectname,['chronic' sessid],'targetbreak',filesep);
analysispathOG.fixbreak=fullfile(OGdir,'analysis',subjectname,['chronic' sessid],'fixbreak',filesep);


%NEW PATHS in SEPARATE directory in home directory
analysispath.home=fullfile(homedir,'analysis',subjectname,['chronic' sessid],filesep);
if ~isfolder(analysispath.home)
    mkdir(analysispath.home);
end
analysispath.bigreward=fullfile(homedir,'analysis',subjectname,['chronic' sessid],'bigreward',filesep);
analysispath.smallreward=fullfile(homedir,'analysis',subjectname,['chronic' sessid],'smallreward',filesep);
analysispath.targetbreak=fullfile(homedir,'analysis',subjectname,['chronic' sessid],'targetbreak',filesep);
analysispath.fixbreak=fullfile(homedir,'analysis',subjectname,['chronic' sessid],'fixbreak',filesep);
trialgroups=fieldnames(analysispath);
trialgroups=trialgroups(2:end);%rm home, only what trial conditions

%%
%Get Data Settings
%Need to get parameters for recording, e.g. sampling rates, NLX CSC site
%labels (e.g. putamen, CN, non-neural phys), etc.
%Not sure if call to getparams REALLY needed, a lot of the code can be
%simplified, most of getparams is for extractsession process
%Merge getparams and getplotsettings and merge param and plotparam vars
%plotparam now part of param as param.plotparam
[param,csc_map,eventcodes]=ef_getparams('patrabipolar','default',ncschannels,'sessnum',sessnum);
%Also not sure getplotsettings needed, just take variables really needed
%and rewrite
plotparam=param.plotParam;  %from ef_getparams
plotparam.pathname=paths{1};        %store data path
ephysIDs=ncschannels(plotparam.lfpid);    %ncschannels from chronicXXchconfigsimple.m loaded above
auxIDs=ncschannels(plotparam.physid);
plotparam.chnums=fscvchs;
plotparam.ephysIDs=ephysIDs;

%Also Modify - (remove) setup rates
plotparam.samplespersec=param.rateFSCV;   %fscv sample rate..

%Shouldn't need this here - set up fscv temporal domain
%This should all be in getplotsettings2 also
%Setup analysis timeframe window
%win=time(1)*param.rateFSCV:1:time(2)*param.rateFSCV;
%plotparam.win=win;  %In samples for fscv, can we convert to Seconds????
%plotparam.interval=2.5;
%Don't define baseline here should be defined during analysis of da
%according to alignment point
%plotparam.baseline='pretarg';           %immediate post-cue/fix appear
%if ~isempty(baseline)
    %user provided
 %   plotparam.baseline=baseline;
%end
plotparam.alignidx=300;         %align to 30s Reward or outcome period, replaced now with trialalnTS=30;
plotparam.trialalnTS=trialalnTS;%timestamp for signal alignment (usually 30s, outcome, as defined in compiletrials as the timepoint around which signals for each trial are aligned to)
plotparam.alnevt=alnevt;        %default alnevt (up above) = targeye

%%
%Load all trial by trial data files, according to user-selected channels
%and trial conditions
%Compile all data, excluding bad trials for each trial condition type
%Currently (previously) loaded trials_chronicXX_bigreward that had da for
%all trials including bad ones
%Clean up the structure and resave with all trial conditions merged

%First load Master trial list (as created by gettrialorder.m) this will form the
%table to list all trials and store all the trial by trial da and ephys
%data extracted below

load([analysispathOG.home 'trlist.mat'],'trlist');  %load trlist from analysis path
%going to incorporate measured data into this master trial list
trlists.trlist=trlist;
maxtrials=length(trlist);
trlists.fscv={};
for ch=1:length(fscvchs)
    %initialize master list so all trials bad as default for each da ch
    ich=fscvchs(ch);
    C=num2cell(false(maxtrials,1)); %how to assign multiple fields to struct array 
    [trlists.fscv{ich}(1:maxtrials).good]=deal(C{:});%set logical false for all channels no good trials initialized
end
%store fscv site names
origpath=fscvpathtemp.bigreward(1:strfind(fscvpathtemp.bigreward,'cvtotxt')-1);
fscvnames=getfscvsites(origpath);%get fscv ch names (e.g. pl1, pl2, etc.)
param.fscvnames=fscvnames;
sites=getsites(sessnum,fscvnames,subjectname,'path',configdir);%get site names based on sess #, actual site in brain (based on depth changes)
trlists.fscvnames=fscvnames;
trlists.fscvsites=sites;

trlists.nlx={};
%Get channel number ids for selected sites or phys signal from NLX
%use ncschannels from chronicXXchconfig or varargin using ephysIDs / auxiDs
trlists.nlxchs=param.NCSchannels;
trlists.nlxnames={ephysIDs{:} auxIDs{:}};
nlxsites=getlfpsites(sessnum,trlists.nlxnames,subjectname,'path',configdir);%get site names based on sess #, actual site in brain (based on depth changes)
trlists.nlxsites=nlxsites;

%Each FSCV channel and Ephys gets its own trial list (redudant) to allow
%for fast data extraction and not too big lists that incorporate all data
%   good: indicating whether good or bad trial for fscv ch (logical)

%%
%Load trial by trial da and ephys data, as stored separately for the
%different task conditions, to compile together regardless of task cond

goodtrials={};          %each cell refers to a condition
ich=1;
for itype=1:length(trialgroups) 
    display(['Current trial condition: ' trialgroups{itype}])
    currpath=getfield(analysispathOG,trialgroups{itype});     %get data path
    plotparam.pathname=currpath;    %store data paths to read into analysis functions
    plotparam.trialtype=trialgroups{itype};

    if strcmp(trialgroups{itype},'fixbreak')
        %shift time to outcome/error, for fix break (no target for alnevt)
        plotparam.alnevt='outcome';        
    end 

    %Load compiled data (from compiletrialsdir.m from extractsession.m)
    compiledname=[currpath 'trials_chronic' sessid '_' trialgroups{itype}];
    load(compiledname,'trialbytrial','samplesfix','samplesfixeye',...
        'samplestarg','samplestargeye','switchtrials','afterfailtrials'); %reloads fscvchs

    numtrials=size(trialbytrial(ich).da,1);
    chbadtrials=[];
    alltrials=1:numtrials;
    goodtrials{itype}{ich}=1:numtrials; %default, all trials goodtrials 5/2019, ch-specific
    %open bad trials if exists in homedir, based on manual fscv data
    %evaluation
    if exist(manlist,'file')
       % disp(['getting good trials from manual list for all fscv chs']);
        goodtrials{itype}=manualtrialsel(manlist,ich,alltrials,trialgroups{itype});          %get recent list of bad trials manually selected in excel sheet
        percgood=size([goodtrials{itype}{:}],2)/maxfscvchs/numtrials*100;
        display(['~' num2str(percgood) '% trials good for analysis for ' trialgroups{itype}]);
    end   
    goodtrialslfp={};
    %open bad trials if exists based on manual lfp data (most for cleo)
    if exist(manlistlfp,'file')
        goodtrialslfp=manualtrialsel(manlistlfp,1,alltrials,trialgroups{itype});          
        goodtrialslfp=goodtrialslfp(1);
    end     

    %remove trials manually if have task event timestamps that don't make
    %sense
    if any(samplesfixeye{1}<100) || any(samplesfixeye{1}<samplesfix{1})
        %if time stampe for cue fixation is less than 100 samples (i.e.
        %before 10 s, then it must be invalid trial because the trial
        %outcom is fixed at 30 s and our task parameter settings were never such
        %that the fixation was such a long time before the outcome
        chbadtrials2=find(samplesfixeye{1}<100);
        chbadtrials2=unique([chbadtrials2 find(samplesfixeye{1}<samplesfix{1})]);
        for iich=1:4
            goodtrials{itype}{iich}=goodtrials{itype}{iich}(find(~ismember(goodtrials{itype}{iich},chbadtrials2)));
        end
        warning([num2str(length(chbadtrials2)) ' initial fix task events outside normal trial time range']);
    end
    if any(samplestarg{1}<samplesfixeye{1}) || any(samplestarg{1}<100)
        chbadtrials2=find(samplestarg{1}<samplesfixeye{1});
        chbadtrials2=unique([chbadtrials2 find(samplestarg{1}<100)]);
        for iich=1:4
        goodtrials{itype}{iich}=goodtrials{itype}{iich}(find(~ismember(goodtrials{itype}{iich},chbadtrials2)));
        end
        warning([num2str(length(chbadtrials2)) ' target task events outside normal trial time range']);
    end
    
    %convert good trial nums for each channel from this task condition 
    % to trial # for all trials conditions from trlist

    %Get indeces of trlist corresponding to current task condition 
    fulltypeids=find(ismember({trlist.type},trialgroups{itype}));   %IDs in full list for the current task condition trials

    %loop through each channel to extract all trial by trial data
    %index into master list, whether good/bad trial, trial #, etc
    %get DA data
    for ch=1:length(fscvchs)
        ich=fscvchs(ch);
        %load([analysispath 'trialtypes']); %good trials load, not useful
        %if going to use the list of all trials irrespective of
        %type/condition...
        alltrialids=alltrials+99;    %convert trial # to trial id (orig trial file id starts at 100)
        chgoodtrialids=[];
        if isempty(goodtrials{itype})
                %if empty, may not have been analyzed yet manually
                disp('Empty good trial list, need to do manual sel');
        else
            chgoodtrialids=goodtrials{itype}{ich}+99;      %convert trial # to trial id within current task condition still
        end
        fulltype=find(ismember({trlist.type},trialgroups{itype}) & ismember([trlist.id],alltrialids));        %IDs in full list corresponding to current task condition 
        fullgood=find(ismember({trlist.type},trialgroups{itype}) & ismember([trlist.id],chgoodtrialids));        %IDs in full list corresponding to current task condition and the good trial #s
        %goodarray=false(1,maxtrials);
        good=true(length(fullgood),1);
        %Convert logical array 2 cell to assign to struct array master list
        C=num2cell(good);
        %Create good/bad trial IDs in master list based on condition (itype) loaded here 
        [trlists.fscv{ich}(fullgood).good]=deal(C{:});    
        %Store trial by trial da signals of this condition into master list
        for it=1:numtrials
            %[trlists.fscv{ich}(fulltype).da]=deal(trialbytrial(ich).da);    
            %C=num2cell(trialbytrial(ich).da);
            trlists.fscv{ich}(fulltype(it)).da=trialbytrial(ich).da(it,:);   
        end
    end

    %Loop through trial by trial ephys data for selected csc channels   
    %Load ephys data, trial by trial (single file each)
    %More complicated, as need to load each trial individually from
    %original directory where trials were extracted
    data={};
    %get reconverted path ie XX_pro directory    
    ephyspath=getfield(fscvpathtemp,trialgroups{itype});     %get data path
    dirlfp=dir(ephyspath);
    %get csc signals defined in NCSchannels
    if isempty(param.NCSchannels)
        error('no ncs channels defined in param');
    end
    for it=1:numtrials
        disp(['Getting Nlx data: trial # ' num2str(it)]);
        tlab=it+100-1;
        %Get row/id in master list for current trial (it)
        trid=find(ismember({trlist.type},trialgroups{itype}) & ismember([trlist.id],tlab)); 
        %Loop through all Nlx channels (ephys + phys)
        for idx=1:length(param.NCSchannels)
            %Load each ch in order of appearing in NCSchannels list
            matchcsc=regexpi({dirlfp.name},['csc' num2str(param.NCSchannels(idx)) '_']);
            idscsc=~cellfun('isempty',matchcsc);       %make array of zeros/1 for matching
            matchtrial=regexpi({dirlfp.name},['_' num2str(tlab)]);
            idstrial=~cellfun('isempty',matchtrial);  
            fileid=find((and(idscsc,idstrial))==1);     %get file id that matches both target trial/cscnum
            tempload=load([ephyspath,dirlfp(fileid).name],'dg_Nlx2Mat_Samples','dg_Nlx2Mat_Timestamps');    
            if idx==1
                %Just the 1st channel, store TS, same for all nlx chs                  
                trlists.trlist(trid).NlxTS=tempload.dg_Nlx2Mat_Timestamps.*1e-6;  %in s
            end
            if it>1 && length(tempload.dg_Nlx2Mat_Samples)~=size(data.lfp{idx},2)
                %if for some reason the length of the data is different
                %from other trials, make it NaN (invalidate)
               % data.lfp{idx}(it,:)=repmat(nan,1,size(data.lfp{idx},2));
                trlists.nlx{idx}(trid).ephys=repmat(nan,1,size(data.lfp{idx},2));
            else
                data.lfp{idx}(it,:)=tempload.dg_Nlx2Mat_Samples;  %store samples
                trlists.nlx{idx}(trid).ephys=tempload.dg_Nlx2Mat_Samples;  %store samples
            end
            %timestamps in tempload are in microseconds and should match tsLFP already loaded 
            %get/store csc name
           % cscnameb=strfind(dirlfp(fileid).name,'_');
            %cscname=[dirlfp(fileid).name(1:cscnameb-1) param.LFPterm];                
            %data.cscNames{idx}=cscname;
        end
    end
end

%Remove data files from trlists to be saved separately (otherwise too big)
savepath=[analysispath.home 'tr' filesep];
if ~isfolder(savepath)
    mkdir(savepath);
end
%Save fscv data
tr_fscv=trlists.fscv;
chname=trlists.nlxnames{idx};
savename=['tr_fscv'];
save([savepath savename],'tr_fscv','-v7.3');
%Save Nlx data, each channel separately
for idx=1:length(param.NCSchannels)
    %Get channel name (e.g. pl1-p5, or pulse)
    chname=trlists.nlxnames{idx};
    savename=['tr_nlx_' chname];
    tr_nlx={trlists.nlx{idx}.ephys}';   %Extract data from trlists for this ch
    save([savepath savename],'tr_nlx','-v7.3');
end
%Save trlists data
trlists=rmfield(trlists,{'nlx','fscv'});%Remove data fields and then save
save([savepath 'trlists'],'trlists')

%Plot testing
%xx=vertcat(tr_nlx{[5:15, 20]});
%figure; plot(xx')
%xx=vertcat(tr_fscv{4}(100:102).da);
%figure; plot(xx')





