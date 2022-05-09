%function updateTRlists(sessnum,varargin)
%05/2022 same as getTRdata, except just updating good/bad trial indeces in
%loaded trlists because did the manual selection of good/bad trials and
%updated the spreadhseets in original direcotry
%assume Patra single subject for now

%%TESTING FUNCTIONS / VARS
sessnum=83; varargin={'fscvchs',1:4};


%default analysis folder in matlab user path
analysisdir=fullfile(userpath,'fscv','analysis',filesep);    %default analysis dir
configdir=fullfile(analysisdir,'config',filesep);   %path with config files for data directory and channels to be analyzed
%Configdir also stores excel sheet for patra

%default fscv channels to analyze (up to 4)
maxfscvchs=4;
fscvchs=1:4;
subjectname='patra';


%REMOVE THESE WHEN IMPROVE to indicate on spreadsheet
homedirsplit=split(userpath,filesep);
homedir=fullfile(homedirsplit{1:end-1},filesep);    %C:\Users....\Documents\"

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

%run config file to get ncschannels & paths
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

%%
%First load Master trial list trlists (generated from getTRdata) containing
%fscv data for each channel, e.g. tr_fscv files
load([analysispath.home 'tr' filesep 'tr_fscv.mat'],'tr_fscv');  %load tr fscv data
load([analysispath.home 'tr' filesep 'trlists.mat'],'trlists');  %load trlists, need trlist for cond type for each trial, since original files separated by cond
trlist=trlists.trlist;
maxtrials=length(trlist);

%%
%Load trial by trial da and ephys data, as stored separately for the
%different task conditions, to compile together regardless of task cond
%   good: indicating whether good or bad trial for fscv ch (logical)


goodtrials={};          %each cell refers to a condition
ich=1;
for itype=1:length(trialgroups) 
    display(['Current trial condition: ' trialgroups{itype}])
    currpath=getfield(analysispath,trialgroups{itype});     %get data path
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

    %index into master list, whether good/bad trial, trial #, etc
    for ch=1:length(fscvchs)
        ich=fscvchs(ch);     
        alltrialids=alltrials+100-1;    %convert trial # to trial id (orig trial file id starts at 100)
        chgoodtrialids=[];
        if isempty(goodtrials{itype})
                %if empty, may not have been analyzed yet manually
                disp('Empty good trial list, need to do manual sel');
        else
            chgoodtrialids=goodtrials{itype}{ich}+100-1;      %convert trial # to trial id within current task condition still
        end
        fulltype=find(ismember({trlist.type},trialgroups{itype}) & ismember([trlist.id],alltrialids));        %IDs in full list corresponding to current task condition 
        fullgood=find(ismember({trlist.type},trialgroups{itype}) & ismember([trlist.id],chgoodtrialids));        %IDs in full list corresponding to current task condition and the good trial #s
        good=true(length(fullgood),1);
        %Convert logical array 2 cell to assign to struct array master list
        C=num2cell(good);
        %Update good/bad trial IDs in master list based on condition (itype) loaded here 
        [tr_fscv{ich}(fullgood).good]=deal(C{:});          
    end

end

%Save updated fscv trial data with updated good/bad trial ids
savepath=[analysispath.home 'tr' filesep];
%Save fscv data
fdata=trlists.fscv;
chname=trlists.nlxnames{idx};
savename=['tr_fscv'];
save([savepath savename],'fdata','-v7.3');





