function analyzeSes(sessnum,varargin)
persistent trlists
%05/2022
%Create workspace for analysis by loading relevant sesion data
%Run trial by trial analysis using trlists and tr_ data stored in 
%\analysis\patra\chronicXX\tr\
%04/30/2022 - attempt to simplify analysis code from extracted trial by
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

%load config file for session to get directory and nlx channels to load
%get paths{1} & paths{2} and ncschannels
run([configdir ['chronic' sessid 'chconfigsimple.m']]);            %get ncschannels
run([configdir 'patra_map_bipolar']);   %run patra_map_bipolar for ch settings


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
            argnum=argnum+1;
            analysisdir=varargin{argnum};        
        case 'nlxch'
            %User provides nlx channels to analyze
            %Overwrite ncschannels from chronicXXchconfig loaded above
            argnum=argnum+1;
            ncschannels=varargin{argnum};    
    end
    argnum=argnum+1;
end


%%%%IMPROVEMENT NEEDED HERE
%get paths{1} just use spreadsheet that has list of sessions recorded and
    %their directories, right now comes from chronicXXchconfig loaded above                                                     

%run config file to get ncschannels & paths
homedirf=strfind(paths{1},'1dr');
subjectdir=paths{1}(1:homedirf-1);

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
param={};
csc_map={};
eventcodes={};
%Not sure if call to getparams REALLY needed, a lot of the code can be
%simplified, most of getparams is for extractsession process
%Merge getparams and getplotsettings and merge param and plotparam vars
%plotparam now part of param as param.plotparam
[param,csc_map,eventcodes]=ef_getparams('patrabipolar','default',ncschannels,'sessnum',sessnum);
%not sure getplotsettings needed, just take variables really needed
plotparam=param.plotParam;  %from ef_getparams
plotparam.pathname=paths{1};        %store data path
ephysIDs=ncschannels(plotparam.lfpid);    %ncschannels from chronicXXchconfigsimple.m loaded above
auxIDs=ncschannels(plotparam.physid);
plotparam.chnums=fscvchs;
plotparam.ephysIDs=ephysIDs;

%Also Modify - (remove) setup rates
ratelfp=param.ratelfp;
rates=[ceil(param.ratelfp) param.rateFSCV];
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
%Load master list and all trial by trial data, with user-selected channels
load([analysispath.home 'tr' filesep 'tr_fscv.mat'],'tr_fscv');  %load tr fscv data
load([analysispath.home 'tr' filesep 'trlists.mat'],'trlists');  %load trlists, need trlist for cond type for each trial, since original files separated by cond
trlist=trlists.trlist;
trlists.param=param;
trlists.fscv=tr_fscv;
maxtrials=length(trlist);
nlxids=find(contains(trlists.nlxnames,ncschannels)); 
trlists.nlxids=nlxids;
trlists.fscvids=fscvchs;
trlists.path=analysispath;
for nlxch=1:length(ncschannels)
    idx=find(contains(trlists.nlxnames,ncschannels{nlxch}));    %get orig ch idx based on all chs saved
    savename=['tr_nlx_' ncschannels{nlxch}];
    load([analysispath.home 'tr' filesep savename],'tr_nlx');
    trlists.nlx{idx}=tr_nlx;%Store nlx data for ch into trlists
end
assignin('base','trlists',trlists)
