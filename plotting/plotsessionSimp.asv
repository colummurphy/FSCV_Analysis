function plotsessionSimp(sessnum,varargin)
%script to run all functions to get compiled trial data lfp/da
%assume patra recording
%opens chconfigsimple to load targeted ncs channels
%many updates following feedback 01/24/2019
homedir=userpath;
homedirsplit=split(userpath,filesep);
homedir=fullfile(homedirsplit{1:end-1},filesep);    %C:\Users....\Documents\"

global plotParam
plotParam={};
argnum=1;
fscvchs=1:4;
targch=2;       %fscv ch used for sorting or getting good trials from list
fbands={};
baseline=[];
fbands2={};         %second band to plot (hf)
fbands3={};         %global beta, need to input argument 'globeta' for only this band
filtgamma=[35 60];
filtgamma=[35 55];
globeta=0;
selbands=[1 2 3];
selbands=[3];            %default 05/2019
chbands;                    %default 05/2019
%plot flags
plotda=0;
plotlfp=0;
plotgam=0;
plotavggam=0;
plotavgda=0;
plotavglfp=0;
plotphys=0;
plotfft=0;
plotcor=0;
plotxvar=0;
scales=[];
replace=0;
subjectname='patra';
            plotxvarall=0;      %selective xvar's
            plotbetaxcov=0;
            plotdaxcov=0;
xchs={};
%plot session flags
plotbig=0;
plotsmall=0;
plottargbreak=0;
plotfixbreak=0;
plotxclust=0;
plotxsess=0;
plotsort=0;
allsess=0;
getinfo=0;
trialtypes={};
%targevents={'intertarg','interfix','targwin'};
targevents={'intertarg','interfix','intertargeye','interfixeye'};
targevents={'interfix','intertarg','intertargeye'};
targevents={'targeye','targ','fix','outcome'};

alnevt='fixeye';
alnevt='targ';
alnevt='fix';
alnevt='targeye';       %aln to when eye on target, base sub before
plotxvarbeh=0;
plotxvarlfp=0;
plotcorall=0;
eventx={};
ncselect={};
lfpflag=0;
daflag=0;
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'fscvchs'
            %user provided fscv selected channels
            argnum=argnum+1;
            fscvchs=varargin{argnum};
            disp(['fscv chs ' num2str(fscvchs)]);
        case 'targch'
            %user provided targ ch
            argnum=argnum+1;
            targch=varargin{argnum};       
        case 'bands'
            %user provides lfp bands to analyze
            argnum=argnum+1;
            fbands=varargin{argnum}; 
        case 'baseline'
            %fscv sub baseline, options, 'pretarg', 'alltrial', 'precue',
            argnum=argnum+1;
            baseline=varargin{argnum};
        case 'bands2'
            %user provides another band of lfps to analyze
            argnum=argnum+1;
            fbands2=varargin{argnum};
        case 'dualbeta'
            %default dual betabands from chbands.m file
            %load fbands and fbands2 from chbands.m
                        %provided groups of beta bands
            %[1 2 3] means all groups low, high, broad
            argnum=argnum+1;
            chbands;
             selbands=varargin{argnum};
        case 'globeta'
            %just plot global 13-28 Hz band of beta defined in chbands
            chbands;            %get fbands3
            globeta=1;
        case 'xchs'
            argnum=argnum+1;
            xchs=varargin{argnum};
        case 'alnevt'
            %event to aln, 'outcome','targeye','fixeye','fix','targ'
            argnum=argnum+1;
            alnevt=varargin{argnum};
        case 'allsessions'
            %plot all sessions, big, small etc.
            plotbig=1;
            plotsmall=1;
            plottargbreak=1;
            plotfixbreak=1;
        case 'plotbig'
            plotbig=1;
        case 'plotsmall'
            plotsmall=1;
        case 'plottargbreak'
            plottargbreak=1;
        case 'plotfixbreak'
            plotfixbreak=1;         
        case 'plotall'
            %plot all types of signals
            plotda=1;
            plotlfp=1;
            plotavgda=1;
            plotavglfp=1;
            plotphys=1;
            plotfft=1;
            plotcor=1;
            plotxvar=1;
            plotxvarall=1;
            plotxvarlfp=1;
            plotxvarbeh=1;
        case 'lfponly'
            plotlfp=1;
            plotavglfp=1;
            plotphys=1;
            plotfft=1;
            lfpflag=1;
            plotda=1;
            fscvchs=[];
        case 'daonly'
            plotda=1;
            plotavgda=1;
            plotphys=1;
            plotxvar=1;
            daflag=1;
        case 'plotda'
            plotda=1;
        case 'plotlfp'
            plotlfp=1;
        case 'plotgam'
            plotgam=1;
        case 'plotavggam'
            plotavggam=1;
        case 'plotavgda'
            plotavgda=1;
        case 'plotavglfp'
            plotavglfp=1;
        case 'plotphys'
            plotphys=1;
        case 'plotfft'
            plotfft=1;
        case 'plotcor'
            plotcor=1;
        case 'plotxvarall'
            plotxvar=1;
            plotxvarall=1;
            plotxvarlfp=1;
            plotxvarbeh=1;
        case 'plotxvarbeh'
            plotxvar=1;
            plotxvarall=1;
            plotxvarlfp=0;
            plotxvarbeh=1;
        case 'nofft'
            plotfft=0;
        case 'plotxclust'
            allsess=1;
            plotxclust=1;       %only after getting all data
            %must supply 'targevents' ie. 'intertarg' after
        case 'plotxsess'
            allsess=1;
            plotxsess=1;       %only after getting all data
        case 'plotsort'
          allsess=1;
            plotsort=1; %also after getting all data after plotxclust
        case 'plotcorall'
            allsess=1;
            plotcorall=1;
        case 'plotpostall'
            %all post analysis
            allsess=1;
            plotcorall=1;
            plotxsess=0;
            plotxclust=1;
        case 'targevents'
            argnum=argnum+1;
            targevents=varargin{argnum};
        case 'getinfo'
            getinfo=1;
        case 'trialtypes'
            argnum=argnum+1; %user provides {'all','phase1', 'phase4'}
            trialtypes.names=varargin{argnum};
        case 'events'
            argnum=argnum+1;    %user selects names which event windows to xcov
            eventx=varargin{argnum};
        case 'replace'
            %replace exported lfp/phys data, so redo settrialax
            replace=1;
        case 'csc'
            %select manual csc channels
            argnum=argnum+1;
            ncselect=varargin{argnum};
        case 'cleo'
            subjectname='cleo';
        case 'daflag'
            daflag=1;
        case 'lfpflag'
            lfpflag=1;
    end
    argnum=argnum+1;
end

if length(selbands)==1 && selbands(1)==3
    globeta=1;
end
sessid='';
if ~isnumeric(sessnum)
    sessid=sessnum;
    sessnum=str2num(sessnum(regexpi(sessnum,'\d')));        %convert to number
else
    sessid=num2str(sessnum);
end

%get dir with config files
pctype=computer;
ispc=strcmpi(pctype,'pcwin64');
%default on putamen pc in lab
configdir='A:\mit\injectrode\experiments\fscv\matlab\analysis\analysis\config\';
if ~ispc
    %chunky dir
    configdir=fullfile(filesep,'home','schwerdt','matlab','analysis','analysis','config',filesep);
end

configdir='C:\Users\putamen\Documents\MATLAB\fscv\analysis\config\';%CHANGED

d=dir(configdir);
filenames={d.name};
targfiles='';
processfiles=[];
targconfigname=[];
if strcmp(subjectname,'patra')
targfiles=strfind(filenames,['chronic' sessid 'chconfigsimple.m']);
processfiles=find(~cellfun(@isempty,targfiles));
targconfigname=filenames{processfiles};
run([configdir targconfigname]);            %get paths{1} & paths{2}
run([configdir 'patra_map_bipolar']);   %run patra_map_bipolar for ch settings
else
    targfiles=strfind(filenames,['cleo_chronic' sessid 'chconfigsimple.m']);    
    processfiles=find(~cellfun(@isempty,targfiles));
    if isempty(processfiles)
        %get chconfig file instead of chconfigsimple (for cleo)
        targfiles=strfind(filenames,['cleo_chronic' sessid '.m']);    
        processfiles=find(~cellfun(@isempty,targfiles));
    end    
    if ~isempty(processfiles)
    targconfigname=filenames{processfiles};
    run([configdir targconfigname]);            %get paths{1} & paths{2}
    run([configdir 'cleo_map_bipolar' sessid]);   %run patra_map_bipolar for ch settings
    else
        %DA only can be processed, no config file, so no LFP recorded
    plotlfp=0;
    plotgam=0;
    plotavggam=0;
    plotavglfp=0;
    plotphys=0;
    plotfft=0;
    plotcor=0;
    plotxvar=0;
    end
end
%run config file to get ncschannels & paths
homedirf=strfind(paths{1},'1dr');
subjectdir=paths{1}(1:homedirf-1);
manlist=[subjectdir 'chronic' sessid '_trialselection.xlsx'];      %list for manually selected bad trials
%assignin('base','fcnStatus',targpathtemp)   %store targpath in workspace
manlistlfp=[subjectdir 'chronic' sessid 'l_trialselection.xlsx'];      %list for manually selected bad trials


%OLD PATHS
fscvpathtemp.bigreward=fullfile(paths{1}, 'matlab','bigreward_pro','analyzed',filesep);
fscvpathtemp.smallreward=fullfile(paths{1}, 'matlab','smallreward_pro','analyzed',filesep);
fscvpathtemp.targetbreak=fullfile(paths{1}, 'matlab','targetbreak_pro','analyzed',filesep);
fscvpathtemp.fixbreak=fullfile(paths{1}, 'matlab','fixbreak_pro','analyzed',filesep);

%NEW PATHS in SEPARATE FOLDER
pctype=computer;
ispc=strcmpi(pctype,'pcwin64');
targpathtemp.bigreward=fullfile(homedir,'analysis',subjectname,['chronic' sessid],'bigreward',filesep);
targpathtemp.smallreward=fullfile(homedir,'analysis',subjectname,['chronic' sessid],'smallreward',filesep);
targpathtemp.targetbreak=fullfile(homedir,'analysis',subjectname,['chronic' sessid],'targetbreak',filesep);
targpathtemp.fixbreak=fullfile(homedir,'analysis',subjectname,['chronic' sessid],'fixbreak',filesep);



targpath={};        %analysis path
fscvpath={};        %original files
if plotbig
    targpath.bigreward=targpathtemp.bigreward;
    fscvpath.bigreward=fscvpathtemp.bigreward;
end
if plotsmall
    targpath.smallreward=targpathtemp.smallreward;
    fscvpath.smallreward=fscvpathtemp.smallreward;
end
if plottargbreak
    targpath.targetbreak=targpathtemp.targetbreak;
    fscvpath.targetbreak=fscvpathtemp.targetbreak;
end
if plotfixbreak
    targpath.fixbreak=targpathtemp.fixbreak;
        fscvpath.fixbreak=fscvpathtemp.fixbreak;
end
trialgroups=fieldnames(targpath);
numtypes=length(trialgroups);

time=[22.5 37.5]; 
if sessnum<66
    time=[25 35]; 
end
if sessnum>=178
    time=[20 35];
end

time=[25 35];                   %MAKE SHORT 01/29/2018 ALLL
param={};
csc_map={};
eventcodes={};
if strcmp(subjectname,'patra')
[param,csc_map,eventcodes]=getparams('patrabipolar','default',ncschannels,'sessnum',sessnum);%ncschannels come from chronicXXchconfig file
else
[param,csc_map,eventcodes]=getparams(['cleo' sessid],'cleo',ncschannels,'sessnum',sessnum);
end
getplotsettings([],ncschannels,[],1000,[]);       
plotparam=plotParam;
lfpchs=ncschannels(plotParam.lfpid);
otherchs=ncschannels(plotParam.physid);

plotparam.chnums=fscvchs;
if isempty(fbands)
    %use default from plotparam, getplotsettings
    fbands=plotparam.fbands;
else
    plotparam.fbands=fbands;
end

fbands1=fbands;

%setup rates
ratelfp=plotparam.ratelfp;
rates=[ceil(plotparam.ratelfp) param.samplerate];
plotparam.samplespersec=param.samplerate;

%set up fscv temporal domain
plotparam.pathname=paths{1};
win=time(1)*param.samplerate:1:time(2)*param.samplerate;
plotparam.win=win;
plotparam.interval=2.5;
plotparam.baseline='pretarg';           %immediate post-cue/fix appear
if ~isempty(baseline)
    %user provided
    plotparam.baseline=baseline;
end
plotparam.alignidx=300;         %align to 30s Reward or outcome period

%get labels for fscv chs from original file names
if ~lfpflag
origpath=paths{1}(1:strfind(paths{1},'cvtotxt')-1);
filesorig=dir(origpath);
filenamesorig={filesorig.name};
targfiles=strfind(filenamesorig,'1dr_');
labeledfiles=find(~cellfun(@isempty,targfiles));
if isempty(labeledfiles)
    
targfiles=strfind(filenamesorig,'1dr2_');
labeledfiles=find(~cellfun(@isempty,targfiles));
end
if isempty(labeledfiles)
    
targfiles=strfind(filenamesorig,'1dr3_');
labeledfiles=find(~cellfun(@isempty,targfiles));
end
filesorig=filenamesorig(labeledfiles);
sitenamessp=strsplit(filesorig{1},'_');
sites=sitenamessp(2:5);
plotparam.sites=sites;
end
plotparam.lfpchs=lfpchs;
%[parameters,csc_map,event_codes]=getparams('patrabipolar',[],[],'sessnum',sessnum);
load([configdir 'colormapparula']);
parula=pmap;        %for chunky where parula color map not available;
plotparam.colormap=parula;
%plotparam.markercolors=[1 1 1; 0 1 1; 1 1 1; 1 0 1; 1 1 0]; %reward color yel
plotparam.markercolors=[1 1 1; 0 1 1; 1 1 1; 1 0 1; 1 0 0];
plotparam.cminshow=0;
plotparam.alnevt=alnevt;
numfbands=1;
allbands={fbands};
if strcmp(subjectname,'cleo')
    fbands3={};             %cleo set default [13 28] for all channels ignore below fbands
for ii=1:length(plotParam.lfpid)
    idxf=ii*2-1;
    fbands3{idxf}=ncschannels{ii};
    fbands3{idxf+1}=[13 28];
end
end
%fbands from opening chbands.m script
if ~isempty(fbands2)
    %2 groups of frequency bands for lfp's to analyze below
    numfbands=2;
    allbands={fbands fbands2};
end
if ~isempty(fbands3)
    %2 groups of frequency bands for lfp's to analyze below
    numfbands=3;
    allbands={fbands fbands2 fbands3};
end
if globeta
    numfbands=1;
    allbands{3}=fbands3;        %get from chbands.m script
end

if ~allsess           
for itype=1:numtypes
    plotparam.baseline='pretarg';  
    %plotparam.alnevt=alnevt;            %swithc =plotparam.alnevt depdnign on targevent
    plotparam.shift=[];
    analysispath=getfield(targpath,trialgroups{itype});
    datapath=getfield(fscvpath,trialgroups{itype});
    disp(analysispath)
    plotparam.pathname=analysispath;
    plotparam.trialtype=trialgroups{itype};
    if ~isdir(analysispath)
        %files/folder not created by sync sigs, skip
        disp(['folder not created, skipping : ' char(10) analysispath]);
        continue
    end
    if strcmp(trialgroups{itype},'targetbreak')
        %shift time to targ
        %plotparam.shift='targ'; 
        %plotparam.baseline='pretarg';  
        plotparam.alnevt='outcome';     %look at outcome period
    end
    if strcmp(trialgroups{itype},'fixbreak')
        %shift time to targ
        %plotparam.shift='fix'; 
        %plotparam.baseline='precue'; 
        plotparam.alnevt='outcome';
        
    end    
    
    lfpexportdir=0;
    lfpexported=0;
    if ~replace
        %only if user does not indicate that data should be replaced
    %check for exported data, so do not need to retrieve raw lfp again
    if ismember(3,selbands) || globeta
        %check if lfp data already exported
        lfpexportdir=[analysispath 'trialbytrial_lfp3' filesep 'all-broadbetalfp.mat'];
        if exist(lfpexportdir)==2
            %file exists
            lfpexported=1;
        end
    end
    if ismember(2,selbands)
        %check if lfp data already exported
        lfpexportdir=[analysispath 'trialbytrial_lfp2' filesep 'all-betalfp2.mat'];
        if exist(lfpexportdir)==2
            %file exists
            lfpexported=1;
        else
            lfpexported=0;
        end
    end
    if ismember(1,selbands)
        %check if lfp data already exported
        lfpexportdir=[analysispath 'trialbytrial_lfp' filesep 'all-betalfp.mat'];
        if exist(lfpexportdir)==2
            %file exists
            lfpexported=1;
        else
            lfpexported=0;
        end
    end
    end
    gamexported=0;
    gamexportdir=[analysispath 'trialbytrial_gam' filesep 'all-gammalfp.mat'];
    if exist(gamexportdir)==2
        %file exists
        gamexported=1;
    else
        gamexported=0;
    end
    %open compiled data
    compiledname=[analysispath 'trials_chronic' sessid '_' trialgroups{itype}];
    fscvchstemp=fscvchs;
    load(compiledname); %reloads fscvchs
    fscvchs=fscvchstemp;
    

    
    ich=targch;      %sort signals by this channel
    numtrials=[];
    if ~lfpflag
    numtrials=size(trialbytrial(ich).da,1);
    else
        ich=1;
        numtrials=length(trialbytrial(ich).alignts);
    end
    chbadtrials=[];
    %goodtrials=1:numtrials;         %default, all trials goodtrials
    alltrials=1:numtrials;
    goodtrials={};
    for iich=1:4
    goodtrials{iich}=1:numtrials; %default, all trials goodtrials 5/2019, ch-specific
    end
    %open bad trials if exists in homedir
    if exist(manlist)>0
        goodtrials=manualtrialsel(manlist,targch,alltrials,trialgroups{itype});          %get recent list of bad trials manually selected in excel sheet
    end   
    goodtrialslfp={};
    if exist(manlistlfp)>0
        goodtrialslfp=manualtrialsel(manlistlfp,1,alltrials,trialgroups{itype});          %get recent list of bad trials manually selected in excel sheet
        goodtrialslfp=goodtrialslfp(1);
    end     
    if any(samplesfixeye{1}<100) || any(samplesfixeye{1}<samplesfix{1})
        chbadtrials2=find(samplesfixeye{1}<100);
        chbadtrials2=unique([chbadtrials2 find(samplesfixeye{1}<samplesfix{1})]);
        for iich=1:4
            goodtrials{iich}=goodtrials{iich}(find(~ismember(goodtrials{iich},chbadtrials2)));
        end
    end
    if any(samplestarg{1}<samplesfixeye{1}) || any(samplestarg{1}<100)
        chbadtrials2=find(samplestarg{1}<samplesfixeye{1});
        chbadtrials2=unique([chbadtrials2 find(samplestarg{1}<100)]);
        for iich=1:4
        goodtrials{iich}=goodtrials{iich}(find(~ismember(goodtrials{iich},chbadtrials2)));
        end
    end
    if ~isfield(trialbytrial,'lfp') && ((plotlfp  && ~lfpexported) || plotfft || (plotgam && ~gamexported) || plotphys)
        %not merged file, get lfps separately based on trial num 100 etc.
        %from original xx_pro directory
        data={};
        %get reconverted path ie XX_pro directory
        pathid=strfind(datapath,'_pro');
        pathlfp=datapath(1:pathid+4);
        dirlfp=dir(pathlfp);
        %get csc signals defined in NCSchannels
        if isempty(param.NCSchannels)
            error('no ncs channels defined in param');
        end
        for it=1:numtrials
            disp(['trial # ' num2str(it)]);
            tlab=it+100-1;
            for idx=1:length(param.NCSchannels)
                %load in order of appearing in NCSchannels list
                matchcsc=regexpi({dirlfp.name},['csc' num2str(param.NCSchannels(idx)) '_']);
                idscsc=~cellfun('isempty',matchcsc);       %make array of zeros/1 for matching
                matchtrial=regexpi({dirlfp.name},['_' num2str(tlab)]);
                idstrial=~cellfun('isempty',matchtrial);  
                fileid=find((and(idscsc,idstrial))==1);     %get file id that matches both target trial/cscnum
                tempload=load([pathlfp,dirlfp(fileid).name]);
                if it>1 && length(tempload.dg_Nlx2Mat_Samples)~=size(data.lfp{idx},2)
                    data.lfp{idx}(it,:)=repmat(nan,1,size(data.lfp{idx},2));
                else
                data.lfp{idx}(it,:)=tempload.dg_Nlx2Mat_Samples;  %store samples, 
                end
                %timestamps in tempload are in microseconds and should match tsLFP already loaded 
                %get/store csc name
                cscnameb=strfind(dirlfp(fileid).name,'_');
                cscname=[dirlfp(fileid).name(1:cscnameb-1) param.LFPterm];                
                data.cscNames{idx}=cscname;
            end          
        end
        trialbytrial(1).lfp=data.lfp;
        %trialbytrial(1).cscNames=data.cscNames;
        trialbytrial(1).cscNames=plotparam.cscNames;
    end
    
    trialbytrial(1).cscNames=plotparam.cscNames;
    %group trial types
    trialbytrial(1).samplesfixeye=samplesfixeye;
    trialbytrial(1).samplesfix=samplesfix;
    trialbytrial(1).samplestarg=samplestarg;
    trialbytrial(1).samplestargeye=samplestargeye;
   % trialtypes=grouptrials(goodtrials,trialbytrial(1),plotparam);
        trialtypes=grouptrials(goodtrials,trialbytrial(1),plotparam,'prev');        %since rerunning this after getting good trials that do not have sw/pos-sw, these types will be empty
        trialtypeslfp={};
        if ~isempty(goodtrialslfp)
        trialtypeslfp=grouptrials(goodtrialslfp,trialbytrial(1),plotparam,'prev');        %since rerunning this after getting good trials that do not have sw/pos-sw, these types will be empty
        trialtypeslfp(1).site='';
        end
    save([analysispath 'trialtypes'], 'trialtypes','trialtypeslfp');
    
    
    %paths to save exported signals from settrialax
    pathda=[analysispath 'trialbytrial_da' filesep];
    if ~isdir(pathda)
        mkdir(pathda)
    end
    pathlfp={};
    pathlfp{1}=[analysispath 'trialbytrial_lfp' filesep];
    pathlfp{2}=[pathlfp{1}(1:end-1) num2str(2) filesep];
    pathlfp{3}=[pathlfp{1}(1:end-1) num2str(3) filesep];
    pathgam=[analysispath 'trialbytrial_gam' filesep];
    if ~isdir(pathlfp{1})
        mkdir(pathlfp{1})
            mkdir(pathlfp{2});
    mkdir(pathlfp{3});
    mkdir(pathgam);
    end


    plotparam.samplesfixeye=samplesfixeye;
    plotparam.samplesfix=samplesfix;
    plotparam.samplestarg=samplestarg;
    plotparam.samplestargeye=samplestargeye;
    plotparam.tfixeye=trialbytrial(1).tfixeye;
    plotparam.tfix=trialbytrial(1).tfix;
    plotparam.ttarg=trialbytrial(1).ttarg;
    plotparam.ttargeye=trialbytrial(1).ttargeye;
    
    plottrials=unique([goodtrials{:}]);
    numtrials=size(plottrials,2);
    plotparam.cscale=[nan 6];    %ch4 chronic 30
    plotparam.vertplotwidth=100;
    plotparam.vertplotwidth2=30;
    plotparam.numtrials=numtrials;
    numcolor=length(fscvchs)+2;
    numbehav=8;     %% behavior plots
    plotparam.figpos=[200 200 1800 900];
    if ~ispc
        plotparam.figpos=[0 0 1000 750];
    end
    figpos=plotparam.figpos;
    plotparam.colorsize=[figpos(3)/(numcolor+2) figpos(4)*3/4];
    plotparam.margins=20;
    %initialize storage variables data
    da={};
    datm={};
    betatm={};      %task modulated signal averages
    gammatm={};
    gammalfp={};
    betalfp={};
    betalfp2={};        %high beta band
    betalfp3={};        %broad beta
    eyedata=[];
    eyexdata=[];
    pulsedata=[];
    lickdata=[];
    eyedistdata=[];
    eyevdata=[];
    blinkdata=[];
    plotparam.chnum=targch;
    eyextm=[];
    
if getinfo
    %just save trialtypes and continue loop
    continue;
end

for ialn=1:length(targevents)
    %cycle different alignment events, fix, targ, targeye, 1 /2019
    switch targevents{ialn}
        case 'fix'
            plotparam.alnevt='fix';
        case 'fixeye'
            plotparam.alnevt='fixeye';
        case 'targ'
            plotparam.alnevt='targ';
        case 'targeye'
            plotparam.alnevt='targeye';
        case 'outcome'
            plotparam.alnevt='outcome';
    end
%%
%plot DA trial by trial side by side relevant lfp channels
%& side by side eyex, eyed, lick

if plotda
    disp(['plotting da chs ' num2str(fscvchs)]);
fig1=figure('visible','off');
if ispc
    fig1=figure('visible','on');
end
if ~isfield(scales,'da')
    for ii=1:length(fscvchs)
    ich=fscvchs(ii);    
    scales.da{ich}=[];
    end
    scales.eye=[];
    scales.lick=[];
    scales.blink=[];
    scales.eyex=[];
    scales.eyev=[];
    scales.eyedist=[];
    if exist([analysispath, 'scales.mat'])==2
        load([analysispath 'scales'], 'scales');
    end
end
    pulsetm=[];
eyetm=[];
eyextm=[];
eyedisttm=[];
licktm=[];
blinktm=[];
eyevtm=[];
hfig=setupTrialPlots(fig1, plotparam,numcolor,numbehav);
pathdaaln=[analysispath 'trialbytrial_da_' targevents{ialn} filesep];
if ~isdir(pathdaaln)
    mkdir(pathdaaln)
end
for it=1:length(trialtypes(1).names)

if ~isempty(fscvchs)
for ii=1:length(fscvchs)
    ich=fscvchs(ii);    
    plotparam.chnum=ich;
        plottrials=trialtypes(ich).nums{it};
    trialname=trialtypes(ich).names{it};
    plotparam.triallabel=trialname;   
    if ~isempty(plottrials)    
    if ii==1
        [datm{ich},da{ich},info]=setTrialAx(trialbytrial(ich),plotparam,hfig{ii},plottrials,...
            'plotnum',1,'sitename',sites(ich),'scales',scales.da{ich});        
    else
        [datm{ich},da{ich},info]=setTrialAx(trialbytrial(ich),plotparam,hfig{ii},plottrials,...
            'sitename',sites(ich),'scales',scales.da{ich});
    end
    if it==1 && ialn==1       
        if itype==1
            scales.da{ich}=info.scales;
        end
        %save trial by trial data with reference ts's for nlx to sort
        %with small rew trials converted at another time, usually outcome
        %ts
        if ii==length(fscvchs)
            alignts=[];
            if isfield(trialbytrial(1),'alignts')
            alignts=trialbytrial(1).alignts;
            end
            trialids=plottrials;
            save([pathda trialname '-data'],'da',...
                'trialids','alignts','trialtypes','info');
            %when reusing data, for specific types of trials 
            %need to index by find(ismember(trialids,phase2trials)==1)
        end
    end
    end
end
end
%for behavior get all good trials across all channels
plottrials=[];
for iip=1:length(trialtypes)
plottrials=sort(unique([plottrials [trialtypes(iip).nums{it}]]));
end
    trialname=trialtypes(1).names{it};
    plotparam.triallabel=trialname;   
rts=setTrialBehavAx(trialbytrial(1),plotparam,hfig,plottrials,length(fscvchs)+3);
if ~isempty(plotparam.cscNames)
if lfpexported 
    if isempty(eyedata)
    %load already saved rclfp data
    disp(['loading: ' analysispath 'trialbytrial_da' filesep 'all-data-beh']);
    load([analysispath 'trialbytrial_da' filesep 'all-data-beh'])
    end
    eyetm=setTrialAx(trialbytrial(1),plotparam,hfig{length(fscvchs)+1},...
        plottrials,'eyed','scales',scales.eye,'rclfp',eyedata,'info',eyeinfo);
    licktm=[];
    if strcmp(subjectname,'patra') || (strcmp(subjectname,'cleo') && sessnum>24)
    licktm=setTrialAx(trialbytrial(1),plotparam,hfig{length(fscvchs)+2},...
        plottrials, 'lickx','scales',scales.lick,'rclfp',lickdata,'info',lickinfo);
    end
   % rts=setTrialBehavAx(trialbytrial(1),plotparam,hfig,plottrials,length(fscvchs)+3);
        blinktm=setTrialAx(trialbytrial(1),plotparam,[],...
                plottrials, 'blink','noplot','rclfp',blinkdata,'info',blinkinfo);
            if ~isempty(eyexdata)
    eyextm=setTrialAx(trialbytrial(1),plotparam,[],...
        plottrials, 'eyex','noplot','rclfp',eyexdata,'info',eyexinfo);
            end
    eyevtm=setTrialAx(trialbytrial(1),plotparam,[],...
        plottrials, 'eyev','noplot','rclfp',eyevdata,'info',eyevinfo);
    eyedisttm=setTrialAx(trialbytrial(1),plotparam,[],...
        plottrials, 'eyedist','noplot','rclfp',eyedistdata,'info',eyedistinfo);
    %calc pulse rate
    pulsetm=[];
    if sessnum>37
    pulsetm=setTrialAx(trialbytrial(1),plotparam,[],...
        plottrials, 'pulse','noplot','rclfp',pulsedata,'info',pulseinfo);
    end
else
    %data not yet exported, get /save data
    [eyetm,eyedatatemp,eyeinfo]=setTrialAx(trialbytrial(1),plotparam,hfig{length(fscvchs)+1},...
        plottrials,'eyed','scales',scales.eye);
        licktm=[];
        lickdatatemp=[];
        lickinfo=[];
        lickinfo.scales=[];
    if strcmp(subjectname,'patra') || (strcmp(subjectname,'cleo') && sessnum>24)
    [licktm,lickdatatemp,lickinfo]=setTrialAx(trialbytrial(1),plotparam,hfig{length(fscvchs)+2},...
        plottrials, 'lickx','scales',scales.lick);
    end
    [blinktm,blinkdatatemp,blinkinfo]=setTrialAx(trialbytrial(1),plotparam,[],...
        plottrials, 'blink','noplot');
    [eyextm,eyexdatatemp,eyexinfo]=setTrialAx(trialbytrial(1),plotparam,[],...
        plottrials, 'eyex','noplot');
    [eyevtm,eyevdatatemp,eyevinfo]=setTrialAx(trialbytrial(1),plotparam,[],...
        plottrials, 'eyev','noplot');
    [eyedisttm,eyedistdatatemp,eyedistinfo]=setTrialAx(trialbytrial(1),plotparam,[],...
        plottrials, 'eyedist','noplot');
    if itype==1 && it==1 && ialn==1
        scales.eye=eyeinfo.scales;
        scales.lick=lickinfo.scales;
        scales.blink=blinkinfo.scales;
        scales.eyex=eyexinfo.scales;
        scales.eyev=eyevinfo.scales;
        scales.eyedist=eyedistinfo.scales;
    end
    %calc pulse rate
    pulsetm=[];
    pulsedatatemp=[];
    pulseinfo=[];
    if sessnum>37
    [pulsetm,pulsedatatemp,pulseinfo]=setTrialAx(trialbytrial(1),plotparam,[],...
        plottrials, 'pulse','noplot');
    end
    if it==1 && ialn==1
        %save when plotting all trials for certain itype trial type
        %these data signals are same for all align evts since unshifted data
        pulsedata=pulsedatatemp;
        eyedistdata=eyedistdatatemp;
        eyevdata=eyevdatatemp;
        eyexdata=eyexdatatemp;
        blinkdata=blinkdatatemp;
        lickdata=lickdatatemp;
        eyedata=eyedatatemp;
        %save data in data path (ie nonshifted data)
        save([pathda trialname '-data-beh'],'eyedata','eyexdata','lickdata',...
        'blinkdata','pulsedata','eyevdata','eyedistdata','eyeinfo',...
        'lickinfo','blinkinfo','eyexinfo','eyevinfo','eyedistinfo',...
        'pulseinfo');
    end
end
end
savefig(fig1,[pathdaaln trialname '.fig']);
saveas(fig1,[pathdaaln trialname '.jpg'],'jpg')
%task modulated signal should be independent of aln evt, save once in data
%path
if ialn==1
save([pathda trialname],'datm','eyetm','licktm',...
    'blinktm','pulsetm','eyextm','eyevtm','eyedisttm','rts');
end
end
close all;
end
%%
%plot lfp's trial by trial
plotparam.vertplotwidth=90;
plotparam.vertplotwidth2=10;
if ~ispc
  plotparam.vertplotwidth=50;
plotparam.vertplotwidth2=10;
end
plotparam.margins=10;
plotparam.colorsize=[figpos(3)/(length(lfpchs)+1) figpos(4)*3/4];
if plotlfp

fig2=figure('visible','off');
if ispc
    fig2=figure('visible','on');
end
hfig2=setupTrialPlots(fig2, plotparam,length(lfpchs),0);
for ifreq=1:length(selbands)
    fbands=allbands{selbands(ifreq)};
    plotparam.fbands=fbands;
    pathlfpaln=[];
    pathlfpdata=[];
    if selbands(ifreq)==1
        if ~isfield(scales,'lfp')
            for ii=1:length(lfpchs)
                scales.lfp{ii}=[];
            end
        end
        pathlfpaln=[analysispath 'trialbytrial_lfp_' targevents{ialn} filesep];
        if ~isdir(pathlfpaln)
            mkdir(pathlfpaln)
        end
    end
    if selbands(ifreq)==2
        if ~isfield(scales,'lfp2')
            for ii=1:length(lfpchs)
                scales.lfp2{ii}=[];
            end
        end
        pathlfpaln=[analysispath 'trialbytrial_lfp2_' targevents{ialn} filesep];
        if ~isdir(pathlfpaln)
            mkdir(pathlfpaln)
        end
    end
    if selbands(ifreq)==3 || globeta
        if ~isfield(scales,'lfp3')
            for ii=1:length(lfpchs)
                scales.lfp3{ii}=[];
            end
        end
        pathlfpaln=[analysispath 'trialbytrial_lfp3_' targevents{ialn} filesep];
        if ~isdir(pathlfpaln)
            mkdir(pathlfpaln)
        end
        if lfpexported
            %load data
            disp(['loading: ' pathlfp{3} 'all-broadbetalfp'])
            load([pathlfp{3} 'all-broadbetalfp'],'betalfp3','info');
            if ~isfield(scales,'lfp')
                load([analysispath 'scales'], 'scales');
            end
        end
    end
    for it=1:length(trialtypes(1).names)
        %plottrials=trialtypes.nums{it};
        plottrials=[];
        for iip=1:length(trialtypes)
        plottrials=sort(unique([plottrials [trialtypes(iip).nums{it}]]));
        end
        if ~isempty(goodtrialslfp)
            plottrials=sort(unique(trialtypeslfp(1).nums{it}));
        end
        if ~isempty(plottrials)    
        trialname=trialtypes(1).names{it};
        plotparam.triallabel=trialname;
        for ii=1:length(lfpchs)
            if it==1 && ialn==1 && ~lfpexported
            [betatm{ii},downlfp,info]=setTrialAx(trialbytrial(1),...
                plotparam,hfig2{ii},plottrials,lfpchs{ii},'fbands',fbands,'getfft',...
                'tapers',[3 5],'smoothwin',.25,'plotnum',ii);
                if itype==1
                    if selbands(ifreq)==1
                        scales.lfp{ii}=info.scales;
                    end
                    if selbands(ifreq)==2
                        scales.lfp{ii}=info.scales;
                    end
                    if selbands(ifreq)==3 || globeta
                        scales.lfp{ii}=info.scales;
                    end
                end
                if selbands(ifreq)==2
                    betalfp2{ii}=downlfp;
                end
                if selbands(ifreq)==1
                    betalfp{ii}=downlfp;
                end
                if selbands(ifreq)==3 || globeta
                    betalfp3{ii}=downlfp;
                end
                if ii==length(lfpchs) 
                    %save betalfp data
                    alignts=[];
                    if isfield(trialbytrial(1),'alignts')
                        alignts=trialbytrial(1).alignts;
                    end
                    trialids=plottrials;
                    if selbands(ifreq)==1
                    save([pathlfp{1} trialname '-betalfp'],'betalfp',...
                        'trialids','alignts','fbands','trialtypes','info');
                    end
                    if selbands(ifreq)==2
                    save([pathlfp{2} trialname '-betalfp2'],'betalfp2',...
                        'trialids','alignts','fbands','trialtypes','info');
                    end
                    if selbands(ifreq)==3 || globeta
                    save([pathlfp{3} trialname '-broadbetalfp'],'betalfp3',...
                        'trialids','alignts','fbands','trialtypes','info','-v7.3');
                    end                   
                end

            else
            rclfp=[];
            if selbands(ifreq)==1
                rclfp=betalfp{ii};
            end
            if selbands(ifreq)==2
                rclfp=betalfp2{ii};
            end
            if selbands(ifreq)==3 || globeta
                rclfp=betalfp3{ii};
            end               
            %already got all betalfp data in first 'all trial' it
            betatm{ii}=setTrialAx(trialbytrial(1),...
                plotparam,hfig2{ii},plottrials,lfpchs{ii},'rclfp',rclfp,'fbands',fbands,...
                'tapers',[3 5],'smoothwin',0.25,'plotnum',ii,'scales',scales.lfp{ii},'info',info);
            end
        end
        end
        if globeta==0
            %savefig(fig2,[pathlfp{selbands(ifreq)} trialname]);
            %saveas(fig2,[pathlfp{selbands(ifreq)} trialname '.eps'],'epsc')
            saveas(fig2,[pathlfpaln trialname '.jpg'],'jpg')
            if ialn==1
            save([pathlfp{selbands(ifreq)} trialname],'betatm');
            end
        else
            %broad beta saved in selbands(ifreq) 3
            savefig(fig2,[pathlfpaln trialname]);
            %saveas(fig2,[pathlfp{3} trialname '.eps'],'epsc')
            saveas(fig2,[pathlfpaln trialname '.jpg'],'jpg')
            if ialn==1
            save([pathlfp{3} trialname],'betatm');
            end
        end
   
    end
end
close all;
end
if isfield(scales,'lfp') && isfield(scales,'da')
    %got both lfp & da scales save them
    save([analysispath 'scales'], 'scales');
end
%plot gamma
%%
if plotgam
fig2=figure('visible','off');
if ispc
    fig2=figure('visible','on');
end
hfig2=setupTrialPlots(fig2, plotparam,length(lfpchs),0);
pathgamaln=[];
if ~isfield(scales,'gam')
    for ii=1:length(lfpchs)
        scales.gam{ii}=[];
    end
end
pathgamaln=[analysispath 'trialbytrial_gam_' targevents{ialn} filesep];
if ~isdir(pathgamaln)
    mkdir(pathgamaln)
end
if gamexported
    %load data, previously exported
    disp(['loading: ' pathgam 'all-gammalfp'])
    load([pathgam 'all-gammalfp'],'gammalfp','info');
    if ~isfield(scales,'gam')
        load([analysispath 'scales'], 'scales');
    end
end
if ~isdir(pathgam)
    mkdir(pathgam)
end
for it=1:length(trialtypes(1).names)
   % plottrials=trialtypes.nums{it};
    plottrials=[];
    for iip=1:length(trialtypes)
    plottrials=sort(unique([plottrials [trialtypes(iip).nums{it}]]));
    end
    trialname=trialtypes(1).names{it};
    plotparam.triallabel=trialname;
    for ii=1:length(lfpchs)
        if it==1 && ialn==1 && ~gamexported
        [gammatm{ii},downlfp,info]=setTrialAx(trialbytrial(1),...
            plotparam,hfig2{ii},plottrials,lfpchs{ii},'bandfilt',filtgamma,'getfft',...
            'tapers',[3 5],'smoothwin',.05,'plotnum',ii);
            if itype==1                
                scales.gam{ii}=info.scales;                
            end
            gammalfp{ii}=downlfp;            
            if ii==length(lfpchs) 
                %save betalfp data
                alignts=[];
                if isfield(trialbytrial(1),'alignts')
                    alignts=trialbytrial(1).alignts;
                end
                trialids=plottrials;
                save([pathgam trialname '-gammalfp'],'gammalfp',...
                    'trialids','alignts','filtgamma','trialtypes','info');                           
            end
        else
            %already gamma exported & reloaded
            rclfp=[];
            rclfp=gammalfp{ii};         
            %already got all betalfp data in first 'all trial' it
            gammatm{ii}=setTrialAx(trialbytrial(1),...
                plotparam,hfig2{ii},plottrials,lfpchs{ii},'rclfp',rclfp,'bandfilt',filtgamma,...
                'tapers',[3 5],'smoothwin',0.05,'plotnum',ii,'scales',scales.gam{ii},'info',info);
        end
    end
    savefig(fig2,[pathgamaln trialname]);
    saveas(fig2,[pathgamaln trialname '.jpg'],'jpg')
    if ialn==1
        save([pathgam trialname],'gammatm');
    end
end
save([analysispath 'scales-gamma'], 'scales');
close all;
end
%%
%plot averages da
if plotavgda
figsigrep=figure('visible','off'); 
if ispc
    figsigrep=figure('visible','on');
end
set(figsigrep, 'Color', [1 1 1],'position',[100 50 1000 950]);
set(0,'CurrentFigure',figsigrep); 
axsigrep={};
cscaleda=[-2 4];
cscaleda=[-4 4];
for ii=1:8
    axsigrep{ii}=subplot(4,2,ii);
end
pathavg=[analysispath 'avg_da_' targevents{ialn} filesep];
if ~isdir(pathavg)
    mkdir(pathavg)
end
if isempty(da)
    %load da data, 'da' & scales
    if exist([pathda 'all-data.mat'])>0
        load([pathda 'all-data'],'da','info');
    end
end
if ~isempty(da)
if ~isfield(scales,'da')
    load([analysispath 'scales'], 'scales');
end
for it=1:length(trialtypes(1).names)
   % plottrials=trialtypes.nums{it};   
    trialname=trialtypes(1).names{it};
    plotparam.triallabel=trialname;
    countp=0;
    for ii=1:length(fscvchs)
        countp=countp+1;
        ich=fscvchs(ii);
        plotparam.chnum=ich;
                plottrials=trialtypes(ich).nums{it};
                if ~isempty(plottrials)
        setTrialAxAvg(da{ich},plotparam,axsigrep{countp},plottrials, 'sitename',sites(ich),'plotnum',countp);
        countp=countp+1;
        setTrialAxAvg(da{ich},plotparam,axsigrep{countp},plottrials, 'avg','scales',scales.da{ich},'sitename',sites(ich),'plotnum',countp);
                end
    end
    savefig(figsigrep,[pathavg trialname]);
    saveas(figsigrep,[pathavg trialname '.jpg'],'jpg')
end
end
end
%%
%plot averages phys
if plotavglfp
figsigrep2=figure('visible','off'); 
if ispc
    figsigrep2=figure('visible','on');
end
set(figsigrep2, 'Color', [1 1 1],'position',[100 50 1500 900]);
set(0,'CurrentFigure',figsigrep2); 
axsigrep2={};
if length(lfpchs)<=5
    for ii=1:length(lfpchs)*2
        axsigrep2{ii}=subplot(length(lfpchs),2,ii);
    end
else
    set(figsigrep2, 'Color', [1 1 1],'position',[100 50 1500 1000]);
     for ii=1:length(lfpchs)*2
        axsigrep2{ii}=subplot(ceil(length(lfpchs)/2),4,ii);
     end
end
pathavglfp={};
pathavglfp{1}=[analysispath 'avg_lfp_' targevents{ialn} filesep];
pathavglfp{2}=[analysispath 'avg_lfp2_' targevents{ialn}  filesep];
pathavglfp{3}=[analysispath 'avg_lfp3_' targevents{ialn} filesep];
if ~isfield(scales,'lfp')
    load([analysispath 'scales'], 'scales');
end
if ismember(1,selbands)
if ~isdir(pathavglfp{1})
    mkdir(pathavglfp{1})
end
end
if ismember(2,selbands)
if ~isdir(pathavglfp{2})
    mkdir(pathavglfp{2})
end
end
if ismember(3,selbands) || globeta
if ~isdir(pathavglfp{3})
    mkdir(pathavglfp{3})
end
end
   
for ifreq=1:length(selbands)
    fbands=allbands{selbands(ifreq)};
    plotparam.fbands=fbands;
    if ismember(1,selbands) && ~globeta
        load([analysispath 'trialbytrial_lfp' filesep 'all-betalfp'],'info')
        %load infos wi sample rate etc
        %need to load
        if isempty(betalfp)
        load([analysispath 'trialbytrial_lfp' filesep 'all-betalfp'],'betalfp','info')
        end
    end
    if ismember(2,selbands) && ~globeta
        load([analysispath 'trialbytrial_lfp2' filesep 'all-betalfp2'],'info')
         if isempty(betalfp2) 
            load([analysispath 'trialbytrial_lfp2' filesep 'all-betalfp2'],'betalfp2','info')
         end
    end
    if ismember(3,selbands) || globeta
        load([analysispath 'trialbytrial_lfp3' filesep 'all-broadbetalfp'],'info')
        if isempty(betalfp3)
        load([analysispath 'trialbytrial_lfp3' filesep 'all-broadbetalfp'],'betalfp3','info')
        end
    end
    for it=1:length(trialtypes(1).names)
        %plottrials=trialtypes.nums{it};
        trialname=trialtypes(1).names{it};
        plotparam.triallabel=trialname;
            plottrials=[];
    for iip=1:length(trialtypes)
    plottrials=sort(unique([plottrials [trialtypes(iip).nums{it}]]));
    end
    if ~isempty(goodtrialslfp)
        plottrials=sort(unique(trialtypeslfp(1).nums{it}));
    end
    if ~isempty(plottrials)
        for ii=1:length(lfpchs)
            if selbands(ifreq)==1
                downlfp=betalfp{ii};
            end
            if selbands(ifreq)==2
                downlfp=betalfp2{ii};
            end                
            if selbands(ifreq)==3 || globeta
                downlfp=betalfp3{ii};
            end
            setTrialAxAvg(downlfp,plotparam,axsigrep2{ii*2-1},...
                plottrials, ...
                lfpchs{ii},'fbands',fbands,'info',info,'scales',scales.lfp{ii});
            setTrialAxAvg(downlfp,plotparam,axsigrep2{ii*2},...
                plottrials, ...
                lfpchs{ii},'avg','fbands',fbands,'info',info,'scales',scales.lfp{ii});
        end
        if globeta==0
           % savefig(figsigrep2,[pathavglfp{selbands(ifreq)} trialname]);
            saveas(figsigrep2,[pathavglfp{selbands(ifreq)} trialname '.jpg'],'jpg')
        else
            savefig(figsigrep2,[pathavglfp{3} trialname]);
            saveas(figsigrep2,[pathavglfp{3} trialname '.jpg'],'jpg')
        end
    end
    end
end
end
%%
%plot averages gamma
if plotavggam
figsigrep2=figure('visible','off'); 
if ispc
    figsigrep2=figure('visible','on');
end
set(figsigrep2, 'Color', [1 1 1],'position',[100 50 1500 900]);
set(0,'CurrentFigure',figsigrep2); 
axsigrep2={};
if length(lfpchs)<=5
    for ii=1:length(lfpchs)*2
        axsigrep2{ii}=subplot(length(lfpchs),2,ii);
    end
else
    set(figsigrep2, 'Color', [1 1 1],'position',[100 50 1500 1000]);
     for ii=1:length(lfpchs)*2
        axsigrep2{ii}=subplot(ceil(length(lfpchs)/2),4,ii);
     end
end
pathavggam=[analysispath 'avg_gam_' targevents{ialn} filesep];
scalestemp=scales;
if ~isfield(scales,'gam')
    load([analysispath 'scales-gamma'], 'scales');
    scalestemp.gam=scales.gam;
    scales=scalestemp;
end
if ~isdir(pathavggam)
mkdir(pathavggam);
end
%need to load
load([analysispath 'trialbytrial_gam' filesep 'all-gammalfp'],'info')
if isempty(gammalfp)
    load([analysispath 'trialbytrial_gam' filesep 'all-gammalfp'],'gammalfp','info')
end
for it=1:length(trialtypes(1).names)
        plottrials=[];
    for iip=1:length(trialtypes)
    plottrials=sort(unique([plottrials [trialtypes(iip).nums{it}]]));
    end
%plottrials=trialtypes.nums{it};
trialname=trialtypes(1).names{it};
plotparam.triallabel=trialname;
for ii=1:length(lfpchs)
   downlfp=gammalfp{ii};    
   if ~isfield(scales,'gam')
    setTrialAxAvg(downlfp,plotparam,axsigrep2{ii*2-1},...
        plottrials, ...
        lfpchs{ii},'bandfilt',filtgamma,'info',info);
    setTrialAxAvg(downlfp,plotparam,axsigrep2{ii*2},...
        plottrials, ...
        lfpchs{ii},'avg','bandfilt',filtgamma,'info',info);       
   else       
    setTrialAxAvg(downlfp,plotparam,axsigrep2{ii*2-1},...
        plottrials, ...
        lfpchs{ii},'bandfilt',filtgamma,'info',info,'scales',scales.gam{ii});
    setTrialAxAvg(downlfp,plotparam,axsigrep2{ii*2},...
        plottrials, ...
        lfpchs{ii},'avg','bandfilt',filtgamma,'info',info,'scales',scales.gam{ii});
   end
end
savefig(figsigrep2,[pathavggam trialname]);
saveas(figsigrep2,[pathavggam trialname '.jpg'],'jpg')
end
end

%%
%plot averages physiological response
if plotphys
figsigrep3=figure('visible','off'); 
figsigrep4=figure('visible','off'); 
if ispc
    set(figsigrep3,'visible','on'); 
    set(figsigrep4,'visible','on'); 
end
set(figsigrep3, 'Color', [1 1 1],'position',[100 100 1000 1000]);
set(figsigrep4, 'Color', [1 1 1],'position',[100 100 1000 1000]);
if ~ispc
set(figsigrep3, 'Color', [1 1 1],'position',[100 100 1000 750]);
set(figsigrep4, 'Color', [1 1 1],'position',[100 100 1000 750]);
end
axsigrep3={};
axsigrep4={};
for ii=1:8
    set(0,'CurrentFigure',figsigrep3); 
    axsigrep3{ii}=subplot(4,2,ii);
    set(0,'CurrentFigure',figsigrep4); 
    axsigrep4{ii}=subplot(4,2,ii);
end

pathavgb=[analysispath 'avg_beh_' targevents{ialn} filesep];
pathavgbeye=[analysispath 'avg_beheye_' targevents{ialn} filesep];
if ~isdir(pathavgb)
    mkdir(pathavgb)
    mkdir(pathavgbeye)
end
if ~isfield(scales,'eye')
    load([analysispath 'scales'], 'scales');
end
%load behavioral data acquired above, emptied for every itype run
if isempty(eyedata) && exist([analysispath 'trialbytrial_da' filesep 'all-data-beh.mat'])>0
load([analysispath 'trialbytrial_da' filesep 'all-data-beh'])
end
if ~isempty(eyedata) || ~isempty(lickdata) || ~isempty(pulsedata)
%loads eyedata, lickdata, blinkdata, eyexdata, pulsedata, etc.
for it=1:length(trialtypes(1).names)
    %plottrials=trialtypes.nums{it};
        plottrials=[];
    for iip=1:length(trialtypes)
    plottrials=sort(unique([plottrials [trialtypes(iip).nums{it}]]));
    end
    trialname=trialtypes(1).names{it};
    plotparam.triallabel=trialname;
if ~isempty(plottrials)
    if ~isempty(lickdata)
        %if strcmp(subjectname,'patra') || (strcmp(subjectname,'cleo') && sessnum>24)
    setTrialAxAvg(lickdata,plotparam,axsigrep3{3},plottrials, 'info',lickinfo,'lickx');
    setTrialAxAvg(lickdata,plotparam,axsigrep3{4},plottrials, 'info',lickinfo,'lickx','avg','scales',scales.lick);
        %end
    end
    if ~isempty(pulsedata)
    %if sessnum>37
    setTrialAxAvg(pulsedata,plotparam,axsigrep3{5},plottrials, 'info',pulseinfo,'pulse');
    setTrialAxAvg(pulsedata,plotparam,axsigrep3{6},plottrials, 'info',pulseinfo,'pulse','avg');
    %end
    end
        if ~isempty(eyedata)
    setTrialAxAvg(eyedata,plotparam,axsigrep3{1},plottrials, 'info',eyeinfo,'eyed');
    setTrialAxAvg(eyedata,plotparam,axsigrep3{2},plottrials, 'info',eyeinfo,'eyed','avg','scales',scales.eye);
    setTrialAxAvg(blinkdata,plotparam,axsigrep3{7},plottrials, 'info',blinkinfo,'blink');
    setTrialAxAvg(blinkdata,plotparam,axsigrep3{8},plottrials, 'info',blinkinfo,'blink','avg','scales',scales.blink);
    setTrialAxAvg(eyexdata,plotparam,axsigrep4{1},plottrials, 'info',eyexinfo,'eyex');
    setTrialAxAvg(eyexdata,plotparam,axsigrep4{2},plottrials, 'info',eyexinfo,'eyex','avg','scales',scales.eyex);
    setTrialAxAvg(eyevdata,plotparam,axsigrep4{3},plottrials, 'info',eyevinfo,'eyevel');
    setTrialAxAvg(eyevdata,plotparam,axsigrep4{4},plottrials, 'info',eyevinfo,'eyevel','avg','scales',scales.eyev);
    setTrialAxAvg(eyedistdata,plotparam,axsigrep4{5},plottrials, 'info',eyedistinfo,'eyedist');
    setTrialAxAvg(eyedistdata,plotparam,axsigrep4{6},plottrials, 'info',eyedistinfo,'eyedist','avg','scales',scales.eyedist);
        end

        %savefig(figsigrep3,[pathavgb trialname ])
    saveas(figsigrep3,[pathavgb trialname '.jpg'],'jpg')
           % savefig(figsigrep4,[pathavgbeye trialname ])
    saveas(figsigrep4,[pathavgbeye trialname '.jpg'],'jpg')
end

end
end

if ialn==1
%reaction times only for first alignment event needed, since just targ rt
binwidth=25;              %nm
binmax=600;
binedges=0:binwidth:binmax;
pathrts=[analysispath 'rts' filesep];
if ~isdir(pathrts)
    mkdir(pathrts)
end
for it=1:length(trialtypes(1).names)
    trialname=trialtypes(1).names{it};
    load([pathda trialname],'rts');
    binsrts=zeros(1,length(binedges));
    binsrtsright=zeros(1,length(binedges));
    if isfield(rts,'left')
    [binsrts,bb]=histc(rts.left.*1000,binedges);
    end
    if isfield(rts,'right')
    [binsrtsright,bb]=histc(rts.right.*1000,binedges);
    end
    fighist=figure('visible','off');
    set(0,'CurrentFigure',fighist); 
    xlims=[0 binmax];
    subhist(1)=subplot(4,2,1);
    set(fighist, 'Color', [1 1 1]);
    set(fighist,'Position',[50,50,700,250]);
    countsub=1;
    subhist(1)=subplot(1,2,1);
    subhist(2)=subplot(1,2,2);
    bar(subhist(countsub),binedges,binsrts,'FaceColor',[0 0 0],'BarWidth',1)
    xlabel(subhist(countsub),'reaction time (left) (ms)');
    set(subhist(countsub),'xtick',binedges(1:10:end));
    set(subhist(countsub),'xlim',xlims,'box','off');
    ylabel(subhist(countsub), '# trials counted')
    countsub=countsub+1;
    bar(subhist(countsub),binedges,binsrtsright,'FaceColor',[0 0 0],'BarWidth',1);
    xlabel(subhist(countsub),'reaction time (right) (ms)');
    set(subhist(countsub),'xtick',binedges(1:10:end));
    set(subhist(countsub),'xlim',xlims,'box','off');
    ylim2=get(subhist(countsub),'ylim');
    set(subhist(countsub),'xlim',xlims,'box','off');
    set(subhist(countsub-1),'ylim',ylim2); set(subhist(countsub),'ylim',ylim2);
    set(subhist(countsub),'yticklabel','')
    savefig(fighist,[pathrts trialname]);
    %saveas(figsigrep3,[PathName 'tracesnlxbehav.eps'],'epsc')
    saveas(fighist,[pathrts trialname],'tif')
end
end
close all;
end
%%
%print average FFT spectrograms
if plotfft
pathffts=[analysispath 'avgffts_' targevents{ialn} filesep];
if ~isdir(pathffts)
    mkdir(pathffts)
end
fs=1000;        %sampler ate nlx
figffts=figure('visible','off');
set(figffts, 'Color', [1 1 1]);
set(figffts,'Position',[300,150,1200,500]);
if ispc
    set(figffts,'visible','on');
end
if ~ispc
set(figffts, 'position',[0 0 1000 500]);
end

axfft=axes;
hold(axfft,'on'); 
%figff=figure('visible','off');
%set(figff, 'Color', [1 1 1]);
%set(figff,'Position',[50,50,700,700]);
%axff=axes;
evffts={};
%hold(axff,'on');
fftlims=[5 50];
plotparam.ffttapers=[3 5];      %nw (time-bandiwth rpdocut, tapers)
%plotparam.ffttapers=[1.8 1];      %nw (time-bandiwth rpdocut, tapers)
plotparam.fftwin=[.75 0.15];
for it=1:length(trialtypes(1).names)
    %plottrials=trialtypes.nums{it};
        plottrials=[];
    for iip=1:length(trialtypes)
    plottrials=sort(unique([plottrials [trialtypes(iip).nums{it}]]));
    end
    if ~isempty(goodtrialslfp)
        plottrials=sort(unique(trialtypeslfp(1).nums{it}));
    end
    trialname=trialtypes(1).names{it};
    plotparam.triallabel=trialname;
    if ~isempty(plottrials)
    for ii=1:length(lfpchs)
    cla(axfft);
    %cla(axff);
    cscid=find(ismember(trialbytrial(1).cscNames,lfpchs{ii}));
    if length(cscid)>1
        cscid=cscid(1);
    end
    
    evffts{ii}=setavgfft(trialbytrial(1).lfp{cscid},fs,plotparam,axfft,plottrials,...
        'freq',fftlims,'sitename',lfpchs{ii});
    
    save([pathffts trialname],'evffts','plotparam');
    pathsfftchs{ii}=[pathffts lfpchs{ii} filesep];
    mkdir(pathsfftchs{ii})
    
    %saveas(figffts,[pathsfftchs{ii} trialname '.eps'],'epsc')
    saveas(figffts,[pathsfftchs{ii} trialname '.jpg'],'jpeg')
    savefig(figffts,[pathsfftchs{ii} trialname])
    %{
    evffts{ii}=setavgfft(trialbytrial(1).lfp{cscid},fs,plotparam,axfft,plottrials,...
        'events','fix',samplesfix{1},...
        'events','target',samplestarg{1},...
        'events','fixeye',samplesfixeye{1},...    
        'freq',fftlims,'sitename',lfpchs{ii},'align','target');
    saveas(figffts,[pathsfftchs{ii} trialname '-targaln.jpg'],'jpeg')
    savefig(figffts,[pathsfftchs{ii} trialname '-targaln'])
    %}
    end
    end
end
end
%%
%calculate trial by trial correlation matrix
if plotcor && ialn==1
cortypes={'fixwin','fixpeak','targwin','targpeak',...
    'targimwin','targimpeak',...
    'rewimpeak','rewshortwin','rewshortpeak','rewwin','rewpeak'};
if strcmp(trialgroups{itype},'fixbreak')
cortypes={'fixwin','fixpeak',...
    'rewimpeak','rewshortwin','rewshortpeak','rewwin','rewpeak'};
end
pathcor={};
pathcor{1}=[analysispath 'corr' filesep];
pathcor{2}=[analysispath 'corr2' filesep];      %HF beta
pathcor{3}=[analysispath 'corr3' filesep];      %broad beta

if ~isdir(pathcor{1})
mkdir(pathcor{1});
end
if ~isdir(pathcor{2})
mkdir(pathcor{2});
end
if ~isdir(pathcor{3})
mkdir(pathcor{3});
end
for ifreq=1:length(selbands)
    fbands=allbands{selbands(ifreq)};
    plotparam.fbands=fbands;
    for it=1:length(trialtypes(1).names)
        trialname=trialtypes(1).names{it};
        betatm={};
        datm={};
        if globeta==0
            if selbands(ifreq)==1
                load([pathlfp{1} trialname],'betatm');
            end
            if selbands(ifreq)==2
                load([pathlfp{2} trialname],'betatm');
            end   
        end
        if selbands(ifreq)==3 || globeta
            if exist([pathlfp{3} trialname '.mat'])>0
            load([pathlfp{3} trialname],'betatm');
            end
        end
        load([pathda trialname],'datm');
        if ~isempty(betatm) && ~isempty(datm)
        corrmat=cormatrix(datm,betatm,cortypes,'sitelabels','sametype');
        pathsave=pathcor{1};
        if globeta==0
            if selbands(ifreq)==2
                pathsave=pathcor{2};
            end
            if selbands(ifreq)==3
                pathsave=pathcor{3};
            end
        end
        if globeta
            pathsave=pathcor{3};
        end
        save([pathsave trialname '-sametype'],'corrmat');
        corrmat=cormatrix(datm,betatm,cortypes,'sitelabels');
        save([pathsave trialname],'corrmat');
        if it==1 && ifreq==1
            %get behavior correlations for da, just for "all" trials
            load([pathda trialname]);
            behtm={eyetm,pulsetm,licktm,blinktm,eyevtm,eyedisttm};            
            if sessnum<=37
                behtm={eyetm,blinktm,eyevtm,eyedisttm};
            end
            if ~isempty(eyetm.trialnums)
                %invalid eye data chronic 16/14 cleo
            corrmat=cormatrix(datm,behtm,cortypes,'sitelabels','sametype');
            save([pathsave trialname '-beh-sametype'],'corrmat');
            end
        end
        end
    end
end

%%
for ifreq=1:length(selbands)
    fbands=allbands{selbands(ifreq)};
    plotparam.fbands=fbands;
    pathsave={};
    for it=1:length(trialtypes(1).names)
        pathsave=pathcor{selbands(ifreq)};
        if globeta
            pathsave=pathcor{3};
        end
        trialname=trialtypes(1).names{it};
        if exist([pathsave trialname '.mat'])>0
        load([pathsave trialname],'corrmat');
        cor=corrmat;
        figcorr=figure('visible','off');
        if ~ispc
            set(figcorr,'position',[0 0 1000 750],'color',[1 1 1]);
        else
        set(figcorr,'position',[100 100 1400 900],'color',[1 1 1]);
        end
        set(0,'CurrentFigure',figcorr);    %set figure handle to current figure

        %figure(figcorr,'visible','off')
        %axcorr=axes(figcorr);
        axcorr={};
        for ii=1:length(fscvchs)
            axcorr{ii}=subplot(2,2,ii);
        end
        range=1:length(cortypes);
        xticks=range;
        xlabels=cortypes;
        cmapp=[1 .5 0;1 .55 0; 1 .6 0; 1 .65 0; 1 .7 0; 1 .75 0; 1 .8 0; 1 .85 0; 1 .9 0];
        cmapc=[1 .5 0;1 .55 0; 1 .6 0; 1 .65 0; 1 .7 0; 1 .75 0; 1 .8 0; 1 .85 0; 1 .9 0];
                    ranjitter=-rand(1,length(lfpchs))*.3-.15;

        sitep=find(contains(lfpchs,'p')==1);
        sitec=find(contains(lfpchs,'c')==1);
        sitecolors=zeros(length(lfpchs),3);
        cmapp1=cool; 
        cmapp=flipud(cmapp1(1:12:end,:))-.35;
        if length(sitep)>6
            cmapp=flipud(cmapp1(1:6:end,:))-.35;
        end
        cmapc1=winter; 
        cmapc=flipud(cmapc1(1:12:end,:))-.35;
        if length(sitec)>6
            cmapc=flipud(cmapc1(1:6:end,:))-.35;
        end
        cmapp(cmapp<0)=0;
        cmapc(cmapc<0)=0;
        if ~isempty(sitep)            
            sitecolors(sitep,:)=cmapp(1:length(sitep),:);
        end
        if ~isempty(sitec)    
            sitecolors(sitec,:)=cmapc(1:length(sitec),:);
        end
        for ii=1:length(fscvchs)
            fregion=plotparam.sites{fscvchs(ii)}(1);    %get letter c or p to know region
            fc=cmapp(1,:);
            if strcmp(fregion,'c')
                fc=cmapc(1,:);
            end
            axis(axcorr{ii},'square')
            set(axcorr{ii},'xlim',[0 max(xticks)])
            set(axcorr{ii},'XTick',xticks)
            set(axcorr{ii},'xticklabel',xlabels,'xticklabelrotation',270)
            set(axcorr{ii},'ylim',[0 max(xticks)])
            set(axcorr{ii},'yTick',xticks)
            set(axcorr{ii},'yticklabel',xlabels)    
            ylabel(axcorr{ii},plotparam.sites{fscvchs(ii)})
            title(axcorr{ii},[trialname ' | comod | da-' plotparam.sites{fscvchs(ii)} ', beta lfp'],...
                'color',fc);
            hold(axcorr{ii},'on');
            if ~isfield(cor{fscvchs(ii)},'corrdata')
                continue
            end
            curdata=cor{fscvchs(ii)}.corrdata;
                    divisions{1}=find(contains(cortypes,'fix')==1);
                divisions{2}=find(contains(cortypes,'targ')==1);
                divisions{2}=[divisions{2} find(contains(cortypes,'rewpre')==1)];
                if ~isempty(divisions{2})
                divisions{3}=find(contains(cortypes,'rew')==1 & ...
                    ~contains(cortypes,'rewpre'));
                aa=plot(axcorr{ii},[divisions{2}(1)-1 divisions{2}(1)-1], ...
                    [0 max(xticks)],'linewidth',1,'color',[0 0 0]);
                aa.Color(4)=.2;
                aa=plot(axcorr{ii},[divisions{3}(1)-1 divisions{3}(1)-1], ...
                    [0 max(xticks)],'linewidth',1,'color',[0 0 0]);
                aa.Color(4)=.2;
                aa=plot(axcorr{ii},[0 max(xticks)],...
                    [divisions{2}(1)-1 divisions{2}(1)-1], ...
                    'linewidth',1,'color',[0 0 0]);
                aa.Color(4)=.2;
                aa=plot(axcorr{ii},[0 max(xticks)],...
                    [divisions{3}(1)-1 divisions{3}(1)-1], ...
                    'linewidth',1,'color',[0 0 0]);
                aa.Color(4)=.2;        
                end

            for irow=1:size(curdata,2)
                if contains(curdata(irow).type1,'task') ||...
                    contains(curdata(irow).type2,'task')
                continue
                end
                %get x,y coord to plot based on type labels
                ycord=find(ismember(cortypes,curdata(irow).type1)==1);
                xcord=find(ismember(cortypes,curdata(irow).type2)==1)-0.25;
                eregion=curdata(irow).d2site(1);
                m='+';
                %c=[1 0 0];
                siteidx=strcmp(lfpchs,curdata(irow).d2site);
                pval=curdata(irow).p;
                c=sitecolors(siteidx,:);
                fsize=50*log10(1/pval);
                if curdata(irow).r<0
                    m='o';
                    %c=[0 0 1];
                end
                %{
                if ~strcmp(fregion,eregion)
                    %different regions, offset x slightly & color
                    %c(c~=0)=c(c~=0)-.65;
                    xcord=xcord+.2;
                    ycord=ycord+.2;
                    %fsize=13;
                end
                %}
                %a=text(axcorr{ii},xcord,ycord,m,'fontsize',fsize,'color',c);
        a=scatter(axcorr{ii},xcord+ranjitter(siteidx),ycord+ranjitter(siteidx),fsize,m,'markeredgecolor',c,'MarkerEdgeAlpha',.5,'linewidth',1.5);

            end
            if ii==length(fscvchs)
                yy1=13;
                for ilfp=1:length(lfpchs)
                    c=sitecolors(ilfp,:);
                    text(axcorr{ii},13,yy1,lfpchs{ilfp},'fontsize',10,'color',c);
                    yy1=yy1-1;
                end
            end
            hold(axcorr{ii},'off');
        end
        savefig(figcorr,[pathsave trialname]);
        saveas(figcorr,[pathsave trialname],'jpg')
        end
    end
end

%plot figure behavioral correlations to da
fbands=allbands{selbands(ifreq)};
plotparam.fbands=fbands;
pathsave={};
pathsave=pathcor{selbands(ifreq)};
if globeta
    pathsave=pathcor{3};
end
trialname=trialtypes(1).names{1};
cor=[];
if exist([pathsave trialname '-beh-sametype.mat'])>0
load([pathsave trialname '-beh-sametype']);
cor=corrmat;
end
figcorr=figure('visible','off');
if ~ispc
    set(figcorr,'position',[0 0 1000 750],'color',[1 1 1]);
else
set(figcorr,'position',[100 100 1400 900],'color',[1 1 1]);
end
set(0,'CurrentFigure',figcorr);    %set figure handle to current figure        
axcorr={};
for ii=1:length(fscvchs)
    axcorr{ii}=subplot(2,2,ii);
end

%reaction times / da correlations
rside=find(trialbytrial(1).rewardside==1);
lside=find(trialbytrial(1).rewardside==0);
trialname=trialtypes(1).names{1};
load([pathda trialname]);      %load datm
rs={}; ps={}; rts={};
for ich=1:length(fscvchs)
    plottrials=sort(unique(([trialtypes(fscvchs(ich)).nums{:}])));      %all good trials
rsel=intersect(plottrials,rside);
lsel=intersect(plottrials,lside);
target_rrt=trialbytrial(1).target_rt(rsel);
target_lrt=trialbytrial(1).target_rt(lsel);
fix_rt=trialbytrial(1).fix_rt(plottrials);

    if ~isempty(rsel)
        tnums=intersect(datm{fscvchs(ich)}.trialnums,rsel);
        daids=find(ismember(datm{fscvchs(ich)}.trialnums,tnums)==1);     
[rs(fscvchs(ich)).targwinr,ps(fscvchs(ich)).targwinr]=corr(datm{fscvchs(ich)}.targwin(daids)',trialbytrial(1).target_rt(tnums)','rows','complete');
[rs(fscvchs(ich)).targimwinr,ps(fscvchs(ich)).targimwinr]=corr(datm{fscvchs(ich)}.targimwin(daids)',trialbytrial(1).target_rt(tnums)','rows','complete');
[rs(fscvchs(ich)).targpeakr,ps(fscvchs(ich)).targpeakr]=corr(datm{fscvchs(ich)}.targpeak(daids)',trialbytrial(1).target_rt(tnums)','rows','complete');
[rs(fscvchs(ich)).targimpeakr,ps(fscvchs(ich)).targimpeakr]=corr(datm{fscvchs(ich)}.targimpeak(daids)',trialbytrial(1).target_rt(tnums)','rows','complete');
    end
    if ~isempty(lsel)
                tnums=intersect(datm{fscvchs(ich)}.trialnums,lsel);
        daids=find(ismember(datm{fscvchs(ich)}.trialnums,tnums)==1);   
[rs(fscvchs(ich)).targwinl,ps(fscvchs(ich)).targwinl]=corr(datm{fscvchs(ich)}.targwin(daids)',trialbytrial(1).target_rt(tnums)','rows','complete');
[rs(fscvchs(ich)).targimwinl,ps(fscvchs(ich)).targimwinl]=corr(datm{fscvchs(ich)}.targimwin(daids)',trialbytrial(1).target_rt(tnums)','rows','complete');
[rs(fscvchs(ich)).targpeakl,ps(fscvchs(ich)).targpeakl]=corr(datm{fscvchs(ich)}.targpeak(daids)',trialbytrial(1).target_rt(tnums)','rows','complete');
[rs(fscvchs(ich)).targimpeakl,ps(fscvchs(ich)).targimpeakl]=corr(datm{fscvchs(ich)}.targimpeak(daids)',trialbytrial(1).target_rt(tnums)','rows','complete');
[rs(fscvchs(ich)).fixwin,ps(fscvchs(ich)).fixwin]=corr(datm{fscvchs(ich)}.fixwin(daids)',trialbytrial(1).target_rt(tnums)','rows','complete');
[rs(fscvchs(ich)).fixpeak,ps(fscvchs(ich)).fixpeak]=corr(datm{fscvchs(ich)}.fixpeak(daids)',trialbytrial(1).target_rt(tnums)','rows','complete');
    end
end
rts.rs=rs;
rts.ps=ps;
save([pathsave 'rts_corr'],'rts');


range=1:length(cortypes);
xticks=range;
xlabels=cortypes;
cmapp=[1 .3 0;.4 1 0; .3 .6 .2; 0 .65 1; 1 .7 0; 1 .75 0; 1 .8 0; 1 .85 0; 1 .9 0];
if ~isempty(cor)
for ii=1:length(fscvchs)
    btypes=unique({cor{fscvchs(ii)}.corrdata(1:end).d2site});
    if ~isempty(btypes)
    btypes=btypes(~contains(btypes,'na'));
cmapp=cool; 
    numcol=length(btypes(~contains(btypes,'na')));
    if ~isempty(btypes) 
numdiv=ceil(size(cmapp,1)/numcol);
cmapp=flipud(cmapp(1:numdiv:end,:))-.25;
cmapp(cmapp<0)=0;
sitecolors=cmapp;
    axis(axcorr{ii},'square')
    set(axcorr{ii},'xlim',[0 max(xticks)])
    set(axcorr{ii},'XTick',xticks)
    set(axcorr{ii},'xticklabel',xlabels,'xticklabelrotation',270)
    set(axcorr{ii},'ylim',[0 max(xticks)])
    set(axcorr{ii},'yTick',xticks)
    set(axcorr{ii},'yticklabel',xlabels)    
    ylabel(axcorr{ii},plotparam.sites{fscvchs(ii)})
    title(axcorr{ii},[trialname ' | comod | da-' plotparam.sites{fscvchs(ii)} ', beh'],...
        'color',[0 0 0]);
    hold(axcorr{ii},'on');
    if ~isfield(cor{fscvchs(ii)},'corrdata')
        continue
    end
    curdata=cor{fscvchs(ii)}.corrdata;
    divisions{1}=find(contains(cortypes,'fix')==1);
    divisions{2}=find(contains(cortypes,'targ')==1);
    divisions{2}=[divisions{2} find(contains(cortypes,'rewpre')==1)];
    if ~isempty(divisions{2})
    divisions{3}=find(contains(cortypes,'rew')==1 & ...
        ~contains(cortypes,'rewpre'));
    aa=plot(axcorr{ii},[divisions{2}(1) divisions{2}(1)], ...
        [0 max(xticks)],'linewidth',1,'color',[0 0 0]);
    aa.Color(4)=.2;
    aa=plot(axcorr{ii},[divisions{3}(1) divisions{3}(1)], ...
        [0 max(xticks)],'linewidth',1,'color',[0 0 0]);
    aa.Color(4)=.2;
    aa=plot(axcorr{ii},[0 max(xticks)],...
        [divisions{2}(1)-.5 divisions{2}(1)-.5], ...
        'linewidth',1,'color',[0 0 0]);
    aa.Color(4)=.2;
    aa=plot(axcorr{ii},[0 max(xticks)],...
        [divisions{3}(1)-.5 divisions{3}(1)-.5], ...
        'linewidth',1,'color',[0 0 0]);
    aa.Color(4)=.2;        
    end
    vrows=find((~contains({curdata(1:end).type1},'task') &  ~contains({curdata(1:end).type2},'task'))==1);
    ranjitter=rand(1,length(btypes))*.25-.125;
    for irow=1:length(vrows)
        %get x,y coord to plot based on type labels
        ycord=find(ismember(cortypes,curdata(vrows(irow)).type1)==1);
        xcord=find(ismember(cortypes,curdata(vrows(irow)).type2)==1)-0.25;
        eregion=curdata(vrows(irow)).d2site(1);
        m='+';
        siteidx=strcmp(btypes,curdata(vrows(irow)).d2site);
        pval=curdata(vrows(irow)).p;
        c=sitecolors(siteidx,:);
        fsize=40*log10(1/pval);
        if curdata(vrows(irow)).r<0
            m='o';
        end
        %a=text(axcorr{ii},xcord,ycord,m,'fontsize',fsize,'color',c);
        a=scatter(axcorr{ii},xcord+ranjitter(siteidx),ycord+ranjitter(siteidx),fsize,m,'markeredgecolor',c,'MarkerEdgeAlpha',.5,'linewidth',1.5);
    end
    if ii==length(fscvchs)
        yy1=13;
        for ilfp=1:length(btypes)
            c=sitecolors(ilfp,:);
            text(axcorr{ii},13,yy1,btypes{ilfp},'fontsize',10,'color',c);
            yy1=yy1-2;
        end
    end
    hold(axcorr{ii},'off');
end
    end
end
savefig(figcorr,[pathsave trialname '-beh']);
saveas(figcorr,[pathsave trialname '-beh'],'jpg')
end
end
%%
%get cross covariance just selected trials above
if plotxvar && ialn==1
    %set up event windows to analyze around task for xvar    
    %load eyedata,lickdata,blinkdata,pulsedata
    fscvchstemp=fscvchs;
if isempty(da)
    %need to load
    if exist([pathda 'all-data.mat'])>0
    load([analysispath 'trialbytrial_da' filesep 'all-data'],'da')
    end
end

datax.da=da;    

plotparam.xrates=rates;
datax.rates=rates;

eventmarks={};
eventmarks{1}.da=samplesfix{1};
eventmarks{2}.da=samplestarg{1};
eventmarks{3}.da=samplestargeye{1};
eventmarks{4}.da=repmat(plotparam.alignidx,1,length(samplesfix{1}));
%just provide middle of win rather than range
%let xvardata/plotxvartrials determine win around align idx
eventmarks{2}.name='targ';
eventmarks{1}.name='fix'; 
eventmarks{3}.name='targeye';       
eventmarks{4}.name='outcome';       

eventmarks{1}.tcsc=trialbytrial(1).tfix;
eventmarks{2}.tcsc=trialbytrial(1).ttarg;
eventmarks{3}.tcsc=trialbytrial(1).ttargeye;
eventmarks{4}.tcsc=eventmarks{4}.da./param.samplerate;

if ~isempty(eventx)
    %user selected events to use
    %replace default eventnames and eventmarks
    chosevents=contains(eventnames,eventx);
    eventntemp=eventnames(chosevents);
    eventmtemp=eventmarks(chosevents);
    eventnames=eventntemp;
    eventmarks=eventmtemp;        
end    
if strcmp(trialgroups{itype},'fixbreak')
    eventmarks={};
    eventnames={};
    eventmarks{1}.da=samplesfix{1};
    eventmarks{1}.name='fix';   
    eventmarks{2}.da=samplesfixeye{1};
    eventnames{2}.name='interfixeye';       %1/10/19 added aligned to eye start
end
%extended win for targ
eventmarks3={};
eventmarks3{1}.da=samplestarg{1};
eventmarks3{1}.name='targ_ext';
eventmarks3{1}.tcsc=trialbytrial(1).ttarg;
extpad=2;
%extended win for outcome
eventmarks2={};
eventmarks2{1}.da=repmat(plotparam.alignidx,1,length(samplesfix{1}));
eventmarks2{1}.name='outcome_ext';
eventmarks2{1}.tcsc=eventmarks2{1}.da./param.samplerate;
if sessnum<67
    extpad=1.5;
end

if plotxvarbeh
behname=[analysispath 'trialbytrial_da' filesep 'all-data-beh.mat'];
if exist(behname)>0
load(behname);
trialbytrial(1).eyedata=eyedata;
trialbytrial(1).lickdata=lickdata;
trialbytrial(1).pulsedata=pulsedata;
trialbytrial(1).eyexdata=eyexdata;
trialbytrial(1).eyedistdata=eyedistdata;
trialbytrial(1).eyevdata=eyevdata;
else
    plotxvarbeh=0;
end
end

if plotxvarlfp && ~daflag
if ismember(1,selbands) && isempty(betalfp)
    load([analysispath 'trialbytrial_lfp' filesep 'all-betalfp'],'betalfp','info')
end
if ismember(2,selbands) && isempty(betalfp2)
    load([analysispath 'trialbytrial_lfp2' filesep 'all-betalfp2'],'betalfp2','info')
end
if (ismember(3,selbands) || globeta) && isempty(betalfp3)
    load([analysispath 'trialbytrial_lfp3' filesep 'all-broadbetalfp'],'betalfp3','info')
end
datax.lfp=betalfp;
if globeta
    %broadband beta only
    datax.lfp=betalfp3;
end
%assume only globeta now..
selbands(1)=3;
fbands=allbands{selbands(1)};
plotparam.fbands=fbands;
end
if plotxvarlfp && daflag
    datax.lfp=[];           %only get da signals, no xvar
end

fscvchs=fscvchstemp;
it=1;
%plottrials=trialtypes.nums{it};
plottrials=[];
trialtypesall=trialtypes;
if ~isempty(trialtypeslfp)
trialtypesall=[trialtypes trialtypeslfp];
end
for iip=1:length(trialtypesall)
plottrials=sort(unique([plottrials sort(unique([trialtypesall(iip).nums{:}]))]));
datax.trials{iip}=sort(unique([trialtypesall(iip).nums{:}]));
end
trialname=trialtypesall(1).names{it};
trialname=[];

if plotxvarlfp
trialname=[trialname '3'];
plotparam.triallabel=trialname;
if ~isempty(datax.da)
xcovdata=plotxvartrials(datax,eventmarks,plotparam,'type','none',...
'trials',plottrials,'label',trialname,'cluster');
else
    xcovdata=plotxvartrials(datax,eventmarks,plotparam,'type','none',...
'trials',plottrials,'label',trialname,'onlylfp');
end
%{
xcovdata=plotxvartrials(datax,eventmarks2,plotparam,'type','none',...
'trials',plottrials,'label',trialname,'cluster','xwin',[-.5 3]);
if ~strcmp(trialgroups{itype},'fixbreak')
xcovdata=plotxvartrials(datax,eventmarks3,plotparam,'type','none',...
'trials',plottrials,'label',trialname,'cluster','xwin',[-.5 extpad]);
end
%}
end

if plotxvarbeh
plotparam.triallabel=[];
if ~isempty(lickdata)
    %only plot for sessions with lick data, that means includes proper eye
    %data
if strcmp(subjectname,'patra')
    if sessnum>37
    datax.lfp={eyedata,lickdata,pulsedata,eyexdata,eyevdata,eyedistdata};
    datax.sitesbeh={'eye','lick','pulse','eyex','eyev','eyedist'};
    datax.ratesbeh=[eyeinfo.samplespersec,lickinfo.samplespersec,pulseinfo.samplespersec,...
    eyexinfo.samplespersec,eyevinfo.samplespersec,eyedistinfo.samplespersec];
    else
     datax.lfp={eyedata,lickdata,eyexdata,eyevdata,eyedistdata};
    datax.sitesbeh={'eye','lick','eyex','eyev','eyedist'};
    datax.ratesbeh=[eyeinfo.samplespersec,lickinfo.samplespersec,...
    eyexinfo.samplespersec,eyevinfo.samplespersec,eyedistinfo.samplespersec];
    end
else
    datax.lfp={eyedata,eyexdata,eyevdata,eyedistdata};
    datax.sitesbeh={'eye','eyex','eyev','eyedist'};
    datax.ratesbeh=[eyeinfo.samplespersec,...
    eyexinfo.samplespersec,eyevinfo.samplespersec,eyedistinfo.samplespersec];
end
xcovdata=plotxvartrials(datax,eventmarks,plotparam,'type','none',...
'trials',plottrials,'label',trialname,'cluster','beh');    
end
end
    
    
    
%split trial types as needed post-pro
end
%end session types
end
%end if plotxclust==0
end
end
%END ~allsess
%END INDIVIDUAL SESSION TYPES (BIG,SMALL, ETC) ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ANALYSIS ACROSS SESSION TYPES

if plotcorall
analysispath=getfield(targpath,trialgroups{1});
fend=strfind(analysispath,filesep);
savepath=[analysispath(1:fend(end-1)) 'all' filesep 'corr' filesep];
plotparam.savepath=savepath;
if ~isdir(savepath)
    mkdir(savepath);
end
cortypes={'fixwin','fixpeak','targwin','targpeak',...
'targimwin','targimpeak','targpeakabs',...
'rewimpeak','rewshortwin','rewshortpeak','rewwin','rewpeak'};
cortypesb={'fixwin','fixpeak','targwin','targpeak',...
'targimwin','targimpeak',...
'rewimpeak','rewshortwin','rewshortpeak','rewwin','rewpeak'};
datmall=[];
betatmall=[];
behtmall=[];
datms={};
betatms={};
behtms={};
tlist={};
trials={};
getvars={cortypes{:}, 'site','trialnums'};
getvarsb={cortypesb{:},'site','trialnums'};
for itype=1:2
    %just get big and small
    analysispath=getfield(targpath,trialgroups{itype});
    plotparam.pathname=analysispath;
    plotparam.trialtype=trialgroups{itype};
    load([analysispath 'trialtypes']); %good trials load
    load([analysispath 'trials_chronic' sessid '_' trialgroups{itype}],'trialbytrial'); %load trial data for given session type
            goodtrials=[];
            swtrials=[];
            postswtrials=[];
                            trialtypesall=trialtypes;
            if exist('trialtypeslfp')
                trialtypesall=[trialtypes trialtypeslfp];
            end
    for iip=1:length(trialtypesall)
        %integrate all channels good trials
    goodtrials=sort(unique([goodtrials [trialtypesall(iip).nums{1}]]));
        swtrials=sort(unique([swtrials trialtypesall(iip).nums{find(strcmp(trialtypesall(1).names,'switch'))}]));          %retain orig sw/post-sw trials
        postswtrials=sort(unique([swtrials trialtypesall(iip).nums{find(strcmp(trialtypesall(1).names,'postswitch'))}]));          %retain orig sw/post-sw trials
    end
    %goodtrials=trialtypes.nums{1};
    %postswtrials=trialtypes.nums{find(strcmp(trialtypes.names,'postswitch'))};
    %    save([analysispath 'trialtypes'], 'trialtypes');

   % trialtypes=grouptrials(goodtrials,trialbytrial(1),plotparam,'prev');        %since rerunning this after getting good trials that do not have sw/pos-sw, these types will be empty
    %trialtypes.nums{find(strcmp(trialtypes.names,'switch'))}=swtrials;
    %trialtypes.nums{find(strcmp(trialtypes.names,'postswitch'))}=postswtrials;
    %goodtrials=trialtypes.nums{1};

    %paths to save exported signals from settrialax
    pathda=[analysispath 'trialbytrial_da' filesep];
    pathlfp=[analysispath 'trialbytrial_lfp3' filesep];    
    %load task-modulated signals for compiling data big/small together
    existda=0;
    if exist([pathda 'all.mat'])>0
        existda=1;
        load([pathda 'all']);   %loads datm, blinktm, eyetm, etc. for good all trials   
                behtm={eyetm,pulsetm,licktm,blinktm,eyevtm,eyedisttm};
                datms{itype}=datm;
    behtms{itype}=behtm;
    end
    existlfp=0;
    if exist([pathlfp 'all.mat'])>0
        existlfp=1;
        load([pathlfp 'all'],'betatm');
            betatms{itype}=betatm;
    end
        
   % if strcmp(subjectname,'cleo')
   % behtm={eyetm,blinktm,eyevtm,eyedisttm};
   % end        
   % if strcmp(subjectname,'cleo') && sessnum<=24
    
for ib=1:length(behtm)
    if itype==1
        behtmall{ib}=[];        %initialize
    end
    if ~isempty(behtm{ib})
        if ~isempty(behtm{ib}.trialnums)
    for ivar=1:length(getvars)
        trialids=find(ismember(getfield(behtm{ib},'trialnums'),goodtrials)); %make sure to retrieve only good trial data
        if itype==1 || strcmp(getvars{ivar},'site')
            curval=[]; %initialize
        else
            curval=getfield(behtmall{ib},getvars{ivar});
        end
        newval=getfield(behtm{ib},getvars{ivar});
        if ~strcmp(getvars{ivar},'site')
            %get only targeted trial data
            newval=newval(trialids);
        end
        putval=[curval newval];
        behtmall{ib}=setfield(behtmall{ib},getvars{ivar},putval);        
    end   
        end
    end 
end

if existda
for ich=1:length(fscvchs)
    ida=fscvchs(ich);
    datrials=trialtypesall(ida).nums{1};
    if itype==1
       datmall{ida}=[];        %initialize
    end
    if ~isempty(datm{ida})
    for ivar=1:length(getvars)
        trialids=find(ismember(getfield(datm{ida},'trialnums'),datrials)); %make sure to retrieve only good trial data
        if itype==1 || strcmp(getvars{ivar},'site')
            curval=[]; %initialize
        else
            %concatenating to last type, get last type data
            curval=getfield(datmall{ida},getvars{ivar});
        end
        newval=getfield(datm{ida},getvars{ivar});
        if ~strcmp(getvars{ivar},'site')
            %get only targeted trial data
            newval=newval(trialids);
        end
        putval=[curval newval];
        datmall{ida}=setfield(datmall{ida},getvars{ivar},putval);        
    end        
    end
    tlist(ida).sizecurr=[];
tlist(ida).sizeprev=[];
tlist(ida).outcomeprev=[];
tlist(ida).sidecurr=[];
    trials.groups{itype}=trialtypesall;       %all good (datm indexed) trials categorized 
%list of previous trial big/sm/non rewarded
lastendidx=length(tlist(ida).sizecurr);      %previous last index of trial ids
tlist(ida).sizecurr=[tlist(ida).sizecurr; repmat({trialgroups{itype}(1:3)},length(trialtypesall(ida).nums{1}),1)];
tlist(ida).sizeprev=[tlist(ida).sizeprev; repmat({'non'},length(trialtypesall(ida).nums{1}),1)];
tlist(ida).sidecurr=[tlist(ida).sidecurr; repmat({'nann'},length(trialtypesall(ida).nums{1}),1)];
afterbigids=find(ismember(trialtypesall(ida).nums{1},trialtypesall(ida).nums{find(contains(trialtypesall(ida).names,'afterbig'))}))+lastendidx;
tlist(ida).sizeprev(afterbigids)={'big'};
aftersmids=find(ismember(trialtypesall(ida).nums{1},trialtypesall(ida).nums{find(contains(trialtypesall(ida).names,'aftersm'))}))+lastendidx;
tlist(ida).sizeprev(aftersmids)={'sma'};
tlist(ida).outcomeprev=[tlist(ida).outcomeprev; repmat({'nann'},length(trialtypesall(ida).nums{1}),1)];
aftersuccids=find(ismember(trialtypesall(ida).nums{1},trialtypesall(ida).nums{find(contains(trialtypesall(ida).names,'aftersuccess'))}))+lastendidx;
afterfailids=find(ismember(trialtypesall(ida).nums{1},trialtypesall(ida).nums{find(contains(trialtypesall(ida).names,'afterfail'))}))+lastendidx;
tlist(ida).outcomeprev(aftersuccids)={'succ'};
tlist(ida).outcomeprev(afterfailids)={'fail'};
left=find(ismember(trialtypesall(ida).nums{1},trialtypesall(ida).nums{find(contains(trialtypesall(ida).names,'left'))}))+lastendidx;
right=find(ismember(trialtypesall(ida).nums{1},trialtypesall(ida).nums{find(contains(trialtypesall(ida).names,'right'))}))+lastendidx;
tlist(ida).sidecurr(left)={'left'};
tlist(ida).sidecurr(right)={'right'};

end
end
if existlfp
for ib=1:length(betatm)
    if itype==1
       betatmall{ib}=[];        %initialize
    end
    if ~isempty(betatm{ib})
    for ivar=1:length(getvarsb)
        trialids=find(ismember(getfield(betatm{ib},'trialnums'),goodtrials)); %make sure to retrieve only good trial data
        if itype==1 || strcmp(getvarsb{ivar},'site')
            curval=[]; %initialize
        else
            curval=getfield(betatmall{ib},getvarsb{ivar});
        end
        newval=getfield(betatm{ib},getvarsb{ivar});
        if ~strcmp(getvarsb{ivar},'site')
            %get only targeted trial data
            newval=newval(trialids);
        end
        putval=[curval newval];
        betatmall{ib}=setfield(betatmall{ib},getvarsb{ivar},putval);        
    end        
    end
end 
end
if itype==1
trials.type{itype}='big';
end
if itype==2
trials.type{itype}='small';
end
end

filename='corr-bigsmall';
corrmat={};
%{
if ~isempty(datmall) && ~isempty(betatmall)
corrmat=cormatrix(datmall,betatmall,cortypesb,'sitelabels');
end
%}
%anova analysis
%p = anovan(datmall{4}.targpeakabs',{tlist.sizecurr,tlist.sizeprev},'interaction');
%ab=plothistogroups(datmall{4}.targpeakabs,{tlist.sizecurr,tlist.sizeprev},'varnames',{'current','previous'},'responsename','targpeakabs','plot');
%ab=plothistogroups(datmall{4}.targpeak,{tlist.sizecurr,tlist.sizeprev},'varnames',{'current','previous'},'responsename','targwin','plot');
%p = anovan(datmall{4}.targpeakabs,{tlist.sizecurr,tlist.sizeprev,tlist.outcomeprev},3,3);
%convert cell entries of site labels to char for indexing;
%fixed cormatrix 2/2091
%{
for ib=1:length(corrmat)
    if ~isempty(corrmat{ib})
        corrtemp=corrmat{ib}.corrdata;
        for ix=1:length(corrmat{ib}.corrdata)
            dsit=getfield(corrmat{ib}.corrdata,{ix},'d1site');
            if iscell(dsit)
                dsit=dsit{:};
            end
            corrtemp=setfield(corrtemp,{ix},'d1site',dsit);
        end
        corrmat{ib}.corrdata=corrtemp;
    end
end
%}
save([savepath filename],'corrmat','datmall','betatmall','trials','tlist');
rewexp={};
sideexp={};

if ~isempty(corrmat)
corplot(corrmat,cortypesb,plotparam,'label',filename);

fscvchs=find(~cellfun(@isempty,corrmat));
fsize=20;
fontsize=12;
ccol=[0 0 0; 1 0 0; 0 0 1];
m={'.','o','x'};
figcorr=figure('visible','off');
if ~ispc
    set(figcorr,'position',[0 0 1000 750],'color',[1 1 1]);
else
set(figcorr,'position',[100 100 1400 900],'color',[1 1 1]);
set(figcorr,'visible','on');
end
axcorr={};
axsiz=[300 250];
for ii=1:max(fscvchs)
    axcorr{ii}=subplot(2,2,ii);
    set(axcorr{ii},'units','pixels');
end
for ii=1:length(fscvchs)
    ida=fscvchs(ii);
    hold(axcorr{ida},'on');
    dasite=plotparam.sites{ida};
    
    %ANOVA tests & plots
    %p = anovan(datmall{4}.targpeakabs',{tlist.sizecurr,tlist.sizeprev},'interaction');
    rewexp{ida}=plothistogroups(datmall{ida},...
        {tlist(ida).sizecurr,tlist(ida).sizeprev},'varnames',{'current','previous'},...
        'response','targpeakabs','plot','savename',[savepath 'rewexphisto_' dasite]);
    sideexp{ida}=plothistogroups(datmall{ida},...
            {tlist(ida).sizecurr,tlist(ida).sidecurr},'varnames',{'current','current'},...
            'response','targpeakabs','plot','savename',[savepath 'sideexphisto_' dasite]);
    
[~,lfptarg]=getsitepairs(dasite);
targrows=find(contains({corrmat{ida}.corrdata.type1},'targ') & ...
    contains({corrmat{ida}.corrdata.type2},'targ') & ...
    strcmp({corrmat{ida}.corrdata.d1site},dasite) & ...
    contains({corrmat{ida}.corrdata.d2site},lfptarg));
if ~isempty(targrows)
targrow=targrows(1);
datype=corrmat{ida}.corrdata(targrow(1)).type1;
dadata=getfield(datmall{ida},datype);
[sortedda,sortda]=sort(dadata);
lfptype=corrmat{ida}.corrdata(targrow(1)).type2;
lfpch=corrmat{ida}.corrdata(targrow(1)).d2ch;
lfpdata=getfield(betatmall{lfpch},lfptype);
lfpsite=corrmat{ida}.corrdata(targrow(1)).d2site;
rval=corrmat{ida}.corrdata(targrow(1)).r;
pval=corrmat{ida}.corrdata(targrow(1)).p;
xdata=dadata(sortda)';
xnan=find(isnan(xdata)==1);
ydata=lfpdata(sortda)';
if ~isempty(xnan)
xdata=xdata(1:xnan(1)-1);
ydata=ydata(1:xnan(1)-1);
end
ynan=isnan(ydata);
if ~isempty(ynan)
ydata=ydata(~ynan);
xdata=xdata(~ynan);
end
scatter(axcorr{ida},xdata,ydata,fsize,m{1},'markeredgecolor',ccol(1,:),'MarkerEdgeAlpha',.5,'linewidth',1);
%linear regression model
b1=xdata\ydata;
[p,s]=polyfit(xdata,ydata,1);
slope=p(1);
intercep=p(2);
yfit=slope*xdata+intercep;
yresid=ydata-yfit;
ssresid=sum(yresid.^2);
sstotal=(length(ydata)-1)*var(ydata);
rsq=1-ssresid/sstotal;
fitline=plot(axcorr{ida},xdata,yfit,'-','color',ccol(3,:));
fitline.Color(4)=0.5;
xsiz=get(axcorr{ida},'position');
ylims=get(axcorr{ida},'ylim');
text(axcorr{ida},xsiz(3)-120,xsiz(4)-50,...
    {['slope: ' num2str(slope)],...
    ['intercept: ' num2str(intercep)],['r: ' num2str(rval)],...
    ['p: ' num2str(pval)]},'color',ccol(3,:),'units','pixels');
xlabel(axcorr{ida},['\DeltaDA (nM) | ' datype ' | ' dasite])
ylabel(axcorr{ida},['beta-lfp (\muV^2) | ' lfptype ' | ' lfpsite])
set(findall(axcorr{ida},'-property','FontSize'),'FontSize',fontsize)
end
end

savefig(figcorr,[savepath filename '-scatter']);
saveas(figcorr,[savepath filename '-scatter'],'tif')
end
%delete(findall(figcorr,'type','text')) 
corrmat={};
%{
if ~isempty(datmall) && ~isempty(betatmall)
corrmat=cormatrix(datmall,betatmall,cortypesb,'sitelabels','sametype');
end
%}

%ignore corrmat > 05/15/2019 especially for "all" variables with mixed-type
%trials... not accurate with cormatrix
save([savepath 'corr-bigsmall-sametype'],'corrmat','datmall','datms','betatmall','betatms','trials','tlist','rewexp','sideexp');

%get behavior correlations 
%never use behtms in plotmultiple....ignore for now
%{
if ~isempty(datmall) && ~isempty([behtmall{:}])
corrmat=cormatrix(datmall,behtmall,cortypesb,'sitelabels','sametype');
save([savepath 'corr-bigsmall-beh-sametype'],'corrmat',...
    'behtmall','behtms','datmall','datms','trials');
corrmat=cormatrix(datmall,behtmall,cortypesb,'sitelabels');
save([savepath 'corr-bigsmall-beh'],'corrmat',...
    'behtmall','datmall','trials');
corplot(corrmat,cortypesb,plotparam,'label',[filename '-beh'],'bplot');
end
%}

close all;

end

%%%%%X Var analyses on already processed xinfo data across session types
%{

if plotxclust==1
%%
%post processing all sessions comparison of xvar's
typecount=0;
bcount=0;   %behavior data types
%targevents={'interfix','intertarg','intertargeye','interfixeye'};
fscvchstemp=fscvchs;        %store because when load trialbytrial, loads all again replacing fscvchs
%get lfp ch groups p & c
pgroup=find(contains(plotparam.lfpchs,'p')==1);
cgroup=find(contains(plotparam.lfpchs,'c')==1);
%analysispath=getfield(targpath,trialgroups{1});
%group trial types
%numb=length(trialgroups)*length(targevents);
%binfo(numb)=struct();
trialinfo=[];
ifreq=1;
    fbands=allbands{selbands(ifreq)};     %only low band modulated seems
    plotparam.fbands=fbands;  
    xinfo=[];
    binfo=[];
    xbinfo=[];
analysispath=getfield(targpath,trialgroups{1});
fend=strfind(analysispath,filesep);
savepath=[analysispath(1:fend(end-1)) 'all' filesep 'xsess' filesep];
plotparam.savepath=savepath;
if ~isdir(savepath)
    mkdir(savepath);
end

for itype=1:2
    %big, small, targbreak    
analysispath=getfield(targpath,trialgroups{itype});
plotparam.pathname=analysispath;
plotparam.trialtype=trialgroups{itype};
load([analysispath 'trialtypes']);
disp(analysispath)
compiledname=[analysispath 'trials_chronic' sessid '_' trialgroups{itype}];
load(compiledname);
fscvchs=fscvchstemp;
%load eyedata,lickdata,blinkdata,pulsedata
 %goodtrials=trialtypes(1).nums{1};
 %{
    if exist(manlist)>0
    goodtrials=manualtrialsel(manlist,2,goodtrials,trialgroups{itype});          %get recent list of bad trials manually selected in excel sheet
    end  
trialtypes=grouptrials(goodtrials,trialbytrial(1),plotparam,'prev');
save([analysispath 'trialtypes'], 'trialtypes');
 %}
behname=[analysispath 'trialbytrial_da' filesep 'all-data-beh.mat'];
trialbytrial(1).eyedata=[];
trialbytrial(1).lickdata=[];
trialbytrial(1).pulsedata=[];
trialbytrial(1).eyexdata=[];
trialbytrial(1).eyedistdata=[];
if exist(behname)>0
load(behname);
trialbytrial(1).eyedata=eyedata;
trialbytrial(1).lickdata=lickdata;
trialbytrial(1).pulsedata=pulsedata;
trialbytrial(1).eyexdata=eyexdata;
trialbytrial(1).eyedistdata=eyedistdata;
end

vartype='none';
num3='3';
plotparam.chnums=fscvchs;    
plotparam.triallabel=[];
for ii=1:length(fscvchs)
    ich=fscvchs(ii);    
    for ievent=1:length(targevents)        
        pathxvar=[analysispath 'xvar_' targevents{ievent} '_' ...
            vartype num3 filesep];
        if exist([pathxvar 'xvar_' sites{ich} '.mat'])==2
            %file exists
            load([pathxvar 'xvar_' sites{ich}],'xcovdata'); 
            %organize behavioral parameters for trial group (big/small/etc)
            %for given event (interfix/targ/etc.)
             %GET ALIGNED BEHAVIOR TRACES WITH XCOV--MODIFY AS INPUT
           
            if ii==1 
                bcount=bcount+1;
                seltrials=[];
                winids=[];
                if daflag
                   seltrials=xcovdata.seltrials;
                winids=xcovdata.winids;  
                else
                seltrials=xcovdata{1}.seltrials;
                winids=xcovdata{1}.winids;     
                end
                binfotemp=orgbeh(trialbytrial(1),plotparam,...
                    targevents{ievent},seltrials,winids,'win',[-2 2]);
                binfotemp.sessiontype=trialgroups{itype};
                binfotemp.event=targevents{ievent};
                binfo=[binfo binfotemp];
                %get trial types again
                %trialtypes={};      %reset, redefined below in loop
                %trialtypes=grouptrials(goodtrials,trialbytrial(1),plotparam,'prev');

               % save([analysispath 'trialtypes'], 'trialtypes');
                plotparam.trialtypes=trialtypes;    %needed for trialinfo /xplot
            end
            xinfotemp=[];
            if daflag
                xinfotemp=xclust(xcovdata,plotparam,'event',targevents{ievent},'nolfp');    
            else
                xinfotemp=xclust(xcovdata,plotparam,'event',targevents{ievent});    
            end
            %xinfotemp=xclust(xcovdata,plotparam);
            xinfo=[xinfo xinfotemp];        %concatenate all info      
        end
       
        pathxvarb=[analysispath 'xvar_beh_' targevents{ievent} '_' ...
            vartype num3 filesep];
        if exist([pathxvarb 'xvar_' sites{ich} '.mat'])==2
        load([pathxvarb 'xvar_' sites{ich}],'xcovdata'); 
        xbinfotemp=xclust(xcovdata,plotparam,'event',targevents{ievent},'beh');                
        xbinfo=[xbinfo xbinfotemp]; 
        end
    end
end

    %save trial ids for other trial classifications in group
    %same across all event types
    trialinfo(itype).trialtypes=plotparam.trialtype;
    trialinfo(itype).sessiontype=trialgroups{itype};
        
end
save([savepath 'xinfo' ], 'xinfo','trialinfo');
save([savepath 'binfo' ], 'binfo','trialinfo');
save([savepath 'xbinfo' ], 'xbinfo','trialinfo');

plotparam.savepath=savepath;        %needed for xplot

    
end

if plotxsess && ~daflag
    %plot multiple types of sess (big,small/targ...et)
%post processing all sessions comparison after all data saved from xclust
%FIX 2/17/2019, trialtypes were loading only big reward for all trial types
%so did not make correct plotxtracesavgov (plotxallsess
pgroup=find(contains(plotparam.lfpchs,'p')==1);
cgroup=find(contains(plotparam.lfpchs,'c')==1);
ifreq=1;
    fbands=allbands{selbands(ifreq)};     %only low band modulated seems
    plotparam.fbands=fbands;  
    xinfo=[];
analysispath=getfield(targpath,trialgroups{1});
fend=strfind(analysispath,filesep);
savepath=[analysispath(1:fend(end-1)) 'all' filesep 'xsess' filesep];
plotparam.savepath=savepath;
if ~isdir(savepath)
    mkdir(savepath);
end

limpath=strfind(analysispath,filesep);
xnam='xsess';      %NOW DEFAULT
load([savepath 'xinfo' ]);
%{
for itype=1:3
analysispath=getfield(targpath,trialgroups{itype});
load([analysispath 'trialtypes'], 'trialtypes');
trialinfo(itype).trialtypes=trialtypes;
end
%}

trialinfo=[];
trialgrps=[];
for itype=1:2
analysispath=getfield(targpath,trialgroups{itype});
load([analysispath 'trialtypes'], 'trialtypes');
%trialgrps(1).trialinfo(itype)=trialtypes;
%goodtrials=trialtypes.nums{1};
load([analysispath 'trials_chronic' sessid '_' trialgroups{itype}],'trialbytrial'); %load trial data for given session type
%trialtypes=grouptrials(goodtrials,trialbytrial(1),plotparam,'prev');
trialgrps(1).trialinfo(itype)=trialtypes;                           %REDO TRIAL TYPES TO INCLUDE PREVIOSU BIG/SM
trt=trialbytrial(1).target_rt;
frt=trialbytrial(1).fix_rt;
targrows=find(contains({xinfo.sessiontype},trialgroups{itype}));
xinfotemp=xinfo;
for ix=targrows
xinfotemp=setfield(xinfotemp,{ix},'trt',trt);
xinfotemp=setfield(xinfotemp,{ix},'frt',frt);
end
xinfo=xinfotemp;
end



datypes={'dapos','daall'};
wins=[-1 4];
if sessnum<67
    wins=[-1 3];
end
[~,lfptarg]=getsitepairs(plotparam.sites(plotparam.chnums));

trialinfo=trialgrps(1).trialinfo;

for datype=1:length(datypes)
    if datype==1
        plotxtraces(xinfo,trialinfo,plotparam,datypes{datype},'win',wins,'plotz','evtype',{'targ'},'lfps',lfptarg,'sesstype',{'reward'},'sort','damaxts');
        plotxtraces(xinfo,trialinfo,plotparam,datypes{datype},'win',wins,'evtype',{'targ'},'lfps',lfptarg,'sesstype',{'reward'},'sort','damaxts');
                plotxtraces(xinfo,trialinfo,plotparam,datypes{datype},'win',wins,'evtype',{'targeye'},'lfps',lfptarg,'sesstype',{'reward'},'sort','damaxts');
                
plotxtraces(xinfo,trialinfo,plotparam,datypes{datype},'win',wins,'plotz','evtype',{'targeye'},'sesstype',{'reward'},'lfps',lfptarg,'sort','rts');
plotxtraces(xinfo,trialinfo,plotparam,'daall','win',wins,'plotz','evtype',{'targeye'},'sesstype',{'reward'},'lfps',lfptarg,'sort','fixrts');
plotxtraces(xinfo,trialinfo,plotparam,'daall','win',wins,'evtype',{'targeye'},'sesstype',{'reward'},'lfps',lfptarg,'sort','rts');
plotxtraces(xinfo,trialinfo,plotparam,'daall','win',wins,'plotz','evtype',{'targeye'},'sesstype',{'reward'},'lfps',lfptarg,'sort','damaxts','plotmarks','lfpmints');
plotxtraces(xinfo,trialinfo,plotparam,'daall','win',wins,'condition',{'left','right'},'evtype',{'targeye'},'sesstype',{'big'},'lfps',lfptarg,'sort','damax');
plotxtraces(xinfo,trialinfo,plotparam,'daall','win',wins,'condition',{'left','right'},'evtype',{'targeye'},'sesstype',{'small'},'lfps',lfptarg,'sort','damax');
plotxtraces(xinfo,trialinfo,plotparam,'daall','win',wins,'condition',{'left','right'},'evtype',{'targeye'},'sesstype',{'reward'},'lfps',lfptarg,'sort','damax');
    end
    %temporary comment, already finished for some sessions, dnr
    plotxtraces(xinfo,trialinfo,plotparam,datypes{datype},'win',wins);
    
    %plot averages overlay comparison diff conditions
    plotxtracesavgov(xinfo,trialinfo,plotparam,'evtype','targ','sesstypes',...
        {'bigreward'},'trialtypes',{'sleft','sright'},datypes{datype},'win',wins);
    plotxtracesavgov(xinfo,trialinfo,plotparam,'evtype','targ','sesstypes',...
        {'bigreward','smallreward'},'trialtypes',{'sleft'},datypes{datype},'win',wins);
    plotxtracesavgov(xinfo,trialinfo,plotparam,'evtype','targeye','sesstypes',...
        {'bigreward','smallreward'},'trialtypes',{'sleft'},datypes{datype},'win',wins);
    %plotxtracesavgov(xinfo,trialinfo,plotparam,'evtype','targeye','sesstypes',...
    %    {'bigreward','smallreward'},'trialtypes',{'sright'},datypes{datype},'win',wins);
    
    plotxtracesavgov(xinfo,trialinfo,plotparam,'evtype','targ','sesstypes',...
        {'bigreward','smallreward'},datypes{datype},'win',wins);
end
datype=2;
        plotxtraces(xinfo,trialinfo,plotparam,datypes{datype},'sort','damax');

end

if plotsort
    %plot xvar signals sorted side bys ide   
analysispath=getfield(targpath,trialgroups{1});
load([analysispath 'trialtypes'], 'trialtypes');
limpath=strfind(analysispath,filesep);
xnam='xsess3';      %NOW DEFAULT
loadpath=[analysispath(1:limpath(end-1)) xnam filesep];
load([loadpath 'xinfo' ], 'xinfo','binfo','trialinfo');
plotparam.trialtypes=trialtypes;  
targevents={'interfix','intertarg','intertargeye','interfixeye'};
plotparam.savepath=loadpath;        %needed for xplot

for datype=1:length(datypes)
        plotxtraces(xinfo,trialinfo,plotparam,datypes{datype},'sort','damax');
end
end

%}
end
