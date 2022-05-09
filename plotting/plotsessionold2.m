function plotsession(sessnum,varargin)
%script to run all functions to get compiled trial data lfp/da
%assume patra recording
%opens chconfigsimple to load targeted ncs channels
global plotParam
plotParam={};
argnum=1;
fscvchs=1:4;
targch=2;       %fscv ch used for sorting or getting good trials from list
fbands={};
baseline=[];
fbands2={};         %second band to plot (hf)
fbands3={};         %global beta, need to input argument 'globeta' for only this band
globeta=0;
selbands=[1 2 3];
%plot flags
plotda=0;
plotlfp=0;
plotavgda=0;
plotavglfp=0;
plotphys=0;
plotfft=0;
plotcor=0;
plotxvar=0;
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
getinfo=0;
trialtypes={};
%targevents={'intertarg','interfix','targwin'};
targevents={'intertarg','interfix','intertargeye','interfixeye'};

eventx={};
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
            chbands;
            globeta=1;
        case 'xchs'
            argnum=argnum+1;
            xchs=varargin{argnum};
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
            plotbetaxcov=1;
            plotdaxcov=1;
        case 'plotda'
            plotda=1;
        case 'plotlfp'
            plotlfp=1;
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
        case 'plotxvar'
            plotxvar=1;
        case 'plotxvarall'
            plotxvar=1;
            plotxvarall=1;
            plotdaxcov=1;
            plotbetaxcov=1;
        case 'plotdaxcov'
            %da inc/dec xcov's
            plotxvar=1;
            plotdaxcov=1;
        case 'plotbetaxcov'
            %beta burst xcov's
            plotxvar=1;
            plotbetaxcov=1;
        case 'nofft'
            plotfft=0;
        case 'plotxclust'
            plotxclust=1;       %only after getting all data
            %must supply 'targevents' ie. 'intertarg' after
        case 'plotxsess'
            plotxsess=1;       %only after getting all data
        case 'plotsort'
            plotsort=1; %also after getting all data after plotxclust
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
    end
    argnum=argnum+1;
end
xchs{1}=fscvchs(1:2);
if length(fscvchs)>2
xchs{2}=fscvchs(3:end);
end

sessid=num2str(sessnum);
%get dir with config files
pctype=computer;
ispc=strcmpi(pctype,'pcwin64');
%default on putamen pc in lab
configdir='A:\mit\injectrode\experiments\fscv\matlab\analysis\analysis\config\';
if ~ispc
    %chunky dir
    configdir=fullfile(filesep,'home','schwerdt','matlab','analysis','analysis','config',filesep);
end

d=dir(configdir);
filenames={d.name};
targfiles=strfind(filenames,['chronic' sessid 'chconfigsimple.m']);
processfiles=find(~cellfun(@isempty,targfiles));
targconfigname=filenames{processfiles};
run([configdir 'patra_map_bipolar']);   %run patra_map_bipolar for ch settings
run([configdir targconfigname]);        
%run config file to get ncschannels & paths

targpathtemp.bigreward=fullfile(paths{1}, 'matlab','bigreward_pro','analyzed',filesep);
targpathtemp.smallreward=fullfile(paths{1}, 'matlab','smallreward_pro','analyzed',filesep);
targpathtemp.targetbreak=fullfile(paths{1}, 'matlab','targetbreak_pro','analyzed',filesep);
targpathtemp.fixbreak=fullfile(paths{1}, 'matlab','fixbreak_pro','analyzed',filesep);
homedirf=strfind(paths{1},'1dr');
homedir=paths{1}(1:homedirf-1);
manlist=[homedir 'chronic' sessid '_trialselection.xlsx'];
assignin('base','fcnStatus',targpathtemp)   %store targpath in workspace

targpath={};
if plotbig
    targpath.bigreward=targpathtemp.bigreward;
end
if plotsmall
    targpath.smallreward=targpathtemp.smallreward;
end
if plottargbreak
    targpath.targetbreak=targpathtemp.targetbreak;
end
if plotfixbreak
    targpath.fixbreak=targpathtemp.fixbreak;
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

[param,csc_map,eventcodes]=getparams('patrabipolar','default',ncschannels);
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
rates=[plotparam.ratelfp param.samplerate];
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
origpath=paths{1}(1:strfind(paths{1},'cvtotxt')-1);
filesorig=dir(origpath);
filenamesorig={filesorig.name};
targfiles=strfind(filenamesorig,'1dr_');
labeledfiles=find(~cellfun(@isempty,targfiles));
filesorig=filenamesorig(labeledfiles);
sitenamessp=strsplit(filesorig{1},'_');
sites=sitenamessp(2:5);
plotparam.sites=sites;
plotparam.lfpchs=lfpchs;
[parameters,csc_map,event_codes]=getparams('patrabipolar');
load([configdir 'colormapparula']);
parula=pmap;        %for chunky where parula color map not available;
plotparam.colormap=parula;
plotparam.markercolors=[1 1 1; 1 0 1; 0 1 1];
plotparam.cminshow=0;

numfbands=1;
allbands={fbands};
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
    allbands={fbands3};
end

if plotxclust==0 && plotxsess==0 && ~plotsort
for itype=1:numtypes
    plotparam.baseline='pretarg';  
    plotparam.shift=[];
    currpath=getfield(targpath,trialgroups{itype});
    disp(currpath)
    plotparam.pathname=currpath;
    plotparam.trialtype=trialgroups{itype};
    if ~isdir(currpath)
        %files/folder not created by sync sigs, skip
        disp(['folder not created, skipping : ' char(10) currpath]);
        continue
    end
    if strcmp(trialgroups{itype},'targetbreak')
        %shift time to targ
        plotparam.shift='targ'; 
        plotparam.baseline='pretarg';  
    end
    if strcmp(trialgroups{itype},'fixbreak')
        %shift time to targ
        plotparam.shift='fix'; 
        plotparam.baseline='precue'; 
    end    
    
    %open compiled data
    compiledname=[currpath 'trials_chronic' sessid '_' trialgroups{itype} '_pro'];
    fscvchstemp=fscvchs;
    load(compiledname); %reloads fscvchs
    fscvchs=fscvchstemp;
    
    %set up event windows to analyze around task for xvar
    eventsrew=repmat(plotparam.alignidx,1,length(samplesfix{1}));
    rewwin=5;           %5s post rew
    blackoutwin=.5;               %.5s away from event marker
    fixwin=1;
    intergap=2;
    targgap=1;
    iti1=1;
    iti3=3;
    iti5=5;
    eventmarks={};
  %  eventmarks{2}=[samplestarg{1}-targgap*param.samplerate; samplestarg{1}+targgap*param.samplerate;].*rates(1)./rates(2);
   % eventmarks{3}=[samplesfix{1}-param.samplerate*blackoutwin; samplestarg{1}].*rates(1)./rates(2);
  %  eventmarks{3}=[samplestarg{1}; samplestarg{1}+intergap*param.samplerate].*rates(1)./rates(2);
    %eventmarks{5}=[eventsrew; eventsrew+param.samplerate*rewwin].*rates(1)./rates(2);
    %eventmarks{1}=[samplesfix{1}-param.samplerate*fixwin; eventsrew-blackoutwin*param.samplerate].*rates(1)./rates(2);
    %eventmarks{1}=[samplesfix{1}-param.samplerate*fixwin; samplesfix{1}+param.samplerate*fixwin].*rates(1)./rates(2);
   % eventmarks{4}=[eventsrew+param.samplerate*1; eventsrew+param.samplerate*3].*rates(1)./rates(2);
   % eventmarks{5}=[eventsrew+param.samplerate*2; eventsrew+param.samplerate*4].*rates(1)./rates(2);
    %eventmarks{6}=[eventsrew+param.samplerate*4; eventsrew+param.samplerate*6].*rates(1)./rates(2);
   % eventmarks{3}=[samplestargeye{1}-targgap*param.samplerate; samplestargeye{1}+targgap*param.samplerate;].*rates(1)./rates(2);
    eventmarks{1}=samplesfix{1};
    eventmarks{2}=samplestarg{1};
    eventmarks{3}=samplestargeye{1};
        eventmarks{4}=samplesfixeye{1};

    %just provide middle of win rather than range
    %let xvardata/plotxvartrials determine win around align idx
    % eventnames{1}='prerew';
    eventnames{2}='intertarg';
   %eventnames{3}='fixtarg';
    %eventnames{3}='targwin';            %change 12/02/2018 so not variable across trials, but fixed 2 s
    %eventnames{5}='rewwin';
    eventnames{1}='interfix'; 
    %eventnames{4}='ititwo';
    %eventnames{5}='itithree';
    %eventnames{6}='itifive';
    eventnames{3}='intertargeye';       %1/10/19 added aligned to eye start
      eventnames{4}='interfixeye';       %1/10/19 added aligned to eye start
  
    if ~isempty(eventx)
        %user selected events to use
        %replace default eventnames and eventmarks
        chosevents=contains(eventnames,eventx);
        eventntemp=eventnames(chosevents);
        eventmtemp=eventmarks(chosevents);
        eventnames=eventntemp;
        eventmarks=eventmtemp;        
    end
    
    if isfield(plotparam,'shift')
        if strcmp(plotparam.shift,'fix')
            eventmarks={};
            eventnames={};
            eventmarks{1}=samplesfix{1};
           % eventmarks{1}=[samplesfix{1}-param.samplerate*fixwin; samplesfix{1}+param.samplerate*intergap].*rates(1)./rates(2);
            eventnames{1}='fix';   
        end
    end
    
    ich=targch;      %sort signals by this channel
    numtrials=size(trialbytrial(ich).da,1);
    chbadtrials=[];
    
    %open bad trials if exists in homedir
    if exist(manlist)>0
        %read bad trials from excel spreadsheet
        sheet = 1;      %sheet 1 default big
        if ~isempty(strfind(trialgroups{itype},'small'))
            sheet=2;    %small reward sheet is #2
        elseif ~isempty(strfind(trialgroups{itype},'target'))
            sheet=3;    %target error sheet
        elseif ~isempty(strfind(trialgroups{itype},'fix'))
            sheet=4;    %target error sheet
        end
        [~,sheetsav] = xlsfinfo(manlist); %find # sheets in worksheet
        if sheet<=length(sheetsav)
            xlsdata = xlsread(manlist,sheet);
            chbadtrials=xlsdata-100;        % 101 first trial
            chbadtrials=xlsdata-99;        %default, 100 first trial
        end
    end
    goodtrials=1:numtrials;         %default, all trials goodtrials
    alltrials=1:numtrials;
    if ~isempty(chbadtrials)        
        %use manual selected trials
        if size(chbadtrials,2)>1
            goodtrials=find(~ismember(alltrials,chbadtrials(:,ich))==1);
        else
            %only single ch in excel sheet
            goodtrials=find(~ismember(alltrials,chbadtrials)==1);
        end        
    end
    if any(samplesfixeye{1}<100) || any(samplesfixeye{1}<samplesfix{1})
        chbadtrials2=find(samplesfixeye{1}<100);
        chbadtrials2=unique([chbadtrials2 find(samplesfixeye{1}<samplesfix{1})]);
        goodtrials=goodtrials(find(~ismember(goodtrials,chbadtrials2)==1));
    end
    if any(samplestarg{1}<samplesfixeye{1}) || any(samplestarg{1}<100)
        chbadtrials2=find(samplestarg{1}<samplesfixeye{1});
        chbadtrials2=unique([chbadtrials2 find(samplestarg{1}<100)]);
        goodtrials=goodtrials(find(~ismember(goodtrials,chbadtrials2)==1));
    end
    if ~isfield(trialbytrial,'lfp') && (plotlfp || plotavglfp)
        %not merged file, get lfps separately based on trial num 100 etc.
        %from original xx_pro directory
        data={};
        %get reconverted path ie XX_pro directory
        pathid=strfind(currpath,'_pro');
        pathlfp=currpath(1:pathid+4);
        dirlfp=dir(pathlfp);
        %get csc signals defined in NCSchannels
        if isempty(param.NCSchannels)
            error('no ncs channels defined in param');
        end
        for it=1:numtrials
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
    
    %group trial types
    trialbytrial(1).samplesfixeye=samplesfixeye;
    trialbytrial(1).samplesfix=samplesfix;
    trialbytrial(1).samplestarg=samplestarg;
    trialbytrial(1).samplestargeye=samplestargeye;
    trialtypes=grouptrials(goodtrials,trialbytrial(1),plotparam);
    save([currpath 'trialtypes'], 'trialtypes');
    
pathda=[currpath 'trialbytrial_da' filesep];
pathlfp={};
pathlfp{1}=[currpath 'trialbytrial_lfp' filesep];
pathlfp{2}=[pathlfp{1}(1:end-1) num2str(2) filesep];
pathlfp{3}=[pathlfp{1}(1:end-1) num2str(3) filesep];

if ~isdir(pathlfp{1})
    mkdir(pathlfp{1})
end
mkdir(pathlfp{2});
mkdir(pathlfp{3});


    plottrials=goodtrials;
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
    betalfp={};
    betalfp2={};        %high beta band
    betalfp3={};        %broad beta
    plotparam.chnum=targch;
    
if getinfo
    %just save trialtypes and continue loop
    continue;
end
%%
%plot DA trial by trial

if plotda
    disp(['plotting da chs ' num2str(fscvchs)]);
fig1=figure('visible','off');
cscales=[-3 8; -3 8; -3 8; -3 8]; 
eyedscale=[-2.5 .5].*1e-3;   
lickscale=[1.3e-4 4.8e-4];
hfig=setupTrialPlots(fig1, plotparam,numcolor,numbehav);
pathda=[currpath 'trialbytrial_da' filesep];
if ~isdir(pathda)
    mkdir(pathda)
end
for it=1:length(trialtypes.names)
    plottrials=trialtypes.nums{it};
    trialname=trialtypes.names{it};
    plotparam.triallabel=trialname;
for ii=1:length(fscvchs)
    ich=fscvchs(ii);    
    plotparam.chnum=ich;
    if ii==1
        [datm{ich},da{ich}]=setTrialAx(trialbytrial(ich),plotparam,hfig{ii},plottrials,...
            samplesfix{1},samplestarg{1},samplesfixeye{1},'plotnum',1,'sitename',sites(ich));
        if ~isempty(cscales)
        [datm{ich},da{ich}]=setTrialAx(trialbytrial(ich),plotparam,hfig{ii},plottrials,...
            samplesfix{1},samplestarg{1},samplesfixeye{1},'plotnum',1,'cscale',cscales(ii,:),'sitename',sites(ich));
        end
    else
        [datm{ich},da{ich}]=setTrialAx(trialbytrial(ich),plotparam,hfig{ii},plottrials,...
            samplesfix{1},samplestarg{1},samplesfixeye{1},'sitename',sites(ich));
        if ~isempty(cscales)
        [datm{ich},da{ich}]=setTrialAx(trialbytrial(ich),plotparam,hfig{ii},plottrials,...
            samplesfix{1},samplestarg{1},samplesfixeye{1},'cscale',cscales(ii,:),'sitename',sites(ich));
        end
    end
    if it==1        
        %save trial by trial data with reference ts's for nlx to sort
        %with small rew trials converted at another time
        %dadata=da{ich};
        if ii==length(fscvchs)
            alignts=[];
            if isfield(trialbytrial(1),'alignts')
            alignts=trialbytrial(1).alignts;
            end
            trialids=plottrials;
            save([pathda trialname '-data-ch' num2str(ich)],'da',...
                'trialids','alignts','trialtypes');
            %when reusing data, for specific types of trials 
            %need to index by find(ismember(trialids,phase2trials)==1)
        end
    end
end
[eyetm,eyedata]=setTrialAx(trialbytrial(1),plotparam,hfig{length(fscvchs)+1},...
    plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},'eyed',...
    'cscale',eyedscale);
[licktm,lickdata]=setTrialAx(trialbytrial(1),plotparam,hfig{length(fscvchs)+2},...
    plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},'lickx',...
    'cscale',lickscale);
rts=setTrialBehavAx(trialbytrial(1),plotparam,hfig,plottrials,length(fscvchs)+3);
[blinktm,blinkdata]=setTrialAx(trialbytrial(1),plotparam,[],...
    plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},'blink','noplot');
[eyextm,eyexdata]=setTrialAx(trialbytrial(1),plotparam,[],...
    plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},'eyex','noplot');
%calc pulse rate
pulsetm=[];
pulsedata=[];
if sessnum>37
[pulsetm,pulsedata]=setTrialAx(trialbytrial(1),plotparam,[],...
    plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},'pulse','noplot');
end
if it==1
    save([pathda trialname '-data-beh'],'eyedata','eyexdata','lickdata',...
    'blinkdata','pulsedata');
end
savefig(fig1,[pathda trialname '.fig']);
saveas(fig1,[pathda trialname '.tif'],'tif')
save([pathda trialname],'datm','eyetm','licktm',...
    'blinktm','pulsetm','eyextm','rts');
end
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
hfig2=setupTrialPlots(fig2, plotparam,length(lfpchs),0);
for ifreq=1:length(selbands)
    fbands=allbands{selbands(ifreq)};
    plotparam.fbands=fbands;
    for it=1:length(trialtypes.names)
        plottrials=trialtypes.nums{it};
        trialname=trialtypes.names{it};
        plotparam.triallabel=trialname;
        for ii=1:length(lfpchs)
            if it==1
            [betatm{ii},downlfp]=setTrialAx(trialbytrial(1),...
                plotparam,hfig2{ii},plottrials,samplesfix{1},samplestarg{1},...
                samplesfixeye{1},lfpchs{ii},'fbands',fbands,'getfft',...
                'tapers',[3 5],'smoothwin',0.3,'plotnum',ii);
                if globeta==0
                    if selbands(ifreq)==2
                        betalfp2{ii}=downlfp;
                    end
                    if selbands(ifreq)==1
                        betalfp{ii}=downlfp;
                    end
                    if selbands(ifreq)==3
                        betalfp3{ii}=downlfp;
                    end
                else
                    %only broad beta
                    betalfp3{ii}=downlfp;
                end
                if ii==length(lfpchs)
                    %save betalfp data
                    alignts=[];
                    if isfield(trialbytrial(1),'alignts')
                        alignts=trialbytrial(1).alignts;
                    end
                    trialids=plottrials;
                    if globeta==0
                        if selbands(ifreq)==1
                        save([pathlfp{1} trialname '-betalfp'],'betalfp',...
                            'trialids','alignts','trialtypes');
                        end
                        if selbands(ifreq)==2
                        save([pathlfp{2} trialname '-betalfp2'],'betalfp2',...
                            'trialids','alignts','trialtypes');
                        end
                        if selbands(ifreq)==3
                        save([pathlfp{3} trialname '-broadbetalfp'],'betalfp3',...
                            'trialids','alignts','trialtypes');
                        end
                    end
                    if globeta
                        %only broad beta
                         save([pathlfp{3} trialname '-broadbetalfp'],'betalfp3',...
                        'trialids','alignts','trialtypes');
                    end
                end

            else
                rclfp=[];
                if globeta==0
                    if selbands(ifreq)==1
                        rclfp=betalfp{ii};
                    end
                    if selbands(ifreq)==2
                        rclfp=betalfp2{ii};
                    end
                    if selbands(ifreq)==3
                        rclfp=betalfp3{ii};
                    end
                end
                if globeta
                    rclfp=betalfp3{ii};
                end
                  
            %already got all betalfp data in first 'all trial' it
            betatm{ii}=setTrialAx(trialbytrial(1),...
                plotparam,hfig2{ii},plottrials,samplesfix{1},samplestarg{1},...
                samplesfixeye{1},lfpchs{ii},'rclfp',rclfp,'fbands',fbands,'getfft',...
                'tapers',[3 5],'smoothwin',0.3,'plotnum',ii);
            end
        end
        if globeta==0
            %savefig(fig2,[pathlfp{selbands(ifreq)} trialname]);
            saveas(fig2,[pathlfp{selbands(ifreq)} trialname '.eps'],'epsc')
            saveas(fig2,[pathlfp{selbands(ifreq)} trialname '.tif'],'tif')
            save([pathlfp{selbands(ifreq)} trialname],'betatm');
        else
            %broad beta saved in selbands(ifreq) 3
            %savefig(fig2,[pathlfp{3} trialname]);
            saveas(fig2,[pathlfp{3} trialname '.eps'],'epsc')
            saveas(fig2,[pathlfp{3} trialname '.tif'],'tif')
            save([pathlfp{3} trialname],'betatm');
        end
   
    end
end
end
%%
%plot averages da
if plotavgda
figsigrep=figure('visible','off'); 
set(figsigrep, 'Color', [1 1 1],'position',[100 50 1000 950]);
set(0,'CurrentFigure',figsigrep); 
axsigrep={};
cscaleda=[-2 4];
cscaleda=[-4 4];
for ii=1:8
    axsigrep{ii}=subplot(4,2,ii);
end
pathavg=[currpath 'avg_da' filesep];
if ~isdir(pathavg)
    mkdir(pathavg)
end
for it=1:length(trialtypes.names)
    plottrials=trialtypes.nums{it};
    trialname=trialtypes.names{it};
    plotparam.triallabel=trialname;
    countp=0;
    for ii=1:length(fscvchs)
        countp=countp+1;
        ich=fscvchs(ii);
        plotparam.chnum=ich;
        setTrialAxAvg(trialbytrial(ich),plotparam,axsigrep{countp},plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},'sitename',sites(ich),'plotnum',countp);
        countp=countp+1;
        setTrialAxAvg(trialbytrial(ich),plotparam,axsigrep{countp},plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},'avg','cscale',cscaleda,'sitename',sites(ich),'plotnum',countp);
    end
    if it==1
        savefig(figsigrep,[pathavg trialname]);
    end
    saveas(figsigrep,[pathavg trialname '.tif'],'tif')
end
end
%%
%plot averages phys
if plotavglfp
figsigrep2=figure('visible','off'); 
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
pathavglfp{1}=[currpath 'avg_lfp' filesep];
pathavglfp{2}=[pathavglfp{1}(1:end-1) '2' filesep];
pathavglfp{3}=[pathavglfp{1}(1:end-1) '3' filesep];

if ~isdir(pathavglfp{1})
    mkdir(pathavglfp{1})
end
    %store high beta separate folder 2
    mkdir(pathavglfp{2});

    %store broad beta separate folder 3
    mkdir(pathavglfp{3});

for ifreq=1:length(selbands)
    fbands=allbands{selbands(ifreq)};
    plotparam.fbands=fbands;
    for it=1:length(trialtypes.names)
        plottrials=trialtypes.nums{it};
        trialname=trialtypes.names{it};
        plotparam.triallabel=trialname;
        if isempty(betalfp)
            %need to load
            if globeta==0
                load([currpath 'trialbytrial_lfp' filesep 'all-betalfp'],'betalfp')
                if ismember(2,selbands)>=2
                    load([currpath 'trialbytrial_lfp2' filesep 'all-betalfp2'],'betalfp2')
                end
            end
            if ismember(3,selbands)==3 || globeta
                load([currpath 'trialbytrial_lfp3' filesep 'all-broadbetalfp'],'betalfp3')
            end
        end
        for ii=1:length(lfpchs)
            if globeta==0
                downlfp=betalfp{ii};
                if selbands(ifreq)==2
                    downlfp=betalfp2{ii};
                end
            end
            if selbands(ifreq)==3 || globeta
                downlfp=betalfp3{ii};
            end
            setTrialAxAvg(trialbytrial(1),plotparam,axsigrep2{ii*2-1},...
                plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},...
                lfpchs{ii},'fbands',fbands,'rclfp',downlfp);
            setTrialAxAvg(trialbytrial(1),plotparam,axsigrep2{ii*2},...
                plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},...
                lfpchs{ii},'avg','fbands',fbands,'rclfp',downlfp);
        end
        if globeta==0
           % savefig(figsigrep2,[pathavglfp{selbands(ifreq)} trialname]);
            saveas(figsigrep2,[pathavglfp{selbands(ifreq)} trialname '.tif'],'tif')
        else
           % savefig(figsigrep2,[pathavglfp{3} trialname]);
            saveas(figsigrep2,[pathavglfp{3} trialname '.tif'],'tif')
        end
    end
end
end
%%
%plot averages physiological response
if plotphys
figsigrep3=figure('visible','off'); 
set(figsigrep3, 'Color', [1 1 1],'position',[100 100 1000 1000]);
if ~ispc
set(figsigrep3, 'Color', [1 1 1],'position',[100 100 1000 750]);
end
set(0,'CurrentFigure',figsigrep3); 
axsigrep3={};
for ii=1:8
    axsigrep3{ii}=subplot(4,2,ii);
end
pathavgb=[currpath 'avg_beh' filesep];
if ~isdir(pathavgb)
    mkdir(pathavgb)
end
for it=1:length(trialtypes.names)
    plottrials=trialtypes.nums{it};
    trialname=trialtypes.names{it};
    plotparam.triallabel=trialname;
    setTrialAxAvg(trialbytrial(1),plotparam,axsigrep3{1},plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},'eyed');
    setTrialAxAvg(trialbytrial(1),plotparam,axsigrep3{2},plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},'eyed','avg');
    setTrialAxAvg(trialbytrial(1),plotparam,axsigrep3{3},plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},'lickx');
    setTrialAxAvg(trialbytrial(1),plotparam,axsigrep3{4},plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},'lickx','avg');
    setTrialAxAvg(trialbytrial(1),plotparam,axsigrep3{5},plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},'pulse');
    setTrialAxAvg(trialbytrial(1),plotparam,axsigrep3{6},plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},'pulse','avg');
    setTrialAxAvg(trialbytrial(1),plotparam,axsigrep3{7},plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},'blink');
    setTrialAxAvg(trialbytrial(1),plotparam,axsigrep3{8},plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},'blink','avg');
    %savefig(figsigrep3,[pathavgb trialname])
    saveas(figsigrep3,[pathavgb trialname '.tif'],'tif')
end
%reaction times
binwidth=25;              %nm
binmax=600;
pathrts=[currpath 'rts' filesep];
if ~isdir(pathrts)
    mkdir(pathrts)
end
for it=1:length(trialtypes.names)
    trialname=trialtypes.names{it};
    load([pathda trialname],'rts');
    [binsrts,bb]=histc(rts.left.*1000,[0:binwidth:binmax]);
    [binsrtsright,bb]=histc(rts.right.*1000,[0:binwidth:binmax]);
    fighist=figure('visible','off');
    set(0,'CurrentFigure',fighist); 
    xlims=[0 binmax];
    subhist(1)=subplot(4,2,1);
    set(fighist, 'Color', [1 1 1]);
    set(fighist,'Position',[50,50,700,250]);
    countsub=1;
    binedges=0:binwidth:binmax;
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
%%
%print average FFT spectrograms
if plotfft
pathffts=[currpath 'avgffts' filesep];
if ~isdir(pathffts)
    mkdir(pathffts)
end
fs=1000;        %sampler ate nlx
figffts=figure('visible','off');
set(figffts, 'Color', [1 1 1]);
set(figffts,'Position',[300,150,1200,500]);
if ~ispc
set(figffts, 'position',[0 0 1000 500]);
end

axfft=axes;
hold(axfft,'on'); 
figff=figure('visible','off');
set(figff, 'Color', [1 1 1]);
set(figff,'Position',[50,50,700,700]);
axff=axes;
hold(axff,'on');
fftlims=[5 50];
plotparam.ffttapers=[3 5];      %nw (time-bandiwth rpdocut, tapers)
%plotparam.ffttapers=[1.8 1];      %nw (time-bandiwth rpdocut, tapers)
plotparam.fftwin=[.75 0.15];
for it=1:length(trialtypes.names)
    plottrials=trialtypes.nums{it};
    trialname=trialtypes.names{it};
    plotparam.triallabel=trialname;
    for ii=1:length(lfpchs)
    cla(axfft);
    cla(axff);
    cscid=find(ismember(trialbytrial(1).cscNames,lfpchs{ii}));
    if length(cscid)>1
        cscid=cscid(1);
    end
    evffts{ii}=setavgfft(trialbytrial(1).lfp{cscid},fs,plotparam,axfft,plottrials,...
        'events','fix',samplesfix{1},...
        'events','target',samplestarg{1},...
        'events','fixeye',samplesfixeye{1},...
        'alignrelavgwin',{'outcome','target','target','target','fix'},...
        [1 0.5; 2 .5; .3 .2; -0.3 0.2;-1 .25],...
        'freq',fftlims,'sitename',lfpchs{ii},'fax',axff);

    save([pathffts trialname],'evffts','plotparam');
    pathsfftchs{ii}=[pathffts lfpchs{ii} filesep];
    mkdir(pathsfftchs{ii})
    
    %saveas(figffts,[pathsfftchs{ii} trialname '.eps'],'epsc')
    saveas(figffts,[pathsfftchs{ii} trialname '.jpg'],'jpeg')
    savefig(figffts,[pathsfftchs{ii} trialname])
    evffts{ii}=setavgfft(trialbytrial(1).lfp{cscid},fs,plotparam,axfft,plottrials,...
        'events','fix',samplesfix{1},...
        'events','target',samplestarg{1},...
        'events','fixeye',samplesfixeye{1},...    
        'freq',fftlims,'sitename',lfpchs{ii},'align','target');
    saveas(figffts,[pathsfftchs{ii} trialname '-targaln.jpg'],'jpeg')
    savefig(figffts,[pathsfftchs{ii} trialname '-targaln'])
    end
end
end
%%
%calculate trial by trial correlation matrix
if plotcor
cortypes={'fixwin','fixpeak','targwin','targpeakvprefix','targpeak',...
    'targpeakdiff','targimwin','targimpeak','targimpeakdiff',...
    'rewprewin','rewprepeak','rewprepeakdiff','rewimpeak',...
    'rewimpeakdiff','rewimpeakvprefix','rewshortpeakvprefix',...
    'rewshortwin','rewshortpeak','rewshortpeakdiff','rewwin',...
    'rewpeak','rewpeakvprefix','rewpeakdiff'};
if isfield(plotparam,'shift')
    if strcmp(plotparam.shift,'fix')
        cortypes={'fixwin','fixpeak','rewimpeak'};  
    end
end
pathcor={};
pathcor{1}=[currpath 'corr' filesep];
pathcor{2}=[currpath 'corr2' filesep];      %HF beta
pathcor{3}=[currpath 'corr3' filesep];      %broad beta

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
    for it=1:length(trialtypes.names)
        trialname=trialtypes.names{it};
        if globeta==0
            if selbands(ifreq)==1
                load([pathlfp{1} trialname],'betatm');
            end
            if selbands(ifreq)==2
                load([pathlfp{2} trialname],'betatm');
            end   
        end
        if selbands(ifreq)==3 || globeta
            load([pathlfp{3} trialname],'betatm');
        end
        load([pathda trialname],'datm');
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
        if it==1 && selbands(ifreq)==1
            %get behavior correlations for da, just for "all" trials
            load([pathda trialname]);
            behtm={eyetm,pulsetm,licktm,blinktm};
            corrmat=cormatrix(datm,behtm,cortypes,'sitelabels','sametype');
            save([pathsave trialname '-beh-sametype'],'corrmat');
        end
    end
end

%%
for ifreq=1:length(selbands)
    fbands=allbands{selbands(ifreq)};
    plotparam.fbands=fbands;
    pathsave={};
    for it=1:length(trialtypes.names)
        pathsave=pathcor{selbands(ifreq)};
        if globeta
            pathsave=pathcor{3};
        end
        trialname=trialtypes.names{it};
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
                    ranjitter=rand(1,length(lfpchs))*.3-.15;

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
                yy1=15;
                for ilfp=1:length(lfpchs)
                    c=sitecolors(ilfp,:);
                    text(axcorr{ii},range(end)+10,yy1,lfpchs{ilfp},'fontsize',10,'color',c);
                    yy1=yy1-2;
                end
            end
            hold(axcorr{ii},'off');
        end
        savefig(figcorr,[pathsave trialname]);
        saveas(figcorr,[pathsave trialname],'tif')
    end
end

%plot figure behavioral correlations to da
fbands=allbands{1};
plotparam.fbands=fbands;
pathsave=pathcor{1};
trialname=trialtypes.names{1};
load([pathsave trialname '-beh-sametype']);
cor=corrmat;
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
plottrials=trialtypes.nums{1};      %all good trials
rsel=intersect(plottrials,rside);
lsel=intersect(plottrials,lside);
target_rrt=trialbytrial(1).target_rt(rsel);
target_lrt=trialbytrial(1).target_rt(lsel);
fix_rt=trialbytrial(1).fix_rt(plottrials);
trialname=trialtypes.names{1};
load([pathda trialname]);      %load datm
rs={}; ps={}; rts={};
for ich=1:length(fscvchs)
[rs(fscvchs(ich)).targwinr,ps(fscvchs(ich)).targwinr]=corr(datm{fscvchs(ich)}.targwin(find(ismember(datm{fscvchs(ich)}.trialnums,rsel)==1))',target_rrt','rows','complete');
[rs(fscvchs(ich)).targimwinr,ps(fscvchs(ich)).targimwinr]=corr(datm{fscvchs(ich)}.targimwin(find(ismember(datm{fscvchs(ich)}.trialnums,rsel)==1))',target_rrt','rows','complete');
[rs(fscvchs(ich)).targpeakr,ps(fscvchs(ich)).targpeakr]=corr(datm{fscvchs(ich)}.targpeak(find(ismember(datm{fscvchs(ich)}.trialnums,rsel)==1))',target_rrt','rows','complete');
[rs(fscvchs(ich)).targimpeakr,ps(fscvchs(ich)).targimpeakr]=corr(datm{fscvchs(ich)}.targimpeak(find(ismember(datm{fscvchs(ich)}.trialnums,rsel)==1))',target_rrt','rows','complete');
[rs(fscvchs(ich)).targwinl,ps(fscvchs(ich)).targwinl]=corr(datm{fscvchs(ich)}.targwin(find(ismember(datm{fscvchs(ich)}.trialnums,lsel)==1))',target_lrt','rows','complete');
[rs(fscvchs(ich)).targimwinl,ps(fscvchs(ich)).targimwinl]=corr(datm{fscvchs(ich)}.targimwin(find(ismember(datm{fscvchs(ich)}.trialnums,lsel)==1))',target_lrt','rows','complete');
[rs(fscvchs(ich)).targpeakl,ps(fscvchs(ich)).targpeakl]=corr(datm{fscvchs(ich)}.targpeak(find(ismember(datm{fscvchs(ich)}.trialnums,lsel)==1))',target_lrt','rows','complete');
[rs(fscvchs(ich)).targimpeakl,ps(fscvchs(ich)).targimpeakl]=corr(datm{fscvchs(ich)}.targimpeak(find(ismember(datm{fscvchs(ich)}.trialnums,lsel)==1))',target_lrt','rows','complete');
[rs(fscvchs(ich)).fixwin,ps(fscvchs(ich)).fixwin]=corr(datm{fscvchs(ich)}.fixwin(find(ismember(datm{fscvchs(ich)}.trialnums,lsel)==1))',target_lrt','rows','complete');
[rs(fscvchs(ich)).fixpeak,ps(fscvchs(ich)).fixpeak]=corr(datm{fscvchs(ich)}.fixpeak(find(ismember(datm{fscvchs(ich)}.trialnums,lsel)==1))',target_lrt','rows','complete');
end
rts.rs=rs;
rts.ps=ps;
save([pathsave 'rts_corr'],'rts');


range=1:length(cortypes);
xticks=range;
xlabels=cortypes;
cmapp=[1 .3 0;.4 1 0; .3 .6 .2; 0 .65 1; 1 .7 0; 1 .75 0; 1 .8 0; 1 .85 0; 1 .9 0];

for ii=1:length(fscvchs)
    btypes=unique({cor{fscvchs(ii)}.corrdata(1:end).d2site});
    btypes=btypes(~contains(btypes,'na'));
cmapp=cool; 
    numcol=length(btypes(~contains(btypes,'na')));
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
        yy1=15;
        for ilfp=1:length(btypes)
            c=sitecolors(ilfp,:);
            text(axcorr{ii},range(end)+10,yy1,btypes{ilfp},'fontsize',10,'color',c);
            yy1=yy1-2;
        end
    end
    hold(axcorr{ii},'off');
end
savefig(figcorr,[pathsave trialname '-beh']);
saveas(figcorr,[pathsave trialname '-beh'],'tif')
   
end
%%
%get cross covariance just selected trials above
if plotxvar
    
plotparam.xrates=rates;
if ismember(1,selbands) && isempty(betalfp)
    load([currpath 'trialbytrial_lfp' filesep 'all-betalfp'],'betalfp')
end
if ismember(2,selbands) && isempty(betalfp2)
    load([currpath 'trialbytrial_lfp2' filesep 'all-betalfp2'],'betalfp2')
end
if ismember(3,selbands) || globeta
    load([currpath 'trialbytrial_lfp3' filesep 'all-broadbetalfp'],'betalfp3')
end
if isempty(da)
    %need to load
    load([currpath 'trialbytrial_da' filesep 'all-data-ch' num2str(max(fscvchs))],'da')
end
datax.lfp=betalfp;
if globeta
    %broadband beta only
    datax.lfp=betalfp3;
end
datax.da=da;
%for idagroup=1:length(xchs)
%plotparam.chnums=xchs{idagroup};
%for ifreq=1:length(allbands)
%burst align event is the middle of evch event window
burstaln={};
for ievent=1:length(eventmarks)
    %NO NOT MIDDLE OF WIN NECESSARY ALIGN IDX
   % burstaln{ievent}=mean(eventmarks{ievent},1);
   %changed 1/10/19 just align idx
        burstaln{ievent}=eventmarks{ievent};

end
for ifreq=1:length(selbands)
    fbands=allbands{selbands(ifreq)};
    plotparam.fbands=fbands;
    if selbands(ifreq)==2
        datax.lfp=betalfp2;
    end
    if selbands(ifreq)==3
        datax.lfp=betalfp3;
    end
    if plotxvarall
        for it=1
            plottrials=trialtypes.nums{it};
            trialname=trialtypes.names{it};
            trialname=[];
            if selbands(ifreq)==2
                trialname=[trialname '2'];
            end
            if selbands(ifreq)==3 || globeta
                trialname=[trialname '3'];
            end
            plotparam.triallabel=trialname;
            xcovdata=plotxvartrials(datax,eventmarks,eventnames,plotparam,'type','none',...
            'trials',plottrials,'label',trialname,'cluster');
            %xcovdata=plotxvartrials(datax,eventmarks,eventnames,plotparam,'type','none',...
           % 'trials',plottrials,'label',trialname,'sortxvar');
        %split trial types as needed post-pro
        end  
    end
    %{
    bursttrialtypes=contains(trialtypes.names,{'all','phase1','phase4'});
    btrialnames=trialtypes.names(bursttrialtypes);
    btrialnums=trialtypes.nums(bursttrialtypes);
    for it=1
        %just do all, select phases/groups of trials post later..more effi
        %burst trial types & da increases/decreases groups after event
        plottrials=btrialnums{it};
        trialname=btrialnames{it};
        if selbands(ifreq)==2
            trialname=[trialname '2'];
        end
        if selbands(ifreq)==3 || globeta
            trialname=[trialname '3'];
        end
        plotparam.triallabel=trialname;
        if plotbetaxcov
        xcovdata=plotxvartrials(datax,eventmarks,eventnames,plotparam,'type','none',...
        'trials',plottrials,'label',trialname,'burstgroup',burstaln,[]);
        end
        if plotdaxcov
        xcovdata=plotxvartrials(datax,eventmarks,eventnames,plotparam,'type','none',...
        'trials',plottrials,'label',trialname,'dasignaltypes');
        end
    end
    %}
end


end
%end session types
end
%end if plotxclust==0
end
if plotxclust==1
%%
%post processing all sessions comparison of xvar's
typecount=0;
bcount=0;   %behavior data types
targevents={'interfix','intertarg','intertargeye'};
fscvchstemp=fscvchs;        %store because when load trialbytrial, loads all again replacing fscvchs
%get lfp ch groups p & c
pgroup=find(contains(plotparam.lfpchs,'p')==1);
cgroup=find(contains(plotparam.lfpchs,'c')==1);
%currpath=getfield(targpath,trialgroups{1});



%group trial types
xinfo=[];
%numb=length(trialgroups)*length(targevents);
%binfo(numb)=struct();
binfo=[];
trialinfo=[];
for ifreq=1:length(selbands)
    fbands=allbands{selbands(ifreq)};     %only low band modulated seems
    plotparam.fbands=fbands;  
    xinfo=[];

for itype=1:3
    %big, small, targbreak
currpath=getfield(targpath,trialgroups{itype});
load([currpath 'trialtypes']);
disp(currpath)
%get behavioral properties trialgroup
compiledname=[currpath 'trials_chronic' sessid '_' trialgroups{itype} '_pro'];
load(compiledname);
%load eyedata,lickdata,blinkdata,pulsedata
behname=[currpath 'trialbytrial_da' filesep 'all-data-beh.mat'];
load(behname);
trialbytrial(1).eyedata=eyedata;
trialbytrial(1).lickdata=lickdata;
trialbytrial(1).pulsedata=pulsedata;
plotparam.pathname=currpath;
plotparam.trialtype=trialgroups{itype};

if ~isfield(trialbytrial(1),'eyexdata')
    %previously not saved, get
    goodtrials=trialtypes.nums{1};
    chids=param.NCSchannels;
    selpath=currpath;
    eyedatatemp=gettrialdata(chids,currpath,'selch',34);
    eyedatatemp.ratelfp=trialbytrial(1).ratelfp;
    trialbytrial(1).lfp=eyedatatemp.lfp;
    trialbytrial(1).cscNames=plotparam.cscNames;
    eyedatatemp.cscNames={'eyex'};
    eyedatatemp.relts=trialbytrial(1).relts;
    eyedatatemp.tscsc=trialbytrial(1).tscsc;
    eyedatatemp.alignts=trialbytrial(1).alignts;
    [eyextm,eyexdata]=setTrialAx(eyedatatemp,plotparam,[],...
        goodtrials,samplesfix{1},samplestarg{1},samplesfixeye{1},'eyex','noplot');   
    trialbytrial(1).eyexdata=eyexdata;
end
if ~isdir(currpath)
    %files/folder not created by sync sigs, skip
    disp(['folder not created, skipping : ' char(10) currpath]);
    continue
end
fscvchs=fscvchstemp;
plotparam.pathname=currpath;
plotparam.trialtype=trialgroups{itype};
if ~isdir(currpath)
    %files/folder not created by sync sigs, skip
    disp(['folder not created, skipping : ' char(10) currpath]);
    continue
end
savepath=[currpath(1:strfind(currpath,'matlab')+6) 'xsess' filesep];
if ~isdir(savepath)
    mkdir(savepath);
end
if selbands(ifreq)==2
    savepath=[currpath(1:strfind(currpath,'matlab')+6) 'xsess2' filesep];
    mkdir(savepath);
end
if selbands(ifreq)==3 || globeta
    savepath=[currpath(1:strfind(currpath,'matlab')+6) 'xsess3' filesep];
    mkdir(savepath);
end

    typecount=typecount+1;    
    vartype='none';
    plotparam.chnums=fscvchs;    
        addname=[];
        if selbands(ifreq)==2
            addname=[addname '2'];
        end
        if selbands(ifreq)==3 || globeta
            addname=[addname '3'];
        end
        plotparam.triallabel=addname;
        for ii=1:length(fscvchs)
            ich=fscvchs(ii);    
            for ievent=1:length(targevents)
                pathxvar=[currpath 'xvar_' targevents{ievent} '_' ...
                    vartype addname filesep];
                load([pathxvar 'xvar_' sites{ich}],'xcovdata'); 
                %organize behavioral parameters for trial group (big/small/etc)
                %for given event (interfix/targ/etc.)
                if ii==1 && ifreq==1
                    bcount=bcount+1;
                    seltrials=xcovdata{1}.seltrials;
                    winids=xcovdata{1}.winids;
                    binfotemp=orgbeh(trialbytrial(1),plotparam,...
                        targevents{ievent},seltrials,winids,'win',[-2 2]);
                    binfotemp.sessiontype=trialgroups{itype};
                    binfotemp.event=targevents{ievent};
                    binfo=[binfo binfotemp];
                    %get trial types again
                    trialtypes={};      %reset, redefined below in loop
                    trialtypes=grouptrials(seltrials,trialbytrial(1),plotparam);
                    save([currpath 'trialtypes'], 'trialtypes');
                    plotparam.trialtypes=trialtypes;    %needed for trialinfo /xplot
                end
                
                xinfotemp=xclust(xcovdata,plotparam,'event',targevents{ievent});                
                %xinfotemp=xclust(xcovdata,plotparam);
                xinfo=[xinfo xinfotemp];        %concatenate all info              

            end
        end
     
        %save trial ids for other trial classifications in group
        %same across all event types
        if ifreq==1
        %(phase1/4..etc)
            trialinfo(itype).trialtypes=plotparam.trialtypes;
            trialinfo(itype).sessiontype=trialgroups{itype};
        end
        
end
save([savepath 'xinfo' ], 'xinfo','binfo','trialinfo');
plotparam.savepath=savepath;        %needed for xplot

end
    
end

if plotxsess
    %plot multiple types of sess (big,small/targ...et)
%post processing all sessions comparison after all data saved from xclust
currpath=getfield(targpath,trialgroups{1});
load([currpath 'trialtypes'], 'trialtypes');
loadpath=[currpath(1:strfind(currpath,'matlab')+6) 'xsess3' filesep];
load([loadpath 'xinfo' ], 'xinfo','binfo','trialinfo');
plotparam.trialtypes=trialtypes;  
targevents={'interfix','intertarg','intertargeye'};
datypes={'dapos','daall','daneg'};
plotparam.savepath=loadpath;        %needed for xplot

for datype=1:length(datypes)
    
    %temporary comment, already finished for some sessions, dnr
    plotxtraces(xinfo,trialinfo,plotparam,datypes{datype});
    plotx(xinfo,trialinfo,plotparam,datypes{datype});
    plotxallsess(xinfo,trialinfo,plotparam,datypes{datype});
    %plot subplots big vs small vs targ break, averages
    plotxtracesavg(xinfo,trialinfo,plotparam,datypes{datype});
    
    %plot averages overlay comparison diff conditions
    plotxtracesavgov(xinfo,trialinfo,plotparam,'evtype','intertarg','sesstypes',...
        {'bigreward'},'trialtypes',{'sleft','sright'},datypes{datype});
    plotxtracesavgov(xinfo,trialinfo,plotparam,'evtype','intertarg','sesstypes',...
        {'bigreward','smallreward'},'trialtypes',{'sleft'},datypes{datype});
    plotxtracesavgov(xinfo,trialinfo,plotparam,'evtype','intertarg','sesstypes',...
        {'bigreward','smallreward'},'trialtypes',{'sright'},datypes{datype});
     plotxtracesavgov(xinfo,trialinfo,plotparam,'evtype','intertarg','sesstypes',...
        {'smallreward'},'trialtypes',{'sleft','sright'},datypes{datype});
    
     plotxtracesavgov(xinfo,trialinfo,plotparam,'evtype','intertargeye','sesstypes',...
        {'bigreward'},'trialtypes',{'sleft','sright'},datypes{datype});
    plotxtracesavgov(xinfo,trialinfo,plotparam,'evtype','intertargeye','sesstypes',...
        {'bigreward','smallreward'},'trialtypes',{'sleft'},datypes{datype});
    plotxtracesavgov(xinfo,trialinfo,plotparam,'evtype','intertargeye','sesstypes',...
        {'bigreward','smallreward'},'trialtypes',{'sright'},datypes{datype});
    plotxtracesavgov(xinfo,trialinfo,plotparam,'evtype','intertargeye','sesstypes',...
        {'smallreward'},'trialtypes',{'sleft','sright'},datypes{datype});
    
    plotxtracesavgov(xinfo,trialinfo,plotparam,'evtype','intertarg','sesstypes',...
        {'bigreward','smallreward'},datypes{datype});
    plotxtracesavgov(xinfo,trialinfo,plotparam,'evtype','interfix','sesstypes',...
        {'bigreward','targetbreak'},datypes{datype});
    %plot eye movements wi da
    %NEED WIDER WINDOW SAVED FOR EYE TRACES!
    
    plotxtracesavgov(xinfo,trialinfo,plotparam,'evtype','intertarg','sesstypes',...
        {'bigreward','smallreward'},'trialtypes',{'sright'},'binfo',binfo,...
        datypes{datype});
    plotxtracesavgov(xinfo,trialinfo,plotparam,'evtype','intertarg','sesstypes',...
        {'bigreward','smallreward'},'trialtypes',{'sleft'},'binfo',binfo,...
        datypes{datype});
    plotxtracesavgov(xinfo,trialinfo,plotparam,'evtype','interfix','sesstypes',...
        {'bigreward','targetbreak'},'binfo',binfo,...
        'win',[-1 1],datypes{datype});
    plotxtracesavgov(xinfo,trialinfo,plotparam,'evtype','intertarg','sesstypes',...
        {'bigreward'},'trialtypes',{'sright','sleft'},'binfo',binfo,...
        datypes{datype});
    plotxtracesavgov(xinfo,trialinfo,plotparam,'evtype','intertarg','sesstypes',...
        {'smallreward'},'trialtypes',{'sright','sleft'},'binfo',binfo,...
        datypes{datype});

       plotxtracesavgov(xinfo,trialinfo,plotparam,'evtype','intertargeye','sesstypes',...
        {'bigreward','smallreward'},'trialtypes',{'sright'},'binfo',binfo,...
        datypes{datype});
    plotxtracesavgov(xinfo,trialinfo,plotparam,'evtype','intertargeye','sesstypes',...
        {'bigreward','smallreward'},'trialtypes',{'sleft'},'binfo',binfo,...
        datypes{datype});
    plotxtracesavgov(xinfo,trialinfo,plotparam,'evtype','intertargeye','sesstypes',...
        {'bigreward'},'trialtypes',{'sright','sleft'},'binfo',binfo,...
        datypes{datype});
    plotxtracesavgov(xinfo,trialinfo,plotparam,'evtype','intertargeye','sesstypes',...
        {'smallreward'},'trialtypes',{'sright','sleft'},'binfo',binfo,...
        datypes{datype}); 
end

end

if plotsort
    %plot signals sorted side bys ide
currpath=getfield(targpath,trialgroups{1});
load([currpath 'trialtypes'], 'trialtypes');
loadpath=[currpath(1:strfind(currpath,'matlab')+6) 'xsess3' filesep];
load([loadpath 'xinfo' ], 'xinfo','binfo','trialinfo');
plotparam.trialtypes=trialtypes;  
targevents={'interfix','intertarg','intertargeye'};
datypes={'daall'};
plotparam.savepath=loadpath;        %needed for xplot

for datype=1:length(datypes)
        plotxtraces(xinfo,trialinfo,plotparam,datypes{datype},'sortda');
end
end

end
