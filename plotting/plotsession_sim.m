function plotsession_sim(sessnum,varargin)
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
fbands2={};
%plot flags
plotda=0;
plotlfp=0;
plotavgda=0;
plotavglfp=0;
plotphys=0;
plotfft=0;
plotcor=0;
plotxvar=0;
xchs={};
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'fscvchs'
            %user provided fscv selected channels
            argnum=argnum+1;
            fscvchs=varargin{argnum};
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
            chbands;
        case 'xchs'
            argnum=argnum+1;
            xchs=varargin{argnum};
        case 'plotall'
            plotda=1;
            plotlfp=1;
            plotavgda=1;
            plotavglfp=1;
            plotphys=1;
            plotfft=1;
            plotcor=1;
            plotxvar=1;
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

targpath.bigreward=fullfile(paths{1}, 'matlab','bigreward_pro','analyzed',filesep);
targpath.smallreward=fullfile(paths{1}, 'matlab','smallreward_pro','analyzed',filesep);
targpath.targetbreak=fullfile(paths{1}, 'matlab','targetbreak_pro','analyzed',filesep);
targpath.fixbreak=fullfile(paths{1}, 'matlab','fixbreak_pro','analyzed',filesep);
homedirf=strfind(paths{1},'1dr');
homedir=paths{1}(1:homedirf-1);
manlist=[homedir 'chronic' sessid '_trialselection.xlsx'];
assignin('base','fcnStatus',targpath)   %store targpath in workspace

trialgroups=fieldnames(targpath);
numtypes=length(trialgroups);

time=[22.5 37.5]; 
if sessnum<66
    time=[25 35]; 
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
if ~isempty(fbands2)
    %2 groups of frequency bands for lfp's to analyze below
    numfbands=2;
    allbands={fbands fbands2};
end


for itype=1:numtypes
    plotparam.baseline='pretarg';  
    currpath=getfield(targpath,trialgroups{itype});
    disp(currpath)
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
    load(compiledname);
    
    %set up event windows to analyze around task
    eventsrew=repmat(plotparam.alignidx,1,length(samplesfix{1}));
    rewwin=5;           %5s post rew
    blackoutwin=.5;               %.5s away from event marker
    fixwin=1;
    intergap=2;
    targgap=1;
    eventmarks={};
    eventmarks{2}=[samplestarg{1}-targgap*param.samplerate; samplestarg{1}+targgap*param.samplerate;].*rates(1)./rates(2);
    eventmarks{3}=[samplesfix{1}-param.samplerate*blackoutwin; samplestarg{1}].*rates(1)./rates(2);
    eventmarks{4}=[samplestarg{1}; eventsrew-blackoutwin*param.samplerate].*rates(1)./rates(2);
    eventmarks{5}=[eventsrew; eventsrew+param.samplerate*rewwin].*rates(1)./rates(2);
    eventmarks{1}=[samplesfix{1}-param.samplerate*fixwin; eventsrew-blackoutwin*param.samplerate].*rates(1)./rates(2);
    eventmarks{6}=[samplesfix{1}-param.samplerate*blackoutwin; eventsrew-blackoutwin*param.samplerate].*rates(1)./rates(2);
    eventnames{1}='prerew';
    eventnames{2}='intertarg';
    eventnames{3}='fixtarg';
    eventnames{4}='targwin';
    eventnames{5}='rewwin';
    eventnames{6}='fixrew';    
    
    if isfield(plotparam,'shift')
        if strcmp(plotparam.shift,'fix')
            eventmarks={};
            eventnames={};
            eventmarks{1}=[samplesfix{1}-param.samplerate*fixwin; samplesfix{1}+param.samplerate*intergap].*rates(1)./rates(2);
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
        [~,sheetsav,~] = xlsfinfo(manlist); %find # sheets in worksheet
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
    if any(samplesfixeye{1}<100)
        chbadtrials2=find(samplesfixeye{1}<100);
        chbadtrials2=unique([chbadtrials2 find(samplesfixeye{1}<samplesfix{1})]);
        goodtrials=find(~ismember(goodtrials,chbadtrials2)==1);
    end
    if ~isfield(trialbytrial,'lfp')
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
        trialbytrial(1).cscNames=data.cscNames;
        trialbytrial(1).cscNames=plotparam.cscNames;
    end
    
    %group trial types
    trialtypes=grouptrials(goodtrials,trialbytrial(1),plotparam);

    plottrials=goodtrials;
    numtrials=size(plottrials,2);
    plotparam.cscale=[nan 6];    %ch4 chronic 30
    plotparam.vertplotwidth=100;
    plotparam.vertplotwidth2=30;
    plotparam.numtrials=numtrials;
    numcolor=length(fscvchs)+2;
    numbehav=8;     %% behavior plots
    plotparam.figpos=[200 200 1800 900];
    figpos=plotparam.figpos;
    plotparam.colorsize=[figpos(3)/(numcolor+2) figpos(4)*3/4];
    plotparam.margins=10;

    %initialize storage variables data
    da={};
    datm={};
    betatm={};      %task modulated signal averages
    betalfp={};
    betalfp2={};        %high beta band
    plotparam.chnum=targch;
%%
%plot DA trial by trial

if plotda
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
        if isfield(trialbytrial(1),'alignts')
            %save trial by trial data with reference ts's for nlx to sort
            %with small rew trials converted at another time
            dadata=da{ich};
            alignts=trialbytrial(1).alignts;
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
%[blinktm,blinkdata]=setTrialAx(trialbytrial(1),plotparam,[],...
%    plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},'blink','noplot');
%calc pulse rate
%[pulsetm,pulsedata]=setTrialAx(trialbytrial(1),plotparam,[],...
%    plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},'pulse','noplot');
if it==1
    save([pathda trialname '-data-beh'],'eyedata','lickdata');
end
saveas(fig1,[pathda trialname '.fig']);
saveas(fig1,[pathda trialname '.tif'],'tif')
save([pathda trialname],'datm','eyetm','licktm','rts');
end
end
%%
%plot lfp's trial by trial
plotparam.vertplotwidth=90;
plotparam.vertplotwidth2=10;
plotparam.margins=10;
plotparam.colorsize=[figpos(3)/(length(lfpchs)+1) figpos(4)*3/4];
if plotlfp
fig2=figure('visible','off');
hfig2=setupTrialPlots(fig2, plotparam,length(lfpchs),0);
pathlfp=[currpath 'trialbytrial_lfp' filesep];
pathlfp2=[pathlfp(1:end-1) num2str(2) filesep];
if ~isdir(pathlfp)
    mkdir(pathlfp)
end
if length(allbands)>1
    mkdir(pathlfp2)
end
for ifreq=1:length(allbands)
    fbands=allbands{ifreq};
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
                if ifreq>1
                    betalfp2{ii}=downlfp;
                else
                    betalfp{ii}=downlfp;
                end
                if ii==length(lfpchs)
                    %save betalfp data
                    alignts=[];
                    if isfield(trialbytrial(1),'alignts')
                        alignts=trialbytrial(1).alignts;
                    end
                    trialids=plottrials;
                    if ifreq==1
                    save([pathlfp trialname '-betalfp'],'betalfp',...
                        'trialids','alignts','trialtypes');
                    else
                    save([pathlfp2 trialname '-betalfp2'],'betalfp2',...
                        'trialids','alignts','trialtypes');
                    end
                end

            else
            %already got all betalfp data in first 'all trial' it
            betatm{ii}=setTrialAx(trialbytrial(1),...
                plotparam,hfig2{ii},plottrials,samplesfix{1},samplestarg{1},...
                samplesfixeye{1},lfpchs{ii},'rclfp',betalfp{ii},'fbands',fbands,'getfft',...
                'tapers',[3 5],'smoothwin',0.3,'plotnum',ii);
            end
        end
        if ifreq==1
            saveas(fig2,[pathlfp trialname '.fig']);
            saveas(fig2,[pathlfp trialname '.eps'],'epsc')
            saveas(fig2,[pathlfp trialname '.tif'],'tif')
            save([pathlfp trialname],'betatm');
        else
            %store high beta separate folder 2
            saveas(fig2,[pathlfp2 trialname '.fig']);
            saveas(fig2,[pathlfp2 trialname '.eps'],'epsc')
            saveas(fig2,[pathlfp2 trialname '.tif'],'tif')
            save([pathlfp2 trialname],'betatm');
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
cscaleda=[-4 2];
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
        setTrialAxAvg(trialbytrial(ich),plotparam,axsigrep{countp},plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},'sitename',sites(ich));
        countp=countp+1;
        setTrialAxAvg(trialbytrial(ich),plotparam,axsigrep{countp},plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},'avg','cscale',cscaleda,'sitename',sites(ich));
    end
    saveas(figsigrep,[pathavg trialname '.fig'])
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
pathavglfp=[currpath 'avg_lfp' filesep];
pathavglfp2=[pathavglfp(1:end-1) '2' filesep];
if ~isdir(pathavglfp)
    mkdir(pathavglfp)
end
if length(allbands)>1
    %store high beta separate folder 2
    mkdir(pathavglfp2);
end
for ifreq=1:length(allbands)
    fbands=allbands{ifreq};
    plotparam.fbands=fbands;
    for it=1:length(trialtypes.names)
        plottrials=trialtypes.nums{it};
        trialname=trialtypes.names{it};
        plotparam.triallabel=trialname;
        if isempty(betalfp)
            %need to load
            load([currpath 'trialbytrial_lfp' filesep 'all-betalfp'],'betalfp')
            if length(allbands)>1
                load([currpath 'trialbytrial_lfp2' filesep 'all-betalfp2'],'betalfp2')
            end
        end
        for ii=1:length(lfpchs)
            downlfp=betalfp{ii};
            if ifreq>1
                downlfp=betalfp2{ii};
            end
            setTrialAxAvg(trialbytrial(1),plotparam,axsigrep2{ii*2-1},...
                plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},...
                lfpchs{ii},'fbands',fbands,'rclfp',downlfp);
            setTrialAxAvg(trialbytrial(1),plotparam,axsigrep2{ii*2},...
                plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},...
                lfpchs{ii},'avg','fbands',fbands,'rclfp',downlfp);
        end
        if ifreq==1
            saveas(figsigrep2,[pathavglfp trialname '.fig']);
            saveas(figsigrep2,[pathavglfp trialname '.tif'],'tif')
        else
            saveas(figsigrep2,[pathavglfp2 trialname '.fig']);
            saveas(figsigrep2,[pathavglfp2 trialname '.tif'],'tif')
        end
    end
end
end
%%
%plot averages physiological response
if plotphys
figsigrep3=figure('visible','off'); 
set(figsigrep3, 'Color', [1 1 1],'position',[100 100 1000 1000]);
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
   % setTrialAxAvg(trialbytrial(1),plotparam,axsigrep3{5},plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},'pulse');
   % setTrialAxAvg(trialbytrial(1),plotparam,axsigrep3{6},plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},'pulse','avg');
%    setTrialAxAvg(trialbytrial(1),plotparam,axsigrep3{7},plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},'blink');
 %   setTrialAxAvg(trialbytrial(1),plotparam,axsigrep3{8},plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},'blink','avg');
    saveas(figsigrep3,[pathavgb trialname '.fig'])
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
    set(fighist,'Position',[300,50,700,250]);
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
    saveas(fighist,[pathrts trialname '.fig']);
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
axfft=axes;
hold(axfft,'on'); 
figff=figure('visible','off');
set(figff, 'Color', [1 1 1]);
set(figff,'Position',[300,150,700,700]);
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
    saveas(figffts,[pathsfftchs{ii} trialname '.fig'])
    evffts{ii}=setavgfft(trialbytrial(1).lfp{cscid},fs,plotparam,axfft,plottrials,...
        'events','fix',samplesfix{1},...
        'events','target',samplestarg{1},...
        'events','fixeye',samplesfixeye{1},...    
        'freq',fftlims,'sitename',lfpchs{ii},'align','target');
    saveas(figffts,[pathsfftchs{ii} trialname '-targaln.jpg'],'jpeg')
    saveas(figffts,[pathsfftchs{ii} trialname '-targaln' '.fig'])
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
pathcor=[currpath 'corr' filesep];
pathcor2=[currpath 'corr2' filesep];
if ~isdir(pathcor)
mkdir(pathcor);
end
if ~isdir(pathcor2)
mkdir(pathcor2);
end
for ifreq=1:length(allbands)
    fbands=allbands{ifreq};
    plotparam.fbands=fbands;
    for it=1:length(trialtypes.names)
        trialname=trialtypes.names{it};
        if ifreq==1
            load([pathlfp trialname],'betatm');
        else
            load([pathlfp2 trialname],'betatm');
        end
        
        load([pathda trialname],'datm');
        corrmat=cormatrix(datm,betatm,cortypes,'sitelabels','sametype');
        pathsave=pathcor;
        if ifreq>1
            pathsave=pathcor2;
        end
        save([pathsave trialname '-sametype'],'corrmat');
        corrmat=cormatrix(datm,betatm,cortypes,'sitelabels');
        save([pathsave trialname],'corrmat');
    end
end

%%
for ifreq=1:length(allbands)
    fbands=allbands{ifreq};
    plotparam.fbands=fbands;
    pathsave=pathcor;
    if ifreq>1
        pathsave=pathcor2;
    end
    for it=1:length(trialtypes.names)
        trialname=trialtypes.names{it};
        load([pathsave trialname],'corrmat');
        cor=corrmat;
        figcorr=figure('position',[100 100 1400 900],'color',[1 1 1],'visible','off');
        figure(figcorr)
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
        cmapp=cool; cmapp=flipud(cmapp(1:12:end,:))-.35;
        cmapc=winter; cmapc=flipud(cmapc(1:12:end,:))-.35;
        cmapp(cmapp<0)=0;
        cmapc(cmapc<0)=0;
        sitep=find(contains(lfpchs,'p')==1);
        sitec=find(contains(lfpchs,'c')==1);
        sitecolors=zeros(length(lfpchs),3);
        sitecolors(sitep,:)=cmapp(1:length(sitep),:);
        sitecolors(sitec,:)=cmapc(1:length(sitec),:);
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
                fsize=10*log10(1/pval);
                if curdata(irow).r<0
                    m='o';
                    %c=[0 0 1];
                end
                if ~strcmp(fregion,eregion)
                    %different regions, offset x slightly & color
                    %c(c~=0)=c(c~=0)-.65;
                    xcord=xcord+.2;
                    ycord=ycord+.2;
                    %fsize=13;
                end
                a=text(axcorr{ii},xcord,ycord,m,'fontsize',fsize,'color',c);
               % alpha(a,'.5');

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
        saveas(figcorr,[pathsave trialname '.fig']);
        saveas(figcorr,[pathsave trialname],'tif')
    end
end
end
%%
%get cross covariance just selected trials above
if plotxvar
plotparam.xrates=rates;
if isempty(betalfp)
    %need to load
    load([currpath 'trialbytrial_lfp' filesep 'all-betalfp'],'betalfp')
    if length(allbands)>1
        load([currpath 'trialbytrial_lfp2' filesep 'all-betalfp2'],'betalfp2')
    end
end
if isempty(da)
    %need to load
    load([currpath 'trialbytrial_da' filesep 'all-data-ch' num2str(max(fscvchs))],'da')
end
datax.lfp=betalfp;
datax.da=da;
plotparam.chnums=xchs{1};
for ifreq=1:length(allbands)
    fbands=allbands{ifreq};
    plotparam.fbands=fbands;
    datax.lfp=betalfp;
    if ifreq>1
        datax.lfp=betalfp2;
    end
    for it=1:3
        plottrials=trialtypes.nums{it};
        trialname=trialtypes.names{it};
        if ifreq>1
            trialname=[trialname '2'];
        end
        plotparam.triallabel=trialname;
        xcovdata=plotxvartrials(datax,eventmarks,eventnames,plotparam,'type','none',...
        'trials',plottrials,'barlags','label',trialname);
    end
end
plotparam.chnums=xchs{2};
for ifreq=1:length(allbands)
    fbands=allbands{ifreq};
    plotparam.fbands=fbands;
    datax.lfp=betalfp;
    if ifreq>1
        datax.lfp=betalfp2;
    end
    for it=1:3
        plottrials=trialtypes.nums{it};
        trialname=trialtypes.names{it};
        if ifreq>1
            trialname=[trialname '2'];
        end
        plotparam.triallabel=trialname;
        xcovdata=plotxvartrials(datax,eventmarks,eventnames,plotparam,'type','none',...
        'trials',plottrials,'barlags','label',trialname);
    end
end
end

end


