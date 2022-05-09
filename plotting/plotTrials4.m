%Plot trial by trial data as loaded in trials_chronicXX_XX_processed.mat
%generated in compileTrialsPatrDir (multiple channels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting functions
time=[22.5 37.5]; 
%time=[25 35]; 
time=[25 35]; 

%lfpchs={'p2-p3','p3-p5','pl1-pl2','cl1-cl4','cl4-cl5'};         %chronic 38
lfpchs={'p1-p3','pl2-p1','pl3-p3','p1-pl3','cl1-cl4','cl4-cl5','cl1-cl5'}; %chronic58
lfpchs={'p1-p3','pl2-p1','p1-pl3','cl1-cl4','cl4-cl5','cl1-cl5'}; %chronic58
%lfpchs={'p1-p2','p1-p3','p2-p3','pl2-p1','pl2-pl3','p1-pl3',...
%    'cl1-cl4','cl4-cl5','cl1-cl5'};     %chronic 113
lfpchs={'pl1-p5','p1-p5','p2-p5','cl1-cl4','cl3-cl4','cl4-cl6'};   %chronic65/67
%frequency pass bands for lfp sites as saved in evffts from setavgfft code
%if bandwidth <5 Hz, then extend from center +/-3 hz
global plotParam
plotParam={};
otherchs={'eyed','lickx','pulse'};
ncschs=[lfpchs otherchs];
[param,csc_map,eventcodes]=getparams('patrabipolar','default',ncschs);
getplotsettings([],ncschs,[],1000,[]);       
plotparam=plotParam;
%plotparam.fbands=fbands;
plotparam.shift='targ';     %shift data to target time point for all trials in settrialax and settrialaxavg
%plotparam.shift='fix';     %shift data to target time point for all trials in settrialax and settrialaxavg

[FileName,PathName] = uigetfile('*.*','Select file');
load([PathName,FileName]);
dadata=trialbytrial.da;     %store original da data as variable
%sort signals
fscvchs=[  2  4];
plotparam.chnums=fscvchs;

ich=2;      %sort signals by this channel
numtrials=size(trialbytrial(ich).da,1);

if ~isfield(trialbytrial,'lfp')
    %not merged file, get lfps separately based on trial num 100 etc.
    %from original xx_pro directory
    data={};
    %get reconverted path ie XX_pro directory
    pathid=strfind(PathName,'_pro');
    pathlfp=PathName(1:pathid+4);
    namx=strfind(FileName,'.mat');
    trialnum=FileName(namx-3:namx-1);      %get trialnum to find corresponding files    
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

%%
%loads also fscvchs by default
namy=strfind(FileName,'_');
trialtype=FileName(namy(2)+1:namy(3)-1);
plotparam.trialtype=trialtype;
plotparam.samplespersec=10;
plotparam.pathname=PathName;
win=time(1)*plotparam.samplespersec:1:time(2)*plotparam.samplespersec;
plotparam.win=win;
plotparam.interval=2.5;

%read bad trials from excel spreadsheet
sheet = 1;
if ~isempty(strfind(FileName,'small'))
    sheet=2;    %small reward sheet is #2
elseif ~isempty(strfind(FileName,'big'))
    sheet=1;
elseif ~isempty(strfind(FileName,'error'))
    sheet=3;    %target error sheet
end
%read bad trials from excel spreadsheet
[filexls,pathxls] = uigetfile('*.*','Select spreadsheet file with bad trials');
chbadtrials=[];
if filexls~=0
    %xlRange = 'B2:C3';
    %subsetA = xlsread(filename,sheet,xlRange);
    xlsdata = xlsread([pathxls filexls],sheet);
    chbadtrials=xlsdata-100;        % 101 first trial
    chbadtrials=xlsdata-99;        %default, 100 first trial
end

%get labels for fscv chs from original file names
origpath=PathName(1:strfind(PathName,'cvtotxt')-1);
filesorig=dir(origpath);
filenamesorig={filesorig.name};
targfiles=strfind(filenamesorig,'1dr_');
labeledfiles=find(~cellfun(@isempty,targfiles));
%contains function does not work for 2013 for cells 
filesorig=filenamesorig(labeledfiles);
sitenamessp=strsplit(filesorig{1},'_');
sites=sitenamessp(2:5);
plotparam.sites=sites;
plotparam.lfpchs=lfpchs;
[parameters,csc_map,event_codes]=getparams('patrabipolar');

plotparam.colormap=parula;
plotparam.markercolors=[1 1 1; 1 0 1; 0 1 1];
plotparam.baseline='postreward';        %look at intertrial change win
plotparam.baseline='fix';           %default from fix appear to fix start
plotparam.baseline='fixeye';           %immediate post-fix
plotparam.baseline='postcue';           %immediate post-cue/fix appear
plotparam.baseline='precue';           %immediate pre fix appear
plotparam.baseline='prereward';
plotparam.baseline='alltrial';      %baseline is entire duration of trial
plotparam.baseline='pretarg';           %immediate post-cue/fix appear

plotparam.alignidx=300;         %align to 30s Reward period
plotparam.cminshow=0;

if strcmp(plotparam.baseline,'postreward')
    %win needs to be changed for inter trial
    time=[30 40];
    win=time(1)*plotparam.samplespersec:1:time(2)*plotparam.samplespersec;
    plotparam.win=win;
    plotparam.interval=2.5;
end

ratelfp=round(trialbytrial(1).ratelfp);
rates=[ratelfp plotparam.samplespersec];
eventsrew=repmat(plotparam.alignidx,1,length(samplesfix{1}));
rewwin=5;           %5s post rew
longwin=3.5;          %3.5s long window
blackoutwin=.5;               %.5s away from event marker
fixwin=1;
posttargwin=2;
intergap=2;
interwin=3;
intertarg=2;
targgap=1;
eventmarks={};
eventmarks{2}=[samplestarg{1}-targgap*plotparam.samplespersec; samplestarg{1}+targgap*plotparam.samplespersec;].*rates(1)./rates(2);
eventmarks{3}=[samplesfix{1}-plotparam.samplespersec*blackoutwin; samplestarg{1}].*rates(1)./rates(2);
eventmarks{4}=[samplestarg{1}; eventsrew-blackoutwin*plotparam.samplespersec].*rates(1)./rates(2);
eventmarks{5}=[eventsrew; eventsrew+plotparam.samplespersec*rewwin].*rates(1)./rates(2);
eventmarks{1}=[samplesfix{1}-plotparam.samplespersec*fixwin; eventsrew-blackoutwin*plotparam.samplespersec].*rates(1)./rates(2);
eventmarks{6}=[samplesfix{1}-plotparam.samplespersec*blackoutwin; eventsrew-blackoutwin*plotparam.samplespersec].*rates(1)./rates(2);

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
eventmarks{1}=[samplesfix{1}-plotparam.samplespersec*fixwin; samplesfix{1}+plotparam.samplespersec*intergap].*rates(1)./rates(2);
eventnames{1}='fix';   
    end
end
%eventnames{5}='targwin1s';
%eventmarks{7}=[eventsrew; eventsrew+plotparam.samplespersec*longwin].*rates(1)./rates(2);
%eventnames{7}='peakpostrew35s';
%eventnames{8}='fixwin2spre15spost';
%eventmarks{8}=[samplesfix{1}-plotparam.samplespersec*2; samplesfix{1}+plotparam.samplespersec*1.5].*rates(1)./rates(2);



goodtrials=1:numtrials;       %ignore trend finding algo, bad
%first 50% of trials
%goodtrials=1:round(numtrials/2);
if exist('chbadtrials')
    numtrials=size(trialbytrial(ich).da,1);
    alltrials=1:numtrials;
    if ~isempty(chbadtrials)
    %use manual selected trials
    if size(chbadtrials,2)>1
        goodtrials=find(~ismember(alltrials,chbadtrials(:,ich))==1);
    else
        goodtrials=find(~ismember(alltrials,chbadtrials)==1);
    end
    end
end

%get only trials where big side switched sides
sideshiftdetect=circshift(diff(trialbytrial(1).rewardside),[0 1]);
sideshiftedtrials=find(abs(sideshiftdetect)==1);
%circ shift to right 1 index since diff is between first point &next

%only trials where big reward side same as current trial side (for error
%trials)
bigrewsideonly=and(trialbytrial(1).rewardside,trialbytrial(1).bigside)+and(~trialbytrial(1).rewardside,~trialbytrial(1).bigside);
bigrewsidetrials=find(bigrewsideonly==1);
%trialbytrial(1).rewardside accurate for what side is big reward for
%calculating switch
findswitch=diff(trialbytrial(1).rewardside);
findswitch=[findswitch(1) findswitch];
switchtrialsall=find(abs(findswitch)>0);
switchtrials=goodtrials(find(ismember(goodtrials,switchtrialsall)));

afterfail=goodtrials(find(ismember(goodtrials,afterfailtrials{1})));
aftersuccess=goodtrials(find(~ismember(goodtrials,afterfailtrials{1})));

nextfailtrials=find(trialbytrial(1).nextsuccess==0);
beforesuccess=goodtrials(find(~ismember(goodtrials,nextfailtrials)));
beforefail=goodtrials(find(ismember(goodtrials,nextfailtrials)));

shiftside=goodtrials(find(ismember(goodtrials,sideshiftedtrials)==1));
rightsidetrials=find(trialbytrial(1).rewardside==1);
lsidetrials=find(trialbytrial(1).rewardside==0);
goodrightsidetrials=goodtrials(find(ismember(goodtrials,rightsidetrials)==1));
goodleftsidetrials=goodtrials(find(ismember(goodtrials,lsidetrials)==1));

%sort reaction times
[aa, bb]=sort(trialbytrial(ich).target_rt);
rtsortedtrials=bb;
rtsortedsel=rtsortedtrials(find(ismember(rtsortedtrials,goodtrials)==1));
quart=floor(length(goodtrials)/4);
quarttrials=1:quart:length(goodtrials);
phase1trials=goodtrials(1:quarttrials(2)-1);
phase2trials=goodtrials(quarttrials(2):quarttrials(3)-1);
phase3trials=goodtrials(quarttrials(3):quarttrials(4)-1);
phase4trials=goodtrials(quarttrials(4):end);
phase1rtrials=intersect(goodrightsidetrials,phase1trials);
phase1ltrials=intersect(goodleftsidetrials,phase1trials);
phase4rtrials=intersect(goodrightsidetrials,phase4trials);
phase4ltrials=intersect(goodleftsidetrials,phase4trials);

%trialtypes.names={'all','switch','left','right',...
 %   'phase1','phase2','phase1right','phase1left','phase2right','phase2left'};
trialtypes.names={'all','switch','sleft','sright','aftersuccess','afterfail',...
    'beforesuccess','beforefail','shiftside','phase1','phase2','phase3','phase4',...
    'srightphase1','sleftphase1','srightphase4','sleftphase4'};
trialtypes.nums={goodtrials,switchtrials,goodleftsidetrials,goodrightsidetrials,...
    aftersuccess,afterfail,beforesuccess,beforefail,shiftside,phase1trials,phase2trials,phase3trials,phase4trials,...
    phase1rtrials,phase1ltrials,...
    phase4rtrials,phase4ltrials};
if strcmp(plotparam.trialtype,'fixbreak')
    %just one type of trial, all
    trialtypes={};
    trialtypes.names{1}='all';
    trialtypes.nums{1}=goodtrials;
    nanrts=find(isnan(trialbytrial(1).fix_rt)==1);
    numnan=length(nanrts);
    if numnan>2
        %use estimate provided by markers instead
        trialbytrial(1).fix_rt=(samplesfixeye{1}-samplesfix{1})./plotparam.samplespersec;
    end
    %make trial types based on rt's, engaged vs not engaged
    fastrttrials=find(trialbytrial(1).fix_rt<=0.2);
    slowrttrials=find(trialbytrial(1).fix_rt>0.2);
    trialtypes.names={'all','fastrt','slowrt'};
    trialtypes.nums={goodtrials,fastrttrials,slowrttrials};
end
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

%high beta freq
plotparam.hfbands={                 
    'p2-p5',    [22 28]...
    'p3-p5',    [22 28]...
    'pl1-p1',   [22 28]...
    'p1-p5',    [22 28]...      
    'pl2-p1',   [22 28]...   
    'pl1-p5',   [22 28]...
    'pl2-pl3',  [22 28]...    
    'cl1-cl4',  [22 28]...    
    'cl1-cl5',  [22 28]...
    'cl3-cl4',  [22 28]...
    'cl4-cl5',  [22 28]...    
    'cl4-cl6',  [22 28]...
    'cl3-cl6',  [22 28]...
    's6-s5',    [22 28]...
    's4-s3',    [22 28]...
    's3-s2',    [22 28]...
    };
fbands=plotparam.fbands;
if isfield(plotparam,'hfbands')
    fbands=plotparam.hfbands;
end

plotparam.lfbands={                 
    'p2-p5',    [12 17]...
    'p3-p5',    [13 17]...
    'pl1-p1',   [13 17]...
    'p1-p5',    [13 17]...      
    'pl2-p1',   [13 17]...   
    'pl1-p5',   [13 17]...
    'pl2-pl3',  [13 17]...    
    'cl1-cl4',  [13 22]...    
    'cl1-cl5',  [13 22]...
    'cl3-cl4',  [13 22]...
    'cl4-cl5',  [13 22]...    
    'cl4-cl6',  [13 22]...
    'cl3-cl6',  [13 22]...
    's6-s5',    [13 22]...
    's4-s3',    [13 22]...
    's3-s2',    [13 22]...
    };
fbands=plotparam.fbands;
if isfield(plotparam,'lfbands')
    fbands=plotparam.lfbands;
end

%%
%plot DA trial by trial

fig1=figure;
cscales=[];
cscales=[-3 8; -3 8; -3 8; -3 8]; 
eyedscale=[-2.5 .5].*1e-3;   
%eyedscale=[5 20].*1e-5;   %chronic 93
lickscale=[1.3e-4 4.8e-4];
%eyedscale=[];
%lickscale=[];
hfig=setupTrialPlots(fig1, plotparam,numcolor,numbehav);
averages={};
da={};
datm={};
pathda=[PathName 'trialbytrial_da' filesep];
if ~isdir(pathda)
    mkdir(pathda)
end
for it=1:length(trialtypes.names)
    %it=2;
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
[blinktm,blinkdata]=setTrialAx(trialbytrial(1),plotparam,[],...
    plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},'blink','noplot');
%calc pulse rate
[pulsetm,pulsedata]=setTrialAx(trialbytrial(1),plotparam,[],...
    plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},'pulse','noplot');
if it==1
    save([pathda trialname '-data-beh'],'eyedata','lickdata',...
    'blinkdata','pulsedata');
end
%savefig(fig1,[pathda trialname '.fig']);
saveas(fig1,[pathda trialname '.eps'],'epsc')
saveas(fig1,[pathda trialname '.tif'],'tif')
save([pathda trialname],'datm','eyetm','licktm',...
    'blinktm','pulsetm','rts');

end

%%
%plot lfp's trial by trial
fig2=figure;
plotparam.vertplotwidth=90;
plotparam.vertplotwidth2=10;
plotparam.margins=10;
plotparam.colorsize=[figpos(3)/(length(lfpchs)+1) figpos(4)*3/4];
hfig2=setupTrialPlots(fig2, plotparam,length(lfpchs),0);
betatm={};      %task modulated signal averages
betalfp={};
pathlfp=[PathName 'trialbytrial_lfp' filesep];
if ~isdir(pathlfp)
    mkdir(pathlfp)
end
for it=1:length(trialtypes.names)
    plottrials=trialtypes.nums{it};
    trialname=trialtypes.names{it};
    plotparam.triallabel=trialname;
%plotting trial by trial
for ii=1:length(lfpchs)
 %   [betaavgs{ii}]=setTrialAx(trialbytrial(1),plotparam,hfig2{ii},plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},'lfp2',lfpchs{ii});
        if it==1
    [betatm{ii},betalfp{ii}]=setTrialAx(trialbytrial(1),...
        plotparam,hfig2{ii},plottrials,samplesfix{1},samplestarg{1},...
        samplesfixeye{1},lfpchs{ii},'fbands',fbands,'getfft',...
        'tapers',[3 5],'smoothwin',0.3,'plotnum',ii);
        else
            %already got all betalfp data in first 'all trial' it
       betatm{ii}=setTrialAx(trialbytrial(1),...
        plotparam,hfig2{ii},plottrials,samplesfix{1},samplestarg{1},...
        samplesfixeye{1},lfpchs{ii},'rclfp',betalfp{ii},'fbands',fbands,'getfft',...
        'tapers',[3 5],'smoothwin',0.3,'plotnum',ii);
        end
end
%savefig(fig2,[PathName 'trialbytriallfpsbeta.fig']);
saveas(fig2,[pathlfp trialname '.eps'],'epsc')
saveas(fig2,[pathlfp trialname '.tif'],'tif')
        save([pathlfp trialname],'betatm');
end

%plot averages da
%%
figsigrep=figure; set(figsigrep, 'Color', [1 1 1],'position',[100 50 1000 950]);
set(0,'CurrentFigure',figsigrep); 
axsigrep={};
cscaleda=[-2 4];
cscaleda=[-4 2];
for ii=1:8
    axsigrep{ii}=subplot(4,2,ii);
end
pathavg=[PathName 'avg_da' filesep];
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
saveas(figsigrep,[pathavg trialname '.eps'],'epsc')
saveas(figsigrep,[pathavg trialname '.tif'],'tif')
%save([pathavg trialname],'betatm');
end

%plot averages phys
%%
figsigrep2=figure; set(figsigrep2, 'Color', [1 1 1],'position',[100 50 1500 900]);
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

pathavglfp=[PathName 'avg_lfp' filesep];
if ~isdir(pathavglfp)
    mkdir(pathavglfp)
end
for it=1:length(trialtypes.names)
    plottrials=trialtypes.nums{it};
    trialname=trialtypes.names{it};
    plotparam.triallabel=trialname;
for ii=1:length(lfpchs)
setTrialAxAvg(trialbytrial(1),plotparam,axsigrep2{ii*2-1},plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},lfpchs{ii},'fbands',fbands,'rclfp',betalfp{ii});
setTrialAxAvg(trialbytrial(1),plotparam,axsigrep2{ii*2},plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},lfpchs{ii},'avg','fbands',fbands,'rclfp',betalfp{ii});
end
saveas(figsigrep2,[pathavglfp trialname '.eps'],'epsc')
saveas(figsigrep2,[pathavglfp trialname '.tif'],'tif')
%save([pathavg trialname],'betatm');
end
%%
%plot averages physiological response

figsigrep3=figure; set(figsigrep3, 'Color', [1 1 1],'position',[100 100 1000 1000]);
set(0,'CurrentFigure',figsigrep3); 
axsigrep3={};
for ii=1:8
    axsigrep3{ii}=subplot(4,2,ii);
end

pathavgb=[PathName 'avg_beh' filesep];
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

saveas(figsigrep3,[pathavgb trialname '.eps'],'epsc')
saveas(figsigrep3,[pathavgb trialname '.tif'],'tif')
%save([pathavg trialname],'betatm');
end


%reaction times
binwidth=25;              %nm
binmax=600;
pathrts=[PathName 'rts' filesep];
if ~isdir(pathrts)
    mkdir(pathrts)
end
for it=1:length(trialtypes.names)
    trialname=trialtypes.names{it};
    load([pathda trialname],'rts');
    [binsrts,bb]=histc(rts.left.*1000,[0:binwidth:binmax]);
    [binsrtsright,bb]=histc(rts.right.*1000,[0:binwidth:binmax]);
    fighist=figure;
    xlims=[0 binmax];
    subhist(1)=subplot(4,2,1);
    set(fighist, 'Color', [1 1 1]);
    %set(fig1, 'Position', [100, 50, 1450, 650]);
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
    savefig(fighist,[pathrts trialname]);
    %saveas(figsigrep3,[PathName 'tracesnlxbehav.eps'],'epsc')
    saveas(fighist,[pathrts trialname],'tif')
end
close all;
%%
%print average FFT spectrograms
pathffts=[PathName 'avgffts\'];
mkdir(pathffts)
fs=1000;        %sampler ate nlx
figffts=figure;
set(figffts, 'Color', [1 1 1]);
set(figffts,'Position',[300,150,1200,500]);
axfft=axes;
hold(axfft,'on'); 
figff=figure;
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

saveas(figffts,[pathsfftchs{ii} trialname '.eps'],'epsc')
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
%save([pathavg trialname],'betatm');
end
close all;

%%
%calculate trial by trial correlation matrix
%Are trial by trial task-modulated signals changing in the same manner?
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
pathcor=[PathName 'corr' filesep];
if ~isdir(pathcor)
mkdir(pathcor);
end
for it=1:length(trialtypes.names)
    trialname=trialtypes.names{it};
    load([pathlfp trialname],'betatm');
    load([pathda trialname],'datm');

corrmat=cormatrix(datm,betatm,cortypes,'sitelabels','sametype');
save([pathcor trialname '-sametype'],'corrmat');
corrmat=cormatrix(datm,betatm,cortypes,'sitelabels');
save([pathcor trialname],'corrmat');
end
close all;

%%
for it=1:length(trialtypes.names)
trialname=trialtypes.names{it};
load([pathcor trialname],'corrmat');
cor=corrmat;
figcorr=figure('position',[100 100 1400 900],'color',[1 1 1]);
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
        %text(axcorr{ii},xcord,ycord,m,'fontsize',18,'fontweight','bold','color',c);
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
savefig(figcorr,[pathcor trialname]);
%saveas(figsigrep3,[PathName 'tracesnlxbehav.eps'],'epsc')
saveas(figcorr,[pathcor trialname],'tif')
end
close all;
%%
%get cross covariance just selected trials above
plotparam.xrates=rates;
datax.lfp=betalfp;
datax.da=da;
plotparam.chnums=[2 4];
for it=1:4
    plottrials=trialtypes.nums{it};
    trialname=trialtypes.names{it};
    plotparam.triallabel=trialname;
xcovdata=plotxvartrials(datax,eventmarks,eventnames,plotparam,'type','coeff',...
'trials',plottrials,'barlags','label',trialname);
%xcovdata=plotxvartrials(datax,eventmarks,eventnames,plotparam,'type','none',...
%'trials',plottrials,'barlags','label',trialname,'cluster','pos');
close all;
end

%%
%plot relevant lfps/beta side by side
time=[25 30]; 
trialname='';
plotparam.triallabel=trialname;
plotparam.win=time(1)*plotparam.samplespersec:1:time(2)*plotparam.samplespersec;
ii=4; 
ich=fscvchs(ii);
sorttype='targwin';
[aa,bb]=sort(getfield(datm{ich},sorttype));
%[aa,bb]=sort(getfield(betatm{ilfp},sorttype));

plottrials=datm{ich}.trialnums(bb);
%[aa, rtsrt]=sort(trialbytrial(1).target_rt);
%plottrials=rtsrt(find(ismember(rtsrt,plottrials)==1));

%ilfp=find(ismember(lfpchs,lfptargs{2})==1);
%[aa,bb]=sort(getfield(betatm{ilfp},sorttype));
%plottrials=betatm{ilfp}.trialnums(bb);

numtrials=length(plottrials);
plotparam.chnum=ich;
lfptargs={'p1-p5','pl1-p5'};
%lfptargs={'cl4-cl5','cl1-cl4'};
%lfptargs={'p2-p5'};
lfptargs={'cl4-cl5','cl1-cl5'};
%lfptargs={'cl4-cl5','cl1-cl5'};
%lfptargs={'cl3-cl4'};

plotparam.baseline='precue';           %immediate post-cue/fix appear
%plotparam.baseline='prereward';           %immediate post-cue/fix appear
plotparam.baseline='pretarg';           %immediate post-cue/fix appear

plotparam.cscale=[nan 6];    %ch4 chronic 30
plotparam.vertplotwidth=100;
plotparam.vertplotwidth2=30;
plotparam.numtrials=numtrials;
numcolor=length(lfptargs)+3;
numbehav=8;     %% behavior plots
plotparam.figpos=[200 200 1800 900]; figpos=plotparam.figpos;
plotparam.colorsize=[figpos(3)/(numcolor+2) figpos(4)*3/4];
plotparam.margins=10;
fig1=figure;
cscales=[-3 8; -3 8; -3 8; -3 8]; 
eyedscale=[-2.5 .5].*1e-3;   
lickscale=[1.3e-4 4.8e-4];
hfig=setupTrialPlots(fig1, plotparam,numcolor,numbehav);

setTrialAx(trialbytrial(ich),plotparam,hfig{1},plottrials,...
    samplesfix{1},samplestarg{1},samplesfixeye{1},'plotnum',1,'cscale',cscales(1,:),'sitename',sites(ich));
for ii=1:length(lfptargs)
    ilfp=find(ismember(lfpchs,lfptargs{ii})==1);
    setTrialAx(trialbytrial(1),...
    plotparam,hfig{ii+1},plottrials,samplesfix{1},samplestarg{1},...
    samplesfixeye{1},lfpchs{ilfp},'fbands',fbands,'smoothwin',.6);
end
setTrialAx(trialbytrial(1),plotparam,hfig{length(lfptargs)+2},...
   plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},'eyed',...
    'cscale',eyedscale);
setTrialAx(trialbytrial(1),plotparam,hfig{length(lfptargs)+3},...
    plottrials,samplesfix{1},samplestarg{1},samplesfixeye{1},'lickx',...
    'cscale',lickscale);
setTrialBehavAx(trialbytrial(1),plotparam,hfig,plottrials,numcolor+1);
%plot averages phys
%%
figsigrep2=figure; set(figsigrep2, 'Color', [1 1 1],'position',[100 50 1000 1000]);
set(0,'CurrentFigure',figsigrep2); 
axsigrep2={};
numt=16;
        ilfp=find(ismember(lfpchs,lfptargs{1})==1);

for ii=1:numt
        axsigrep2{ii}=subplot(numt/4,4,ii);
end
n=25;
for ii=1:8
    %ptrials=plottrials(length(plottrials)-100:length(plottrials));
       ptrials=plottrials((ii-1)*n+1:ii*(n));
       ptrials=ptrials(ptrials>0 & ptrials<=length(goodtrials));
setTrialAxAvg(trialbytrial(1),plotparam,axsigrep2{ii*2-1},ptrials,samplesfix{1},samplestarg{1},samplesfixeye{1},lfpchs{ilfp},'fbands',fbands,'rclfp',betalfp{ilfp});
%setTrialAxAvg(trialbytrial(1),plotparam,axsigrep2{ii*2},ptrials,samplesfix{1},samplestarg{1},samplesfixeye{1},lfpchs{ilfp},'avg','fbands',plotparam.fbands,'rclfp',betalfp{ilfp});
%setTrialAxAvg(trialbytrial(1),plotparam,axsigrep2{ii*2},ptrials,samplesfix{1},samplestarg{1},samplesfixeye{1},lfpchs{ilfp},'avg','fbands',plotparam.fbands,'rclfp',betalfp{ilfp});
    %countp=countp+1;
    %ich=fscvchs(ii);
plotparam.chnum=ich;
setTrialAxAvg(trialbytrial(ich),plotparam,axsigrep2{ii*2},ptrials,samplesfix{1},samplestarg{1},samplesfixeye{1});

end
