function plotbehcorsummary(xinfos,binfos,plotparam,varargin)
%2/17/2019 pre, now updated xinfos so not arrayed by sess num, all in one
%chart
argnum=1;
datype='dapos';
evts={};
sestypes={};
grpreg=0;

while argnum<=length(varargin)
    switch varargin{argnum}
        case 'daneg'
            datype='daneg';
        case 'dareb'
            datype='dareb';
        case 'daall'
            datype='daall';
        case 'evts'
            argnum=argnum+1;
            evts=varargin{argnum};
        case 'sestypes'
            %if "reward', then big and small reward
            argnum=argnum+1;
            sestypes=varargin{argnum};  
        case 'grpreg'
            grpreg=1;       %group regions            
        case 'ttypes'
            argnum=argnum+1;
            ttypes=varargin{argnum};        %condition types grouped eg. {{'big','left'},{'small','left','aftersm'}}
    end
    argnum=argnum+1;
end
fontsize=9;
savepath=plotparam.savepath;
savepath=[savepath(1:end-1) 'beh' filesep];
        if ~isdir(savepath)
            mkdir(savepath);
        end
sessnums=plotparam.sessnums;
trialgrps=plotparam.trialgrps;
trialinfo=trialgrps(1).trialinfo;
lfpchs=plotparam.lfpchs;
%get lfp ch groups p & c
pgroup=find(contains(lfpchs,'p')==1);
cgroup=find(contains(lfpchs,'c')==1);
targdasites=plotparam.dasites;
figpos=[50,50,1400,900];
figsess=figure('position',figpos,'color',[1 1 1]);
set(0,'CurrentFigure',figsess);    %set figure handle to current figure
axa={};

[dapair,lfppair]=getsitepairs(targdasites);
targdasites=plotparam.dasites;

    sites=getsites(sessnums,targdasites);
uniquesites=unique({sites(1:end).site});
cnsites=uniquesites(contains(uniquesites,'c'));
psites=uniquesites(contains(uniquesites,'p'));

binfo=binfos;
xinfo=xinfos;
eventrows={xinfo(1:end).event};
eventtypes=unique(eventrows);
if ~isempty(evts)
    eventtypes=evts;
end
sessionrows={xinfo(1:end).sessiontype};
sessiontypes=unique({xinfo(1:end).sessiontype});
if ~isempty(sestypes)
    sessiontypes=sestypes;
    if strcmp(sestypes{1},'reward')
        savepath=[savepath(1:end-1) 'rew' filesep];
        if ~isdir(savepath)
            mkdir(savepath);
        end
    end
end
params={'damax','damaxts','damin','damints','dafallts','darisets',...
    'lfpmax','lfpmaxts','lfpmin','lfpmints','lfprisets','lfpfallts',...
    'lfppostmax','lfppostmaxts'};
params_mag={'damax','damin',...
    'lfpmax','lfpmin','lfppostmax','zlagcoef'};
params_xcov={'maxprelagts',...
    'minprelagts','maxpostlagts','minpostlagts','maxprecoef',...
    'minprecoef','maxpostcoef','minpostcoef'};
params_timing_damax={'delt_lfpmin_damax','delt_lfprise_damax',...
    'delt_lfppostmax_damax','lfppostmaxts','damaxts','lfpmints'}; 

params=[params_mag params_timing_damax];   
ygrp={};
ynam={};
ygrp{1}=params_mag;
ygrp{2}=params_xcov;
ygrp{3}=params_timing_damax;

ynam{1}='mag';
ynam{2}='xcov';
ynam{3}='timing_damax';
%bparams={'fix_rt','target_lrt','target_rrt','eyed','leyed','reyed',...
 %   'pulse','lickpre','lickpost'};
bparams=unique({binfo(1:end).sitelfp});
bparamtarg='lfpmax';    %get max amp of phys sig within targ period
targevent=evts;
axpos={};
axsiz=[750,320];
axoff=100;
mark=[0 .1 .9; 0 .5 .3; 0 .9 .1];
mark=[.7 0 0; 0.6 0 .7; 0 .7 .1];
mark=[0 0.1 .7; 0.8 0 .5; 0 .7 0;.7 0 0; 0.4 0 .7];

markr=[1 .5 .0; .3 .3 0; .9 .1 0];
typeoff=[-0.3 -0.15 0];
sessmark={'.','+','o','^','s','*','v','p','h','x'};
sessmarktxt={'.','+','o','\delta','sq','*','\nabla','p','h','x'};
sessmark={'o','sq','o','^','+','v','*','p','h','x','d','<','>','o','sq','o','^','+','v','*','p','h','x','d','<','>'};
sessmarktxt={'.',' sq','o','\Delta','+','\nabla','*','p','h','x','d','<','>','.',' sq','o','\Delta','+','\nabla','*','p','h','x','d','<','>'};
marksize=[50 200 150 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 50 200 150 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200];
cmapp=[1 .5 0;1 .55 0; 1 .6 0; 1 .65 0; 1 .7 0; 1 .75 0; 1 .8 0; 1 .85 0; 1 .9 0];
cmapc=[1 .5 0;1 .55 0; 1 .6 0; 1 .65 0; 1 .7 0; 1 .75 0; 1 .8 0; 1 .85 0; 1 .9 0];
cmap2=cool;
cmap1=hot;
markc=cmap1(1:5:end,:);
markp=cmap2(1:5:end,:);
condlabels={};
condfull=[];
for ic=1:length(ttypes)
        if ic==1
            condlabels=ttypes{ic}(1:3);
        else
            condlabels=[condlabels '_' ttypes{ic}(1:3)];
        end
        
end
condfull=condlabels;

for yid=1:length(ygrp)  
    avgda={};
yvals={};
bvals={};
nn=0;
sitenum=0;
countc=0;
countp=0;
cursessids=[];
curda=[];
dasites=[];

    %timing/magnitudes/xcov parameters to plot in diff figs    
    clf(figsess,'reset');
    set(figsess,'color',[1 1 1]);
    axpos={};
    xlines=[2 4 6 8];
    for ip=1:2
        axa{ip}=subplot(1,2,ip);   hold(axa{ip},'on');
        set(axa{ip},'units','pixels');
        axis('square')         
        set(axa{ip},'ytick',1:length(ygrp{yid}));
        set(axa{ip},'yticklabel',ygrp{yid});
        set(axa{ip},'ylim',[0 length(ygrp{yid})+.5]);
        set(axa{ip},'TickLabelInterpreter','none');
        set(axa{ip},'xtick',1:length(bparams));
        set(axa{ip},'xticklabel',bparams);
        set(axa{ip},'xlim',[0 length(bparams)+.5]);
        set(axa{ip},'xTickLabelRotation',90)
        axpos{ip}=get(axa{ip},'position');
        for iline=1:min(length(bparams),length(ygrp{yid}))
        aa=plot(axa{ip},[iline iline], ...
            [0 length(ygrp{yid})],'linewidth',1,'color',[0 0 0]);
        aa.Color(4)=.2;
        aa=plot(axa{ip},[0 length(bparams)], ...
            [iline iline],'linewidth',1,'color',[0 0 0]);
        aa.Color(4)=.2;
        end    
    end    
    labels={};
    titletext=[targevent ' | ' ];  
    set(axa{1},'units','pixels')
    text(axa{1},-30,axpos{1}(4)+30,titletext,'units','pixels','fontweight','bold');
    %set legend
    labelx={};
    labelx={};
ranj=[];        %random jitter x axis for each site
countc=1;
countp=1;
ytxtcor=550;
for is=1:length(uniquesites)
    labelx{is}=[num2str(uniquesites{is})];  
    axpos=get(axa{1},'position');
     daregion=contains(uniquesites(is),'c');     %1=='c'
    if daregion
        %CN
        if ~grpreg
        ax=text(axa{1},is*60-60,ytxtcor,labelx{is},'units','pixels','fontweight','bold','color',markc(countc,:));
        ax=text(axa{1},is*60-60,ytxtcor-25,labelx{is},'units','pixels','fontweight','bold','color',markp(countc,:));
        else
        ax=text(axa{1},is*60-60,ytxtcor,labelx{is},'units','pixels','fontweight','bold','color',markc(countc,:));
        end
        countc=countc+1;
    else
        %put
        if ~grpreg
        ax=text(axa{1},is*60-60,ytxtcor,labelx{is},'units','pixels','fontweight','bold','color',markc(countp,:));
         ax=text(axa{1},is*60-60,ytxtcor-25,labelx{is},'units','pixels','fontweight','bold','color',markp(countp,:));
        else
        ax=text(axa{1},is*60-60,ytxtcor,labelx{is},'units','pixels','fontweight','bold','color',markp(countp,:));
        end
        countp=countp+1;
    end   
end
for is=1:length(sessnums)
    text(axa{1},is*75-60,ytxtcor+50,[num2str(sessnums(is)) ' ' sessmarktxt{is}],'interpreter','tex','units','pixels','fontweight','bold','color',[0 0 0]);
    ranj(is)=rand(1,1)*.8-.5;
end
text(axa{1},-60,ytxtcor-50,[ynam{yid} ' ' datype ' '  targevent],'interpreter','none','units','pixels','fontweight','bold','color',[0 0 0]);

   % for xleg=1:length(sessnums)
   %     labelx{xleg}=[num2str(sessnums(xleg))];  
   %     text(axa{1},100+xleg*50,axpos{1}(4)+30,labelx{xleg},'units','pixels','fontweight','bold','color',mark(xleg,:));
   % end
    curgrp=ygrp{yid};
    count=1;
    cc=1;
    cors={};
    plotdata=[];    
    sesslab=[];
    maxsites=20;
    boxposc=0:1/maxsites:1;
    boxposc=(boxposc+.01).*.85;
    boxposp=0:1/maxsites:1;
    boxposp=(boxposp+.05).*.85;
   boxpos=[];
   countc=1;
countp=1;
cboxc=1;
cboxp=1;

for ise=1:length(sessnums)
trialinfo=trialgrps(ise).trialinfo;
cursessids=find([sites.sessnum]==sessnums(ise));
curda={sites(cursessids).probeid};
dasites={sites(cursessids).site};
for ida=1:length(curda)
%each da separate subplot  
lfpsites={};
daregion=contains(curda(ida),'c');     %1=='c'
dasite=dasites(ida);
sitecolor=[0 0 1];
sitecolorn=[0 1 0];
sitenum=countc+countp+1;
if daregion
    %CN
    countc=countc+1;
    pid=1;
    sitecolor=markc(find(contains(cnsites,dasite)),:);
    sitecolorn=markp(find(contains(cnsites,dasite)),:);       %separate color group for p
    boxpos=boxposc(countc);
    if grpreg
        pid=1;      %group regions together, just separate conditions to separate plots
        sitecolor=markc(find(contains(cnsites,dasite)),:);       %color for positive R
        sitecolorn=markp(find(contains(cnsites,dasite)),:);       %separate color group for negative R
    else
        
        sitenum=countc;
    end
else
    %put
    countp=countp+1;
    pid=2;
    sitecolor=markc(find(contains(psites,dasite)),:);
    sitecolorn=markp(find(contains(psites,dasite)),:);       %separate color group for negative R
        boxpos=boxposp(countp);

    if grpreg
        pid=1;      %group regions together, just separate conditions to separate plots
        sitecolor=markc(find(contains(psites,dasite)),:);       %color for positive R
        sitecolorn=markp(find(contains(psites,dasite)),:);       %separate color group for negative R
    else
        sitenum=countp;
    end
end     
    
%get ch # for current sess and probe id
targsi=find([sites.sessnum]==sessnums(ise) & ...
    strcmp({sites.probeid},curda(ida)));
%if found ch# continue, otherwise skip loop
if ~isempty(targsi)
    ich=sites(targsi(1)).ch;        %ch # defined by fscv recording ch

 daid=find(strcmp(dapair,curda(ida)));
    if isempty(daid)
        disp(['empty ' targdasites(ida) ]);
    end
    lfptarg=lfppair{daid(1)};           %JUST FIRST GOOD LFP PAIR
    
for iy=1:length(curgrp)
%scroll parameters for group (ie. all timing variables)
%curytrials=getfield(xdata,'trials');        %trial nums for cur y variable
%curydata=getfield(xdata,curgrp{iy});       %y data
%for 'reward', both big and small sess types
%separate ax da reg
%text(axa{icond},-100,axpos{itype}(4)+20,condlabels{icond},'units','pixels','fontweight','bold','color',[0 0 0]);
%get trials corresponding to groups of trial conditions in current group ttypes{icond}
%get trial types (big/small/break) for current group
targtrialtypes=[any(strcmp('big',ttypes))  ...
    any(strcmp('small',ttypes)) ...
    any(strcmp('targetbreak',ttypes)) ...
    any(strcmp('fixbreak',ttypes)) ]; %get logical array for if want big/small/targ/fix
targtrialtypes=find(targtrialtypes==1);
nottypes=find(~strcmp(ttypes,'big') & ...
    ~strcmp(ttypes,'small') &...
    ~strcmp(ttypes,'targetbreak') &...
   ~strcmp(ttypes,'fixbreak'));
%get trial #'s for current condition, ie. not trial type identifier
%get da values for current group condition
for ib=1:length(bparams)
%scroll b params to get data using trials found above
bvals{pid}{sitenum}{ib}=[];
btrials{ib}=[];
yvals{pid}{sitenum}{ib}=[];
ytrials={};

for itype=targtrialtypes
    trialnumsb=[];
    trialnums=[];
    targrow=find((contains({xinfo.siteda},curda(ida)) & ...
    contains({xinfo.sitelfp},lfptarg) & ...
    strcmp({xinfo.event},targevent)) & ...
    contains({xinfo.sessionid},num2str(sessnums(ise))) &...
    contains({xinfo.sessiontype},sessiontypes(itype))==1);
if ~isempty(targrow)
    targrowb=find(contains({binfo.siteda},curda(ida)) & ...
    strcmp({binfo.sitelfp},bparams{ib}) & ...
    strcmp({binfo.event},targevent) & ...
    contains({binfo.sessionid},num2str(sessnums(ise))) &...
    contains({binfo.sessiontype},sessiontypes(itype))==1);
for tt=nottypes
    curcond=ttypes{tt};
    targt=find(contains(trialinfo(itype).names,curcond)==1);
    trialnums=[trialnums trialinfo(itype).nums{targt}]; 
end
trialnums=unique(trialnums);        %trial nums for current trial type big/small
%get b data for all bparams
trialnumsb=trialnums;          %intersect with those in all bparams to get trials exist in all groups
%get trials that exist in both bparams and da/lfp selected trials
bdata=getfield(binfo(targrowb),datype);    %'dapos' or 'daneg' types
trialnumsb=intersect(bdata.trials,trialnumsb);
%get da/lfp data, same for all bparams so only get once
if ~isempty(targrow)
    xdata=getfield(xinfo(targrow),datype);    %'dapos' or 'daneg' types
    %trialids=find(ismember(xdata.trials,trialnumsb)==1); %trial ids for targeted metric
    %ytrials=[ytrials intersect(xdata.trials,trialnums)];
    ytrials{itype}=intersect(xdata.trials,trialnumsb);
    trialids=find(ismember(xdata.trials,ytrials{itype})==1); %trial ids for targeted metric
    davals=getfield(xdata,curgrp{iy});    %all da values for good trials for sess/type
    yvals{pid}{sitenum}{ib}=[yvals{pid}{sitenum}{ib} davals(trialids)];         %da values for specific group of conditions (ie. trial nums)
end
if length(targrowb)>1
    warning('more than one targ row b')
end
if ~isempty(targrowb)
    bdata=getfield(binfo(targrowb),datype);    %'dapos' or 'daneg' types
    trialidsb=find(ismember(bdata.trials,ytrials{itype})==1); %trial ids for targeted metric
tempvals=getfield(bdata,bparamtarg);    %all da values for good trials for sess/type
bvals{pid}{sitenum}{ib}=[bvals{pid}{sitenum}{ib} tempvals(trialidsb)];         %da values for specific group of conditions (ie. trial nums)
end
end
end

%cor curr bparam for all condition variables with yvar data
ycordata=yvals{pid}{sitenum}{ib};
bcordata=bvals{pid}{sitenum}{ib};
cors=setfield(cors,{cc},'ydata',ycordata);
cors=setfield(cors,{cc},'yvar',curgrp{iy});
cors=setfield(cors,{cc},bparams{ib},bcordata);
cors=setfield(cors,{cc},'event',targevent);
cors=setfield(cors,{cc},'session',sessnums(ise));
cors=setfield(cors,{cc},'seltrials',ytrials);        
rcor=[];
pcor=[];
if ~isempty(ycordata) && ~isempty(bcordata)
    [rcor,pcor]=corr(bcordata',ycordata','rows','complete');
end
cors=setfield(cors,{cc},'rcor',rcor);
cors=setfield(cors,{cc},'pcor',pcor);
%plot if sig (p<=0.05)
if pcor<=0.05
    m=sessmark{ise};
    c=sitecolor;
    fsize=40*log10(1/pcor);
    if rcor<0
       c=sitecolorn;
    end
    a=scatter(axa{pid},ib-boxpos,iy-boxpos,...
        fsize,m,'markeredgecolor',c,'MarkerEdgeAlpha',.5,'linewidth',1.5);
end
cc=cc+1;
end
end
end
end
end

save([savepath 'behcor_' ynam{yid} '_' datype '_'  targevent '_' condfull...
    ],'cors');
savefig(figsess,[savepath 'behcor_' ynam{yid} '_' datype '_' targevent '_' condfull]);
saveas(figsess,[savepath 'behcor_' ynam{yid} '_' datype '_' targevent '_' condfull],'tif')
delete(findall(figsess,'type','text')) 
   

end
end

