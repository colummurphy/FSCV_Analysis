function plotxtracesavgov(xinfo,trialinfo,plotparam,varargin)
%2/17/2019 found something wrong with plotting left/right conditions, seems
%to just plto BOTH --this was fixed in plottracesummary.m function
%plot bar scatter plots of properties of cross-covariance lags/waveforms
%12/31/2018 updated for timing characteristics from new xclust 
%same as plotx, but plot big/small/targ on same plots
%for plotting multiple sessions
%1/3/2018, udpates as in plotx/plotxallsess
%PLOT TRACES for DA GROUPS, mean
%OVERLAY BIG VS SAMLL and defined groups
rate=10;      %10hz default for downsampled da/lfp from xcov
win=[-2 2];     %+/-2 s from aln idx
interval=1;       %in seconds
ratelfp=1000;
argnum=1;
datype='dapos';
evtype={};
sesstypes={};
trtypes={'all'};
binfo=[];
bplot=0;
ntype=[];
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'daneg'
            datype='daneg';
        case 'dareb'
            datype='dareb';
        case 'daall'
            datype='daall';
        case 'evtype'
            %user selects what event to plot 
            %ie interfix or intertarg
            argnum=argnum+1;
            evtype=varargin{argnum};
        case 'sesstypes'
            %user selects sess types for comparison (big v small..)
            argnum=argnum+1;
            sesstypes=varargin{argnum};
        case 'trialtypes'
            %select sleft/sright etc.
            argnum=argnum+1;
            trtypes=varargin{argnum};
        case 'binfo'
            %plot beh eye data
            argnum=argnum+1;
            binfo=varargin{argnum};
            bplot=1;
        case 'nest'
            %select within trialtype defined another nested group
            %ie trialtype of 'sleft', and want only 'phase1' of these
            argnum=argnum+1;
            ntype=varargin{argnum};
        case 'win'
            %user specified window +/- s
            argnum=argnum+1;
            win=varargin{argnum};
            
    end
    argnum=argnum+1;
end
fontsize=10;
savepath=plotparam.savepath;
%{
if ~strcmp(datype,'dapos')
    savepath=[savepath(1:end-1) '_' datype filesep];
    if ~isdir(savepath)
        mkdir(savepath);
    end
end
%}
mark=[0 .2 .7; 0 .5 .3; 0 .7 .1];
markr=[.2 .7 .0; .4 .3 0; .7 .1 0];
mark=[0 0.1 .7; 0.8 0 .5; 0 .7 0;.7 0 0; 0.4 0 .7];
markr=[1 .5 .0; .3 .3 0; .9 .1 0];
mark=[0 0 0; 1 0 0; 0 .7 0;.7 0 0; 0.4 0 .7];

tsx=[win(1):1/rate:win(2)];
tsxlfp=[win(1):1/ratelfp:win(2)];
lfpchs=plotparam.lfpchs;
%get lfp ch groups p & c
pgroup=find(contains(lfpchs,'p')==1);
cgroup=find(contains(lfpchs,'c')==1);
dasites={xinfo(1:end).siteda};
lfprows={xinfo(1:end).sitelfp};
%lfpchs=unique(lfprows);
targdasites=unique(dasites);
figpos=[50,50,700,700];
%figsess=figure('color',[1 1 1]);
%set(figsess,'position',figpos);
figsess=figure('visible','off');     %figure for each channel
if ispc
figsess=figure('visible','on');     %figure for each channel
end
set(figsess,'position',figpos,'color',[1 1 1]);

set(0,'CurrentFigure',figsess);    %set figure handle to current figure
axa={};
axsiz=[500 500];
off=150;
mar=100;
eventrows={xinfo(1:end).event};
eventtypes=unique(eventrows);
if ~isempty(evtype)
    if iscell(evtype)
        %multiple types ie. {'interfix','intertarg'}
        eventtypes=evtype;
    else
        %single string
        eventtypes={};
        eventtypes{1}=evtype;
    end    
end
sessiontypes=unique({xinfo(1:end).sessiontype});
if ~isempty(sesstypes)
    if iscell(sesstypes)
        %multiple types ie. {'bigreward','smallreward'}
        sessiontypes=sesstypes;
    else
        %single string
        sessiontypes={};
        sessiontypes{1}=sesstypes;
    end  
end

for ida=1:length(targdasites)
    %each da separate subplot  
    lfpsites={};
    daregion=contains(targdasites(ida),'c');     %1=='c'
    if daregion==1
        lfpsites=lfpchs(cgroup);
    else
        lfpsites=lfpchs(pgroup);
    end
    if ~isempty(binfo)
        lfpsites={};
    %plot only eye data in comparison with da, not lfp
    lfpsites={'eyex'};    
end

for ilfp=1:length(lfpsites)    
for ievent=1:length(eventtypes)
    scalesda=[];
    scaleslfp=[];
    maxda=[];
    minda=[];
    maxlfp=[];
    minlfp=[];
    clf(figsess,'reset');
    set(figsess,'color',[1 1 1]);
    iplot=1;
    axa{iplot}=subplot(1,length(sessiontypes),iplot);   hold(axa{iplot},'on');
            set(axa{iplot},'units','pixels');
    axpos{iplot}=get(axa{iplot},'position');         
    set(axa{iplot},'Units','Pixels','Position', [axsiz(1)*(iplot-1)+off*(iplot-1)+mar axpos{iplot}(2) axsiz(1) axsiz(2)]);

    lc=1;
    legtext={};
for itype=1:length(sessiontypes)
    
    %for each pair, event, plot da & lfp signals side by side trial by
    %trial for selected pos/neg/reb groups
    targrow=find((contains({xinfo.siteda},targdasites(ida)) & ...
        strcmp({xinfo.event},eventtypes(ievent))) & ...
        contains({xinfo.sessiontype},sessiontypes(itype))==1);
    if ~isempty(targrow)
    targrow=targrow(1);
    if ~bplot
     targrow=find((contains({xinfo.siteda},targdasites(ida)) & ...
        contains({xinfo.sitelfp},lfpsites(ilfp)) & ...
        strcmp({xinfo.event},eventtypes(ievent))) & ...
        contains({xinfo.sessiontype},sessiontypes(itype))==1);
    end
    for itrial=1:length(trtypes)
         trialtypes=trialinfo(itype); %get trial groups for session   
         ttype=find(ismember(trialtypes.names,trtypes{itrial})==1);
        typename=trialtypes.names{ttype};
        %curdata=xinfo(targrow).dapos;
        curdata=getfield(xinfo(targrow),datype);
        seltrials=...
        find(ismember(curdata.trials,trialtypes.nums{ttype})==1);
    %curdata=xinfo(targrow).dapos;
    curdata=getfield(xinfo(targrow),datype);
    alnidx=curdata.mididx;
    wins=[alnidx+win(1)*rate:alnidx+win(2)*rate];       %time window to plot
    da=curdata.datracesaln(seltrials,wins);                             %signals aligned to aln idx
    lfp=[];
    if ~isempty(curdata)
    if bplot
        targbrow=find(strcmp({binfo.event},eventtypes(ievent)) & ...
        strcmp({binfo.sessiontype},sessiontypes(itype)));
        curbdata=[];
        seltrialsb=[];
         trialnumsgrp=intersect(curdata.trials,trialtypes.nums{ttype});
        if strcmp(eventtypes{ievent},'interfix')
             seltrialsb=...
                find(ismember(binfo(targbrow).seltrials,trialnumsgrp)==1);
            curbdata=getfield(binfo(targbrow),'eyextrace');
        elseif strcmp(trtypes{itrial},'sleft')
            seltrialsb=...
                find(ismember(binfo(targbrow).seltrialsl,trialnumsgrp)==1);
            curbdata=getfield(binfo(targbrow),'leyextrace');
        elseif strcmp(trtypes{itrial},'sright')
            seltrialsb=...
                find(ismember(binfo(targbrow).seltrialsr,trialnumsgrp)==1);
            curbdata=getfield(binfo(targbrow),'reyextrace');
        end
        %data should already be aligned
        %need lfp rate since eye data not downcoverted in xcov, not used in
        %xcov
        alignidxb=median(1:size(curbdata,2));
        winsb=[win(1)*ratelfp+alignidxb:win(2)*ratelfp+alignidxb]; 
        lfp=curbdata(seltrialsb,winsb).*1e3;
    else
            lfp=curdata.lfptracesaln(seltrials,wins);                             %signals aligned to aln idx
    end
        
    scaleslfp(itype)=nanmean(nanstd(lfp,[],2));
    scalesda(itype)=nanmean(nanstd(da,[],2));
    daavgs(itype,:)=nanmean(da,1);
    lfpavgs(itype,:)=nanmean(lfp,1);
    dastd=nanstd(da,[],1);
    lfpstd=nanstd(lfp,[],1);
    daci(itype,:)=dastd./sqrt(size(da,1))*1.96;
    lfpci(itype,:)=lfpstd./sqrt(size(lfp,1))*1.96;
    
    maxda(itype)=max(daavgs(itype,:))+daci(itype,(find(daavgs(itype,:)==max(daavgs(itype,:)))));    
    minda(itype)=min(daavgs(itype,:))-daci(itype,(find(daavgs(itype,:)==min(daavgs(itype,:)))));
    maxlfp(itype)=max(lfpavgs(itype,:))+lfpci(itype,(find(lfpavgs(itype,:)==max(lfpavgs(itype,:)))));    
    minlfp(itype)=min(lfpavgs(itype,:))-lfpci(itype,(find(lfpavgs(itype,:)==min(lfpavgs(itype,:)))));
    end
    end

end 
end
%set scales uniform for all types big,small/targ
limsmaxd=round(nanmean(scalesda)*15)/10;
limsmind=-round(nanmean(scalesda)*15)/10;
if min(minda)<limsmind
    limsmind=floor(min(minda)*10)/10;
end
if max(maxda)>limsmaxd
    limsmaxd=ceil(max(maxda)*10)/10;
end 
limsmaxl=round(nanmean(scaleslfp)*250)/100;
limsminl=round(nanmean(scaleslfp)*100)/100;
if bplot
    limsminl=-round(nanmean(scaleslfp)*250)/100;
    if max(maxlfp)<0
        limsmaxl=max(maxlfp);
    end
end
if min(minlfp)<limsminl
    limsminl=floor(min(minlfp)*100)/100;
end
if max(maxlfp)>limsmaxl
    limsmaxl=ceil(max(maxlfp)*100)/100;
end    
savelab=[];
cc=1;
for itype=1:length(sessiontypes)   
    %for each pair, event, plot da & lfp signals side by side trial by
    %trial for selected pos/neg/reb groups
    targrow=find((contains({xinfo.siteda},targdasites(ida)) & ...
        strcmp({xinfo.event},eventtypes(ievent))) & ...
        contains({xinfo.sessiontype},sessiontypes(itype))==1);
    if ~isempty(targrow)
    targrow=targrow(1);
    if ~bplot
     targrow=find((contains({xinfo.siteda},targdasites(ida)) & ...
        contains({xinfo.sitelfp},lfpsites(ilfp)) & ...
        strcmp({xinfo.event},eventtypes(ievent))) & ...
        contains({xinfo.sessiontype},sessiontypes(itype))==1);
    end
for itrial=1:length(trtypes)
         trialtypes=trialinfo(itype); %get trial groups for session   
         ttype=find(ismember(trialtypes.names,trtypes{itrial})==1);
        typename=trialtypes.names{ttype};

    %curdata=xinfo(targrow).dapos;
    curdata=getfield(xinfo(targrow),datype);
    seltrials=...
    find(ismember(curdata.trials,trialtypes.nums{ttype})==1);
    alnidx=curdata.mididx;
    wins=[alnidx+win(1)*rate:alnidx+win(2)*rate];       %time window to plot
    da=curdata.datracesaln(seltrials,wins);                             %signals aligned to aln idx
    lfp=[];
    tsplot=tsx;
    tsplotr=tsx;
    if bplot
        targbrow=find(strcmp({binfo.event},eventtypes(ievent)) & ...
        contains({binfo.sessiontype},sessiontypes(itype))==1);
        curbdata=[];
       seltrialsb=[];
       trialnumsgrp=intersect(curdata.trials,trialtypes.nums{ttype});
        if strcmp(eventtypes{ievent},'interfix')
             seltrialsb=...
    find(ismember(binfo(targbrow).seltrials,trialnumsgrp)==1);
            curbdata=getfield(binfo(targbrow),'eyextrace');
        elseif strcmp(trtypes{itrial},'sleft')
            seltrialsb=...
    find(ismember(binfo(targbrow).seltrialsl,trialnumsgrp)==1);
            curbdata=getfield(binfo(targbrow),'leyextrace');
        elseif strcmp(trtypes{itrial},'sright')
            seltrialsb=...
    find(ismember(binfo(targbrow).seltrialsr,trialnumsgrp)==1);
            curbdata=getfield(binfo(targbrow),'reyextrace');
        end
        %winsb=[win(1)*ratelfp:win(2)*ratelfp]; 
        alignidxb=median(1:size(curbdata,2));
        winsb=[win(1)*ratelfp+alignidxb:win(2)*ratelfp+alignidxb]; 
        lfp=curbdata(seltrialsb,winsb).*1e3;
        tsplotr=tsxlfp;
    else
            lfp=curdata.lfptracesaln(seltrials,wins);                             %signals aligned to aln idx
    end    
    daavg=nanmean(da,1);
    lfpavg=nanmean(lfp,1);
    dastd=nanstd(da,[],1);
    lfpstd=nanstd(lfp,[],1);
    daci=dastd./sqrt(size(da,1))*1.96;
    lfpci=lfpstd./sqrt(size(lfp,1))*1.96;
    mididx=median(1:length(wins));
    freqband=getfield(xinfo(targrow),'freq'); 
   
    %plot da left
    yyaxis(axa{1},'left');
    aa=plot(axa{1},tsplot,([-daci; daci]+daavg)','-','linewidth',1,'color',mark(cc,:));
    aa(1).Color(4)=.2;      %1st line
    aa(2).Color(4)=.2;      %2nd line
    am=plot(axa{1},tsplot,daavg','-','linewidth',1.5,'color',mark(cc,:));
    am.Color(4)=.6;
     legtext{lc}=[sessiontypes{itype} ' ' trtypes{itrial} ' da'];
     lc=lc+1;
    %plot lfp right
    yyaxis(axa{1},'right');    
    aa=plot(axa{1},tsplotr,([-lfpci; lfpci]+lfpavg)','--','linewidth',1,'color',mark(cc,:));
    aa(1).Color(4)=.2;      %1st line
    aa(2).Color(4)=.2;      %2nd line
    am=plot(axa{1},tsplotr,lfpavg','--','linewidth',1.5,'color',mark(cc,:));
    am.Color(4)=.6;
    legtext{lc}=[sessiontypes{itype} ' ' trtypes{itrial} ' lfp (---)'];
    lc=lc+1;    
%    curtext=legtext{ileg};
if ~isempty(trtypes)
    savelab=[savelab sessiontypes{itype}(1:3) '_' trtypes{itrial}(1:3) '_'];
else
    savelab=[savelab sessiontypes{itype}(1:3) '_' ];
end
    cc=cc+1;
end
    end  
end
 yyaxis(axa{1},'left');
    set(axa{1},'ycolor',mark(1,:));
    clabel='DA';            
    ylabel(axa{1},clabel)
        set(axa{1},'ylim',[limsmind limsmaxd])
xticklabels=min(win):interval:max(win);
    xticklabels=round(xticklabels.*rate)./rate;
    xticklabels=num2str(xticklabels');
    xticks=1:round(interval*rate):size(da,2);
  %  set(axa{1},'xtick',xticks,'xticklabel',xticklabels);        
    xlabel(axa{1},'time (s)')
   % set(axa{1},'xlim',[win(1) win(2)]);     

yyaxis(axa{1},'right');
set(axa{1},'ycolor',mark(1,:));
    clabel_lfp=['beta-LFP ' num2str(freqband(1)) ...
                 '-' num2str(freqband(2)) ' hz'];               
    hy=ylabel(axa{1},clabel_lfp);
    %rotate right yaxis
    set(hy,'rotation',270,'units','pixels');
    hypos=get(hy,'position');
    set(hy,'position',[hypos(1)+25,hypos(2),hypos(3)]);  
    set(axa{1},'ylim',[limsminl limsmaxl])
    set(findall(axa{1},'-property','FontSize'),'FontSize',fontsize)
%legend text
for ileg=1:length(legtext)
    legc=[];
    if mod(ileg,2)
        legc=mark(ceil(ileg/2),:);
    else
        legc=mark(ceil(ileg/2),:);
    end
text(axa{1},400,axpos{1}(4)+50-25*ileg,legtext{ileg},'units','pixels',...
    'fontsize',fontsize+1,'fontweight','bold','color',legc);
end


    %title text
titletext=[datype ' | ' ...
        eventtypes{ievent} ' | da ' targdasites{ida} ...
        ' | lfp beta ' lfpsites{ilfp}];  
text(axa{1},0,axpos{1}(4)+10,titletext,'units','pixels',...
    'fontsize',fontsize+2,'fontweight','bold');

savedir=[savepath 'allsess_tracesavg_ov'  filesep ...
    targdasites{ida} 'x' lfpsites{ilfp} filesep];
if ~isdir(savedir)
    mkdir(savedir);
end
savename=[savedir  datype '_' eventtypes{ievent} ...
    '_' savelab targdasites{ida} 'x' lfpsites{ilfp}];

savefig(figsess,savename);
saveas(figsess,savename,'tif')
delete(findall(figsess,'type','text')) 


end
end
end

close all;

end