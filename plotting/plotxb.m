function plotdata=plotxb(xinfo,plotparam,varargin)
%plot bar scatter plots of properties of cross-covariance lags/waveforms
%according to burst groups, only look at "all" trials around specified
%events, include behavior
label='xinfo_';
argnum=1;
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'label'
            argnum=argnum+1;
            label=varargin{argnum};
    end
    argnum=argnum+1;
end
axplot={};
figsess={};
fontsize=8;
maxvar=0.02;
maxvar2=0.02;
savepath=plotparam.savepath;
betaband=1;
if str2num(savepath(end-1))==2
    %beta band 2
    betaband=2;
end
if str2num(savepath(end-1))==3
    %beta band 2
    betaband=3;
end
lfpchs=plotparam.lfpchs;
%get lfp ch groups p & c
pgroup=find(contains(lfpchs,'p')==1);
cgroup=find(contains(lfpchs,'c')==1);
maxplots=max(length(pgroup),length(cgroup));
dasites={xinfo(1:end).siteda};
lfprows={xinfo(1:end).sitelfp};
%lfpchs=unique(lfprows);
targdasites=unique(dasites);
figpos=[50,50,1200,800];
figsess=figure('position',figpos,'color',[1 1 1]);
set(0,'CurrentFigure',figsess);    %set figure handle to current figure
axa={};
axb={};

eventrows={xinfo(1:end).event};
eventtypes=unique(eventrows);

sessionrows={xinfo(1:end).sessiontype};
sessiontypes=unique({xinfo(1:end).sessiontype});
sessperiods=unique({xinfo(1:end).sessionperiod});       %burst vs no burst
sessperiods=plotparam.trialtypes.names;

sessperiodrows={xinfo(1:end).sessionperiod};
containsnumrows=regexp(sessperiodrows{1},'\d');
if betaband>=2 && ~isempty(containsnumrows)
    for ib=1:length(sessperiods)
    sessperiods{ib}=[sessperiods{ib} num2str(betaband)];
    end
end

posids={xinfo(1:end).posids};       %positive cluster trials
negids={xinfo(1:end).negids};       %positive cluster trials


target_lrt={xinfo(1:end).target_lrt};
target_rrt={xinfo(1:end).target_rrt};
fix_rt={xinfo(1:end).fix_rt};
leyed={xinfo(1:end).leyed};
reyed={xinfo(1:end).reyed};
eyed={xinfo(1:end).eyed};
pulse={xinfo(1:end).pulse};
lickpre={xinfo(1:end).lickpre};
lickpost={xinfo(1:end).lickpost};
    
    
mark=[0 .2 .5; .2 .5 0; .5 0 .2];

for ida=1:length(targdasites)
    %each da separate subplot  
    lfpsites={};
    daregion=contains(targdasites(ida),'c');     %1=='c'
    if daregion==1
        lfpsites=lfpchs(cgroup);
    else
        lfpsites=lfpchs(pgroup);
    end
for ilfp=1:length(lfpsites)    
for ievent=1:length(eventtypes)
    %each da-lfp pair separate figure
    %& each window event type separate figure
   % axa=axplot{ilfp,ievent}.axa;
    %axb=axplot{ilfp,ievent}.axb;     
    targids=find((contains({xinfo.siteda},targdasites(ida)) & ...
        contains({xinfo.sitelfp},lfpsites(ilfp)) & ...
        contains({xinfo.event},eventtypes(ievent)))==1);
    aid=1;  
    clf(figsess,'reset');
    set(figsess,'color',[1 1 1]);
    for ip=1:2
        ip2=(ip-1)*2+1;
        axa{ip}=subplot(2,2,ip2);   hold(axa{ip},'on');
        axb{ip}=subplot(2,2,ip2+1); hold(axb{ip},'on');
    end
    numsitepairs=length(sessiontypes)*length(sessperiods);
    labels={};
    maxvar=0.02;
    maxvar2=0.02;
    titletext=[label ' | ' eventtypes{ievent} ' | ' targdasites{ida} ' x ' lfpsites{ilfp}];  
    %parfig=get(axa,'parent');
    %parpos=get(parfig,'position');
    %if ievent==1
    set(axa{1},'units','pixels')
    axpos=get(axa{1},'position');
    text(axa{1},-30,axpos(4)+30,titletext,'units','pixels','fontweight','bold');
    
    lrt=[];
    rrt=[];
    if strcmp(eventtypes{ievent},'interfix')
        leyed={xinfo(1:end).eyed};
        reyed=[];
        lrt=fix_rt;
        rrt=[];
    end
    if strcmp(eventtypes{ievent},'intertarg')
        lrt=target_lrt;
        rrt=target_rrt;
        leyed={xinfo(1:end).leyed};
        reyed={xinfo(1:end).reyed};
    end
        
%end
for isess=1:length(sessperiods)  
    sessids=find(strcmp(sessperiodrows,sessperiods{isess})==1);
    sesstargids=intersect(targids,sessids);
    for itype=1:length(sessiontypes)
        typeids=find(contains(sessionrows,sessiontypes{itype})==1);
        currid=intersect(sesstargids,typeids);    
        %eye rt
        if length(currid)>1
            currid=currid(1);
        end
        if isempty(currid)
            aid=aid+1;
            continue
        end
        %yyaxis(axa{1},'left');
        xdata=repmat(aid,1,length(lrt{currid}));
        ranjitter=rand(1,length(lrt{currid}))*.4;
        xa=xdata-ranjitter;  
        cirt=nanstd(lrt{currid})./sqrt(length(lrt{currid}))*1.96;
        ciline=[nanmean(lrt{currid})-cirt  nanmean(lrt{currid})+cirt];
        plot(axa{1},[aid-.25 aid-.25],ciline,'-k','linewidth',1,'color',[0.3 0.3 0.3]);
        scatter(axa{1},xa,lrt{currid},15,'o','markeredgecolor',mark(1,:),'MarkerEdgeAlpha',.3,'linewidth',1.5);
        if ~isempty(rrt)
            xdata=repmat(aid,1,length(rrt{currid}));
            ranjitter=rand(1,length(rrt{currid}))*.4;
            xaa=xdata+ranjitter;
            cirt=nanstd(rrt{currid})./sqrt(length(rrt{currid}))*1.96;
            ciline=[nanmean(rrt{currid})-cirt  nanmean(rrt{currid})+cirt];
            plot(axa{1},[aid+.25 aid+.25],ciline,'-k','linewidth',1,'color',[0.6 0.6 0.6]);
            scatter(axa{1},xaa,rrt{currid},15,'x','markeredgecolor',mark(1,:),'MarkerEdgeAlpha',.3,'linewidth',1.5);
        end
        
        %pupil diam
        %yyaxis(axb{1},'left');
        xdata=repmat(aid,1,length(leyed{currid}));
        ranjitter=rand(1,length(leyed{currid}))*.4;
        xa=xdata-ranjitter;  
        cirt=nanstd(leyed{currid})./sqrt(length(leyed{currid}))*1.96;
        ciline=[nanmean(leyed{currid})-cirt  nanmean(leyed{currid})+cirt];
        plot(axb{1},[aid-.25 aid-.25],ciline,'-k','linewidth',1,'color',[0.3 0.3 0.3]);
        scatter(axb{1},xa,leyed{currid},15,'o','markeredgecolor',mark(1,:),'MarkerEdgeAlpha',.3,'linewidth',1.5);
        if ~isempty(reyed)
            xdata=repmat(aid,1,length(reyed{currid}));
            ranjitter=rand(1,length(reyed{currid}))*.4;
            xaa=xdata+ranjitter;
            cirt=nanstd(reyed{currid})./sqrt(length(reyed{currid}))*1.96;
            ciline=[mean(reyed{currid})-cirt  mean(reyed{currid})+cirt];
            plot(axb{1},[aid+.25 aid+.25],ciline,'-k','linewidth',1,'color',[0.6 0.6 0.6]);
            scatter(axb{1},xaa,reyed{currid},15,'x','markeredgecolor',mark(1,:),'MarkerEdgeAlpha',.3,'linewidth',1.5);
        end
        
        %pulse
        %yyaxis(axa{2},'left');
        xdata=repmat(aid,1,length(pulse{currid}));
        ranjitter=rand(1,length(pulse{currid}))*.3;
        xa=xdata-ranjitter;  
        cirt=std(pulse{currid})./sqrt(length(pulse{currid}))*1.96;
        ciline=[mean(pulse{currid})-cirt  mean(pulse{currid})+cirt];
        plot(axa{2},[aid-.25 aid-.25],ciline,'-k','linewidth',1,'color',[0.3 0.3 0.3]);
        scatter(axa{2},xa,pulse{currid},15,'o','markeredgecolor',mark(1,:),'MarkerEdgeAlpha',.3,'linewidth',1.5);
        
        %lick
        %yyaxis(axb{2},'left');
        xdata=repmat(aid,1,length(lickpre{currid}));
        ranjitter=rand(1,length(lickpre{currid}))*.4;
        xa=xdata-ranjitter;  
        cirt=std(lickpre{currid})./sqrt(length(lickpre{currid}))*1.96;
        ciline=[mean(lickpre{currid})-cirt  mean(lickpre{currid})+cirt];
        plot(axb{2},[aid-.25 aid-.25],ciline,'-k','linewidth',1,'color',[0.3 0.3 0.3]);
        scatter(axb{2},xa,lickpre{currid},15,'o','markeredgecolor',mark(1,:),'MarkerEdgeAlpha',.3,'linewidth',1.5);
        xdata=repmat(aid,1,length(lickpost{currid}));
        ranjitter=rand(1,length(lickpost{currid}))*.4;
        xaa=xdata+ranjitter;
        cirt=std(lickpost{currid})./sqrt(length(lickpost{currid}))*1.96;
        ciline=[mean(lickpost{currid})-cirt  mean(lickpost{currid})+cirt];
        plot(axb{2},[aid+.25 aid+.25],ciline,'-k','linewidth',1,'color',[0.6 0.6 0.6]);
        scatter(axb{2},xaa,lickpost{currid},15,'x','markeredgecolor',mark(1,:),'MarkerEdgeAlpha',.3,'linewidth',1.5);
        
        labels{aid}=[sessiontypes{itype}(1:4) ' | ' sessperiods{isess}];          
       
        aid=aid+1;
    end    
end
    ylabel(axa{1},'reaction time left/right (s)')
    ylabel(axb{1},'pupil diameter left/right')
    ylabel(axa{2},'pulse')
    ylabel(axb{2},'lick pre/post')
for ip=1:2
    set(axa{ip},'xtick',1:numsitepairs,'xticklabel',labels);
    set(axa{ip},'xlim',[0 numsitepairs+1]);
    set(axa{ip},'xTickLabelRotation',90)
    set(axb{ip},'xtick',1:numsitepairs,'xticklabel',labels);
    set(axb{ip},'xlim',[0 numsitepairs+1]);    
    set(axb{ip},'xTickLabelRotation',90)
    set(findall(axa{ip},'-property','FontSize'),'FontSize',fontsize)
    set(findall(axb{ip},'-property','FontSize'),'FontSize',fontsize)
end
savefig(figsess,[savepath label eventtypes{ievent} '_beh_' targdasites{ida} 'x' lfpsites{ilfp}]);
saveas(figsess,[savepath label eventtypes{ievent} '_beh_' targdasites{ida} 'x' lfpsites{ilfp}],'tif')
delete(findall(figsess,'type','text')) 
        
end
end
end


end