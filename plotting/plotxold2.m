function plotdata=plotx(xinfo,plotparam,varargin)
%plot bar scatter plots of properties of cross-covariance lags/waveforms
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
    %beta band 3 broad
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
sessperiods=unique({xinfo(1:end).sessionperiod});
sessperiods=plotparam.trialtypes.names;
sessperiodrows={xinfo(1:end).sessionperiod};
containsnumrows=regexp(sessperiodrows{1},'\d');
if betaband>=2 && ~isempty(containsnumrows)
    for ib=1:length(sessperiods)
    sessperiods{ib}=[sessperiods{ib} num2str(betaband)];
    end
end

xcovdashiftpos={xinfo(1:end).xcovdashiftpos};           %positive clus wfs'
xcovdashiftneg={xinfo(1:end).xcovdashiftneg};
posids={xinfo(1:end).posids};       %positive cluster trials
negids={xinfo(1:end).negids};       %positive cluster trials

lagsdefneg={xinfo(1:end).lagsdefneg};
lagsnegminids=1e9;      %very large number that data never reaches
for ii=1:size(lagsdefneg,2)
    if size(lagsdefneg{ii},1)<lagsnegminids
        lagsnegminids=size(lagsdefneg{ii},1);
    end
end
        
lagsdefpos={xinfo(1:end).lagsdefpos};
lagsneg={xinfo(1:end).lagsneg};  
lagspos={xinfo(1:end).lagspos};
    
ratiopos={xinfo(1:end).ratiopos};
rationeg={xinfo(1:end).rationeg};
varwf={xinfo(1:end).varwf};
varwfposaln={xinfo(1:end).varwfposaln};
varwfnegaln={xinfo(1:end).varwfnegaln};
    mark=[0 .2 .5; .2 .5 0; .5 0 .2];

cmark=[0 0 0];
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
    %end
for isess=1:length(sessperiods)  
    sessids=find(strcmp(sessperiodrows,sessperiods{isess})==1);
    sesstargids=intersect(targids,sessids);
    for itype=1:length(sessiontypes)
        typeids=find(contains(sessionrows,sessiontypes{itype})==1);
        currid=intersect(sesstargids,typeids);    
        if length(currid)>1
            %NEED TO FIX ISSUE
            currid=currid(1);
        end
        if isempty(currid)
            aid=aid+1;
            continue
        end
        yyaxis(axa{1},'left');
        xdata=repmat(aid,1,length(lagsdefpos{currid}(:,1)));
        ranjitter=rand(1,length(lagsdefpos{currid}(:,1)))*.5;
        xa=xdata-ranjitter;                
        cilags=std(lagsdefpos{currid}(:,1))./sqrt(length(lagsdefpos{currid}(:,1)))*1.96;
        ciline=[mean(lagsdefpos{currid}(:,1))-cilags  mean(lagsdefpos{currid}(:,1))+cilags];
        plot(axa{1},[aid aid],ciline,'-k','linewidth',1);
        scatter(axa{1},xa,lagsdefpos{currid}(:,1),15,'.','markeredgecolor',mark(1,:),'MarkerEdgeAlpha',.3,'linewidth',1.5);
        yyaxis(axa{1},'right');
        plot(axa{1},aid,var(lagsdefpos{currid}(:,1)),'o','color',mark(2,:));    
        yyaxis(axb{1},'left');
        %variances of only aligned pos xvar within positive cluster (posids)
        varpos=var(xcovdashiftpos{currid}(posids{currid},:),0,1);
        plot(axb{1},aid-.25,mean(varpos),'o','markersize',3,'color',mark(1,:));
        plot(axb{1},aid-.25,max(varpos),'^','markersize',5,'color',mark(1,:));
        %all wf's pos/neg included non clustered
        plot(axb{1},aid-.25,mean(varwf{currid}),'+','markersize',6,'color',mark(1,:));
        %NO DIFFERENCES RATIO TAKE OUT ratiopos{currid}
        yyaxis(axb{1},'right');                
        %NO DIFFERENCES IN COEFF TAKE OUT mean(lagsdefpos{currid}(:,2))
        klagspos=kurtosis(lagspos{currid}(:,1));
        scatter(axb{1},aid+.25,klagspos,30,'o','markeredgecolor',mark(2,:),'MarkerEdgeAlpha',.3,'linewidth',1.5);
         klagsposdef=kurtosis(lagsdefpos{currid}(:,1));
        scatter(axb{1},aid+.25,klagsposdef,30,'x','markeredgecolor',mark(2,:),'MarkerEdgeAlpha',.3,'linewidth',1.5);
       
        if max(varwfposaln{currid})>maxvar
            maxvar=max(varwfposaln{currid});
        end
        
        labels{aid}=[sessiontypes{itype}(1:4) ' | ' sessperiods{isess}];  
        
        %plot negs xvars on {2}
        yyaxis(axb{2},'left');
        yyaxis(axa{2},'left');
        xdata=repmat(aid,1,length(lagsdefneg{currid}(:,1)));
        ranjitter=rand(1,length(lagsdefneg{currid}(:,1)))*.5;
        cilags=std(lagsdefneg{currid}(:,1))./sqrt(length(lagsdefneg{currid}(:,1)))*1.96;
        ciline=[mean(lagsdefneg{currid}(:,1))-cilags  mean(lagsdefneg{currid}(:,1))+cilags];
        plot(axa{2},[aid aid],ciline,'-k','linewidth',1);
        scatter(axa{2},xdata-ranjitter,lagsdefneg{currid}(:,1),10,'sq','markeredgecolor',mark(1,:),'MarkerEdgeAlpha',.3,'linewidth',1.5);
        yyaxis(axa{2},'right');
        plot(axa{2},aid,var(lagsdefneg{currid}(:,1)),'d','color',mark(2,:));    
        %confidence intervals of xvar aligend to pos or neg peak instead
        %NO CI makes variances smaller misleadingly for higher ## trials
        varneg=var(xcovdashiftneg{currid}(negids{currid},:),0,1);
        plot(axb{2},aid-.25,mean(varneg),'sq','markersize',3,'color',mark(1,:))
        plot(axb{2},aid-.25,max(varneg),'^','markersize',5,'color',mark(1,:))
        yyaxis(axb{2},'right');     
        klagsneg=kurtosis(lagsneg{currid}(:,1));
        scatter(axb{2},aid+.25,klagsneg,30,'o','markeredgecolor',mark(2,:),'MarkerEdgeAlpha',.3,'linewidth',1.5);
        klagsnegdef=kurtosis(lagsdefneg{currid}(:,1));
        scatter(axb{2},aid+.25,klagsnegdef,30,'x','markeredgecolor',mark(2,:),'MarkerEdgeAlpha',.3,'linewidth',1.5);

        
        if max(varwfnegaln{currid})>maxvar2
            maxvar2=max(varwfnegaln{currid});
        end        
        aid=aid+1;
    end    
end
for ip=1:2
    set(axa{ip},'xtick',1:numsitepairs,'xticklabel',labels);
    set(axa{ip},'xlim',[0 numsitepairs+1]);
    set(axa{ip},'xTickLabelRotation',90)
    %title(axa{ip},eventtypes{ievent}) 
    yyaxis(axb{ip},'left');
    set(axb{ip},'ycolor',mark(1,:))
    set(axb{ip},'xtick',1:numsitepairs,'xticklabel',labels);
        set(axb{ip},'xlim',[0 numsitepairs+1]);
    if ip==1
    ylabel(axb{ip},'mean/max variance + defined & alnd xcov wfs')
    else
    ylabel(axb{ip},'mean/max variance - defined & alnd xcov wfs')
    end
    set(axb{ip},'xTickLabelRotation',90)
    yyaxis(axb{ip},'right');
    set(axb{ip},'ycolor',mark(2,:));
    %set(axb{1},'ylim',[0 maxvar]);
        %set(axb{1},'ylim',[0 maxvar2]);
    if ip==1
        hh=ylabel(axb{ip},'kurtosis + xcov peak lags all/defined');
    else
        hh=ylabel(axb{ip},'kurtosis - xcov peak lags all/defined');
    end
    set(hh,'rotation',270,'units','pixels');
    hhpos=get(hh,'position');
    set(hh,'position',[hhpos(1)+10,hhpos(2),hhpos(3)]);        
    yyaxis(axa{ip},'left');
    set(axa{ip},'ycolor',mark(1,:))
    set(axa{ip},'ylim',[-1 1])
    ylabel(axa{ip},'peak defined lags (s)')
    yyaxis(axa{ip},'right');
    set(axa{ip},'ycolor',mark(2,:))
    set(axa{ip},'ylim',[0 .1])
    hh=ylabel(axa{ip},'variance peak lag (s^2) defined');
    set(hh,'rotation',270,'units','pixels');
    hhpos=get(hh,'position');
    set(hh,'position',[hhpos(1)+10,hhpos(2),hhpos(3)]);
set(findall(axa{ip},'-property','FontSize'),'FontSize',fontsize)
set(findall(axb{ip},'-property','FontSize'),'FontSize',fontsize)
end
savefig(figsess,[savepath label eventtypes{ievent} '_' targdasites{ida} 'x' lfpsites{ilfp}]);
saveas(figsess,[savepath label eventtypes{ievent} '_' targdasites{ida} 'x' lfpsites{ilfp}],'tif')
delete(findall(figsess,'type','text')) 
        
end
end
end


end