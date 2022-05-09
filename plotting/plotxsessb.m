function plotdata=plotxsessb(xinfos,plotparam)
axplot={};
figsess={};
fontsize=8;
maxvar=0.02;
maxvar2=0.02;
savepath=plotparam.savepath;
lfpchs=plotparam.lfpchs;
targdasites=plotparam.dasites;
sessnums=plotparam.sessnums;
%get lfp ch groups p & c
pgroup=find(contains(lfpchs,'p')==1);
cgroup=find(contains(lfpchs,'c')==1);
maxplots=max(length(pgroup),length(cgroup));

figpos=[50,50,1400,800];
figsess{1}=figure('position',figpos,'color',[1 1 1]);

set(0,'CurrentFigure',figsess{1});    %set figure handle to current figure
axa={};
axb={};

xinfo=xinfos{1};
    eventrows={xinfo(1:end).event};
    eventtypes=unique(eventrows);
sessionrows={xinfo(1:end).sessiontype};
    sessiontypes=unique({xinfo(1:end).sessiontype});        
    sessperiods=plotparam.trialtypes.names;
    sessperiodrows={xinfo(1:end).sessionperiod};
    
    numdiv=50/length(xinfos);
    mark=cool;
    mark=mark-.4;
    mark(mark<0)=0;
    mark=mark(1:round(numdiv):end,:);
    markr=mark;
    sesslab=[];

for xx=1:length(sessnums)
    sesslab=[sesslab '_' num2str(sessnums(xx))];
end
    ranjitter=rand(length(sessnums),1)*.4-.2;

   %mark=[0 .2 .5; .2 .5 0; .5 0 .2];    
    
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
        labeled=0;

    %each da-lfp pair separate figure
    %& each window event type separate figure
   % axa=axplot{ilfp,ievent}.axa;
    %axb=axplot{ilfp,ievent}.axb; 
    clf(figsess{1},'reset');
    set(figsess{1},'color',[1 1 1]);   
    set(0,'CurrentFigure',figsess{1});
    for ip=1:2
        ip2=(ip-1)*2+1;
        axa{ip}=subplot(2,2,ip2);   hold(axa{ip},'on');
        axb{ip}=subplot(2,2,ip2+1); hold(axb{ip},'on');
    end
      
    numsitepairs=length(sessiontypes)*length(sessperiods);
    labels={};
    maxvar=0.02;
    maxvar2=0.02;
    titletext=[eventtypes{ievent} ' | ' targdasites{ida} ' x ' lfpsites{ilfp}];  
    set(axa{1},'units','pixels')
    axpos=get(axa{1},'position');
    text(axa{1},-30,axpos(4)+30,titletext,'units','pixels','fontweight','bold');
for xsess=1:length(xinfos)
    aid=1;  
    xinfo=xinfos{xsess};  
    posids={xinfo(1:end).posids};       %positive cluster trials
    negids={xinfo(1:end).negids};       %positive cluster trials   
    sessionrows={xinfo(1:end).sessiontype};
    sessiontypes=unique({xinfo(1:end).sessiontype});        
    sessperiods=plotparam.trialtypes.names;
    sessperiodrows={xinfo(1:end).sessionperiod};
target_lrt={xinfo(1:end).target_lrt};
target_rrt={xinfo(1:end).target_rrt};
fix_rt={xinfo(1:end).fix_rt};
leyed={xinfo(1:end).leyed};
reyed={xinfo(1:end).reyed};
eyed={xinfo(1:end).eyed};
pulse={xinfo(1:end).pulse};
lickpre={xinfo(1:end).lickpre};
lickpost={xinfo(1:end).lickpost};
    
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
    targids=find((contains({xinfo.siteda},targdasites(ida)) & ...
        contains({xinfo.sitelfp},lfpsites(ilfp)) & ...
        contains({xinfo.event},eventtypes(ievent)))==1);
    
for isess=1:length(sessperiods)  
    sessids=find(strcmp(sessperiodrows,sessperiods{isess})==1);
    sesstargids=intersect(targids,sessids);
    aids=(isess-1)*length(sessiontypes)+1:(isess)*length(sessiontypes);
    assignin('base','aids',aids);
    for itype=1:length(sessiontypes)
        typeids=find(contains(sessionrows,sessiontypes{itype})==1);
        currid=intersect(sesstargids,typeids);  
        %eye rt
         if xsess==1
            labels{aid}=[sessiontypes{itype}(1:4) ' | ' sessperiods{isess}];  
        end
        if isempty(currid)
            aid=aid+1;
            continue
        end
        if length(currid)>1
            aid=aid+1;
            currid=currid(1);
        end
       if isempty(lrt{currid}) 
           aid=aid+1;
            continue
        end
        if isnan(nanmean(lrt{currid})) 
            aid=aid+1;
            continue
        end
        scatter(axa{1},aid+ranjitter(xsess),nanmean(lrt{currid}),30,'o','markeredgecolor',mark(xsess,:),'MarkerEdgeAlpha',.7);
        if ~isempty(rrt)            
            scatter(axa{1},aid+ranjitter(xsess),nanmean(rrt{currid}),15,'x','markeredgecolor',mark(xsess,:),'MarkerEdgeAlpha',.3,'linewidth',1.5);
        end
        
        %pupil diam
        scatter(axb{1},aid+ranjitter(xsess),nanmean(leyed{currid}),30,'o','markeredgecolor',mark(xsess,:),'MarkerEdgeAlpha',.7);
        if ~isempty(reyed)            
            scatter(axb{1},aid+ranjitter(xsess),nanmean(reyed{currid}),15,'x','markeredgecolor',mark(xsess,:),'MarkerEdgeAlpha',.3,'linewidth',1.5);
        end       
                  

        
        %pulse
        scatter(axa{2},aid+ranjitter(xsess),nanmean(pulse{currid}),30,'o','markeredgecolor',mark(xsess,:),'MarkerEdgeAlpha',.7);

        %lick        
        scatter(axb{2},aid+ranjitter(xsess),nanmean(lickpre{currid}),30,'o','markeredgecolor',mark(xsess,:),'MarkerEdgeAlpha',.7);
        scatter(axb{2},aid+ranjitter(xsess),nanmean(lickpost{currid}),15,'x','markeredgecolor',mark(xsess,:),'MarkerEdgeAlpha',.7);

        aid=aid+1;
    end    
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
    set(axa{ip},'xgrid','on');
    set(axb{ip},'xtick',1:numsitepairs,'xticklabel',labels);
        set(axb{ip},'xlim',[0 numsitepairs+1]);
        set(axb{ip},'xgrid','on');  
    set(axb{ip},'xTickLabelRotation',90)   
set(findall(axa{ip},'-property','FontSize'),'FontSize',fontsize)
set(findall(axb{ip},'-property','FontSize'),'FontSize',fontsize)

end
%set legend
for xleg=1:length(xinfos)
    label{xleg}=[num2str(sessnums(xleg))];  
    axpos=get(axa{1},'position');
    text(axa{1},100+xleg*50,axpos(4)+30,label{xleg},'units','pixels','fontweight','bold','color',mark(xleg,:));
%    axpos=get(axc{1},'position');
   % text(axc{1},100+xleg*50,axpos(4)+30,label{xleg},'units','pixels','fontweight','bold','color',mark(xleg,:));

end
for ip=1

savefig(figsess{ip},[savepath 'multisessx_burst_beh_' sesslab '_' eventtypes{ievent} '_' targdasites{ida} 'x' lfpsites{ilfp} ]);
saveas(figsess{ip},[savepath 'multisessx_burst_beh_' sesslab '_' eventtypes{ievent} '_' targdasites{ida} 'x' lfpsites{ilfp} ],'tif')

delete(findall(figsess{ip},'type','text')) 
end

end
end
end


end