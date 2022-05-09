function plotx(axplot,xinfo,plotparam)
fontsize=8;
maxvar=0.1;
lfpchs=plotparam.lfpchs;
%get lfp ch groups p & c
pgroup=find(contains(lfpchs,'p')==1);
cgroup=find(contains(lfpchs,'c')==1);
maxplots=max(length(pgroup),length(cgroup));

dasites={xinfo(1:end).siteda};
lfprows={xinfo(1:end).sitelfp};
%lfpchs=unique(lfprows);
%targdasites=plotparam.chnums;
%if isempty(targdasites)
    targdasites=unique(dasites);
%end

eventrows={xinfo(1:end).event};
eventtypes=unique(eventrows);

sessionrows={xinfo(1:end).sessiontype};
sessiontypes=unique({xinfo(1:end).sessiontype});
sessperiods=unique({xinfo(1:end).sessionperiod});
sessperiodrows={xinfo(1:end).sessionperiod};

lagsdefneg={xinfo(1:end).lagsdefneg};
lagsdefpos={xinfo(1:end).lagsdefpos};
ratiopos={xinfo(1:end).ratiopos};
rationeg={xinfo(1:end).rationeg};
varwf={xinfo(1:end).varwf};
varwfposaln={xinfo(1:end).varwfposaln};
varwfnegaln={xinfo(1:end).varwfnegaln};
    mark=[0 .2 .5; .2 .5 0; .5 0 .2];

cmark=[0 0 0];
for isess=1:length(sessperiods)
    %each session period (ie, all, phase 1, within separate figure
    %each da-lfp pair separate figure
    axa=axplot{isess}.axa;
    axb=axplot{isess}.axb;
for ii=1:length(targdasites)
    sessdasites=find(contains(sessperiodrows,sessperiods{isess})==1);
    targdarows=find(contains(dasites,targdasites(ii))==1);
    sessdarows=intersect(sessdasites,targdarows);
    daregion=contains(targdasites(ii),'c');     %1=='c'
    lfpsites={};
    if daregion==1
    lfpsites=lfpchs(cgroup);
    end
    if daregion==0
        lfpsites=lfpchs(pgroup);
    end
    colcount=1;
    bid=1;
    plotid=ii;
    numsitepairs=length(lfpsites)*length(sessiontypes)*length(sessperiods)*2;
    labels={};

    for ievent=1:length(eventtypes)
        maxvar=0.1;

        hold(axa{ievent,plotid},'on')
        hold(axb{ievent,plotid},'on')
        aid=1;    
        targevent=eventtypes{ievent};
        targeventrows=find(contains(eventrows,targevent)==1);
        for ilfp=1:length(lfpsites)     
            targlfprow=find(contains(lfprows,lfpsites{ilfp})==1);
            dalfprows=intersect(targlfprow,sessdarows);
            targrows=intersect(dalfprows,targeventrows);
            group=1;
            for isess=1:length(sessiontypes)
                yyaxis(axa{ievent,plotid},'left');
                targrow=intersect(targrows,find(contains(sessionrows,sessiontypes{isess})==1));
                xdata=repmat(aid,1,length(lagsdefpos{targrow}(:,1)));
                ranjitter=rand(1,length(lagsdefpos{targrow}(:,1)))*.5;
                xa=xdata-ranjitter;                
                cilags=std(lagsdefpos{targrow}(:,1))./sqrt(length(lagsdefpos{targrow}(:,1)))*1.96;
                ciline=[mean(lagsdefpos{targrow}(:,1))-cilags  mean(lagsdefpos{targrow}(:,1))+cilags];
                plot(axa{ievent,plotid},[aid aid],ciline,'-k','linewidth',1);
                scatter(axa{ievent,plotid},xa,lagsdefpos{targrow}(:,1),15,'.','markeredgecolor',mark(1,:),'MarkerEdgeAlpha',.3,'linewidth',1.5);
                yyaxis(axa{ievent,plotid},'right');
                plot(axa{ievent,plotid},aid,var(lagsdefpos{targrow}(:,1)),'o','color',mark(2,:));    
                yyaxis(axb{ievent,plotid},'left');
                plot(axb{ievent,plotid},aid,ratiopos{targrow},'o','markersize',5,'color',mark(1,:));
                yyaxis(axb{ievent,plotid},'right');                
                plot(axb{ievent,plotid},aid+.25,mean(varwf{targrow}),'+','markersize',6,'color',mark(2,:));
                plot(axb{ievent,plotid},aid+.25,max(varwf{targrow}),'^','markersize',6,'color',mark(2,:));
                plot(axb{ievent,plotid},aid+.25,mean(varwfposaln{targrow}),'o','markersize',4,'color',mark(2,:));
                plot(axb{ievent,plotid},aid+.25,max(varwfposaln{targrow}),'o','markersize',6,'color',mark(2,:));
                if max(varwfposaln{targrow})>maxvar
                maxvar=max(varwfposaln{targrow});
                end
                yyaxis(axb{ievent,plotid},'left');
                aid=aid+1;
                if ievent==1
                    labels{colcount}=[sessiontypes{isess}(1:4) ' | ' targdasites{ii} ' x ' lfpsites{ilfp}];  
                    if group~=1
                        %don't repeat unnecessary lfp labels
                     labels{colcount}=[sessiontypes{isess}(1:4) ];  
                    end
                    colcount=colcount+1;
                end
                group=group+1;                
            end
            for isess=1:length(sessiontypes)
                yyaxis(axa{ievent,plotid},'left');
                targrow=intersect(targrows,find(contains(sessionrows,sessiontypes{isess})==1));
                xdata=repmat(aid,1,length(lagsdefneg{targrow}(:,1)));
                ranjitter=rand(1,length(lagsdefneg{targrow}(:,1)))*.5;
                cilags=std(lagsdefneg{targrow}(:,1))./sqrt(length(lagsdefneg{targrow}(:,1)))*1.96;
                ciline=[mean(lagsdefneg{targrow}(:,1))-cilags  mean(lagsdefneg{targrow}(:,1))+cilags];
                plot(axa{ievent,plotid},[aid aid],ciline,'-k','linewidth',1);
                scatter(axa{ievent,plotid},xdata-ranjitter,lagsdefneg{targrow}(:,1),10,'sq','markeredgecolor',mark(1,:),'MarkerEdgeAlpha',.3,'linewidth',1.5);
                yyaxis(axa{ievent,plotid},'right');
                plot(axa{ievent,plotid},aid,var(lagsdefneg{targrow}(:,1)),'d','color',mark(2,:));    
                plot(axb{ievent,plotid},aid,rationeg{targrow},'sq','markersize',5,'color',mark(1,:))
               yyaxis(axb{ievent,plotid},'right');                
               plot(axb{ievent,plotid},aid+.25,mean(varwfnegaln{targrow}),'sq','markersize',4,'color',mark(2,:));
                plot(axb{ievent,plotid},aid+.25,max(varwfnegaln{targrow}),'sq','markersize',6,'color',mark(2,:));
                if max(varwfnegaln{targrow})>maxvar
                maxvar=max(varwfnegaln{targrow});
                end
                aid=aid+1;
                if ievent==1
                labels{colcount}=[sessiontypes{isess}(1:4) ' | ' targdasites{ii} ' x ' lfpsites{ilfp}];  
                     if group~=1
                        %don't repeat unnecessary lfp labels
                     labels{colcount}=[sessiontypes{isess}(1:4) ];  
                    end
                colcount=colcount+1;
                end
                group=group+1;
            end
        end
    end
    for ievent=1:length(eventtypes)
        set(axa{ievent,plotid},'xtick',1:numsitepairs,'xticklabel',labels);
        set(axa{ievent,plotid},'xTickLabelRotation',90)
        title(axa{ievent,plotid},eventtypes{ievent})        
        
        yyaxis(axb{ievent,plotid},'left');
        set(axb{ievent,plotid},'ycolor',mark(1,:))
          set(axb{ievent,plotid},'xtick',1:numsitepairs,'xticklabel',labels);
          ylabel(axb{ievent,plotid},'ratio def lags')
        set(axb{ievent,plotid},'xTickLabelRotation',90)
                yyaxis(axb{ievent,plotid},'right');
                set(axb{ievent,plotid},'ycolor',mark(2,:));
                set(axb{ievent,plotid},'ylim',[0 maxvar]);
        hh=ylabel(axb{ievent,plotid},'mean variance waveform');
        set(hh,'rotation',270,'units','pixels');
        hhpos=get(hh,'position');
        set(hh,'position',[hhpos(1)+10,hhpos(2),hhpos(3)]);
        
        yyaxis(axa{ievent,plotid},'left');
        set(axa{ievent,plotid},'ycolor',mark(1,:))
        set(axa{ievent,plotid},'ylim',[-1 1])
        ylabel(axa{ievent,plotid},'def lags (s)')
        yyaxis(axa{ievent,plotid},'right');
        set(axa{ievent,plotid},'ycolor',mark(2,:))
        set(axa{ievent,plotid},'ylim',[0 .1])
        hh=ylabel(axa{ievent,plotid},'variance peak lag (s^2)');
        set(hh,'rotation',270,'units','pixels');
                hhpos=get(hh,'position');

        set(hh,'position',[hhpos(1)+10,hhpos(2),hhpos(3)]);

    set(findall(axa{ievent,plotid},'-property','FontSize'),'FontSize',fontsize)
    set(findall(axb{ievent,plotid},'-property','FontSize'),'FontSize',fontsize)

    end
    
        
end
end

end