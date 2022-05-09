function plotdata=plotxsess(xinfos,plotparam,varargin)
label='multisessx_';
argnum=1;
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'label'
            argnum=argnum+1;
            label=[label varargin{argnum}];
    end
    argnum=argnum+1;
end
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
figsess{2}=figure('position',figpos,'color',[1 1 1]);

set(0,'CurrentFigure',figsess{1});    %set figure handle to current figure
axa={};
axb={};
axc={};

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
    %each da-lfp pair separate figure
    %& each window event type separate figure
   % axa=axplot{ilfp,ievent}.axa;
    %axb=axplot{ilfp,ievent}.axb; 
    clf(figsess{1},'reset');
    set(figsess{1},'color',[1 1 1]);
    clf(figsess{2},'reset');
    set(figsess{2},'color',[1 1 1]);    
    set(0,'CurrentFigure',figsess{1});
    for ip=1:2
        ip2=(ip-1)*2+1;
        axa{ip}=subplot(2,2,ip2);   hold(axa{ip},'on');
        axb{ip}=subplot(2,2,ip2+1); hold(axb{ip},'on');
    end
        set(0,'CurrentFigure',figsess{2});
    for ip=1:2
        axc{ip}=subplot(2,2,ip);   hold(axc{ip},'on');
    end
    numsitepairs=length(sessiontypes)*length(sessperiods);
    labels={};
    maxvar=0.02;
    maxvar2=0.02;
    titletext=[eventtypes{ievent} ' | ' targdasites{ida} ' x ' lfpsites{ilfp}];  
    set(axa{1},'units','pixels')
    axpos=get(axa{1},'position');
    text(axa{1},-30,axpos(4)+30,titletext,'units','pixels','fontweight','bold');
    set(axc{1},'units','pixels')
    axpos=get(axc{1},'position');
    text(axc{1},-30,axpos(4)+30,titletext,'units','pixels','fontweight','bold');

for xsess=1:length(xinfos)
    aid=1;  
    xinfo=xinfos{xsess};

    xcovdashiftpos={xinfo(1:end).xcovdashiftpos};
    xcovdashiftneg={xinfo(1:end).xcovdashiftneg};
    lagsneg={xinfo(1:end).lagsneg};  
    lagspos={xinfo(1:end).lagspos};
    posids={xinfo(1:end).posids};       %positive cluster trials
    negids={xinfo(1:end).negids};       %positive cluster trials
    lagsdefneg={xinfo(1:end).lagsdefneg};  
    lagsdefpos={xinfo(1:end).lagsdefpos};
    varwf={xinfo(1:end).varwf};
    varwfposaln={xinfo(1:end).varwfposaln};
    varwfnegaln={xinfo(1:end).varwfnegaln};
    sessionrows={xinfo(1:end).sessiontype};
    sessiontypes=unique({xinfo(1:end).sessiontype});        
    sessperiods=plotparam.trialtypes.names;
    sessperiodrows={xinfo(1:end).sessionperiod};

    targids=find((contains({xinfo.siteda},targdasites(ida)) & ...
        contains({xinfo.sitelfp},lfpsites(ilfp)) & ...
        contains({xinfo.event},eventtypes(ievent)))==1);
    
for isess=1:length(sessperiods)  
    sessids=find(contains(sessperiodrows,sessperiods{isess})==1);
    sesstargids=intersect(targids,sessids);
    aids=(isess-1)*length(sessiontypes)+1:(isess)*length(sessiontypes);
    assignin('base','aids',aids);
    for itype=1:length(sessiontypes)
        typeids=find(contains(sessionrows,sessiontypes{itype})==1);
        currid=intersect(sesstargids,typeids);            
        %yyaxis(axa{1},'left');
        %ranjitter=rand(1,length(lagsdefpos{currid}(:,1)))*.5;
                if xsess==1
            labels{aid}=[sessiontypes{itype}(1:4) ' | ' sessperiods{isess}];  
        end
        if length(currid)>1
            currid=currid(1);
        end
        if isempty(currid)
            aid=aid+1;
            continue
        end
        scatter(axa{1},aid+ranjitter(xsess),mean(lagsdefpos{currid}(:,1)),30,'o','markeredgecolor',mark(xsess,:),'MarkerEdgeAlpha',.7);
        yyaxis(axb{1},'left');
        varpos=var(xcovdashiftpos{currid}(posids{currid},:),0,1);
        scatter(axb{1},aid+ranjitter(xsess),mean(varpos),30,'o','markeredgecolor',mark(xsess,:),'MarkerEdgeAlpha',.7);
        %scatter(axb{1},aid-.25,max(varpos),35,'^','markeredgecolor',mark(xsess,:),'MarkerEdgeAlpha',.7);
        %yyaxis(axb{1},'right');                

         %coeffpos=mean(lagsdefpos{currid}(:,2));     %mean peak neg
        %scatter(axb{1},aid+ranjitter(xsess),coeffpos,30,'x','markeredgecolor',mark(xsess,:),'MarkerEdgeAlpha',.7);

        %yyaxis(axa{2},'left');
        scatter(axa{2},aid+ranjitter(xsess),mean(lagsdefneg{currid}(:,1)),30,'o','markeredgecolor',mark(xsess,:),'MarkerEdgeAlpha',.7);
        %yyaxis(axa{2},'right');
       % scatter(axa{2},aid+ranjitter(xsess),var(lagsdefneg{currid}(:,1)),15,'.','markeredgecolor',mark(xsess,:),'MarkerEdgeAlpha',.7);   
        yyaxis(axb{2},'left');    
        varneg=var(xcovdashiftneg{currid}(negids{currid},:),0,1);
        
        scatter(axb{2},aid+ranjitter(xsess),mean(varneg),30,'o','markeredgecolor',mark(xsess,:),'MarkerEdgeAlpha',.7);
        %scatter(axb{2},aid-.25,max(varneg),35,'^','markeredgecolor',mark(xsess,:),'MarkerEdgeAlpha',.7);
        %yyaxis(axb{2},'right');
        %NO DIFFERENCES IN COEFF TAKE OUT 12/08
        %coeffneg=mean(lagsdefneg{currid}(:,2));     %mean peak neg
        %scatter(axb{2},aid+ranjitter(xsess),coeffneg,30,'x','markeredgecolor',mark(xsess,:),'MarkerEdgeAlpha',.7);
         yyaxis(axc{1},'left');  
        scatter(axc{1},aid+ranjitter(xsess),kurtosis(lagspos{currid}(:,1)),30,'o','markeredgecolor',mark(xsess,:),'MarkerEdgeAlpha',.7);
         scatter(axc{1},aid+ranjitter(xsess),kurtosis(lagsdefpos{currid}(:,1)),20,'x','markeredgecolor',mark(xsess,:),'MarkerEdgeAlpha',.7);
      
       % yyaxis(axc{1},'right');  
        %        scatter(axc{1},aid+ranjitter(xsess),mean(varwf{currid}),30,'.','markeredgecolor',mark(xsess,:),'MarkerEdgeAlpha',.7);

        scatter(axc{1},aid+ranjitter(xsess),var(lagsdefpos{currid}(:,1)),30,'x','markeredgecolor',mark(xsess,:),'MarkerEdgeAlpha',.7);  
                 yyaxis(axc{2},'left');  
        scatter(axc{2},aid+ranjitter(xsess),kurtosis(lagsneg{currid}(:,1)),30,'o','markeredgecolor',mark(xsess,:),'MarkerEdgeAlpha',.7);
        scatter(axc{2},aid+ranjitter(xsess),kurtosis(lagsdefneg{currid}(:,1)),20,'x','markeredgecolor',mark(xsess,:),'MarkerEdgeAlpha',.7);
  
        %yyaxis(axc{2},'right');   
       % scatter(axc{2},aid+ranjitter(xsess),var(lagsdefneg{currid}(:,1)),30,'x','markeredgecolor',mark(xsess,:),'MarkerEdgeAlpha',.7);   
        
        aid=aid+1;
    end    
end
end

for ip=1:2
    set(axa{ip},'xtick',1:numsitepairs,'xticklabel',labels);
    set(axa{ip},'xlim',[0 numsitepairs+1]);
    set(axa{ip},'xTickLabelRotation',90)
    set(axa{ip},'xgrid','on');
    %title(axa{ip},eventtypes{ievent}) 
    yyaxis(axb{ip},'left');
    set(axb{ip},'ycolor',mark(xsess,:))
    set(axb{ip},'xtick',1:numsitepairs,'xticklabel',labels);
        set(axb{ip},'xlim',[0 numsitepairs+1]);
        set(axb{ip},'xgrid','on');
        set(axc{ip},'xtick',1:numsitepairs,'xticklabel',labels);
    set(axc{ip},'xlim',[0 numsitepairs+1]);
    set(axc{ip},'xTickLabelRotation',90)
    set(axc{ip},'xgrid','on');
    
    ylabel(axb{ip},'mean/max variance grouped xcov waveforms')
    set(axb{ip},'xTickLabelRotation',90)
    %yyaxis(axb{ip},'right');
   % set(axb{ip},'ycolor',mark(2,:));
    %set(axb{1},'ylim',[0 maxvar]);
        %set(axb{1},'ylim',[0 maxvar2]);

   % hh=ylabel(axb{ip},'mean coeff xcov grouped peaks');
    %set(hh,'rotation',270,'units','pixels');
    %hhpos=get(hh,'position');
    %set(hh,'position',[hhpos(1)+10,hhpos(2),hhpos(3)]);        
    yyaxis(axa{ip},'left');
    set(axa{ip},'ycolor',mark(xsess,:))
    set(axa{ip},'ylim',[-.5 .5])
    ylabel(axa{ip},'def lags (s)')
    
    yyaxis(axc{ip},'left');
    set(axc{ip},'ycolor',mark(xsess,:))
    if ip==1
    ylabel(axc{ip},'kurtosis all + lags)')
    else
        ylabel(axc{ip},'kurtosis all - lags)')
    end
    yyaxis(axc{ip},'right');
    set(axc{ip},'ycolor',mark(xsess,:))
    hh=ylabel(axc{ip},'variance peak lag (s^2)');
    set(hh,'rotation',270,'units','pixels');
    hhpos=get(hh,'position');
    set(hh,'position',[hhpos(1)+10,hhpos(2),hhpos(3)]);
    
set(findall(axa{ip},'-property','FontSize'),'FontSize',fontsize)
set(findall(axb{ip},'-property','FontSize'),'FontSize',fontsize)
set(findall(axc{ip},'-property','FontSize'),'FontSize',fontsize)

end
%set legend
labelx={};
for xleg=1:length(xinfos)
    labelx{xleg}=[num2str(sessnums(xleg))];  
    axpos=get(axa{1},'position');
    text(axa{1},100+xleg*50,axpos(4)+30,labelx{xleg},'units','pixels','fontweight','bold','color',mark(xleg,:));
   % axpos=get(axc{1},'position');
  %  text(axc{1},100+xleg*50,axpos(4)+30,label{xleg},'units','pixels','fontweight','bold','color',mark(xleg,:));

end
for ip=1:2

savefig(figsess{ip},[savepath label sesslab '_' eventtypes{ievent} '_' targdasites{ida} 'x' lfpsites{ilfp} '_' num2str(ip)]);
saveas(figsess{ip},[savepath label sesslab '_' eventtypes{ievent} '_' targdasites{ida} 'x' lfpsites{ilfp} '_' num2str(ip)],'tif')

delete(findall(figsess{ip},'type','text')) 
end

end
end
end


end