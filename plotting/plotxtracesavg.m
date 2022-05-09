function plotxtracesavg(xinfo,trialinfo,plotparam,varargin)
%plot bar scatter plots of properties of cross-covariance lags/waveforms
%12/31/2018 updated for timing characteristics from new xclust 
%same as plotx, but plot big/small/targ on same plots
%for plotting multiple sessions
%1/3/2018, udpates as in plotx/plotxallsess
%PLOT TRACES for DA GROUPS, mean
rate=10;      %10hz default for downsampled da/lfp from xcov
win=[-2 2];     %+/-2 s from aln idx
interval=1;       %in seconds

argnum=1;
datype='dapos';
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'daneg'
            datype='daneg';
        case 'dareb'
            datype='dareb';
        case 'daall'
            datype='daall';
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

lfpchs=plotparam.lfpchs;
%get lfp ch groups p & c
pgroup=find(contains(lfpchs,'p')==1);
cgroup=find(contains(lfpchs,'c')==1);
dasites={xinfo(1:end).siteda};
lfprows={xinfo(1:end).sitelfp};
%lfpchs=unique(lfprows);
targdasites=unique(dasites);
figpos=[50,50,1600,500];
figsess=figure('visible','off');     %figure for each channel
if ispc
figsess=figure('visible','on');     %figure for each channel
end
set(figsess,'position',figpos,'color',[1 1 1]);

set(0,'CurrentFigure',figsess);    %set figure handle to current figure
axa={};
axsiz=[350 350];
off=150;
mar=100;
eventrows={xinfo(1:end).event};
eventtypes=unique(eventrows);

sessiontypes=unique({xinfo(1:end).sessiontype});

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
    scalesda=[];
    scaleslfp=[];
    maxda=[];
    minda=[];
    maxlfp=[];
    minlfp=[];
    clf(figsess,'reset');
    set(figsess,'color',[1 1 1]);
for itype=1:length(sessiontypes)
    %for each pair, event, plot da & lfp signals side by side trial by
    %trial for selected pos/neg/reb groups
     targrow=find((contains({xinfo.siteda},targdasites(ida)) & ...
        contains({xinfo.sitelfp},lfpsites(ilfp)) & ...
        strcmp({xinfo.event},eventtypes(ievent))) & ...
        contains({xinfo.sessiontype},sessiontypes(itype))==1);
    if ~isempty(targrow)
    %curdata=xinfo(targrow).dapos;
    curdata=getfield(xinfo(targrow),datype);
    alnidx=curdata.mididx;
    wins=[alnidx+win(1)*rate:alnidx+win(2)*rate];       %time window to plot
    da=curdata.datracesaln(:,wins);                             %signals aligned to aln idx
    lfp=curdata.lfptracesaln(:,wins);                             %signals aligned to aln idx
    scaleslfp(itype)=nanmean(nanstd(curdata.lfptracesaln(:,wins),[],2));
    scalesda(itype)=nanmean(nanstd(curdata.datracesaln(:,wins),[],2));
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

    axa{itype}=subplot(1,length(sessiontypes),itype);   hold(axa{itype},'on');
            set(axa{itype},'units','pixels');
    axpos{itype}=get(axa{itype},'position');         
    set(axa{itype},'Units','Pixels','Position', [axsiz(1)*(itype-1)+off*(itype-1)+mar axpos{itype}(2) axsiz(1) axsiz(2)]);

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
if min(minlfp)<limsminl
    limsminl=floor(min(minlfp)*100)/100;
end
if max(maxlfp)>limsmaxl
    limsmaxl=ceil(max(maxlfp)*100)/100;
end    

for itype=1:length(sessiontypes)   
    %for each pair, event, plot da & lfp signals side by side trial by
    %trial for selected pos/neg/reb groups
     targrow=find((contains({xinfo.siteda},targdasites(ida)) & ...
        contains({xinfo.sitelfp},lfpsites(ilfp)) & ...
        strcmp({xinfo.event},eventtypes(ievent))) & ...
        contains({xinfo.sessiontype},sessiontypes(itype))==1);
    if ~isempty(targrow)
    %curdata=xinfo(targrow).dapos;
    curdata=getfield(xinfo(targrow),datype);
    alnidx=curdata.mididx;
    wins=[alnidx+win(1)*rate:alnidx+win(2)*rate];       %time window to plot
    da=curdata.datracesaln(:,wins);                             %signals aligned to aln idx
    lfp=curdata.lfptracesaln(:,wins);     
    daavg=nanmean(da,1);
    lfpavg=nanmean(lfp,1);
    dastd=nanstd(da,[],1);
    lfpstd=nanstd(lfp,[],1);
    daci=dastd./sqrt(size(da,1))*1.96;
    lfpci=lfpstd./sqrt(size(lfp,1))*1.96;
    mididx=median(1:length(wins));
    freqband=getfield(xinfo(targrow),'freq'); 
   
    %plot da left
    yyaxis(axa{itype},'left');
    set(axa{itype},'ycolor',mark(1,:));
    title(axa{itype},[sessiontypes{itype}]);
    clabel='DA';            
    ylabel(axa{itype},clabel)
    aa=plot(axa{itype},([-daci; daci]+daavg)','--','linewidth',1,'color',mark(1,:));
    aa(1).Color(4)=.2;      %1st line
    aa(2).Color(4)=.2;      %2nd line
    am=plot(axa{itype},daavg','-','linewidth',1.5,'color',mark(1,:));
    am.Color(4)=.6;
    set(axa{itype},'ylim',[limsmind limsmaxd])
    %plot lfp right
    yyaxis(axa{itype},'right');
    set(axa{itype},'ycolor',markr(1,:));
    clabel_lfp=['beta-LFP ' num2str(freqband(1)) ...
                 '-' num2str(freqband(2)) ' hz'];               
    hy=ylabel(axa{itype},clabel_lfp);
    %rotate right yaxis
    set(hy,'rotation',270,'units','pixels');
    hypos=get(hy,'position');
    set(hy,'position',[hypos(1)+25,hypos(2),hypos(3)]);  
    aa=plot(axa{itype},([-lfpci; lfpci]+lfpavg)','--','linewidth',1,'color',markr(1,:));
    aa(1).Color(4)=.2;      %1st line
    aa(2).Color(4)=.2;      %2nd line
    am=plot(axa{itype},lfpavg','-','linewidth',1.5,'color',markr(1,:));
    am.Color(4)=.6;
    xticklabels=min(win):interval:max(win);
    xticklabels=round(xticklabels.*rate)./rate;
    xticklabels=num2str(xticklabels');
    xticks=1:round(interval*rate):size(da,2);
    set(axa{itype},'xtick',xticks,'xticklabel',xticklabels);        
    xlabel(axa{itype},'time (s)')
    set(axa{itype},'xlim',[0 size(da,2)]);     
    set(axa{itype},'ylim',[limsminl limsmaxl])
    set(findall(axa{itype},'-property','FontSize'),'FontSize',fontsize)
    end
end   

titletext=[datype ' | ' ...
        eventtypes{ievent} ' | da ' targdasites{ida} ...
        ' | lfp beta ' lfpsites{ilfp}];  
text(axa{1},0,axpos{1}(4)+10,titletext,'units','pixels',...
    'fontsize',fontsize+2,'fontweight','bold');

savedir=[savepath 'allsess_tracesavg'  filesep ...
    targdasites{ida} 'x' lfpsites{ilfp} filesep];

if ~isdir(savedir)
    mkdir(savedir);
end

savefig(figsess,[savedir  datype '_' eventtypes{ievent} ...
    '_' targdasites{ida} 'x' lfpsites{ilfp}]);
saveas(figsess,[savedir datype '_' eventtypes{ievent} ...
    '_' targdasites{ida} 'x' lfpsites{ilfp}],'tif')
delete(findall(figsess,'type','text')) 


end
end
end


end