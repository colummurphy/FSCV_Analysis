function plottracesummary(sessids,xinfo,dasite,plotparam,varargin)
%OVERLAY BIG VS SAMLL and defined groups
rate=10;      %10hz default for downsampled da/lfp from xcov
win=[-1 4];     %+/-2 s from aln idx
interval=1;       %in seconds
ratelfp=1000;
argnum=1;
datype='goodtrials';
evtype={};
sesstypes={};
trtypes={'all'};
argnum=1;
fontsize=14;
conditions='left';
event='targ';
tstrace={'damaxts'};
tstype='postrials';
sigtype='da';
ftype='siteda';
plotz=0;
plott=0;
plotci=0;
plotse=0;
binfo={};
bplot=0;
btarg='eye';
targda='';
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'daneg'
            datype='negtrials';
        case 'daall'
            datype='goodtrials';
        case 'dapos'
            datype='postrials';
        case 'event'
            %user selects what event aln
            argnum=argnum+1;
            event=varargin{argnum};
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
        case 'condition'
            argnum=argnum+1;
            conditions=varargin{argnum};
        case 'win'
            %user specified window +/- s
            argnum=argnum+1;
            win=varargin{argnum};  
        case 'lfp'
            %signal type lfp
            sigtype='lfp';
            ftype='sitelfp';
            tstrace={'lfppostmaxts','lfpmints'};
        case 'plotz'
            %zscore
            plotz=1;
        case 'traces'
            %plot average timing overlaid
            plott=1;
        case 'ci'
            %plot ci error bars
            plotci=1;
        case 'se'
            %plot se error bars
            plotse=1;
        case 'bplot'
            %plot b variable, must supply xbinfo instead of xinfo
            bplot=1;
            sigtype='lfp';
        case 'targda'
            %for b or lfp plot, pos/neg da, wan t to supply reference da channel
            argnum=argnum+1;
            targda=varargin{argnum};
    end
    argnum=argnum+1;
end
fontsize=14;
savepath=plotparam.savepath;

mark=[0 0 0; 1 0 0; 0 .7 0;.7 0 0; 0.4 0 .7];
trialinfo=plotparam.trialgrps(contains({plotparam.trialgrps.sessid},sessids)).trialinfo;
tsx=[win(1):1/rate:win(2)];
tsxlfp=[win(1):1/ratelfp:win(2)];
figpos=[50,50,600,500];
figsess=figure('visible','off');     %figure for each channel
if ispc
figsess=figure('visible','on');     %figure for each channel
end
set(figsess,'position',figpos,'color',[1 1 1]);
set(0,'CurrentFigure',figsess);    %set figure handle to current figure
axa=axes;
set(axa,'units','pixels');
axpos=get(axa,'position');
axsiz=[400 300];
off=75;
mar=100;
set(axa,'units','pixels','position',[off off axsiz(1) axsiz(2)])
sessiontypes={'big','small'};
linem={'-','--'};
linem={'-','-'};
hold(axa,'on')
for itype=1:length(sessiontypes)
rowlist={xinfo.siteda};
if strcmp(sigtype,'lfp')
    rowlist={xinfo.sitelfp};
end
targrows=find(strcmp(rowlist,dasite) & ...
    strcmp({xinfo.event},event) & ...
    contains({xinfo.sessionid},sessids) &...
    contains({xinfo.sessiontype},sessiontypes(itype)));
if ~isempty(targda)
   targrows=find(strcmp(rowlist,dasite) & ...
       strcmp({xinfo.siteda},targda) & ...
    strcmp({xinfo.event},event) & ...
    contains({xinfo.sessionid},sessids) &...
    contains({xinfo.sessiontype},sessiontypes(itype)));
end
targrow=targrows(1);
trialtypes=trialinfo(itype); %get trial groups for session   
ttype=find(contains(trialtypes.names,conditions)==1);
typename=trialtypes.names{ttype};
curdata=getfield(xinfo(targrow),'daall');
datatrials=getfield(xinfo(targrow),datype);     %get trial ids for pos/neg/all
    seltrials=find(ismember(curdata.trials,trialtypes.nums{ttype}) & ismember(curdata.trials,datatrials));

%curdata=getfield(xinfo(targrow),datype);
alnidx=curdata.mididx;
wins=[alnidx+win(1)*rate:alnidx+win(2)*rate];       %time window to plot
da=curdata.datracesaln(seltrials,wins);                             %signals aligned to aln idx
if strcmp(sigtype,'lfp')
    da=curdata.lfptracesaln(seltrials,wins);
end
if plotz
    da=zscore(da,[],2);
end
tsplot=tsx;
daavg=nanmean(da,1);    
dastd=nanstd(da,[],1);    
daci=dastd./sqrt(size(da,1))*1.96;
dase=dastd./sqrt(size(da,1));
if plotci 
    plotshaded(axa,tsplot,([-daci; daci]+daavg)',mark(itype,:),.2);
%aa=plot(axa,tsplot,([-daci; daci]+daavg)','linestyle',linem{itype},'linewidth',1,'color',mark(itype,:));
%aa(1).Color(4)=.2;      %1st line
%aa(2).Color(4)=.2;      %2nd line
end
if plotse
    plotshaded(axa,tsplot,([-dase; dase]+daavg)',mark(itype,:),.2);
end    
am=plot(axa,tsplot,daavg','linestyle',linem{itype},'linewidth',1,'color',mark(itype,:));
am.Color(4)=.6;
legtext{itype}=[sessiontypes{itype} ];
text(axa,axsiz(1)+10,axsiz(2)+50-25*itype,legtext{itype},'units','pixels',...
    'fontsize',fontsize,'color',mark(itype,:));

if ~isempty(tstrace) && plott
    %plot overlying ts parameter
    markc=[1 .75 1; 0 0 0 ];
    for it=1:length(tstrace)
    %dataforts=getfield(xinfo(targrow),tstype);
    datafortstemp=getfield(xinfo(targrow),'daall');
    tstrials=getfield(xinfo(targrow),tstype);     %get trial ids for pos/neg/all
    tstrialsids=find(ismember(datafortstemp.trials,tstrials));
    tsdataall=getfield(datafortstemp,tstrace{it})-alnidx;
    tsdata=tsdataall(tstrialsids);
    meants=nanmean(tsdata);
    stdts=nanstd(tsdata);
    tsci=stdts./sqrt(length(tsdata))*1.96;
    trange=round(meants-tsci):round(meants+tsci);
    trange=trange-win(1)*rate;
    tsmean=(meants)./rate;
    ttt=plot(axa,tsplot(trange),daavg(trange)','linestyle','-','linewidth',8,'color',markc(it,:));
    ttt.Color(4)=.3;
    end
end
end
clabel='\Delta[DA] (nM)';      
if strcmp(sigtype,'lfp')
    clabel='\beta-lfp (\muV^2)';
end
if bplot
    clabel=dasite;
end
ylabel(axa,clabel)
xticklabels=min(win):interval:max(win);
xticklabels=round(xticklabels.*rate)./rate;
xticklabels=num2str(xticklabels');
xticks=1:round(interval*rate):size(da,2);
xlabel(axa,'time (s)')
xlim(axa,[min(win) max(win)]);
set(findall(axa,'-property','FontSize'),'FontSize',fontsize)
set(axa,'box','off')

titletext=[conditions ' | ' datype ' | '  event ' | ' sigtype ' | ' dasite ' | sess ' sessids{:}];  
text(axa,-50,axsiz(2)+50,titletext,'units','pixels','fontsize',fontsize+1,'fontweight','bold');

savename=[savepath 'tracesavg_' datype '_' event '_' conditions '_' dasite '_sess' sessids{:}];
if plotz
    savename=[savename '_z'];
end
if plotse
    savename=[savename '_se'];
end
if plotci
    savename=[savename '_ci'];
end
if ~isdir(savepath)
    mkdir(savepath);
end
savefig(figsess,savename);
saveas(figsess,savename,'tif')
print(figsess,savename,'-painters','-depsc');

end
