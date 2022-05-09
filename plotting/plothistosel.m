function plothistosel(sessnum,xinfo,dasite,plotparam,trialnums,varargin)
%plot histogram of supplied and indicated signals
%for trials supplied/ sess id or trials before/after supplied trials based on 
%trlists supplied
figpos=[50,50,950,500];
rate=10;      %10hz default for downsampled da/lfp from xcov
win=[-1 4];     %+/-2 s from aln idx
interval=1;       %in seconds
ratelfp=1000;
argnum=1;
datype='goodtrials';
sesstypes={};
trtypes={'all'};
argnum=1;
fontsize=14;
conditions='all';
targsig='trt';      %target reaction time
event='targ';
tx=trialnums(1).cat;
trtypes=unique({trialnums(1).trialgrps.type});
if length(trialnums)>1
    for iix=2:length({trialnums.cat})
        tx=[tx '_vs_' trialnums(iix).cat];
    end
end
label=['trls_' tx '_' trialnums(1).site '_' [trtypes{:}]];

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
        case 'label'
            %group label, type
            argnum=argnum+1;
            labtx=varargin{argnum};
            labtemp=[labtx{:}];
            for ix=1:length(labtemp)
                labtt=labtemp{ix};
                if ~ischar(labtt)
                    labtt=num2str(labtt);
                end
                if ix<length(labtemp)
                    label=[label labtt '_'];
                else
                    label=[label labtt];
                end
                
            end
    end
    argnum=argnum+1;
end
fontsize=14;
savepath=plotparam.savepath;
targsig='datracesaln';
if strcmp(sigtype,'lfp')
targsig='lfptracesaln';
end
legtext={};
if length(trialnums)>1
    %2 groups to plot
    for iix=1:length(trialnums)
        legtext{iix}=trialnums(iix).cat;
    end
end
mark=[1 0 0; 0 0 0; 0 .7 0;.7 0 0; 0.4 0 .7];
if length(trialnums)>2
mark=linspecer(5,'qualitative');
end
trialinfo=plotparam.trialgrps(contains({plotparam.trialgrps.sessid},num2str(sessnum))).trialinfo;
infonames={'bigreward','smallreward','targetbreak','fixbreak'};
tsx=[win(1):1/rate:win(2)];
tsxlfp=[win(1):1/ratelfp:win(2)];
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
linem={'-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-'};
hold(axa,'on')
%for itype=1:length(sessiontypes)
for igrp=1:length(trialnums)
sigtemp=[];
siterows={xinfo.siteda};
if strcmp(sigtype,'lfp')
    siterows={xinfo.sitelfp};
end
for itype=1:length(trialnums(igrp).trialgrps)
    sesstype=trialnums(igrp).trialgrps(itype).type;
targrows=find(strcmp(siterows,dasite) & ...
    strcmp({xinfo.event},event) & ...
    contains({xinfo.sessionid},num2str(sessnum)) &...
    contains({xinfo.sessiontype},sesstype));
if ~isempty(targda)
   targrows=find(strcmp(siterows,dasite) & ...
       strcmp({xinfo.siteda},targda) & ...
    strcmp({xinfo.event},event) & ...
    contains({xinfo.sessionid},num2str(sessnum)) &...
    contains({xinfo.sessiontype},sesstype));
end
targrow=targrows(1);
trtype=find(contains(infonames,sesstype));
trialtypes=trialinfo(trtype); %get trial groups for session   
ttype=find(contains(trialtypes.names,conditions)==1);
typename=trialtypes.names{ttype};
curdata=getfield(xinfo(targrow),'daall');
datatrials=getfield(xinfo(targrow),datype);     %get trial ids for pos/neg/all
%get trials that match condition (trialtypes), datype (datatrials), and
%supplied trials (trialnums)
seltrials=find(ismember(curdata.trials,trialtypes.nums{ttype}) & ismember(curdata.trials,datatrials) & ismember(curdata.trials,trialnums(igrp).trialgrps(itype).trials));
targdata=getfield(curdata,targsig);
alnidx=curdata.mididx;
wins=[alnidx+win(1)*rate:alnidx+win(2)*rate];       %time window to plot
datawin=targdata(seltrials,wins);
sigtemp=[sigtemp; datawin];
end

da=sigtemp;
if plotz
    da=zscore(da,[],2);
end
tsplot=tsx;
daavg=nanmean(da,1);    
dastd=nanstd(da,[],1);    
daci=dastd./sqrt(size(da,1))*1.96;
dase=dastd./sqrt(size(da,1));
if plotci 
    plotshaded(axa,tsplot,([-daci; daci]+daavg)',mark(igrp,:));
%aa=plot(axa,tsplot,([-daci; daci]+daavg)','linestyle',linem{itype},'linewidth',1,'color',mark(itype,:));
%aa(1).Color(4)=.2;      %1st line
%aa(2).Color(4)=.2;      %2nd line
end
if plotse
    plotshaded(axa,tsplot,([-dase; dase]+daavg)',mark(igrp,:));
end    
am=plot(axa,tsplot,daavg','linestyle',linem{igrp},'linewidth',1,'color',mark(igrp,:));
am.Color(4)=.6;
if ~isempty(legtext)
text(axa,axsiz(1)+10,axsiz(2)+50-25*igrp,legtext{igrp},'units','pixels',...
    'fontsize',fontsize,'color',mark(igrp,:),'interpreter','none');
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

titletext=['sess' num2str(sessnum) ' | ' conditions ' | '  dasite  ];  
text(axa,-50,axsiz(2)+100,titletext,'units','pixels','fontsize',fontsize+1,'fontweight','bold','interpreter','none');
text(axa,-50,axsiz(2)+75,label,'units','pixels','fontsize',fontsize-1,'interpreter','none');

savename=[savepath 'tracesavgsel_' 'sess' num2str(sessnum) '_' dasite '_' conditions '_' label];
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
