function outdata=plottracesel(sessnum,xinfo,site,plotparam,trialnums,varargin)
outdata=[];
%OVERLAY BIG VS SAMLL and defined groups
figpos=[50,50,950,500];
rate=10;      %10hz default for downsampled da/lfp from xcov
win=[-1 4];     %+/-2 s from aln idx
interval=1;       %in seconds
ratelfp=1000;
datype='goodtrials';
argnum=1;
fontsize=14;
conditions='all';
event='targ';
sigtype='da';
plotz=0;
ploth=0;
htype='max';        %plot maxima/minima/mean for histograms in defined win
plott=0;
plotci=0;
plotse=0;
binfo={};
bplot=0;
targda='';
bstatic=0;
binwidth=25;              %ms
binmax=600;
binedges=0:binwidth:binmax;
label='';
targtypes={};
trseq={};
trtable={};
argnum=1;
plotline=0;
alphashade=.4;
condtypes={};
trplot=0;
pulsetype='';
ttrtypes={};
plotibi=0;
defaultcolor=0;
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
        case 'binfo'
            %supply binfo data
            argnum=argnum+1;
            binfo=varargin{argnum};
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
        case 'ibi'
            %interburst interval lfp
            sigtype='lfp';
            ftype='sitelfp';
            plotibi=1;
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
        case 'beh'
            %plot b variable, must supply xbinfo instead of xinfo
            bplot=1;
            sigtype='lfp';
        case 'targda'
            %for b or lfp plot, pos/neg da, wan t to supply reference da channel
            argnum=argnum+1;
            targda=varargin{argnum};
        case 'hist'
            ploth=1;
        case 'trtable'
            trplot=1;           %plot according to desired trial sequence, must supply binfo as well
            argnum=argnum+1;
            trtable=varargin{argnum};
            argnum=argnum+1;
            targtypes=varargin{argnum};     %{'big','small'}
            argnum=argnum+1;
            trseq=varargin{argnum};         %trsequence, as would be mapped onto maketrorg2 variable, eg. {{'-3fix','-2fix','-1fix','-1big','-1small'}}
        case 'plotline'
            plotline=1;
        case 'condtypes'
            %plot conditions against each other
            argnum=argnum+1;
            targtypes=varargin{argnum};     %{'big','small'}
            argnum=argnum+1;
            condtypes=varargin{argnum};   
        case 'trtypes'
            %plot trial types ('big', vs 'small') with each group in identified conditions against each other
            argnum=argnum+1;
            targtypes=varargin{argnum};     %{'big','small'}
            argnum=argnum+1;
            ttrtypes=varargin{argnum};   
        case 'rmssd'
            site='pulse';
            bplot=1;
            pulsetype='rmssd';
            binedges=0:1:20;
        case 'rr'
            site='pulse';
            bplot=1;
            pulsetype='rr';
            binedges=350:10:600;
        case 'rrstd'
            site='pulse';
            bplot=1;
            pulsetype='rrstd';
            binedges=0:1:27;
        case 'hr'
            site='pulse';
            bplot=1;
            pulsetype='hr';
            binedges=90:2.5:160;
        case 'dcol'
            %default color
            defaultcolor=1;
    end
    argnum=argnum+1;
end
fontsize=14;
savepath=plotparam.savepath;

%get trial nums for selected da ch if da ch, otherwise get all trialnums
%for all sites in trialtypes
dasite=[];
sitedetail=getsites(sessnum,{site});
if ~isempty(sitedetail)
    dasite=sitedetail.site;
end
%configure trialtypes according to site specified, if da channel, more selective? NO, because already selective in datm and xinfo variables...    
    %NEVERMIND
    
if trplot
    trialnums=getmulttrials(plotparam,xinfo,binfo,sessnum,site,targtypes,event,'trtable',trtable,trseq);
    tx=trialnums(1).cat;
    trtypes=unique({trialnums(1).trialgrps.type});
    if length(trialnums)>1
        for iix=2:length({trialnums.cat})
            tx=[tx '_vs_' trialnums(iix).cat];
        end
    end
    label=['trls_' tx '_' trialnums(1).site '_' [trtypes{:}]];
    outdata=trialnums;
end

trialinfo=plotparam.trialgrps(contains({plotparam.trialgrps.sessid},num2str(sessnum))).trialinfo;
infonames={'bigreward','smallreward','targetbreak','fixbreak'};
if ~isempty(condtypes)
    trialnums=getmulttrials(plotparam,xinfo,binfo,sessnum,site,targtypes,event,'conditions',condtypes);
    tx=trialnums(1).cat;
    trtypes=unique({trialnums(1).trialgrps.type});
    if length(trialnums)>1
        for iix=2:length({trialnums.cat})
            tx=[tx '_vs_' trialnums(iix).cat];
        end
    end
    label=['trls_' tx '_' trialnums(1).site '_' [trtypes{:}]];
    outdata=trialnums;
end
if ~isempty(ttrtypes)
    trialnums=getmulttrials(plotparam,xinfo,binfo,sessnum,site,targtypes,event,'trialtypeconds',ttrtypes);
    tx=trialnums(1).cat;
    ttrtypes=unique({trialnums(1).trialgrps.type});
    if length(trialnums)>1
        for iix=2:length({trialnums.cat})
            tx=[tx '_vs_' trialnums(iix).cat];
        end
    end
    label=['trls_' tx '_' trialnums(1).site '_' [ttrtypes{:}]];
    outdata=trialnums;
end
%get signal type label to retrieve from supplied xinfos, xbinfos, binfos
targsig='datracesaln';  %da traces aligned to specified event
if strcmp(sigtype,'lfp')    
targsig='lfptracesaln'; %lfp or behav traces aligned to specified event
end
if strcmp(site,'trt')
    %reaction times, not time data
    targsig='target_rts';
    bstatic=1;
end
if strcmp(site,'rt')
    targsig='fix_rt';
    bstatic=1;
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
if length(trialnums)>5
mark=linspecer(8,'qualitative');
end
end


tsx=[win(1):1/rate:win(2)];
figsess=figure('visible','off');     %figure for each channel
if ispc
figsess=figure('visible','on');     %figure for each channel
end
if length(trialnums)==2 || defaultcolor
cmap = get(gca,'ColorOrder');
mark=cmap;
end
clf(figsess);
set(figsess,'position',figpos,'color',[1 1 1]);
set(0,'CurrentFigure',figsess);    %set figure handle to current figure
axa=axes;
set(axa,'units','pixels');
axpos=get(axa,'position');
axsiz=[400 300];
off=75;
set(axa,'units','pixels','position',[off off axsiz(1) axsiz(2)])
linem={'-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-'};
hold(axa,'on')

for igrp=1:length(trialnums)
sigtemp=[];
siterows=[];
if ~bstatic
    siterows={xinfo.siteda};
    if strcmp(sigtype,'lfp')
        siterows={xinfo.sitelfp};
    end
end
targrows=[];
datatrials=[];
curdata=[];
alnidx=[];
for itype=1:length(trialnums(igrp).trialgrps)
    sesstype=trialnums(igrp).trialgrps(itype).type;
    trtype=find(contains(infonames,sesstype));
trialtypes=trialinfo(trtype); %get trial groups for session   
ttype=find(contains(trialtypes.names,conditions)==1);
typename=trialtypes.names{ttype};

if ~bstatic
    %time trace series data
    basedata=[];        %baseline pre-fix for eye diameter
targrows=find(strcmp(siterows,site) & ...
    strcmp({xinfo.event},event) & ...
    contains({xinfo.sessionid},num2str(sessnum)) &...
    contains({xinfo.sessiontype},sesstype));
if ~isempty(targda)
   targrows=find(strcmp(siterows,site) & ...
       strcmp({xinfo.siteda},targda) & ...
    strcmp({xinfo.event},event) & ...
    contains({xinfo.sessionid},num2str(sessnum)) &...
    contains({xinfo.sessiontype},sesstype));
end
targrow=targrows(1);
curdata=getfield(xinfo(targrow),'daall');
datatrials=getfield(xinfo(targrow),datype);     %get trial ids for pos/neg/all
%get trials that match condition (trialtypes), datype (datatrials), and
%supplied trials (trialnums)
seltrials=find(ismember(curdata.trials,trialtypes.nums{ttype}) & ismember(curdata.trials,datatrials) & ismember(curdata.trials,trialnums(igrp).trialgrps(itype).trials));
seltrialnums=intersect(curdata.trials,trialtypes.nums{ttype});
seltrialnums=intersect(seltrialnums,trialnums(igrp).trialgrps(itype).trials);
targdata=getfield(curdata,targsig);
alnidx=curdata.mididx;
wins=[alnidx+win(1)*rate:alnidx+win(2)*rate];       %time window to plot
basewin=[alnidx-3*rate:alnidx];          %fix baseline window for eye below
%plothwin=[alnidx-7*rate:alnidx+4*rate];         %for calculating zscore values, plotwin should be extended for histograms that will take mean over window
glitchwidth=2;
if isempty(pulsetype)
    datawin=targdata(seltrials,wins);
    basedata=[];
    if plotibi
    %interburst interval calculation for lfp
        basedata=targdata(:,basewin);              %baseline from all good trials;
        baseline=nanmean(nanstd(basedata,[],2));
        dataz=(datawin-nanmean(datawin,2))./baseline;       %common baseline uniform for all trials, mean trial selective
        %dataz=(datawin-nanmean(nanmean(datawin,2)))./baseline;       %common baseline uniform for all trials, mean 
        %dataz=(datawin-nanmean(datawin,2))./nanstd(datawin,[],2);   
        thres=nanstd(dataz,[],2)*.5;
        rasterdata=zeros(size(dataz));
        rasterhist=nan(size(dataz));            %timestamps of bursts detected each trial
        ibis=nan(size(dataz));
        for ix=1:size(dataz,1)
            idsbeta=find(dataz(ix,:)>=thres(ix));
            breakpts=find(abs(diff(idsbeta))>1);
            idspeakburst=[];           %center max peak of bursts
            if ~isempty(breakpts)
                for ib=1:length(breakpts)+1
                    maxb=nan;
                    if ib==1
                        [maxamp,maxid]=max(dataz(ix,idsbeta(1:breakpts(ib))));
                        maxb=idsbeta(maxid);
                    elseif ib>1 && ib<=length(breakpts)
                       [maxamp,maxid]=max(dataz(ix,idsbeta(breakpts(ib-1)+1:breakpts(ib))));
                       maxb=idsbeta(maxid + breakpts(ib-1));
                    elseif ib>length(breakpts)
                       [maxamp,maxid]=max(dataz(ix,idsbeta(breakpts(ib-1)+1:end)));
                       maxb=idsbeta(maxid + breakpts(ib-1));
                    end
                    idspeakburst(ib)=maxb;
                end
                rasterdata(ix,idspeakburst)=1;
                rasterhist(ix,idspeakburst)=idspeakburst./rate;
                if length(idspeakburst)>1
                shiftedids=circshift(idspeakburst,-1);       %shift to left;
                ibitemp=(shiftedids(1:end-1)-idspeakburst(1:end-1))./rate;
                ibis(ix,idspeakburst(1:end-1))=ibitemp;
                end
            end
        end
        %get means ibi for window starting from targeted align idx w 200 ms
        %padding
        ibimean=nanmean(ibis(:,-win(1)*rate-1:size(ibis,2)),2);
        datawin=ibimean;
    end
    if strcmp(site,'eye')
        basewin=[alnidx-7*rate:alnidx+4*rate];          %fix baseline window for eye below
        %zscore eye wrt baseline as done in olivier et al 2014 frontier
        %remove blinks
        datawin=-datawin;       %flip so original sign
        datawin=deglitchinterpsimple(datawin,3.5,glitchwidth);
        targrows=find(strcmp(siterows,site) & ...
            strcmp({xinfo.event},'fix') & ...
            contains({xinfo.sessionid},num2str(sessnum)) &...
            contains({xinfo.sessiontype},sesstype));
        if ~isempty(targda)
           targrows=find(strcmp(siterows,site) & ...
               strcmp({xinfo.siteda},targda) & ...
            strcmp({xinfo.event},'fix') & ...
            contains({xinfo.sessionid},num2str(sessnum)) &...
            contains({xinfo.sessiontype},sesstype));
        end
       % targrow=targrows(1);
       % fixdata=getfield(xinfo(targrow),'daall');
       % datatrials=getfield(xinfo(targrow),datype);     %get trial ids for pos/neg/all
       % seltrials=find(ismember(fixdata.trials,seltrialnums));
        %fixdata=getfield(fixdata,targsig);
      %  basedata=fixdata(seltrials,basewin);
            basedata=targdata(seltrials,basewin);   %normalizing to pre-fix and targ window allows normalization trial to trial to look at faster changes within 1 sec 06/04 without phase trend eg. 67/68/69

        basedata=-basedata;
      %   basedata=deglitchinterpsimple(basedata,2,2);
        for ix=1:size(basedata,1)
            [xx, xb]=rmoutliers(basedata(ix,:));
            basedata(ix,xb)=nan;
        end
       % baseline=nanmean(nanstd(basedata(:,1:5*rate),[],2));
       baseline=nanstd(basedata,[],2);          %normalized to each trial so independent of phase changes
               dataz=(datawin-nanmean(basedata,2))./baseline;

        %dataz=(datawin-nanmean(datawin,2))./baseline;
        datawin=dataz;
    end    
    
    sigtemp=[sigtemp; datawin];
    
elseif ~isempty(pulsetype)
    % pulse data
    datawin=targdata(seltrials,:);
    switch pulsetype
        case 'hr'
            datawin=datawin(:,wins);
            datawin=nanmean(datawin,2);
        case 'rr'
           rrtemp=1./datawin.*60.*1e3; 
           datawin=rrtemp(:,wins);
           datawin=nanmean(datawin,2);
        case 'rrstd'
            rrtemp=1./datawin.*60.*1e3; 
            datawin=rrtemp(:,wins);
            datawin=nanstd(datawin,[],2);
        case 'rmssd'
            rrtemp=1./datawin.*60.*1e3; 
            datawin=rms(diff(rrtemp(:,wins),1,2),2);    
    end
        sigtemp=[sigtemp; datawin];
end
if ploth && size(sigtemp,1)>1
    sigtemp=nanmean(sigtemp,2);
end
else
    %reaction time static data
targrows=find(strcmp({xinfo.event},event) & ...
contains({xinfo.sessionid},num2str(sessnum)) &...
contains({xinfo.sessiontype},sesstype));
targrow=targrows(1);
seltrials=find(ismember(xinfo(targrow).seltrials,trialtypes.nums{ttype}) & ismember(xinfo(targrow).seltrials,trialnums(igrp).trialgrps(itype).trials));
targdata=getfield(xinfo(targrow),targsig);
datawin=targdata(seltrials).*1e3;       %in ms
sigtemp=[sigtemp datawin];
end
end

da=sigtemp;
if strcmp(site,'eye')
 %   da=-da;
end
if ~ploth && ~plotibi
    %time-series trace plot
if plotz
    da=zscore(da,[],2);
end
tsplot=tsx;
%if size da trials < 3, all nan
if size(da,1)<3
    da(1:size(da,1),:)=nan;
end
%remove columns (time points) where > 50% of data is nan at that ts
da2=isnan(da);
da3=sum(da2,1);
rempts=find(da3>.5*size(da,1));
%get smoothed trace before removing bad ts's for se
dasmooth=da;
for ix=1:size(da,1)
[dasmooth(ix,:),TF] = fillmissing(da(ix,:),'pchip');
end
da(:,rempts)=nan;
daavg=nanmean(da,1);   
daavgforse=daavg;
dastd=nanstd(da,[],1);    
daci=dastd./sqrt(size(da,1))*1.96;
dase=dastd./sqrt(size(da,1));
if any(isnan(dase))
    dastd=nanstd(dasmooth,[],1);    
daci=dastd./sqrt(size(dasmooth,1))*1.96;
dase=dastd./sqrt(size(dasmooth,1));
daavgforse=nanmean(dasmooth,1);   
end
if plotci 
    plotshaded(axa,tsplot,([-daci; daci]+daavgforse)',mark(igrp,:),alphashade);
end
if plotse
    plotshaded(axa,tsplot,([-dase; dase]+daavgforse)',mark(igrp,:),alphashade);
end    
if (~plotse && ~plotci) || plotline
am=plot(axa,tsplot,daavg','linestyle',linem{igrp},'linewidth',1,'color',mark(igrp,:));
am.Color(4)=.6;
end
else
    %histo plot
    if strcmp(site,'rt') || strcmp(site,'trt')     
        %reaction times
        [bindata,bb]=histc(da,binedges);
        histplot= bar(axa,binedges,bindata,'FaceColor',mark(igrp,:),'BarWidth',1,'FaceAlpha',.4,'linestyle','none');  
                cidas=nanstd(da)./sqrt(length(da))*1.96;
        avgline=plot(axa,[nanmean(da)-cidas nanmean(da)-cidas],get(axa,'ylim'),'linestyle','--','linewidth',1,'color',mark(igrp,:));
        avgline=plot(axa,[nanmean(da)+cidas nanmean(da)+cidas],get(axa,'ylim'),'linestyle','--','linewidth',1,'color',mark(igrp,:));
    elseif strcmp(site,'pulse')
        %reaction times
        [bindata,bb]=histc(da,binedges);
        histplot= bar(axa,binedges,bindata,'FaceColor',mark(igrp,:),'BarWidth',1,'FaceAlpha',.4,'linestyle','none');  
        cidas=nanstd(da)./sqrt(length(da))*1.96;
        avgline=plot(axa,[nanmean(da)-cidas nanmean(da)-cidas],get(axa,'ylim'),'linestyle','--','linewidth',1,'color',mark(igrp,:));
        avgline=plot(axa,[nanmean(da)+cidas nanmean(da)+cidas],get(axa,'ylim'),'linestyle','--','linewidth',1,'color',mark(igrp,:));
    else
        %da/lfp data
        da=rmoutliers(da);
        binlims=[nanmean(da)-nanstd(da)*3 nanmean(da)+nanstd(da)*3];
        if plotibi
            binlims(1)=min(da);
            binlims(2)=max(da);
        end
        binwidth=round(ceil(binlims(2)-binlims(1))./20*100)/100;
        binedges=binlims(1):binwidth:binlims(2);
        if length(binedges)<=1
            binwidth=round(ceil((binlims(2)-binlims(1))*1000)./20*100)/100/1000;
            binedges=binlims(1):binwidth:binlims(2);
        end
        [bindata,bb]=histcounts(da,binedges,'Normalization','probability');
        %[bindata,bb]=histcounts(da,binedges);
        %[bindata,bb]=histc(da,binedges);
        histplot= bar(axa,binedges(1:end-1),bindata,'FaceColor',mark(igrp,:),'BarWidth',1,'FaceAlpha',.4,'linestyle','none');
         cidas=nanstd(da)./sqrt(length(da))*1.96;
        avgline=plot(axa,[nanmean(da)-cidas nanmean(da)-cidas],get(axa,'ylim'),'linestyle','--','linewidth',1,'color',mark(igrp,:));
        avgline=plot(axa,[nanmean(da)+cidas nanmean(da)+cidas],get(axa,'ylim'),'linestyle','--','linewidth',1,'color',mark(igrp,:));
    end
end

%legend for grp
if ~isempty(legtext)
text(axa,axsiz(1)+10,axsiz(2)+50-25*igrp,legtext{igrp},'units','pixels',...
    'fontsize',fontsize,'color',mark(igrp,:),'interpreter','none');
end

end

%finish plotting labels & save
if ~ploth && ~plotibi
    %time series trace plot
clabel='\Delta[DA] (nM)';      
if strcmp(sigtype,'lfp')
    clabel='\beta-lfp (\muV^2)';
end
if bplot
    clabel=site;
end
ylabel(axa,clabel)
xticklabels=min(win):interval:max(win);
xticklabels=round(xticklabels.*rate)./rate;
xticklabels=num2str(xticklabels');
xticks=1:round(interval*rate):size(da,2);
xlabel(axa,'time (s)')
xlim(axa,[min(win) max(win)]);
else
    %histo plot
    if strcmp(site,'rt') || strcmp(site,'trt')     
        %reaction times
        xlabel(axa,'reaction time (ms)');
        set(axa,'xtick',binedges(1:10:end));
        set(axa,'xlim',[min(binedges) max(binedges)],'box','off');
        ylabel(axa, '# trials')   
    end
    if strcmp(site,'pulse')
        %reaction times
        xlabel(axa,pulsetype);
        set(axa,'xtick',binedges(1:10:end));
        set(axa,'xlim',[min(binedges) max(binedges)],'box','off');
        ylabel(axa, '# trials')   
    end
    if plotibi
        clabel='\Delta[DA] (nM)';      
if strcmp(sigtype,'lfp')
    clabel='\beta-lfp (\muV^2)';
end
if bplot
    clabel=site;
end
        xlabel(axa,'ibi width (s)')
ylabel(axa,clabel)
    end
end

%finish plot formatting
set(axa,'box','off')
set(findall(axa,'-property','FontSize'),'FontSize',fontsize)
titletext=['sess' num2str(sessnum) ' | ' conditions ' | '  site ' | ' event ];  
text(axa,-50,axsiz(2)+100,titletext,'units','pixels','fontsize',fontsize+1,'fontweight','bold','interpreter','none');
if ~isempty(label)
text(axa,-50,axsiz(2)+75,label,'units','pixels','fontsize',fontsize-1,'interpreter','none');
savename=[savepath 'tracesavgsel_' 'sess' num2str(sessnum) '_' site '_' event '_' conditions '_' label];
else
    savename=[savepath 'tracesavgsel_' 'sess' num2str(sessnum) '_' site '_' event '_' conditions];
end
if plotz
    savename=[savename '_z'];
end
if plotibi
    savename=[savename '_ibi'];
end
if plotse
    savename=[savename '_se'];
end
if plotci
    savename=[savename '_ci'];
end
if ploth
    savename=[savename '_hist'];
end
if ~isempty(pulsetype)
    savename=[savename '_' pulsetype];
end
if ~isdir(savepath)
    mkdir(savepath);
end
savefig(figsess,savename);
saveas(figsess,savename,'tif')
print(figsess,savename,'-painters','-depsc');
end
