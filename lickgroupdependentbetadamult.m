savepath=fullfile(plotparam.savepath,'dep_corr_');
event='targeye';
targbeh='reyed';
condtype={'right','phase1'};
condtype={'left'};

%condtype={'all'};
%condtype={'phase1'};
discretecutoff=0;
cutoff{1}=260;
cutoff{2}=450;          %changed from 400
cutoffbot{1}=230;       %changed from 240
cutoffbot{2}=400;       %changed form 300
grptypes=0;
binavg=1;
pcthigh=75;
pctlow=20;
mididx=300;
rate=10;
norm=1;
normz=0;
logbeta=0;
beh='eye';
%beh='pulse';
beh='trt';
pulsetype='rrstd';
pad=[0 0];
botgroup=1;
offset=40;
pvalcutoff=0.5;
noncorr=0;
rpos=0;
rtarg=-.1;
targimwinonly=0;
%beh='lick';
%beh='error';
figpos=[50,50,1400,600];
axwid=400;
axhei=350;
if binavg
    figpos(4)=900;
end
figsess=figure('visible','off');     %figure for each channel
if ispc
figsess=figure('visible','on');     %figure for each channel
end
set(figsess,'position',figpos,'color',[1 1 1]);
set(0,'CurrentFigure',figsess);    %set figure handle to current figure
axa={};
axpos={};
ax1=gca;
cmap = get(gca,'ColorOrder');
ctop=cmap(1,:);
cbot=cmap(2,:);
ctops=brewermap(11,'PiYG');
cbots=brewermap(9,'OrRd');
cbots=ctops(end-3:end,:);
ctops=ctops(1:5,:);

das={};
betas={};
count=1;
behs{1}=[];
behs{2}=[];
for ise=1:length(sessnums)
sessnum=sessnums(ise);
%binwidth=10;
sesid=find(strcmp({trialgrps.sessid},num2str(sessnum)));
trialinfo=trialgrps(sesid).trialinfo;
condtrialsid=find(contains(trialinfo(1).names,condtype));
condtrials{1}=trialinfo(1).nums{condtrialsid};
condtrials{2}=trialinfo(2).nums{condtrialsid};
condlabel=condtype{1};
if length(condtrialsid)>1
    for ii=1:length(condtrialsid)
    condtrials{1}=intersect(condtrials{1},trialinfo(1).nums{condtrialsid(ii)});
condtrials{2}=intersect(condtrials{2},trialinfo(2).nums{condtrialsid(ii)});
if ii==1
    condlabel=condtype{ii};
else
condlabel=[condlabel '_' condtype{ii}];
end
    end
end

targdasites=plotparam.dasites;
sites=getsites(sessnum,targdasites,'patra');
uniquesites=unique({sites(1:end).site});
cnsites=uniquesites(contains(uniquesites,'c'));
psites=uniquesites(contains(uniquesites,'p'));
sitesessids=find([sites.sessnum]==sessnum);
dachs={sites.probeid};
dasites={sites.site};
daids=[sites.ch];
[dapair,lfppair]=getsitepairs(dachs,'patra');
big=[];
sma=[];
btrials={};
bids={};
bdata={};
if ~strcmp(beh,'error') && ~strcmp(beh,'trt')
    big=find(strcmp({xbinfos.sessionid},num2str(sessnum)) & ...
        strcmp({xbinfos.event},event) & ...
        contains({xbinfos.sessiontype},'big') & ...
        strcmp({xbinfos.sitelfp},beh));
    big=big(1);
    sma=find(contains({xbinfos.sessionid},num2str(sessnum)) & ...
        strcmp({xbinfos.event},event) & ...
        contains({xbinfos.sessiontype},'small') & ...
        strcmp({xbinfos.sitelfp},beh));
    sma=sma(1);

    btrials{1}=intersect(xbinfos(big).daall.trials, condtrials{1});
    btrials{2}=intersect(xbinfos(sma).daall.trials, condtrials{2});
    bids{1}=find(ismember(xbinfos(big).daall.trials,btrials{1}));
    bids{2}=find(ismember(xbinfos(sma).daall.trials,btrials{2}));
    switch beh
        case 'pulse'
            switch pulsetype
            case 'hr'
               bdata{1}=nanmean(xbinfos(big).daall.lfptracesaln(bids{1},mididx+pad(1):mididx+pad(2)+offset),2);
                bdata{2}=nanmean(xbinfos(sma).daall.lfptracesaln(bids{2},mididx+pad(1):mididx+pad(2)+offset),2);
            case 'rr'
                bdatatemp1=xbinfos(big).daall.lfptracesaln(bids{1},mididx+pad(1):mididx+pad(2)+offset);
                bdatatemp2=xbinfos(sma).daall.lfptracesaln(bids{2},mididx+pad(1):mididx+pad(2)+offset);
                bdata1=1./bdatatemp1.*60.*1e3; 
                bdata2=1./bdatatemp2.*60.*1e3; 
                bdata{1}=nanmean(bdata1,2);
                bdata{2}=nanmean(bdata2,2);
            case 'rrstd'
                bdatatemp1=xbinfos(big).daall.lfptracesaln(bids{1},mididx+pad(1):mididx+pad(2)+offset);
                bdatatemp2=xbinfos(sma).daall.lfptracesaln(bids{2},mididx+pad(1):mididx+pad(2)+offset);
                bdata1=1./bdatatemp1.*60.*1e3; 
                bdata2=1./bdatatemp2.*60.*1e3; 
                bdata{1}=nanstd(bdata1,[],2);
                bdata{2}=nanstd(bdata2,[],2);
            case 'rmssd'
                bdatatemp1=xbinfos(big).daall.lfptracesaln(bids{1},mididx+pad(1):mididx+pad(2)+offset);
                bdatatemp2=xbinfos(sma).daall.lfptracesaln(bids{2},mididx+pad(1):mididx+pad(2)+offset);
                bdata1=1./bdatatemp1.*60.*1e3; 
                bdata2=1./bdatatemp2.*60.*1e3; 
                datawin=rms(diff(rrtemp(:,wins),1,2),2);  
                 bdata{1}=rms(diff(bdata1,1,2),2);
                 bdata{2}=rms(diff(bdata2,1,2),2);
            end
        case 'eye'
            %zscore eye wrt baseline as done in olivier et al 2014 frontier
            if offset>0
                bdatatemp1=xbinfos(big).daall.lfptracesaln(bids{1},:).*100000;
                bdatatemp2=xbinfos(sma).daall.lfptracesaln(bids{2},:).*100000;
                
                bdata{1}=nanmean(xbinfos(big).daall.lfptracesaln(bids{1},mididx:mididx+offset),2).*100000;
                bdata{2}=nanmean(xbinfos(sma).daall.lfptracesaln(bids{2},mididx:mididx+offset),2).*100000;
            else
                bdata{1}=nanmean(xbinfos(big).daall.lfptracesaln(bids{1},mididx+offset:mididx),2).*100000;
                bdata{2}=nanmean(xbinfos(sma).daall.lfptracesaln(bids{2},mididx+offset:mididx),2).*100000;
            end
        otherwise
            if offset>0
                bdata{1}=nanmean(xbinfos(big).daall.lfptracesaln(bids{1},mididx:mididx+offset),2).*100000;
                bdata{2}=nanmean(xbinfos(sma).daall.lfptracesaln(bids{2},mididx:mididx+offset),2).*100000;
            else
                bdata{1}=nanmean(xbinfos(big).daall.lfptracesaln(bids{1},mididx+offset:mididx),2).*100000;
                bdata{2}=nanmean(xbinfos(sma).daall.lfptracesaln(bids{2},mididx+offset:mididx),2).*100000;
            end
    end
    %bdata{1}=nanmean(xbinfos(big).daall.lfptracesaln(:,mididx-offset:mididx),2);
    %bdata{2}=nanmean(xbinfos(sma).daall.lfptracesaln(:,mididx-offset:mididx),2);
elseif strcmp(beh,'trt')
    big=find(strcmp({binfos.sessionid},num2str(sessnum)) & ...
        strcmp({binfos.event},event) & ...
        contains({binfos.sessiontype},'big'));
    big=big(1);
    sma=find(contains({binfos.sessionid},num2str(sessnum)) & ...
        strcmp({binfos.event},event) & ...
        contains({binfos.sessiontype},'small'));
    sma=sma(1);

    btrials{1}=intersect(binfos(big).seltrials, condtrials{1});
    btrials{2}=intersect(binfos(sma).seltrials, condtrials{2});
    bids{1}=find(ismember(binfos(big).seltrials,btrials{1}));
    bids{2}=find(ismember(binfos(sma).seltrials,btrials{2}));

    bdata{1}=binfos(big).target_rts(bids{1}).*1000;
    bdata{2}=binfos(sma).target_rts(bids{2}).*1000;
end

if length(sessnums)<2
    figure; histogram(bdata{1})
hold on; histogram(bdata{2})
%sum(binfos(big).lickpost<370)/length(~isnan(binfos(big).lickpost))
ax1=gca;
cmap = get(gca,'ColorOrder');
ctop=cmap(1,:);
cbot=cmap(2,:);

aplot=plot(ax1,prctile(bdata{1},pcthigh)*[1 1],get(ax1,'YLim'),'b');
aplot=plot(ax1,prctile(bdata{1},pctlow)*[1 1],get(ax1,'YLim'),'b');
aplot=plot(ax1,prctile(bdata{2},pcthigh)*[1 1],get(ax1,'YLim'),'r');
aplot=plot(ax1,prctile(bdata{2},pctlow)*[1 1],get(ax1,'YLim'),'r');
xlabel(ax1,beh);
ylabel(ax1,'counts');
 title(ax1,['big / small  histogram, session ' num2str(sessnum) ]);
savename=[savepath beh '_hist_sess' num2str(sessnum)];
savefig(gcf,savename);
saveas(gcf,savename,'jpg')
close(gcf);
end
%set(axa,'position',[axpos(1) axpos(2) axsiz(3) axsiz(4)])
%axsiz=[350 300];
lfpmetric='targimwin';
for ida=1:length(dasites)
    datarg=dachs{ida};
    datargids=find(strcmp(dapair,datarg));
    datargid=datargids(1);
    %clf(figsess);
    tid=find(strcmp(dapair{datargid},dachs));
    targid=daids(tid);
da{2}=datm{sesid}{2}{targid}.targpeak;
da{1}=datm{sesid}{1}{targid}.targpeak;
 lfptarg=lfppair{datargid};
lfpdet=getlfpsites(sessnum,{lfptarg});
if isempty(lfpdet) && length(datargids)>1
    flag=0;
    for inn=2:length(datargids)
        datargid=datargids(inn);
        lfptargtemp=lfppair{datargid};
        lfpdettemp=getlfpsites(sessnum,{lfptargtemp});
        if ~isempty(lfpdettemp) && flag==0
            lfpdet=lfpdettemp;
            lfptarg=lfptargtemp;
            flag=1;
        end
    end
end
lfpsite=lfpdet.site;
if contains(lfptarg,'p')
        lfpmetric='targwin';
else
    lfpmetric='targimwin';
end
if targimwinonly
    lfpmetric='targimwin';
end
 targidlfp=[];
for iich=1:length(betatm{sesid}{1})
if strcmp(betatm{sesid}{1}{iich}.site,lfptarg)
    targidlfp=iich;
end
end  
if isempty(targidlfp)
    continue
end
beta{1}=getfield(betatm{sesid}{1}{targidlfp},lfpmetric);
beta{2}=getfield(betatm{sesid}{2}{targidlfp},lfpmetric);
if logbeta
    beta{1}=log(beta{1});
    beta{2}=log(beta{2});
end
outlierthres=nanmean(horzcat(beta{:}))+7*nanstd(horzcat(beta{:}));
%{
outliers=find(beta{1}>outlierthres);
if ~isempty(outliers)
beta{1}(outliers)=[];
aid=find(ismember(btrials{1},betatm{sesid}{1}{targidlfp}.trialnums(outliers)));
btrials{1}(outliers)=[];
end
outliers=find(beta{2}>outlierthres);
if ~isempty(outliers)
beta{2}(outliers)=[];
aid=find(ismember(btrials{2},betatm{sesid}{2}{targidlfp}.trialnums(outliers)));
btrials{2}(outliers)=[];
end
%}

if norm
    danormval=abs(nanmedian(horzcat(da{:}))); 
    betanormval=abs(nanmedian(horzcat(beta{:})));        %median of all values for all plot conditions for current probe
    for icond=1:length(da)        
        da{icond}=da{icond}./danormval;
        beta{icond}=beta{icond}./betanormval;
    end
end


%{
for ix=1:2
    axa{ix}=subplot(1,2,ix);
    set(axa{ix},'units','pixels');
    axpos{ix}=get(axa{ix},'position');
    hold(axa{ix},'on');
    ylabel(axa{ix},'\beta LFP, \muV^2','interpreter','tex');
    xlabel(axa{ix},'\Delta[DA], nM','interpreter','tex');
end
%}
for itype=1:2
    seltrialstop{itype}=[];
    seltrialsbot{itype}=[];
    if ~strcmp(beh,'error')
    seltrialstop{itype}=btrials{itype}(bdata{itype}>prctile(bdata{itype},pcthigh));
    seltrialsbot{itype}=btrials{itype}(bdata{itype}<prctile(bdata{itype},pctlow));   
    if discretecutoff
            seltrialstop{itype}=btrials{itype}(bdata{itype}>=cutoff{itype});
    seltrialsbot{itype}=btrials{itype}(bdata{itype}<=cutoffbot{itype});   
    end
    else
        condtrialsid=find(contains(trialinfo(1).names,'fail'));
        condtrialstop{itype}=trialinfo(itype).nums{condtrialsid};
        condtrialsid=find(contains(trialinfo(1).names,'success'));
        condtrialsbot{itype}=trialinfo(itype).nums{condtrialsid};
    seltrialstop{itype}=condtrialstop{itype};     %after error
    seltrialsbot{itype}=condtrialsbot{itype};   %after success
    end        
    seltrialstop{itype}=intersect(seltrialstop{itype},datm{sesid}{itype}{targid}.trialnums);
    seltrialstop{itype}=intersect(seltrialstop{itype},betatm{sesid}{itype}{targidlfp}.trialnums);
    seltrialsbot{itype}=intersect(seltrialsbot{itype},datm{sesid}{itype}{targid}.trialnums);
    seltrialsbot{itype}=intersect(seltrialsbot{itype},betatm{sesid}{itype}{targidlfp}.trialnums);    
   
    datop{itype}=da{itype}(find(ismember(datm{sesid}{itype}{targid}.trialnums,seltrialstop{itype})));
    betatop{itype}=beta{itype}(find(ismember(betatm{sesid}{itype}{targidlfp}.trialnums,seltrialstop{itype})));
    betabot{itype}=beta{itype}(find(ismember(betatm{sesid}{itype}{targidlfp}.trialnums,seltrialsbot{itype})));
    dabot{itype}=da{itype}(find(ismember(datm{sesid}{itype}{targid}.trialnums,seltrialsbot{itype})));
    

    topgood=find(~isnan(datop{itype}) & ~isnan(betatop{itype}));
    botgood=find(~isnan(dabot{itype}) & ~isnan(betabot{itype}));
    datop{itype}=datop{itype}(topgood);
    betatop{itype}=betatop{itype}(topgood);
    dabot{itype}=dabot{itype}(botgood);
    betabot{itype}=betabot{itype}(botgood);
        %{
    hold(axa{itype},'on'); 
    if itype==1
    title(axa{itype},['big ' datarg ' ' lfptarg ]);
    else
            title(axa{itype},['small ' datarg ' ' lfptarg ]);
    end
    scatter(axa{itype},datop{itype},betatop{itype},50,'+','MarkerEdgeAlpha',1,'MarkerFaceAlpha',.5,'MarkerFaceColor',ctop,'markeredgecolor',ctop); 
    scatter(axa{itype},dabot{itype},betabot{itype},50,'sq','MarkerEdgeAlpha',1,'MarkerFaceAlpha',.5,'MarkerFaceColor',cbot,'markeredgecolor',cbot); 
    %find(betabigtop>500)
    %betabigtop(64)=[];
    %dabigtop(64)=[];
    [rtop,ptop]=corr(datop{itype}',betatop{itype}');
    [rbot,pbot]=corr(dabot{itype}',betabot{itype}');   
    [polydata,sdata]=polyfit(datop{itype},betatop{itype},1);
    slope=polydata(1);
    intercep=polydata(2);
    yfit=slope*datop{itype}+intercep;
    yresid=betatop{itype}-yfit;
    ssresid=sum(yresid.^2);
    sstotal=(length(betatop{itype})-1)*var(betatop{itype});
    rsq=1-ssresid/sstotal;
    fitline=plot(axa{itype},datop{itype},yfit,'--','color',ctop,'linewidth',1);
    fitline.Color(4)=0.5;
    text(axa{itype},axpos{itype}(3)-100,axpos{itype}(4)-50,...
    {['top ' num2str(pcthigh) ' percentile'],['r: ' num2str(rtop)],...
    ['p: ' num2str(ptop)]},'color',ctop,'units','pixels');
    [polydata,sdata]=polyfit(dabot{itype},betabot{itype},1);
    slope=polydata(1);
    intercep=polydata(2);
    yfit=slope*dabot{itype}+intercep;
    yresid=betabot{itype}-yfit;
    ssresid=sum(yresid.^2);
    sstotal=(length(betabot{itype})-1)*var(betabot{itype});
    rsq=1-ssresid/sstotal;
    fitline=plot(axa{itype},dabot{itype},yfit,'linestyle','--','color',cbot,'linewidth',1);
    fitline.Color(4)=0.5;
    text(axa{itype},axpos{itype}(3)-100,axpos{itype}(4)-100,...
    {['bot ' num2str(pctlow) ' percentile'],['r: ' num2str(rbot)],...
    ['p: ' num2str(pbot)]},'color',cbot,'units','pixels');
%}

das(count).top{itype}=datop{itype};
das(count).bot{itype}=dabot{itype};
betas(count).top{itype}=betatop{itype};
betas(count).bot{itype}=betabot{itype};
das(count).site=dasites{ida};
das(count).sess=sessnum;
betas(count).site=lfpsite;
%behs{1}=[behs{1} bdata{1}];
%behs{2}=[behs{2} bdata{2}];
if grptypes
    das(count).top{1}=[das(count).top{1} datop{itype}];
das(count).bot{1}=[das(count).bot{1} dabot{itype}];
betas(count).top{1}=[betas(count).top{1} betatop{itype}];
betas(count).bot{1}=[betas(count).bot{1} betabot{itype}];
end
end

count=count+1;
end
end

datop={};
dabot={};
betatop={};
betabot={};
corrsites={};
countcor={};
countcor{1}=1;
countcor{2}=1;
        binavgbetabot={};
        binavgbetatop={};
        binstdbetatop={};
        binstdbetabot={};
        bincibetabot={};
        bincibetatop={};
         binavgdabot={};
        binavgdatop={};
        binstddatop={};
        binstddabot={};
        bincidatop={};
        bincidabot={};
                binedges=[0:.1:1];

for ix=1:2
    axa{ix}=subplot(1,2,ix);
    set(axa{ix},'units','pixels');
    axpos{ix}=get(axa{ix},'position');
    hold(axa{ix},'on');
    ylabel(axa{ix},'\beta LFP, \muV^2','interpreter','tex');
    xlabel(axa{ix},'\Delta[DA], nM','interpreter','tex');
end
if binavg
    for ix=1:4
    axa{ix}=subplot(2,2,ix);
    set(axa{ix},'units','pixels');
    axpos{ix}=get(axa{ix},'position');
    hold(axa{ix},'on');
    ylabel(axa{ix},'\beta LFP norm','interpreter','tex');
    xlabel(axa{ix},'\Delta[DA] norm','interpreter','tex');
    if ix>2
     xlabel(axa{ix},'\beta LFP norm','interpreter','tex');
    ylabel(axa{ix},'\Delta[DA] norm','interpreter','tex');       
    end
    set(axa{ix},'position',[axpos{ix}(1) axpos{ix}(2) axwid axhei]);
    set(axa{ix},'fontsize',14)
    end
end
for itype=1:2
    datop{itype}=[];
    betatop{itype}=[];
    dabot{itype}=[];
    betabot{itype}=[];
    if length(das)>4
    for ix=1:length(das)
        if ~isempty(das(ix).top{itype}) && ~isempty(das(ix).bot{itype})
            [rtop,ptop]=corr(das(ix).top{itype}',betas(ix).top{itype}');
        [rbot,pbot]=corr(das(ix).bot{itype}',betas(ix).bot{itype}');  
        chor=rtop;
        chop=ptop;
        datatemp=das(ix).top{itype};
        if botgroup
            chor=rbot;
            chop=pbot;
            datatemp=das(ix).bot{itype};
        end
        meetsitecriteria=0;
        if ~noncorr
            if ~rpos
            if chor<rtarg && chop<pvalcutoff
                meetsitecriteria=1;
            end
            else
            if chor>-rtarg && chop<pvalcutoff
                meetsitecriteria=1;
            end
            end
        elseif noncorr
            if chor<0.2 && chor>-0.2 && chop>0.1
                meetsitecriteria=1;
            end
        end
        if meetsitecriteria && length(datatemp)>5
            %if site shows neg correlation
        %if contains(das(ix).site,'c54') || contains(das(ix).site,'p2') || contains(das(ix).site,'p33')
        %if contains(betas(ix).site,'c')
        datoptemp=das(ix).top{itype};
        dabottemp=das(ix).bot{itype};
        betatoptemp=betas(ix).top{itype};
        betabottemp=betas(ix).bot{itype};
        
        if normz
         %median of all values for all plot conditions for current probe
                datoptemp=zscore(datoptemp);
                dabottemp=zscore(dabottemp);
                 betatoptemp=zscore(betatoptemp);
                betabottemp=zscore(betabottemp);               
        end       
        if rescale
                datoptemp=(datoptemp-min(datoptemp))./(max(datoptemp)-min(datoptemp));
                                dabottemp=(dabottemp-min(dabottemp))./(max(dabottemp)-min(dabottemp));
                betatoptemp=(betatoptemp-min(betatoptemp))./(max(betatoptemp)-min(betatoptemp));
                                betabottemp=(betabottemp-min(betabottemp))./(max(betabottemp)-min(betabottemp));
        end
        
        datop{itype}=[datop{itype} datoptemp];
        dabot{itype}=[dabot{itype} dabottemp];
        betatop{itype}=[betatop{itype} betatoptemp];
        betabot{itype}=[betabot{itype} betabottemp];  
        corrsites{itype}(countcor{itype}).siteda=das(ix).site;
        corrsites{itype}(countcor{itype}).sitelfp=betas(ix).site;
        corrsites{itype}(countcor{itype}).sess=das(ix).sess;
        corrsites{itype}(countcor{itype}).datop=das(ix).top{itype};
        corrsites{itype}(countcor{itype}).dabot=das(ix).bot{itype};
        corrsites{itype}(countcor{itype}).betatop=betas(ix).top{itype};
        corrsites{itype}(countcor{itype}).betabot=betas(ix).bot{itype};
        corrsites{itype}(countcor{itype}).rtop=rtop;
                corrsites{itype}(countcor{itype}).ptop=ptop;
        corrsites{itype}(countcor{itype}).rbot=rbot;
                corrsites{itype}(countcor{itype}).pbot=pbot;
        countcor{itype}=countcor{itype}+1;
        end
        end
       % end
    end
    cla(axa{itype});
    %betatop{itype}=log(betatop{itype});
    %betabot{itype}=log(betabot{itype});
    if ~binavg
    scatter(axa{itype},datop{itype},betatop{itype},50,'+','MarkerEdgeAlpha',1,'MarkerFaceAlpha',.5,'MarkerFaceColor',ctop,'markeredgecolor',ctop); 
    scatter(axa{itype},dabot{itype},betabot{itype},10,'o','MarkerEdgeAlpha',1,'MarkerFaceAlpha',.5,'MarkerFaceColor',cbot,'markeredgecolor',cbot); 
    else
        [daN,daEDGES,bindabot] = histcounts(dabot{itype},binedges);
        [daN,daEDGES,bindatop] = histcounts(datop{itype},binedges);
        for ib=1:length(binedges)-1
           idsbot=find(bindabot==ib);
              idstop=find(bindatop==ib);
            binavgbetabot{itype}(ib)=nanmean(betabot{itype}(idsbot));
            binstdbetabot{itype}(ib)=nanstd(betabot{itype}(idsbot));
             bincibetabot{itype}(ib)=nanstd(betabot{itype}(idsbot))./sqrt(length(idsbot)).*1.96;
            binavgbetatop{itype}(ib)=nanmean(betatop{itype}(idstop));
            binstdbetatop{itype}(ib)=nanstd(betatop{itype}(idstop));
           bincibetatop{itype}(ib)=nanstd(betatop{itype}(idstop))./sqrt(length(idstop)).*1.96;
        end
        [daN,daEDGES,binbetabot] = histcounts(betabot{itype},binedges);
        [daN,daEDGES,binbetatop] = histcounts(betatop{itype},binedges);
        for ib=1:length(binedges)-1
            idsbot=find(binbetabot==ib);
            idstop=find(binbetatop==ib);
            binavgdabot{itype}(ib)=nanmean(dabot{itype}(idsbot));
            binstddabot{itype}(ib)=nanstd(dabot{itype}(idsbot));
            bincidabot{itype}(ib)=nanstd(dabot{itype}(idsbot))./sqrt(length(idsbot)).*1.96;
            binavgdatop{itype}(ib)=nanmean(datop{itype}(idstop));
            binstddatop{itype}(ib)=nanstd(datop{itype}(idstop));
            bincidatop{itype}(ib)=nanstd(datop{itype}(idstop))./sqrt(length(idstop)).*1.96;
        end
         scatter(axa{itype},binedges(1:end-1)+.05,binavgbetatop{itype},50,'+','MarkerEdgeAlpha',1,'MarkerFaceAlpha',.5,'MarkerFaceColor',ctop,'markeredgecolor',ctop); 
         scatter(axa{itype},binedges(1:end-1)+.05,binavgbetabot{itype},10,'o','MarkerEdgeAlpha',1,'MarkerFaceAlpha',.5,'MarkerFaceColor',cbot,'markeredgecolor',cbot); 
         errorbar(axa{itype},binedges(1:end-1)+.05,binavgbetabot{itype},bincibetabot{itype}./1.96,'marker','none','linestyle','none','linewidth',1,'capsize',1,'color',cbot);
                  errorbar(axa{itype},binedges(1:end-1)+.05,binavgbetatop{itype},bincibetatop{itype}./1.96,'marker','none','linestyle','none','linewidth',1,'capsize',1,'color',ctop);
                   cla(axa{itype+2});
                  scatter(axa{itype+2},binedges(1:end-1)+.05,binavgdatop{itype},50,'+','MarkerEdgeAlpha',1,'MarkerFaceAlpha',.5,'MarkerFaceColor',ctop,'markeredgecolor',ctop); 
         scatter(axa{itype+2},binedges(1:end-1)+.05,binavgdabot{itype},10,'o','MarkerEdgeAlpha',1,'MarkerFaceAlpha',.5,'MarkerFaceColor',cbot,'markeredgecolor',cbot); 
         errorbar(axa{itype+2},binedges(1:end-1)+.05,binavgdabot{itype},bincidabot{itype}./1.96,'marker','none','linestyle','none','linewidth',1,'capsize',1,'color',cbot);
                  errorbar(axa{itype+2},binedges(1:end-1)+.05,binavgdatop{itype},bincidatop{itype}./1.96,'marker','none','linestyle','none','linewidth',1,'capsize',1,'color',ctop);
    end
    %find(betabigtop>500)
    %betabigtop(64)=[];
    %dabigtop(64)=[];
    [rtop,ptop]=corr(datop{itype}',betatop{itype}');
    [rbot,pbot]=corr(dabot{itype}',betabot{itype}');   
    [polydata,sdata]=polyfit(datop{itype},betatop{itype},1);
    slope=polydata(1);
    intercep=polydata(2);
    yfit=slope*datop{itype}+intercep;
    yresid=betatop{itype}-yfit;
    ssresid=sum(yresid.^2);
    sstotal=(length(betatop{itype})-1)*var(betatop{itype});
    rsq=1-ssresid/sstotal;
    fitline=plot(axa{itype},datop{itype},yfit,'--','color',ctop,'linewidth',1);
    fitline.Color(4)=0.5;
    text(axa{itype},axpos{itype}(3)-100,axpos{itype}(4)-0,...
    {['top ' num2str(pcthigh) ' percentile'],['r: ' num2str(rtop)],...
    ['p: ' num2str(ptop)]},'color',ctop,'units','pixels');
    [polydata,sdata]=polyfit(dabot{itype},betabot{itype},1);
    slope=polydata(1);
    intercep=polydata(2);
    yfit=slope*dabot{itype}+intercep;
    yresid=betabot{itype}-yfit;
    ssresid=sum(yresid.^2);
    sstotal=(length(betabot{itype})-1)*var(betabot{itype});
    rsq=1-ssresid/sstotal;
    fitline=plot(axa{itype},dabot{itype},yfit,'linestyle','--','color',cbot,'linewidth',1);
    fitline.Color(4)=0.5;
    text(axa{itype},axpos{itype}(3)-100,axpos{itype}(4)-60,...
    {['bot ' num2str(pctlow) ' percentile'],['r: ' num2str(rbot)],...
    ['p: ' num2str(pbot)]},'color',cbot,'units','pixels');
        if binavg
            %make slopes for da vs bet plot on bottom
            [rtop,ptop]=corr(betatop{itype}',datop{itype}');
            [rbot,pbot]=corr(betabot{itype}',dabot{itype}');   
            [polydata,sdata]=polyfit(betatop{itype},datop{itype},1);
            slope=polydata(1);
            intercep=polydata(2);
            yfit=slope*betatop{itype}+intercep;
            yresid=datop{itype}-yfit;
            ssresid=sum(yresid.^2);
            sstotal=(length(datop{itype})-1)*var(datop{itype});
            rsq=1-ssresid/sstotal;
            fitline=plot(axa{itype+2},betatop{itype},yfit,'--','color',ctop,'linewidth',1);
            fitline.Color(4)=0.5;           
            [polydata,sdata]=polyfit(betabot{itype},dabot{itype},1);
            slope=polydata(1);
            intercep=polydata(2);
            yfit=slope*betabot{itype}+intercep;
            yresid=dabot{itype}-yfit;
            ssresid=sum(yresid.^2);
            sstotal=(length(dabot{itype})-1)*var(dabot{itype});
            rsq=1-ssresid/sstotal;
            fitline=plot(axa{itype+2},betabot{itype},yfit,'linestyle','--','color',cbot,'linewidth',1);
            fitline.Color(4)=0.5;          
        end
    else
        %single session plot
        cla(axa{itype});
        for ix=1:length(das)
    datop{itype}=das(ix).top{itype};
        dabot{itype}=das(ix).bot{itype};
        betatop{itype}=betas(ix).top{itype};
        betabot{itype}=betas(ix).bot{itype};
        dasite=das(ix).site;
                betasite=betas(ix).site;

    scatter(axa{itype},datop{itype},betatop{itype},10,'+','MarkerEdgeAlpha',1,'MarkerFaceAlpha',.4,'MarkerFaceColor',ctops(ix,:),'markeredgecolor',ctops(ix,:)); 
    scatter(axa{itype},dabot{itype},betabot{itype},50,'sq','MarkerEdgeAlpha',1,'MarkerFaceAlpha',.4,'MarkerFaceColor',cbots(ix,:),'markeredgecolor',cbots(ix,:)); 
    
    [rtop,ptop]=corr(datop{itype}',betatop{itype}');
    [rbot,pbot]=corr(dabot{itype}',betabot{itype}');   
    [polydata,sdata]=polyfit(datop{itype},betatop{itype},1);
    slope=polydata(1);
    intercep=polydata(2);
    yfit=slope*datop{itype}+intercep;
    yresid=betatop{itype}-yfit;
    ssresid=sum(yresid.^2);
    sstotal=(length(betatop{itype})-1)*var(betatop{itype});
    rsq=1-ssresid/sstotal;
    fitline=plot(axa{itype},datop{itype},yfit,'--','color',ctops(ix,:),'linewidth',1);
    fitline.Color(4)=0.5;
    text(axa{itype},axpos{itype}(3)-100,axpos{itype}(4)-50-(ix-1)*100,...
    {[dasite ' | ' betasite ' top ' num2str(pcthigh) ' percentile'],['r: ' num2str(rtop)],...
    ['p: ' num2str(ptop)]},'color',ctops(ix,:),'units','pixels');
    [polydata,sdata]=polyfit(dabot{itype},betabot{itype},1);
    slope=polydata(1);
    intercep=polydata(2);
    yfit=slope*dabot{itype}+intercep;
    yresid=betabot{itype}-yfit;
    ssresid=sum(yresid.^2);
    sstotal=(length(betabot{itype})-1)*var(betabot{itype});
    rsq=1-ssresid/sstotal;
    fitline=plot(axa{itype},dabot{itype},yfit,'linestyle','--','color',cbots(ix,:),'linewidth',1);
    fitline.Color(4)=0.5;
    text(axa{itype},axpos{itype}(3)-100,axpos{itype}(4)-100-(ix-1)*100,...
    {['bot ' num2str(pctlow) ' percentile'],['r: ' num2str(rbot)],...
    ['p: ' num2str(pbot)]},'color',cbots(ix,:),'units','pixels');
        end
        
    end
end
%{
figure; histogram(behs{1},[200:5:600])
hold on; histogram(behs{2},[200:5:600])
%}
savename=savepath;
if rpos
    savename=[savename 'poscorrgrp_'];
end
if targimwinonly
    savename=[savename 'targimwinall_'];
end
if noncorr
    savename=[savename 'noncorrgrp_'];
end
if binavg
    savename=[savename 'binavg_'];
end
if discretecutoff
    savename=[savename 'cutoff_' num2str(cutoff{1}) 'top_' num2str(cutoffbot{1}) 'bot_' num2str(cutoff{2}) 'smtop_' num2str(cutoffbot{2}) 'smbot_'];
end
if grptypes
    savename=[savename 'grptypes_'];
end
if norm
    savename=[savename 'norm_'];
end
if rescale
    savename=[savename 'rescale_'];
end
if length(sessnums)>1
    savename=[savename 'multsess_'];    
else
    savename=[savename 'sess' num2str(sessnum) '_'];
end    
if botgroup
    savename=[savename 'botgroup_'];
end
savename=[savename  beh '_' event  '_win_' num2str(round(offset/rate)) 's_' condlabel];
save(savename,'corrsites','das','betas');
savefig(figsess,savename);
saveas(figsess,savename,'jpg')
print(figsess,savename,'-painters','-depsc');

