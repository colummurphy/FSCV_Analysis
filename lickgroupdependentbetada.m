savepath=fullfile(plotparam.savepath,'lickpostdep_corr_');
sessnum=92;
event='targ';
targbeh='lickpost';
pcthigh=75;
pctlow=25;
%binwidth=10;
sesid=find(strcmp({trialgrps.sessid},num2str(sessnum)));
trialinfo=trialgrps(sesid).trialinfo;
big=find(strcmp({binfos.sessionid},num2str(sessnum)) & ...
    strcmp({binfos.event},event) & ...
    contains({binfos.sessiontype},'big'));
sma=find(contains({binfos.sessionid},num2str(sessnum)) & ...
    strcmp({binfos.event},event) & ...
    contains({binfos.sessiontype},'small'));
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

cmap=colormap('parula');
[cmap,num,typ]=brewermap(2,'PuOr');
%ctop=cmap(round(size(cmap,1)*.2),:);
%cbot=cmap(round(size(cmap,1)*.80),:);
minlim=min([binfos(big).lickpost binfos(sma).lickpost]);
maxlim=max([binfos(big).lickpost binfos(sma).lickpost]);
binwid=floor((maxlim-minlim)/50);
figure; histogram(binfos(big).lickpost,[minlim:binwid:maxlim])
hold on; histogram(binfos(sma).lickpost,[minlim:binwid:maxlim])
%sum(binfos(big).lickpost<370)/length(~isnan(binfos(big).lickpost))
ax1=gca;
cmap = get(gca,'ColorOrder');
ctop=cmap(1,:);
cbot=cmap(2,:);

aplot=plot(ax1,prctile(binfos(big).lickpost,75)*[1 1],get(ax1,'YLim'),'b');
aplot=plot(ax1,prctile(binfos(big).lickpost,50)*[1 1],get(ax1,'YLim'),'b');
aplot=plot(ax1,prctile(binfos(sma).lickpost,75)*[1 1],get(ax1,'YLim'),'r');
aplot=plot(ax1,prctile(binfos(sma).lickpost,50)*[1 1],get(ax1,'YLim'),'r');
xlabel(ax1,'lick');
ylabel(ax1,'counts');
 title(ax1,['big / small lick histogram, session ' num2str(sessnum) ]);
savepath=fullfile(plotparam.savepath,'lickdep_hist_');
savename=[savepath 'sess' num2str(sessnum)];
savefig(gcf,savename);
saveas(gcf,savename,'jpg')

figpos=[50,50,1400,600];
figsess=figure('visible','off');     %figure for each channel
if ispc
figsess=figure('visible','on');     %figure for each channel
end
set(figsess,'position',figpos,'color',[1 1 1]);
set(0,'CurrentFigure',figsess);    %set figure handle to current figure
axa={};
axpos={};

%set(axa,'position',[axpos(1) axpos(2) axsiz(3) axsiz(4)])
%axsiz=[350 300];
lfpmetric='targimwin';
for ip=1:length(dapair)
    clf(figsess);
    tid=find(strcmp(dapair{ip},dachs));
    targid=daids(tid);
da{2}=datm{sesid}{2}{targid}.targpeak;
da{1}=datm{sesid}{1}{targid}.targpeak;
if contains(lfppair{ip},'p')
        lfpmetric='targwin';
    else
        lfpmetric='targimwin';
    end
     targidlfp=[];
    for iich=1:length(betatm{sesid}{1})
    if strcmp(betatm{sesid}{1}{iich}.site,lfppair{ip})
        targidlfp=iich;
    end
    end  
    if isempty(targidlfp)
        continue
    end
    beta{1}=getfield(betatm{sesid}{1}{targidlfp},lfpmetric);
    beta{2}=getfield(betatm{sesid}{2}{targidlfp},lfpmetric);
    seltrialsbig=binfos(big).seltrials(binfos(big).lickpost>prctile(binfos(big).lickpost,pcthigh));
    seltrialsbigbot=binfos(big).seltrials(binfos(big).lickpost<prctile(binfos(big).lickpost,pctlow));
    seltrialstop{2}=binfos(sma).seltrials(binfos(sma).lickpost>prctile(binfos(sma).lickpost,pcthigh));
    seltrialsbot{2}=binfos(sma).seltrials(binfos(sma).lickpost<prctile(binfos(sma).lickpost,pctlow));
    bigtrials=intersect(seltrialsbig,datm{sesid}{1}{targid}.trialnums);
    bigtrials=intersect(bigtrials,betatm{sesid}{1}{targidlfp}.trialnums);
    bigtrialsbot=intersect(seltrialsbigbot,datm{sesid}{1}{targid}.trialnums);
    bigtrialsbot=intersect(bigtrialsbot,betatm{sesid}{1}{targidlfp}.trialnums);    
   
    dabigtop=da{1}(find(ismember(datm{sesid}{1}{targid}.trialnums,bigtrials)));
    betabigtop=beta{1}(find(ismember(betatm{sesid}{1}{targidlfp}.trialnums,bigtrials)));
    betabigbot=beta{1}(find(ismember(betatm{sesid}{1}{targidlfp}.trialnums,bigtrialsbot)));
    dabigbot=da{1}(find(ismember(datm{sesid}{1}{targid}.trialnums,bigtrialsbot)));
    for ix=1:2
    axa{ix}=subplot(1,2,ix);
    set(axa{ix},'units','pixels');
    axpos{ix}=get(axa{ix},'position');
    hold(axa{ix},'on');
    ylabel(axa{ix},'\beta LFP, \muV^2','interpreter','tex');
    xlabel(axa{ix},'\Delta[DA], nM','interpreter','tex');
    end
    hold(axa{1},'on'); 
    title(axa{1},['big ' dapair{ip} ' ' lfppair{ip} ]);
    topgood=find(~isnan(dabigtop) & ~isnan(betabigtop));
    botgood=find(~isnan(dabigbot) & ~isnan(betabigbot));
    dabigtop=dabigtop(topgood);
    betabigtop=betabigtop(topgood);
    dabigbot=dabigbot(botgood);
    betabigbot=betabigbot(botgood);
    scatter(axa{1},dabigtop,betabigtop,50,'+','MarkerEdgeAlpha',1,'MarkerFaceAlpha',.5,'MarkerFaceColor',ctop,'markeredgecolor',ctop); 
    scatter(axa{1},dabigbot,betabigbot,50,'sq','MarkerEdgeAlpha',1,'MarkerFaceAlpha',.5,'MarkerFaceColor',cbot,'markeredgecolor',cbot); 
    %find(betabigtop>500)
    %betabigtop(64)=[];
    %dabigtop(64)=[];
    [rtop,ptop]=corr(dabigtop',betabigtop');
    [rbot,pbot]=corr(dabigbot',betabigbot');   
    [polydata,sdata]=polyfit(dabigtop,betabigtop,1);
    slope=polydata(1);
    intercep=polydata(2);
    yfit=slope*dabigtop+intercep;
    yresid=betabigtop-yfit;
    ssresid=sum(yresid.^2);
    sstotal=(length(betabigtop)-1)*var(betabigtop);
    rsq=1-ssresid/sstotal;
    fitline=plot(axa{1},dabigtop,yfit,'--','color',ctop);
    fitline.Color(4)=0.5;
    text(axa{1},axpos{1}(3)-100,axpos{1}(4)-50,...
    {['top 25th percentile'],['r: ' num2str(rtop)],...
    ['p: ' num2str(ptop)]},'color',ctop,'units','pixels');
    [polydata,sdata]=polyfit(dabigbot,betabigbot,1);
    slope=polydata(1);
    intercep=polydata(2);
    yfit=slope*dabigbot+intercep;
    yresid=betabigbot-yfit;
    ssresid=sum(yresid.^2);
    sstotal=(length(betabigbot)-1)*var(betabigbot);
    rsq=1-ssresid/sstotal;
    fitline=plot(axa{1},dabigbot,yfit,'--','color',cbot);
    fitline.Color(4)=0.5;
    text(axa{1},axpos{1}(3)-100,axpos{1}(4)-100,...
    {['bottom 50th percentile'],['r: ' num2str(rbot)],...
    ['p: ' num2str(pbot)]},'color',cbot,'units','pixels');

     hold(axa{2},'on');
         title(axa{2},['small ' dapair{ip} ' ' lfppair{ip} ]);
           seltrialstop{2}=intersect(seltrialstop{2},datm{sesid}{2}{targid}.trialnums);
    seltrialstop{2}=intersect(seltrialstop{2},betatm{sesid}{2}{targidlfp}.trialnums);
    seltrialsbot{2}=intersect(seltrialsbot{2},datm{sesid}{2}{targid}.trialnums);
    seltrialsbot{2}=intersect(seltrialsbot{2},betatm{sesid}{2}{targidlfp}.trialnums);    
    
     datop{2}=da{2}(find(ismember(datm{sesid}{2}{targid}.trialnums,seltrialstop{2})));
    betatop{2}=beta{2}(find(ismember(betatm{sesid}{2}{targidlfp}.trialnums,seltrialstop{2})));
    betabot{2}=beta{2}(find(ismember(betatm{sesid}{2}{targidlfp}.trialnums,seltrialsbot{2})));
    dabot{2}=da{2}(find(ismember(datm{sesid}{2}{targid}.trialnums,seltrialsbot{2})));
     topgood=find(~isnan(datop{2}) & ~isnan(betatop{2}));
    botgood=find(~isnan(dabot{2}) & ~isnan(betabot{2}));
    datop{2}=datop{2}(topgood);
    betatop{2}=betatop{2}(topgood);
    dabot{2}=dabot{2}(botgood);
    betabot{2}=betabot{2}(botgood);
    scatter(axa{2},datop{2},betatop{2},50,'+','MarkerEdgeAlpha',1,'MarkerFaceAlpha',.5,'MarkerFaceColor',ctop,'markeredgecolor',ctop); 
    scatter(axa{2},dabot{2},betabot{2},50,'sq','MarkerEdgeAlpha',1,'MarkerFaceAlpha',.5,'MarkerFaceColor',cbot,'markeredgecolor',cbot); 
    %find(betabigtop>500)
    %betabigtop(64)=[];
    %dabigtop(64)=[];
    [rtop,ptop]=corr(datop{2}',betatop{2}');
    [rbot,pbot]=corr(dabot{2}',betabot{2}');   
    [polydata,sdata]=polyfit(datop{2},betatop{2},1);
    slope=polydata(1);
    intercep=polydata(2);
    yfit=slope*datop{2}+intercep;
    yresid=betatop{2}-yfit;
    ssresid=sum(yresid.^2);
    sstotal=(length(betatop{2})-1)*var(betatop{2});
    rsq=1-ssresid/sstotal;
    fitline=plot(axa{2},datop{2},yfit,'--','color',ctop);
    fitline.Color(4)=0.5;
    text(axa{2},axpos{2}(3)-100,axpos{2}(4)-50,...
    {['top 25th percentile'],['r: ' num2str(rtop)],...
    ['p: ' num2str(ptop)]},'color',ctop,'units','pixels');
    [polydata,sdata]=polyfit(dabot{2},betabot{2},1);
    slope=polydata(1);
    intercep=polydata(2);
    yfit=slope*dabot{2}+intercep;
    yresid=betabot{2}-yfit;
    ssresid=sum(yresid.^2);
    sstotal=(length(betabot{2})-1)*var(betabot{2});
    rsq=1-ssresid/sstotal;
    fitline=plot(axa{2},dabot{2},yfit,'--','color',cbot);
    fitline.Color(4)=0.5;
    text(axa{2},axpos{2}(3)-100,axpos{2}(4)-100,...
    {['bottom 50th percentile'],['r: ' num2str(rbot)],...
    ['p: ' num2str(pbot)]},'color',cbot,'units','pixels');

savename=[savepath 'sess' num2str(sessnum) '_' event '_' dapair{ip} '_' lfppair{ip}];
savefig(figsess,savename);
saveas(figsess,savename,'jpg')
print(figsess,savename,'-painters','-depsc');
end
