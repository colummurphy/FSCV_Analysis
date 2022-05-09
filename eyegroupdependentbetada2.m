savepath=fullfile(plotparam.savepath,'eyerdep_corr_');
sessnum=83;
event='targeye';
targbeh='reyed';
condtype={'left','phase1'};
pcthigh=60;
pctlow=40;
mididx=300;
offset=40;
rate=10;
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

big=find(strcmp({xbinfos.sessionid},num2str(sessnum)) & ...
    strcmp({xbinfos.event},event) & ...
    contains({xbinfos.sessiontype},'big') & ...
    strcmp({xbinfos.sitelfp},'eye'));
big=big(1);
sma=find(contains({xbinfos.sessionid},num2str(sessnum)) & ...
    strcmp({xbinfos.event},event) & ...
    contains({xbinfos.sessiontype},'small') & ...
    strcmp({xbinfos.sitelfp},'eye'));
sma=sma(1);

btrials{1}=intersect(xbinfos(big).daall.trials, condtrials{1});
btrials{2}=intersect(xbinfos(sma).daall.trials, condtrials{2});
bids{1}=find(ismember(xbinfos(big).daall.trials,btrials{1}));
bids{2}=find(ismember(xbinfos(sma).daall.trials,btrials{2}));

bdata{1}=nanmean(xbinfos(big).daall.lfptracesaln(bids{1},mididx:mididx+offset),2).*100000;
bdata{2}=nanmean(xbinfos(sma).daall.lfptracesaln(bids{2},mididx:mididx+offset),2).*100000;
%bdata{1}=nanmean(xbinfos(big).daall.lfptracesaln(:,mididx-offset:mididx),2);
%bdata{2}=nanmean(xbinfos(sma).daall.lfptracesaln(:,mididx-offset:mididx),2);

cmap=colormap('parula');
[cmap,num,typ]=brewermap(2,'PuOr');
%ctop=cmap(round(size(cmap,1)*.2),:);
%cbot=cmap(round(size(cmap,1)*.80),:);
minlim=min([bdata{1}; bdata{2}]);
maxlim=max([bdata{1}; bdata{2}]);
binwid=floor((maxlim-minlim)/50);
    
figure; histogram(bdata{1},[minlim:binwid:maxlim])
hold on; histogram(bdata{2},[minlim:binwid:maxlim])
%sum(binfos(big).lickpost<370)/length(~isnan(binfos(big).lickpost))
ax1=gca;
cmap = get(gca,'ColorOrder');
ctop=cmap(1,:);
cbot=cmap(2,:);

aplot=plot(ax1,prctile(bdata{1},pcthigh)*[1 1],get(ax1,'YLim'),'b');
aplot=plot(ax1,prctile(bdata{1},pctlow)*[1 1],get(ax1,'YLim'),'b');
aplot=plot(ax1,prctile(bdata{2},pcthigh)*[1 1],get(ax1,'YLim'),'r');
aplot=plot(ax1,prctile(bdata{2},pctlow)*[1 1],get(ax1,'YLim'),'r');
xlabel(ax1,'lick');
ylabel(ax1,'counts');
 title(ax1,['big / small eye r histogram, session ' num2str(sessnum) ]);
savepath2=fullfile(plotparam.savepath,'eyerdep_hist_');
savename=[savepath2 'sess' num2str(sessnum) '_win_' num2str(round(offset/rate)) 's'];
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
if norm
    danormval=abs(nanmedian(horzcat(da{:}))); 
    betanormval=abs(nanmedian(horzcat(beta{:})));        %median of all values for all plot conditions for current probe
    for icond=1:length(da)        
        da{icond}=da{icond}./danormval;
        beta{icond}=beta{icond}./betanormval;
    end
end
for ix=1:2
    axa{ix}=subplot(1,2,ix);
    set(axa{ix},'units','pixels');
    axpos{ix}=get(axa{ix},'position');
    hold(axa{ix},'on');
    ylabel(axa{ix},'\beta LFP, \muV^2','interpreter','tex');
    xlabel(axa{ix},'\Delta[DA], nM','interpreter','tex');
end
for itype=1:2
    seltrialstop{itype}=btrials{itype}(bdata{itype}>prctile(bdata{itype},pcthigh));
    seltrialsbot{itype}=btrials{itype}(bdata{itype}<prctile(bdata{itype},pctlow));   
    seltrialstop{itype}=intersect(seltrialstop{itype},datm{sesid}{itype}{targid}.trialnums);
    seltrialstop{itype}=intersect(seltrialstop{itype},betatm{sesid}{itype}{targidlfp}.trialnums);
    seltrialsbot{itype}=intersect(seltrialsbot{itype},datm{sesid}{itype}{targid}.trialnums);
    seltrialsbot{itype}=intersect(seltrialsbot{itype},betatm{sesid}{itype}{targidlfp}.trialnums);    
   
    datop{itype}=da{itype}(find(ismember(datm{sesid}{itype}{targid}.trialnums,seltrialstop{itype})));
    betatop{itype}=beta{itype}(find(ismember(betatm{sesid}{itype}{targidlfp}.trialnums,seltrialstop{itype})));
    betabot{itype}=beta{itype}(find(ismember(betatm{sesid}{itype}{targidlfp}.trialnums,seltrialsbot{itype})));
    dabot{itype}=da{itype}(find(ismember(datm{sesid}{itype}{targid}.trialnums,seltrialsbot{itype})));
    
    hold(axa{itype},'on'); 
    if itype==1
    title(axa{itype},['big ' dapair{ip} ' ' lfppair{ip} ]);
    else
            title(axa{itype},['small ' dapair{ip} ' ' lfppair{ip} ]);
    end
    topgood=find(~isnan(datop{itype}) & ~isnan(betatop{itype}));
    botgood=find(~isnan(dabot{itype}) & ~isnan(betabot{itype}));
    datop{itype}=datop{itype}(topgood);
    betatop{itype}=betatop{itype}(topgood);
    dabot{itype}=dabot{itype}(botgood);
    betabot{itype}=betabot{itype}(botgood);
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
end
if norm
    savepath=[savepath 'norm_'];
end
savename=[savepath 'sess' num2str(sessnum) '_' event '_' dapair{ip} '_' lfppair{ip} '_win_' num2str(round(offset/rate)) 's_' condlabel];
savefig(figsess,savename);
saveas(figsess,savename,'jpg')
print(figsess,savename,'-painters','-depsc');
end
