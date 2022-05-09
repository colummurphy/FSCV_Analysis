sessnum=83;
event='targ';

big=find(contains({binfos.sessionid},num2str(sessnum)) & ...
    strcmp({binfos.event},event) & ...
    contains({binfos.sessiontype},'big'));
sma=find(contains({binfos.sessionid},num2str(sessnum)) & ...
    strcmp({binfos.event},event) & ...
    contains({binfos.sessiontype},'small'));

yy=prctile(binfos(58).lickpost,75)
sum(binfos(58).lickpost<325)/length(~isnan(binfos(58).lickpost))
sum(binfos(62).lickpost<275)/length(~isnan(binfos(62).lickpost))
sum(binfos(62).lickpost<240)/length(~isnan(binfos(62).lickpost))
figure; histogram(binfos(130).lickpost,[100:10:600])
figure; histogram(binfos(130).lickpost)
figure; histogram(binfos(130).lickpost,[100:10:600])
hold on; histogram(binfos(134).lickpost,[100:10:600])
sum(binfos(134).lickpost<370)/length(~isnan(binfos(134).lickpost))
yy=prctile(binfos(134).lickpost,70)
yy=prctile(binfos(130).lickpost,70)
figure; histogram(binfos(146).lickpost)
hold on; histogram(binfos(150).lickpost)
figure; histogram(binfos(146).lickpost,[40:10:200])
hold on; histogram(binfos(150).lickpost,[40:10:200])
figure; histogram(binfos(146).lickpost,[40:2:200])
hold on; histogram(binfos(150).lickpost,[40:2:200])
sum(binfos(146).lickpost<80)/length(~isnan(binfos(146).lickpost))
sum(binfos(150).lickpost<65)/length(~isnan(binfos(150).lickpost))
yy=prctile(binfos(146).lickpost,70);
yy=prctile(binfos(150).lickpost,70);
ax1=gca;
aplot=plot(ax1,prctile(binfos(58).lickpost,75)*[1 1],[0 25],'b')
ax1=gca;
aplot=plot(ax1,prctile(binfos(58).lickpost,75)*[1 1],[0 25],'b')
aplot=plot(ax1,prctile(binfos(58).lickpost,75)*[1 1],[0 100],'b')
aplot=plot(ax1,prctile(binfos(58).lickpost,75)*[1 1],get(ax1,'YLim'),'b')
aplot=plot(ax1,prctile(binfos(58).lickpost,50)*[1 1],get(ax1,'YLim'),'b')
aplot=plot(ax1,prctile(binfos(62).lickpost,75)*[1 1],get(ax1,'YLim'),'r')
aplot=plot(ax1,prctile(binfos(62).lickpost,50)*[1 1],get(ax1,'YLim'),'r')
xlabel('lick')
ylabel('counts')
selids=binfos(58).lickpost<365;

thres=prctile(binfos(58).lickpost,75);
seltrials=binfos(58).seltrials(binfos(58).lickpost>thres);
dasmall=datm{8}{2}{4}.targpeak;
dabig=datm{8}{1}{4}.targpeak;
betabig=betatm{8}{1}{7}.targimwin;
betasmall=betatm{8}{2}{7}.targimwin;
seltrialsbig=seltrials;
seltrialsbigbot=binfos(58).seltrials(binfos(58).lickpost<prctile(binfos(58).lickpost,50));
seltrialssmall=binfos(62).seltrials(binfos(62).lickpost>prctile(binfos(62).lickpost,75));
seltrialssmallbot=binfos(62).seltrials(binfos(62).lickpost<prctile(binfos(62).lickpost,50));
dabigtop=dabig(find(ismember(datm{8}{1}{4}.trialnums,seltrialsbig)));
betabigtop=betabig(find(ismember(betatm{8}{1}{7}.trialnums,seltrialsbig)));
betabigbot=betabig(find(ismember(betatm{8}{1}{7}.trialnums,seltrialsbigbot)));
dabigbot=dabig(find(ismember(datm{8}{1}{4}.trialnums,seltrialsbigbot)));
figure; plot(dabigtop,betabigtop,'*');
hold on; plot(dabigbot,betabigbot,'o');
find(betabigtop>500)
betabigtop(64)=[];
dabigtop(64)=[];
figure; plot(dabigtop,betabigtop,'*');
hold on; plot(dabigbot,betabigbot,'o');
[rtop,ptop]=corr(dabigtop',betabigtop');
xy=find(~isnan(dabigbot));
[rbot,pbot]=corr(dabigbot(xy)',betabigbot(xy)');
xdata=dabigbot(xy);
ydata=betabigbot(xy);
[p,s]=polyfit(xdata,ydata,1);
slope=p(1);
intercep=p(2);
yfit=slope*xdata+intercep;
yresid=ydata-yfit;
ssresid=sum(yresid.^2);
sstotal=(length(ydata)-1)*var(ydata);
rsq=1-ssresid/sstotal;
ax=gca;
fitline=plot(ax,xdata,yfit,'-','color',[1 0 0]);
fitline.Color(4)=0.5;
xsiz=get(ax,'position');
set(ax,'units','pixels');
ylims=get(ax,'ylim');
text(ax,xsiz(3)+100,xsiz(4)+50,...
{['slope: ' num2str(slope)],...
['intercept: ' num2str(intercep)],['r: ' num2str(rbot)],...
['p: ' num2str(pbot)]},'color',[1 0 0],'units','pixels');

[p,s]=polyfit(dabigtop,betabigtop,1);
slope=p(1);
intercep=p(2);
yfit=slope*xdata+intercep;
yresid=ydata-yfit;
ssresid=sum(yresid.^2);
sstotal=(length(ydata)-1)*var(ydata);
rsq=1-ssresid/sstotal;
ax=gca;
fitline=plot(ax,xdata,yfit,'-','color',[0 0 1]);
fitline.Color(4)=0.5;
xsiz=get(ax,'position');
set(ax,'units','pixels');
ylims=get(ax,'ylim');
text(ax,xsiz(3)+200,xsiz(4)+200,...
{['slope: ' num2str(slope)],...
['intercept: ' num2str(intercep)],['r: ' num2str(rtop)],...
['p: ' num2str(ptop)]},'color',[0 0 1],'units','pixels');

