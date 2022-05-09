function plotnormscatter(groupvals,plotparam,varargin)
figpos=[50,50,1100,600];
%figsess=figure('color',[1 1 1]);
%set(figsess,'position',figpos);
figsess=figure('visible','off');     %figure for each channel
if ispc
figsess=figure('visible','on');     %figure for each channel
end
set(figsess,'position',figpos,'color',[1 1 1]);
set(0,'CurrentFigure',figsess);    %set figure handle to current figure
axa={};
axsiz=[350 300];
off=75;
mar=100;
mark=[0 0.1 .7; 0.8 0 .5; 0 .7 0;.7 0 0; 0.4 0 .7];
sessmark={'.','+','o','^','s','*','v','p','h','x'};
sessmarktxt={'.','+','o','\Delta','sq','*','\nabla','p','h','x'};
%mark=[linspecer(10,'qualitative'); linspecer(10,'qualitative')];
sessmark={'sq','o','^','+','v','*','p','h','x','d','<','>','o','^','v','<','>','+','*','p','h','x','d','<','>'};
sessmarktxt={' sq','\o','\Delta','+','\nabla','*','p','h','\times','d','<','>','.','\Delta','\nabla','<','>','+','*','p','h','x','d','<','>'};
marksize=[200 200 150 200 200 200 200 200 200 200 200 200 50 50 50 50 50 50 50 50 50 50 50 200 200 200 200 200 200 200 200 200 200 200 200 200];
cmapp=[1 .5 0;1 .55 0; 1 .6 0; 1 .65 0; 1 .7 0; 1 .75 0; 1 .8 0; 1 .85 0; 1 .9 0];
cmapc=[1 .5 0;1 .55 0; 1 .6 0; 1 .65 0; 1 .7 0; 1 .75 0; 1 .8 0; 1 .85 0; 1 .9 0];
cmap2=cool;
cmap1=hot;
markc=cmap1(1:4:end,:);
markp=cmap2(1:5:end,:);
%markc=[linspecer(10,'qualitative'); linspecer(10,'qualitative')];
%markp=brewermap(13,'RdBu');
markp=[brewermap(8,'Accent'); linspecer(10,'qualitative'); brewermap(8,'RdBu');];
markc=[brewermap(8,'Accent'); linspecer(10,'qualitative'); brewermap(8,'RdBu');];

argnum=1;
targval='targpeak';         %targeted metric in datm metrics to use for
fontsize=14;
trialgrps=plotparam.trialgrps;
savepath=fullfile(plotparam.savepath,'scatternorm_');
sessnums=plotparam.sessnums;
daregions={'c','p'};
dalabels={'cn','put'};
conditions={'left','right'};
condlabels={'L','R'};
ratefscv=10;
plotvar={};
event='targ';

yaxlabel='\Delta[DA] (nM) ';
normda=0;
meanda=0;
grpreg=1;
scattypes=1;
plotlfp=0;
datype='goodtrials';        %all good trials for xinfo data
xflag=0;            %get xinfo data
xinfo={};
sessiontypes={};
tlabel=0;
xcflag=0;
plotv=0;
trtable={};
trlabels={};
simple=1;
plotlines=0;
plotsens=1; 
avgsites=0;
subj='patra';
gain=1;
while argnum<=length(varargin)
    switch varargin{argnum}
        
        case 'cleo'
            subj='cleo';
        case 'gain'
            plotsens=1;
            gain=1;
        case 'avgsites'
            avgsites=1;
        case 'norm'
            normda=1;
            
    end
    argnum=argnum+1;
end
for ip=1
    axa{ip}=subplot(1,1,ip);   hold(axa{ip},'on');
    set(axa{ip},'units','pixels');
    axpos{ip}=get(axa{ip},'position');   
    %set(axa{ip},'xTickLabelRotation',90)
    set(axa{ip},'position',[(ip-1)*axsiz(1)+ip*(mar/2)+off,125,axsiz(1),axsiz(2)])
end
sitevals={};
sitevals{1}=groupvals(1).sitevals;
sitevals{2}=groupvals(2).sitevals;
targval=groupvals.metric;
event=groupvals.event;
signal=groupvals.signal;


%plotting yvals
empties=[];
numplots=0;
empties=cellfun('isempty',{sitevals{1}.vals});
empties2=cellfun('isempty',{sitevals{2}.vals});
emp=unique([empties empties2]);
sitevals{1}(emp)=[];
sitevals{2}(emp)=[];
numplots=length(sitevals{1});
allsites=unique([sitevals{1}.site]);
if avgsites
    numplots=length(allsites);
end
ise=1;
countp=1;
countc=1;
sitecount=1;
avgda={};
diffgrps={};
[cmap,num,typ]=brewermap(5,'RdBu');

for sitenum=1:numplots
avgy{1}=nan(1,2);
ciy{1}=nan(1,2);
avgy{2}=nan(1,2);
ciy{2}=nan(1,2);
dasite='';
sessnum=0;
dapercond={};
plotvals={};
pid=1;
if ~avgsites
    dasite=sitevals{1}(sitenum).site;
    sessnum=sitevals{1}(sitenum).sessnum;
    plotvals{1}=sitevals{1}(sitenum).vals;
    plotvals{2}=sitevals{2}(sitenum).vals;

else
    siteids=find(strcmp(allsites,allsites{sitenum}));
    dasite=sitevals{1}(siteids(1)).site;
    sessnum=sitevals{1}(siteids(1)).sessnum;
    for iis=1:length(siteids)
        for icond=1:length(sitevals{1}(1).vals)
            if iis==1
                plotvals{1}{icond}=sitevals{1}(siteids(iis)).vals{icond};
                                plotvals{2}{icond}=sitevals{2}(siteids(iis)).vals{icond};
            else
                plotvals{1}{icond}=[plotvals{1}{icond} sitevals{1}(siteids(iis)).vals{icond}];
                                plotvals{2}{icond}=[plotvals{2}{icond} sitevals{2}(siteids(iis)).vals{icond}];
            end
        end
    end
end
ise=find(sessnums==sessnum);

if ~isempty(plotvals{1})
    for ix=1:length(plotvals)
for icond=1:length(plotvals{ix})
dapercond{ix}{icond}=plotvals{ix}{icond};
if plotv
    dapercond{ix}{icond}=var(dapercond{ix}{icond},'omitnan');
end
end
    end

    sitecolor=[0 0 0];

%normalize values
danormval=1;
if normda
    for ix=1:length(plotvals)
    danormval=abs(nanmedian([plotvals{ix}{:}]));        %median of all values for all plot conditions for current probe
    for icond=1:length(plotvals{ix})        
        plotvals{ix}{icond}=plotvals{ix}{icond}./danormval;
    end
    end
end
%plot probe values if not plotting mean over all or not scatter
for ix=1:length(plotvals)
for icond=1:length(plotvals{ix})
    curplotval=plotvals{ix}{icond};
    ciy{ix}(icond)=nanstd(curplotval)./sqrt(length(curplotval(~isnan(curplotval))))*1.96;
    avgy{ix}(icond)=nanmean(curplotval);
    avgda{ix}(icond,sitecount)=avgy{ix}(icond);
    %is=find(strcmp(uniquesites,sites(targsi).site)==1); %idx for unique site name
end
diffy=avgy{ix}(1)-avgy{ix}(2);
if gain
    diffy=(avgy{ix}(1)-avgy{ix}(2))/(avgy{ix}(1)+avgy{ix}(2));
end
diffgrps{ix}(sitecount)=diffy;
end
%get signifiance for difference between first pair of coniditon groups
sig{1}=0;
sig{2}=0;
for ix=1:length(plotvals)
[sig{ix},pval]=ttest2(dapercond{ix}{1},dapercond{ix}{2});
if isnan(sig{ix})
    sig{ix}=0;
end
end
oppocolor=[0 0 0];
if contains(dasite,'p')
    sitecolor=cmap(2,:);
    oppocolor=cmap(1,:);
else
    sitecolor=cmap(4,:);
    oppocolor=cmap(5,:);
end
facecolor=sitecolor;
edgecolor=oppocolor;
marksize=40;
alphac=0;
alphaf=0;
lwid=1.5;
if sig{1}
    alphac=.75;
else
    alphac=0.15;
end
if sig{2}
    alphaf=.75;
else
    alphaf=.15;
end  
    scatter(axa{1},diffgrps{2}(sitecount),...
        diffgrps{1}(sitecount),marksize,'o','markerfacecolor',facecolor,'markeredgecolor',edgecolor,...
        'MarkerEdgeAlpha',alphac,'MarkerFaceAlpha',alphaf,'linewidth',lwid);     
    
end
sitecount=sitecount+1;
set(findall(axa{1},'-property','FontSize'),'FontSize',fontsize)
end
lims=[-0.4 0.4];
 maxlim=max(nanmean([diffgrps{:}],2))+3*max(nanstd([diffgrps{:}],[],2));
    minlim=min(nanmean([diffgrps{:}],2))-3*max(nanstd([diffgrps{:}],[],2));
    lims=[minlim maxlim];

ylim(axa{1},lims)
xlim(axa{1},lims)
aline=plot(axa{1},[0 0],[-.4 .4],'--','color',[0 0 0],'linewidth',.5);
aline.Color(4)=.5;
aline=plot(axa{1},[-.4 .4],[0 0],'--','color',[0 0 0],'linewidth',.5);
aline.Color(4)=.5;
labels={};
savelabels={};
for ix=1:length(groupvals)
    common=intersect(groupvals(ix).cond{1},groupvals(ix).cond{2});
    unique1=find(~contains(groupvals(ix).cond{1},common));
    unique2=find(~contains(groupvals(ix).cond{2},common));
    if isempty(unique2)
            unique2=find(~strcmp(groupvals(ix).cond{2},common));
    end
    var1=groupvals(ix).cond{1}{unique1};
    var2=groupvals(ix).cond{2}{unique2};    
  %  labels{ix}=['$\mathsf{' var1(1) '-' var2(1) '\over ' var1(1) '+' var2(1) '}$'];
        labels{ix}=['$\mathrm{' var1(1) '-' var2(1) '\over ' var1(1) '+' var2(1) '}$'];

    savelabels{ix}=[var1 '-' var2];
end
ylabel(axa{1},labels{1},'interpreter','latex','fontsize',20);
xlabel(axa{1},labels{2},'interpreter','latex','fontsize',20);

titletext=[signal ' | ' event ' | ' targval];  
set(axa{1},'units','pixels')
axpos=get(axa{1},'position');
text(axa{1},0,axpos(4)+50,titletext,'units','pixels','fontweight','bold','fontsize',fontsize,'interpreter','none');

text(axa{1},axpos(3)-20,axpos(4)-20,'CN','units','pixels','fontsize',fontsize,'interpreter','none','color',cmap(5,:));
text(axa{1},axpos(3)-20,axpos(4)-40,'put','units','pixels','fontsize',fontsize,'interpreter','none','color',cmap(1,:));


savecond=[savelabels{1} 'vs' savelabels{2}];
savename=[savepath event '_' savecond '_' targval ];
if gain
    savename=[savename '_gain'];
end
if strcmp(signal,'lfp')
    savename=[savename '_lfp'];
end
if plotsens
    savename=[savename '_sens'];
end
if normda
    savename=[savename '_z'];
end
if simple
    savename=[savename '_simp'];
end
if avgsites
    savename=[savename '_avgsites'];
end

savefig(figsess,savename);
saveas(figsess,savename,'jpg')
print(figsess,savename,'-painters','-depsc');
end
