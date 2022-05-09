function plotdasummary(datm,plotparam,varargin)
%sites defined from excel sheet indicating different sites based on sessnum
%& site label from getsites function
%plot bar scatter plots of properties of da peaks or win during
%event/targval types, values from datm variables
figpos=[50,50,1800,600];
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
mark=[linspecer(10,'qualitative'); linspecer(10,'qualitative')];
sessmark={'o','sq','o','^','+','v','*','p','h','x','d','<','>','o','sq','o','^','+','v','*','p','h','x','d','<','>'};
sessmarktxt={'.',' sq','o','\Delta','+','\nabla','*','p','h','x','d','<','>','.',' sq','o','\Delta','+','\nabla','*','p','h','x','d','<','>'};
marksize=[50 200 150 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 50 200 150 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200];
cmapp=[1 .5 0;1 .55 0; 1 .6 0; 1 .65 0; 1 .7 0; 1 .75 0; 1 .8 0; 1 .85 0; 1 .9 0];
cmapc=[1 .5 0;1 .55 0; 1 .6 0; 1 .65 0; 1 .7 0; 1 .75 0; 1 .8 0; 1 .85 0; 1 .9 0];
cmap2=cool;
cmap1=hot;
markc=cmap1(1:5:end,:);
markp=cmap2(1:5:end,:);

argnum=1;
targval='targpeak';         %targeted metric in datm metrics to use for
fontsize=14;
trialgrps=plotparam.trialgrps;
savepath=fullfile(plotparam.savepath,'summary_da_');
sessnums=plotparam.sessnums;
daregions={'c','p'};
dalabels={'cn','put'};
conditions={'left','right'};
condlabels={'L','R'};
plotnums=1:length(daregions)*length(conditions);

plotvar={};
event='targ';
targdasites=plotparam.dasites;
sites=getsites(sessnums,targdasites);
uniquesites=unique({sites(1:end).site});
cnsites=uniquesites(contains(uniquesites,'c'));
psites=uniquesites(contains(uniquesites,'p'));
normda=0;
meanda=0;
grpreg=0;
scattypes=0;
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'event'
            argnum=argnum+1;    
            event=varargin{argnum};         %event period, eg. fix, targ, outcome
        case 'metric'
            argnum=argnum+1;
            targval=varargin{argnum};       %datm variable, eg. targpeak, targwin, etc.
        case 'norm'
            normda=1;
        case 'mean'
            meanda=1;
        case 'grpreg'
            grpreg=1;       %group regions
        case 'scattypes'
            scattypes=1;          
    end
    argnum=argnum+1;
end
if grpreg
    %group regions
    plotnums=1:length(conditions);
end
for ic=1:length(conditions)
    if ~grpreg
        %separate plots  for different conditions & regions
        for ir=1:length(daregions)
            plotvar{ir+length(conditions)*(ic-1)}=[dalabels{ir} '_' condlabels{ic}(1)];
        end
    else  
        %separate plots for different conditions
        plotvar{ic}=[condlabels{ic}(1)];
    end
end
sessiontypes={'big','small'};
%plot da differences (max, mean) targ window big vs small for 
%CN Left, CN right, P left, P right
clf(figsess,'reset');
set(figsess,'color',[1 1 1]);
axpos={};
for ip=plotnums
    axa{ip}=subplot(1,length(plotnums),ip);   hold(axa{ip},'on');
    set(axa{ip},'units','pixels');
    axpos{ip}=get(axa{ip},'position');
    if ~scattypes
    set(axa{ip},'xtick',1:length(sessiontypes),'xticklabel',sessiontypes);
    set(axa{ip},'xlim',[0 length(sessiontypes)+1]);
    end
    %set(axa{ip},'xTickLabelRotation',90)
    set(axa{ip},'position',[(ip-1)*axsiz(1)+ip*(mar/2)+off,75,axsiz(1),axsiz(2)])
    ax=text(axa{ip},10,axsiz(2),plotvar{ip},'units','pixels','fontweight','bold');
end
titletext=[event ' | ' targval];  
set(axa{1},'units','pixels')
axpos=get(axa{1},'position');
text(axa{1},-100,axpos(4)+100,titletext,'units','pixels','fontweight','bold');
labelx={};
ranj=[];        %random jitter x axis for each site
countc=1;
countp=1;
for is=1:length(uniquesites)
    labelx{is}=[num2str(uniquesites{is})];  
    axpos=get(axa{1},'position');
        daregion=contains(uniquesites(is),'c');     %1=='c'
    if daregion
        %CN
        if ~grpreg
        ax=text(axa{1},100+is*60,axpos(4)+50,labelx{is},'units','pixels','fontweight','bold','color',markc(countc,:));
        else
                    ax=text(axa{1},100+is*60,axpos(4)+50,labelx{is},'units','pixels','fontweight','bold','color',markc(countc,:));
        end
        countc=countc+1;
    else
        %put
        if ~grpreg
        ax=text(axa{1},100+is*60,axpos(4)+50,labelx{is},'units','pixels','fontweight','bold','color',markp(countp,:));
        else
                    ax=text(axa{1},100+is*60,axpos(4)+50,labelx{is},'units','pixels','fontweight','bold','color',markp(countp,:));

        end
        countp=countp+1;
    end   
end
for is=1:length(sessnums)
    text(axa{1},100+is*75,axpos(4)+100,[num2str(sessnums(is)) ' ' sessmarktxt{is}],'interpreter','tex','units','pixels','fontweight','bold','color',[0 0 0]);
    ranj(is)=rand(1,1)*.8-.5;
end
    
    avgda={};
for icond=1:length(conditions)
    nc=0;
    np=0;
    nn=0;
for ise=1:length(sessnums)
    trialinfo=trialgrps(ise).trialinfo;
    avgy=nan(1,2);
    ciy=nan(1,2);
    cursessids=find([sites.sessnum]==sessnums(ise));
    curda={sites(cursessids).probeid};
    dasites={sites(cursessids).site};
for ida=1:length(curda)
    %each da separate subplot  
    lfpsites={};
    daregion=contains(curda(ida),'c');     %1=='c'
    dasite=dasites(ida);
    sitecolor=[0 0 1];
    if daregion
        %CN
        pid=1+(icond-1)*2;
        sitecolor=markc(find(contains(cnsites,dasite)),:);
        if grpreg
            pid=icond;      %group regions together, just separate conditions to separate plots
            sitecolor=markc(find(contains(cnsites,dasite)),:);       %separate color group for p
        end
        nc=nc+1;
        nn=nc;
    else
        %put
        pid=2+(icond-1)*2;
        sitecolor=markp(find(contains(psites,dasite)),:);
        if grpreg
            pid=icond;      %group regions together, just separate conditions to separate plots
            sitecolor=markp(find(contains(psites,dasite)),:);       %separate color group for p
        end
        np=np+1;
        nn=np;
    end  
    danormval=1;
    if normda
        dapertype={};
    for itype=1:length(sessiontypes)
        %get median value for both types big/small for normalization factor
        targt=find(contains(trialinfo(itype).names,conditions(icond))==1);
        trialnums=trialinfo(itype).nums{targt};
        targsi=find([sites.sessnum]==sessnums(ise) & ...
            strcmp({sites.probeid},curda(ida)));
        if ~isempty(targsi)
        ich=sites(targsi(1)).ch;
        datarg=datm{ise}{itype}{ich};
        if ~isempty(datarg)
        trialids=find(ismember(datarg.trialnums,trialnums)==1);
        davals=getfield(datarg,targval);    %all da values for good trials for sess/type
        dapertype{itype}=davals(trialids);         %da values for specific condition (le/ri) trials
        end
        end
    end
    if ~isempty(dapertype)
        danormval=nanmedian([dapertype{1} dapertype{2}]);
    end
    end
dvals={};
for itype=1:length(sessiontypes)
    %big vs small
    targt=find(contains(trialinfo(itype).names,conditions(icond))==1);
    trialnums=trialinfo(itype).nums{targt};
    targsi=find([sites.sessnum]==sessnums(ise) & ...
        strcmp({sites.probeid},curda(ida)));
    if ~isempty(targsi)
    ich=sites(targsi(1)).ch;
    datarg=datm{ise}{itype}{ich};
    if ~isempty(datarg)
    trialids=find(ismember(datarg.trialnums,trialnums)==1);
    davals=getfield(datarg,targval);    %all da values for good trials for sess/type
    yvals=davals(trialids);         %da values for specific condition (le/ri) trials
    if normda
        yvals=yvals./danormval;
    end
    %xplot=repmat(itype,1,length(yvals))+ranj(ise);
    ciy(itype)=nanstd(yvals)./sqrt(length(yvals(~isnan(yvals))))*1.96;
    avgy(itype)=nanmean(yvals);
    avgda{pid}(itype,nn)=avgy(itype);
    is=find(strcmp(uniquesites,sites(targsi).site)==1); %idx for unique site name
    if ~meanda && ~scattypes
        %plot ci line for each da value
    aci=plot(axa{pid},[itype itype]+ranj(ise),[avgy(itype)-ciy(itype)  avgy(itype)+ciy(itype)],'-','linewidth',1,'color',sitecolor);  
    aci.Color(4)=.5;
    end
    end
    end
    dvals{itype}=yvals;
end
[sig,pval]=ttest2(dvals{1},dvals{2});
if ~scattypes
    %plot as bar graph with types on x-axis
    if sig
        alphac=1;
        lwid=2;
    else
        alphac=0.3;
        lwid=.75;
    end
     if ~strcmp(sessmark{ise},'.')
scatter(axa{pid},(1:length(sessiontypes))+ranj(ise),...
            avgy,marksize(ise),sessmark{ise},'markeredgecolor',sitecolor,...
            'MarkerEdgeAlpha',alphac,'linewidth',lwid);
     else
            scatter(axa{pid},(1:length(sessiontypes))+ranj(ise),...
            avgy,marksize(ise),sessmark{ise},'markerfacecolor',sitecolor,...
            'MarkerEdgeAlpha',alphac,'linewidth',lwid);
     end
end
if scattypes && length(sessiontypes)==2
    %scatter plot with 1 type on y-axis and 2nd type on x-axis
    if sig
        alphac=1;
        lwid=2;
    else
        alphac=0.3;
        lwid=.75;
    end
    if ~strcmp(sessmark{ise},'.')
    scatter(axa{pid},avgy(2),...
            avgy(1),marksize(ise),sessmark{ise},'markeredgecolor',sitecolor,...
            'MarkerEdgeAlpha',alphac,'linewidth',lwid);
            else
        scatter(axa{pid},avgy(2),...
            avgy(1),marksize(ise),sessmark{ise},'markerfacecolor',sitecolor,...
            'MarkerEdgeAlpha',alphac,'MarkerFaceAlpha',alphac,'linewidth',lwid);
    end
end
if ~meanda && ~scattypes
            %plot line between types for each da value pairs
aline=plot(axa{pid},(1:length(sessiontypes))+ranj(ise),avgy,'-','linewidth',1,'color',sitecolor);      
aline.Color(4)=.3;
end
set(findall(axa{pid},'-property','FontSize'),'FontSize',fontsize)

end
end
end
if meanda && ~scattypes
for ip=plotnums
    %mean changes for all values from one type to another type
    meansda=nanmean(avgda{ip},2);
    stdsda=nanstd(avgda{ip},[],2);
    cisda=stdsda./sqrt(size(avgda{ip},2)).*1.96;
    cilines=meansda-cisda ;
    cilinespos=meansda+cisda;
    aci=plot(axa{ip},[1 1],[cilines(1) cilinespos(1)],'-','linewidth',2,'color',[0 0 0]);  
    aci2=plot(axa{ip},[2 2],[cilines(2) cilinespos(2)],'-','linewidth',2,'color',[0 0 0]);  
    aci.Color(4)=.75;     aci2.Color(4)=.75;
    aline=plot(axa{ip},(1:length(sessiontypes)),meansda,'-','linewidth',1,'color',[0 0 0]);  
    aline.Color(4)=.7;
end
end
ylabel(axa{1},'\Delta[DA] (nM)');
if scattypes
    maxlim=max(nanmean([avgda{:}],2))+3*max(nanstd([avgda{:}],[],2));
    for ip=plotnums
        xlim(axa{ip},[0.5 maxlim]);
        ylim(axa{ip},[0.5 maxlim]);
        %plot unitary line
        aci=plot(axa{ip},[0.5 maxlim],[0.5 maxlim],'--','linewidth',1,'color',[0 0 0]);  
        aci.Color(4)=.25; 
        ylabel(axa{ip},['\Delta[DA] (nM) '  sessiontypes{1}]);
    xlabel(axa{ip},['\Delta[DA] (nM) ' sessiontypes{2}]);
    end
end
savename=[savepath event '_' targval];
if normda
    savename=[savename '_z'];
end
if meanda
    savename=[savename '_avg'];
end
if scattypes
    savename=[savename '_scat'];
end

savefig(figsess,savename);
saveas(figsess,savename,'jpg')
print(figsess,savename,'-painters','-depsc');

end
