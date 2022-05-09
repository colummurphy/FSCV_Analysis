function plotxsummary(xinfos,plotparam,varargin)
%sites defined from excel sheet indicating different sites based on sessnum
%& site label from getsites function
%plot bar scatter plots of properties of da peaks or win during
%event/targval types, values from datm variables
figpos=[50,50,1500,800];
%figsess=figure('color',[1 1 1]);
%set(figsess,'position',figpos);
figsess=figure('visible','off');     %figure for each channel
if ispc
figsess=figure('visible','on');     %figure for each channel
end
set(figsess,'position',figpos,'color',[1 1 1]);
set(0,'CurrentFigure',figsess);    %set figure handle to current figure
axa={};
axsiz=[500 450];
off=100;
mar=150;
cmap=colorcube;         %or lines
cmap=cmap(10:end,:);
mark=linspecer(10);
mark=[linspecer(10,'qualitative'); linspecer(10,'qualitative')];
sessmark={'o','sq','o','^','+','v','*','p','h','x','d','<','>','o','sq','o','^','+','v','*','p','h','x','d','<','>'};
sessmarktxt={'.',' sq','o','\Delta','+','\nabla','*','p','h','x','d','<','>','.',' sq','o','\Delta','+','\nabla','*','p','h','x','d','<','>'};
marksize=[50 200 150 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 50 200 150 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200];
%marksize=[250 200 150 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200];
cmap2=cool;
cmap1=hot;
markc=cmap1(1:5:end,:);
markp=cmap2(1:5:end,:);

argnum=1;
datype='postrials';
targval={'damaxts','lfpmints','lfppostmaxts','delt_lfppostmax_damax','delt_lfpmin_damax'};         %targeted metric in datm metrics to use for
fontsize=18;
trialgrps=plotparam.trialgrps;
savepath=fullfile(plotparam.savepath,'summary_da_');
sessnums=plotparam.sessnums;
sesstype={'big','small','big','small','big','small','big','big','small','small'};
rate=10;
daregions={'c','p'};
dalabels={'cn','put'};
conditions={'all','all','left','left','right','right','left','right','left','right'};
plotvar={};
plotx={};
event='targ';
targdasites=plotparam.dasites;
uniquesites=unique({xinfos.dasite});
cnsites=uniquesites(contains(uniquesites,'c'));
psites=uniquesites(contains(uniquesites,'p'));
plotz=0;
grpreg=0;
scattypes=0;

while argnum<=length(varargin)
    switch varargin{argnum}
        case 'event'
            argnum=argnum+1;
            event=varargin{argnum};
        case 'datype'
            argnum=argnum+1;
            if strcmp(varargin{argnum},'daall')
                datype='goodtrials';
            else
                datype=varargin{argnum};
            end
        case 'sesstypes'
            argnum=argnum+1;        %provide sesstype associated with each condition provided below in pairs
            sesstype=varargin{argnum};
        case 'metric'
            argnum=argnum+1;
            targval=varargin{argnum};
        case 'conditions'
            argnum=argnum+1;
            conditions=varargin{argnum};
        case 'norm'
            plotz=1;
             case 'grpreg'
            grpreg=1;
        case 'scattypes'
            scattypes=1;       
    end
    argnum=argnum+1;
end

%make sure lengths of conidtions/sesstype variablle same
if length(conditions)~=length(sesstype)
    error('conditions and sesstype variable should be same length');
end
plotnums=1;
[dapair,lfppair]=getsitepairs(targdasites);
plotvar={};
if ~grpreg
for ir=1:length(daregions)
    plotvar{ir}=[dalabels{ir}];
end
plotnums=1:length(daregions);
else
    plotvar{1}=[dalabels{1} ', ' dalabels{2}];
end

%plot metrics big vs small different conditions
for ix=1:length(conditions)
    plotx{ix}=[sesstype{ix}(1:3) ' ' conditions{ix}(1:3)];
end

for ival=1:length(targval)
cid=1;
xvar={};
xvarc={};
clf(figsess,'reset');
set(figsess,'color',[1 1 1]);
axpos={};
for ip=plotnums
    axa{ip}=subplot(1,length(plotnums),ip);   hold(axa{ip},'on');
    set(axa{ip},'units','pixels');
    axpos{ip}=get(axa{ip},'position');
    if ~scattypes
    set(axa{ip},'xtick',1:length(plotx),'xticklabel',plotx);
    set(axa{ip},'xlim',[0 length(plotx)+1]);
    set(axa{ip},'xTickLabelRotation',45)
    end
    set(axa{ip},'position',[(ip-1)*axsiz(1)+ip*(mar/2)+off,100,axsiz(1),axsiz(2)])
    ax=text(axa{ip},10,axsiz(2),plotvar{ip},'units','pixels','fontweight','bold');
end
titletext=[event ' | ' targval{ival}];  
set(axa{1},'units','pixels')
axpos=get(axa{1},'position');
text(axa{1},-100,axpos(4)+150,titletext,'units','pixels','fontweight','bold','interpreter','none');
labelx={};
ranj=[];        %random jitter x axis for each site

%site legends
countc=1;
countp=1;
for is=1:length(uniquesites)
    labelx{is}=[num2str(uniquesites{is})];  
    axpos=get(axa{1},'position');
        daregion=contains(uniquesites(is),'c');     %1=='c'
    if daregion
        %CN
        ax=text(axa{1},-150+is*65,axpos(4)+50,labelx{is},'units','pixels','fontweight','bold','color',markc(countc,:));
        countc=countc+1;
    else
        %put
        ax=text(axa{1},-150+is*65,axpos(4)+50,labelx{is},'units','pixels','fontweight','bold','color',markp(countp,:));
        countp=countp+1;
    end   
end
%sess legends
for is=1:length(sessnums)
    text(axa{1},-150+is*90,axpos(4)+100,[num2str(sessnums(is)) ' ' sessmarktxt{is}],'interpreter','tex','units','pixels','fontweight','bold','color',[0 0 0]);
    ranj(is)=rand(1,1)*.65-.35;
end

xplot={};
yplot={};
yplotci={};
yplotraw={};
cplot={};
xvar={};
xvarc={};
nc=0;
np=0;
nn=0;
counts=0;
for ise=1:length(sessnums)
    sessid=num2str(sessnums(ise));
    trialinfo=trialgrps(ise).trialinfo;
    avgy=nan(1,2);
    ciy=nan(1,2);
    targrow=find(contains({xinfos.sessionid},sessid));
    targdasites=unique({xinfos(targrow).siteda});     %da sites for curr sess
    %cursessids=find(contains({xinfos.sessionid},sessid));
    %curda={sites(cursessids).probeid};
    dasites=unique({xinfos(targrow).dasite});
for ida=1:length(dasites)
    %each da separate subplot  
    lfpsites={};
    dasite=dasites(ida);
    daregion=contains(dasite,'c');     %1=='c'
    sitecolor=[0 1 0];
    if daregion
        %CN
        pid=1;
        sitecolor=markc(find(contains(cnsites,dasite)),:);
        nc=nc+1;
        nn=nc;
        if grpreg
            counts=counts+1;
        nn=counts;
        end
    else
        %put
        pid=2;
        sitecolor=markp(find(contains(psites,dasite)),:);
        np=np+1;
        nn=np;
        if grpreg
            pid=1;
                        counts=counts+1;
            nn=counts;
        end
    end 
    daid=find(strcmp(dapair,targdasites(ida)));
    if isempty(daid)
        disp(['empty ' targdasites(ida) ]);
    end
    lfptarg=lfppair(daid);
    avgy=nan(1,2);
    ciy=nan(1,2);
for icond=1:length(conditions)
    c=(icond-1)*2;
for itype=1:length(sesstype)
    %big vs small
    ttype=find(contains({'big','small'},sesstype{itype}));
    targt=find(contains(trialinfo(ttype).names,conditions(icond))==1);
    targrow=find(contains({xinfos.dasite},dasite) & ...
        contains({xinfos.sitelfp},lfptarg) & ...
        strcmp({xinfos.event},event) & ... 
        contains({xinfos.sessionid},sessid) &...
        contains({xinfos.sessiontype},sesstype(itype)));    
if ~isempty(targrow)
    targrow=targrow(1);     %get first corresponding lfp pair or matching row
    curdata=getfield(xinfos(targrow),'daall');    %'dapos' or 'daneg' types
    trialnums=trialinfo(ttype).nums{targt};   
    datatrials=getfield(xinfos(targrow),datype);     %pos/neg/all trials
    trialids=find(ismember(curdata.trials,trialnums) & ismember(curdata.trials,datatrials));
    davals=getfield(curdata,targval{ival});    %all da values for good trials for sess/type
    if contains(targval{ival},'ts')
        %time stamp, get in seconds and subtract alignidx
        davals=(davals-curdata.mididx)./rate;
    end
    if contains(targval{ival},'delt')
        %relative time stamp, get in seconds and subtract alignidx
        davals=(davals)./rate;
    end
    yvals=davals(trialids);         %da values for specific condition (le/ri) trials
    %xplot=repmat(itype,1,length(yvals))+ranj(ise);
    ciy(itype)=nanstd(yvals)./sqrt(length(yvals(~isnan(yvals))))*1.96;
    avgy(itype)=nanmean(yvals);
    is=find(strcmp(uniquesites,dasite)); %idx for unique site name
    if ~plotz && ~scattypes
    aci=plot(axa{pid},[itype itype]+ranj(ise)+c,[avgy(itype)-ciy(itype)  avgy(itype)+ciy(itype)],'-','linewidth',1,'color',sitecolor);  
    aci.Color(4)=.35;
    end
   % xvar{pid}{nc}(itype,icond)=itype;
   % xvarc{pid}{nc}(itype,icond)=icond;
    xplot{pid}{nn}(itype,icond)=ise;
    yplot{pid}{nn}(itype,icond)=avgy(itype);
    yplotci{pid}{nn}(itype,icond)=ciy(itype);
    yplotraw{pid}{nn}{itype,icond}=davals(trialids);
    cid=cid+1;
end
end
[sig,pval]=ttest2(yplotraw{pid}{nn}{1,icond},yplotraw{pid}{nn}{2,icond});
if ~plotz && ~scattypes
scatter(axa{pid},(1:length(sesstype))+ranj(ise)+c,...
            avgy,marksize(ise),sessmark{ise},'markeredgecolor',sitecolor,...
            'MarkerEdgeAlpha',.65,'linewidth',1);
aline=plot(axa{pid},(1:length(sesstype))+ranj(ise)+c,avgy,'--','linewidth',1,'color',sitecolor);      
aline.Color(4)=.3;
end
if scattypes && length(sesstype)==2
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
            'MarkerEdgeAlpha',alphac,'MarkerFaceAlpha',alphac,'linewidth',lwid);
    else
        scatter(axa{pid},avgy(2),...
            avgy(1),marksize(ise),sessmark{ise},'markerfacecolor',sitecolor,...
            'MarkerEdgeAlpha',alphac,'MarkerFaceAlpha',alphac,'linewidth',lwid);
    end
end
cplot{pid}{nn}{icond}=sitecolor;

end
set(findall(axa{pid},'-property','FontSize'),'FontSize',fontsize)
end
end

if plotz && ~scattypes
    %percent change
for pi=plotnums
   % xv=xvar{pi};
   % xvc=xvarc{pi};
    xdata=xplot{pi};
    ydata=yplot{pi};
    ydataci=yplotci{pi};
    ydataraw=yplotraw{pi};
    ynormval=[];
    %get median per probe for all conditions (condition 1 = all)
    for nn=1:length(ydataraw)
        ynormval(nn)=nanmedian([ydataraw{nn}{1,1} ydataraw{nn}{2,1}]);
    end
    danorm={};
    meansperprobe={};
    for icond=1:length(conditions)
        c=(icond-1)*2;
        for itype=1:length(sesstype)
            for nn=1:length(ydataraw)
                danorm{icond,itype}{nn}=ydataraw{nn}{itype,icond}./ynormval(nn);
                meansperprobe{icond,itype}(nn)=nanmean(danorm{icond,itype}{nn});
                scatter(axa{pi},itype+ranj(xplot{pi}{nn}(itype,icond))+c,...
                    meansperprobe{icond,itype}(nn),marksize(xplot{pi}{nn}(itype,icond)),sessmark{xplot{pi}{nn}(itype,icond)},'markeredgecolor',cplot{pi}{nn}{icond},...
                    'MarkerEdgeAlpha',.65,'linewidth',1);
            end
        end
    end
    for icond=1:length(conditions)
        c=(icond-1)*2;
        for itype=1:length(sesstype)
            daavg(itype)=nanmean(meansperprobe{icond,itype});
            dastd=nanstd(meansperprobe{icond,itype});
            daci(itype)=dastd./sqrt(length(meansperprobe{icond,itype}))*1.96;
            aci=plot(axa{pi},[itype itype]+c,[daavg(itype)-daci(itype)  daavg(itype)+daci(itype)],'-','linewidth',2,'color',[ 0 0 0]);  
            aci.Color(4)=.75;
        end
        aline=plot(axa{pi},(1:length(sesstype))+c,daavg,'--','linewidth',1,'color',[ 0 0 0]);      
        aline.Color(4)=.75;
    end
end
end
    
if scattypes
    maxlim=max(nanmean([yplot{:}{:}],2))+3*max(nanstd([yplot{:}{:}],[],2));
    minlim=min(nanmean([yplot{:}{:}],2))-3*max(nanstd([yplot{:}{:}],[],2));
    for ip=plotnums
        xlim(axa{ip},[minlim maxlim]);
        ylim(axa{ip},[minlim maxlim]);
        %plot unitary line
        aci=plot(axa{ip},[minlim maxlim],[minlim maxlim],'--','linewidth',1,'color',[0 0 0]);  
        aci.Color(4)=.25; 
        ylabel(axa{ip},[targval{ival}  ' ' sesstype{1} ' ' conditions{1}],'interpreter','none');
    xlabel(axa{ip},[targval{ival} ' ' sesstype{2} ' ' conditions{2}],'interpreter','none');
    end
else       
ylabel(axa{1},targval{ival},'interpreter','none');
end

savename=[savepath event '_' targval{ival}];
if length(sesstype)==2
    savename=[savepath event '_' targval{ival} '_' sesstype{1} '_' conditions{1}(1) '_' sesstype{2} '_' conditions{2}(1)];
end
if plotz
    savename=[savename '_z'];
end
if scattypes
        savename=[savename '_scat'];
end
savefig(figsess,savename);
saveas(figsess,savename,'jpg')
print(figsess,savename,'-painters','-depsc');

end

end
