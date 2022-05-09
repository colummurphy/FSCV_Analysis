function plotsummarycounts(counts,plotparam,varargin)
%counts positive/negative da peaks for given categories
%sites defined from excel sheet indicating different sites based on sessnum
%& site label from getsites function
%plot bar scatter plots of properties of da peaks or win during
%event/targval types, values from counts variables
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
sessmark={'.','sq','o','^','+','v','*','p','h','x','d','<','>'};
sessmarktxt={'.',' sq','o','\Delta','+','\nabla','*','p','h','x','d','<','>'};
marksize=[250 200 150 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200];
sessmark={'o','sq','o','^','+','v','*','p','h','x','d','<','>','o','sq','o','^','+','v','*','p','h','x','d','<','>'};
sessmarktxt={'.',' sq','o','\Delta','+','\nabla','*','p','h','x','d','<','>','.',' sq','o','\Delta','+','\nabla','*','p','h','x','d','<','>'};
marksize=[50 200 150 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200 50 200 150 200 200 200 200 200 200 200 200 200 200 200 200 200 200 200];
cmap2=cool;
cmap1=hot;
markc=cmap1(1:5:end,:);
markp=cmap2(1:5:end,:);

argnum=1;
targval='ratioposneg';         %targeted metric in counts metrics to use for
fontsize=14;
trialgrps=plotparam.trialgrps;
savepath=fullfile(plotparam.savepath,'summary_da_');
sessnums=plotparam.sessnums;
daregions={'c','p'};
dalabels={'cn','put'};
conditions={'left','right'};
condlabels={'L','R'};
plotnums=1:length(daregions);
plotvar={};
event='targ';
targdasites=plotparam.dasites;
sites=getsites(sessnums,targdasites);
uniquesites=unique({sites(1:end).site});
cnsites=uniquesites(contains(uniquesites,'c'));
psites=uniquesites(contains(uniquesites,'p'));

variables=counts{1}{1}.variable;
origvar=counts{1}{1}.variable;
origgrps={};
for ivar=1:length(origvar)
    origgrps{ivar}=[origvar{ivar}{1} '' origvar{ivar}{2}];
end
normda=0;
meanda=0;
grpreg=0;
scattypes=0;

while argnum<=length(varargin)
    switch varargin{argnum}
        case 'event'
            argnum=argnum+1;
            event=varargin{argnum};
        case 'metric'
            argnum=argnum+1;
            targval=varargin{argnum};
        case 'norm'
            normda=1;
        case 'mean'
            meanda=1;
            case 'grpreg'
            grpreg=1;
        case 'scattypes'
            scattypes=1;   
        case 'cond'
            argnum=argnum+1;
            variables=varargin{argnum};
            
    end
    argnum=argnum+1;
end
vargrps={};
for ivar=1:length(variables)
    vargrps{ivar}=[variables{ivar}{1} '' variables{ivar}{2}];
end
plotvar={};
if ~grpreg
for ir=1:length(daregions)
    plotvar{ir}=[dalabels{ir}];
end
plotnums=1:length(daregions);
else
    plotvar{1}=[dalabels{1} ', ' dalabels{2}];
    plotnums=1;
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
    set(axa{ip},'xtick',1:length(vargrps),'xticklabel',vargrps);
    set(axa{ip},'xlim',[0 length(vargrps)+1]);
    set(axa{ip},'xTickLabelRotation',45)
    end
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
        ax=text(axa{1},100+is*60,axpos(4)+50,labelx{is},'units','pixels','fontweight','bold','color',markc(countc,:));
        countc=countc+1;
    else
        %put
        ax=text(axa{1},100+is*60,axpos(4)+50,labelx{is},'units','pixels','fontweight','bold','color',markp(countp,:));
        countp=countp+1;
    end   
end
for is=1:length(sessnums)
    text(axa{1},100+is*75,axpos(4)+100,[num2str(sessnums(is)) ' ' sessmarktxt{is}],'interpreter','tex','units','pixels','fontweight','bold','color',[0 0 0]);
    ranj(is)=rand(1,1)*.8-.5;
end
    
    ratiodas={};
    nc=0;
    np=0;
    nn=0;
    nall=0;
    countspall={};
    countsnall={};
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
    sitecolor=[1 0 0];
    if daregion
        %CN
        pid=1;
        sitecolor=markc(find(contains(cnsites,dasite)),:);
        nc=nc+1;
        nn=nc;
        if grpreg
            nall=nall+1;
                        nn=nall;
        end
    else
        %put
        pid=2;
        sitecolor=markp(find(contains(psites,dasite)),:);
        np=np+1;
        nn=np;
        if grpreg
            pid=1;
            nall=nall+1;
                        nn=nall;
        end
    end  
    danormval=1;
    if normda
        dapertype={};
    for itype=1:length(variables)
        isel=find(contains(origgrps,vargrps{itype}));
        targsi=find([sites.sessnum]==sessnums(ise) & ...
            strcmp({sites.probeid},curda(ida)));
        if ~isempty(targsi)
            ich=sites(targsi(1)).ch;
            datarg=counts{ise}{ich};
            if ~isempty(datarg)
            davals=getfield(datarg,targval);   
            dapertype{itype}=davals{isel};         
            end
        end
    end
    if ~isempty(dapertype)
        danormval=length([dapertype{:}]);
    end
    end
stattable={};
for itype=1:length(variables)    
     isel=find(contains(origgrps,vargrps{itype}));
    targsi=find([sites.sessnum]==sessnums(ise) & ...
        strcmp({sites.probeid},curda(ida)));
    if ~isempty(targsi)
    ich=sites(targsi(1)).ch;
    datarg=counts{ise}{ich};
    if ~isempty(datarg)
        countspos=0;
        countsneg=0;
        if ~isfield(datarg,targval) && strcmp(targval,'ratioposneg')
            davaltemp=getfield(datarg,'response');
            countspos=length(find(davaltemp{isel}>0));
            countsneg=length(find(davaltemp{isel}<0));
        else
            davals=getfield(datarg,targval);    
            countspos=length(davals{isel}{1});
            countsneg=length(davals{isel}{2});
        end
    if countsneg>0
        yvals=countspos/countsneg;         %ratio da counts > 0 / da counts < 0 for specific condition 
    else
        yvals=nan;
    end
    if normda
        yvals=yvals./danormval;
    end
    ratiodas{pid}(itype,nn)=yvals;
    chitest=dg_chi2test2([countspos countsneg]);
    chisig{pid}(itype,nn)=chitest;
    countspall{pid}(itype,nn)=countspos ;
    countsnall{pid}(itype,nn)=countsneg ;
    is=find(strcmp(uniquesites,sites(targsi).site)==1); %idx for unique site name
    end
    end
end
stattable = array2table([countspall{pid}(:,nn)'; countsnall{pid}(:,nn)'],'VariableNames',vargrps,'RowNames',{'pos','neg'});
    %stattable = table(countspall{pid}', countsnall{pid}','VariableNames',vargrps,'RowNames',{'pos','neg'})
if ~scattypes
        alphac=.5;
        lwid=1;  
scatter(axa{pid},(1:length(variables))+ranj(ise),...
            ratiodas{pid}(:,nn),marksize(ise),sessmark{ise},'markeredgecolor',sitecolor,...
            'MarkerEdgeAlpha',alphac,'linewidth',lwid);
end
if scattypes && length(variables)==2
    %scatter plot with 1 type on y-axis and 2nd type on x-axis
    [sig,p,stats] = fishertest(stattable);  %test stat relation between motivational state (eg. current big, previuos no rew & current sma & previous no rew & polarity of da change

    if sig
        alphac=1;
        lwid=2;
    else
        alphac=0.3;
        lwid=.75;
    end
    scatter(axa{pid},ratiodas{pid}(2,nn),...
            ratiodas{pid}(1,nn),marksize(ise),sessmark{ise},'markeredgecolor',sitecolor,...
            'MarkerEdgeAlpha',alphac,'linewidth',lwid);
end
if ~meanda && ~scattypes
aline=plot(axa{pid},(1:length(variables))+ranj(ise),ratiodas{pid}(:,nn),'-','linewidth',1,'color',sitecolor);      
aline.Color(4)=.3;
end
set(findall(axa{pid},'-property','FontSize'),'FontSize',fontsize)

end
end
if meanda
for ip=plotnums
    meansda=nanmean(ratiodas{ip},2);
    stdsda=nanstd(ratiodas{ip},[],2);
    cisda=stdsda./sqrt(size(ratiodas{ip},2)).*1.96;
    cilines=meansda-cisda ;
    cilinespos=meansda+cisda;
    for ivar=1:length(variables)
    aci=plot(axa{ip},[ivar ivar],[cilines(ivar) cilinespos(ivar)],'-','linewidth',2,'color',[0 0 0]);  
    aci.Color(4)=.75;   
    end
   % aline=plot(axa{ip},(1:length(variables)),meansda,'-','linewidth',1,'color',[0 0 0]);  
    %aline.Color(4)=.7;

end
end
savename=[savepath event '_' targval];
if ~scattypes
    ylabel(axa{1},targval);
else
    for ip=plotnums
    ylabel(axa{ip},[targval ' '  vargrps{1} ]);
    xlabel(axa{ip},[targval ' ' vargrps{2} ]);
    maxlims=max(nanmean(ratiodas{ip},2))+3*max(nanstd(ratiodas{ip},[],2));
    xlim(axa{ip},[0 maxlims]);
    ylim(axa{ip},[0 maxlims]);
    %plot unitary line
        aci=plot(axa{ip},[0 maxlims],[0 maxlims],'--','linewidth',1,'color',[0 0 0]);  
        aci.Color(4)=.25; 
    end
savename=[savename '_' vargrps{1} '_' vargrps{2} '_scat'];
end
savefig(figsess,savename);
saveas(figsess,savename,'jpg')
print(figsess,savename,'-painters','-depsc');

end
