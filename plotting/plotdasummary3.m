function groupvals=plotdasummary2(datm,plotparam,varargin)
%3/25/19 
%plot da targeted values for grouped conditions eg. {{'big','left'},{'small','left','aftersm'}}
%sites defined from excel sheet indicating different sites based on sessnum
%& site label from getsites function
%plot bar scatter plots of properties of da peaks or win during
%event/targval types, values from datm variables
%11/15/2019, change mintrials from <3 to <=5
figpos=[50,50,1500,600];
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
%mark=[linspecer(10,'qualitative'); linspecer(10,'qualitative')];
sessmark={'sq','o','^','+','v','*','p','h','x','d','<','>','o','^','v','sq','o','^','+','v','*','p','h','x','d','<','>','o','^','v','sq','o','^','+','v','*','p','h','x','d','<','>','o','^','v'};
marksize=[200 200 200 200 200 200 200 200 200 200 200 200 50 50 50 50 50 50 50 50 50 50 50 300 300 300 300 300 300 300 300 300 300 300 300 300 200 200 200 200 200 200 200 200 200 200 200 200 200];
cmapp=[1 .5 0;1 .55 0; 1 .6 0; 1 .65 0; 1 .7 0; 1 .75 0; 1 .8 0; 1 .85 0; 1 .9 0];
cmapc=[1 .5 0;1 .55 0; 1 .6 0; 1 .65 0; 1 .7 0; 1 .75 0; 1 .8 0; 1 .85 0; 1 .9 0];
cmap2=cool;
cmap1=hot;
markc=cmap1(1:5:end,:);
markp=cmap2(1:5:end,:);
%markc=[linspecer(10,'qualitative'); linspecer(10,'qualitative')];
%markp=brewermap(13,'RdBu');
%markp=[brewermap(8,'Accent'); linspecer(10,'qualitative'); brewermap(8,'RdBu'); brewermap(8,'Pastel2');];
%markc=[brewermap(8,'Accent'); linspecer(10,'qualitative'); brewermap(8,'RdBu'); brewermap(8,'Pastel2');];
markc=[cmap1(1:8:end-10,:); cmap2(1:10:end,:);  linspecer(10,'qualitative'); brewermap(8,'RdBu');];
markp=[cmap1(1:8:end-10,:); cmap2(1:10:end,:); linspecer(10,'qualitative'); brewermap(8,'RdBu'); brewermap(8,'Pastel2');];
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
ratefscv=10;
plotvar={};
event='targ';

yaxlabel='\Delta[DA] (nM) ';
normda=0;
mintrials=6;
meanda=0;
grpreg=0;
scattypes=0;
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
simple=0;
plotlines=0;
plotsens=0; 
avgsites=0;
subj='patra';
gain=0;
offset=[0.2 -0.2];
xwin=[0 4];
ovrdays=0;
ovrsess=0;
plotsensors=0;
minmaxnorm=0;
cleoonly=0;
patraonly=0;
numpairs=0;
constrainpairs=0;
assignedsites={};
uiponly=0;
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'minmaxnorm'
            minmaxnorm=1;
        case 'sens'
            plotsens=1;         %ie. cond1/cond2 for cn vs put
        case 'avgsites'
            %average per site across sessions
            avgsites=1;
        case 'xinfo'
            xflag=1;
            argnum=argnum+1;
            xinfo=varargin{argnum};
            sessiontypes=unique({xinfo(1:end).sessiontype});
        case 'xcov'
            %also xinfo flag, but more complex retrieval
            xflag=1;
            xcflag=1;
            argnum=argnum+1;
            xinfo=varargin{argnum};
            sessiontypes=unique({xinfo(1:end).sessiontype});    
            %targval options are minlag, maxlag, mincoef, maxcoef
        case 'xwin'
            argnum=argnum+1;
            xwin=varargin{argnum};
        case 'datype'
            argnum=argnum+1;
            datype=varargin{argnum};    %'postrials' or 'negtrials'
        case 'event'
            argnum=argnum+1;    
            event=varargin{argnum};         %event period, eg. fix, targ, outcome
        case 'metric'
            argnum=argnum+1;
            targval=varargin{argnum};       %datm variable, eg. targpeak, targwin, etc.
            if contains(targval,'ts')
                yaxlabel='ts';
                tlabel=1;
            end
        case 'norm'
            normda=1;           %normalize to plot overall changes
        case 'mean'
            meanda=1;       %plot mean error bars
        case 'grpreg'
            grpreg=1;       %group regions
        case 'scattypes'
            scattypes=1;        %scatter plot of data
        case 'ttypes'
            argnum=argnum+1;
            ttypes=varargin{argnum};        %condition types grouped eg. {{'big','left'},{'small','left','aftersm'}}
        case 'lfp'
            plotlfp=1;          %betatm given instead of datm
        case 'plotvar'
            plotv=1;      %plot variances
        case 'trtable'
            %supply trial list table with labels of which groups to plot
            %(as generaetd by gettrtable & maketrorg)
            argnum=argnum+1;
            trtable=varargin{argnum};
            argnum=argnum+1;
            ttypes=varargin{argnum}; %eg. {{'big','-3break'},{'big','-2break'},{'big','-1break'}};
        case 'simple'
            %simple bar plot
            simple=1;
        case 'plotlines'
            plotlines=1;
        case 'cleo'
            subj='cleo';
        case 'gain'
            plotsens=1;
            gain=1;
        case 'overdays'
            %plot over days (get date from xcel)
            ovrdays=1;
        case 'oversess'
            %plot over session # (get from xcel)
            ovrsess=1;
        case 'cleopatra'
            cleopatra=1;        %both subjects;
        case 'sensors'
            plotsensors=1;
       case 'uiponly'
            uiponly=1;
        case 'patra'
            patraonly=1;
            sessnums=sessnums(sessnums>35);
        case 'cleo'
            sessnums=sessnums(sessnums<35);
            cleoonly=1;
        case 'numpairs'
            %limit # lfp pairs per da site
            argnum=argnum+1;
            numpairs=varargin{argnum};
        case 'constrainpairs'
            constrainpairs=1;
            argnum=argnum+1;
            assignedsites=varargin{argnum};
            
    end
    argnum=argnum+1;
end
targdasites=plotparam.dasites;
sites=getsites(sessnums,targdasites,subj);
uniquesites=unique({sites(1:end).site});
cnsites=uniquesites(contains(uniquesites,'c'));
psites=uniquesites(contains(uniquesites,'p'));
if numpairs==0
    [dapair,lfppair]=getsitepairs(targdasites,subj);
else
    %not zero, num pair limit defined
    [dapair,lfppair]=getsitepairs(targdasites,subj,'numpairs',numpairs);
end
    
if length(sessnums)>23
    fontsize=11;
end
if plotlfp && strcmp(subj,'patra')
    targdasites=plotparam.lfpchs;
    sites=getlfpsites(sessnums,targdasites,subj);
    uniquesites=unique({sites(1:end).site});
    cnsites=uniquesites(contains(uniquesites,'c'));
    psites=uniquesites(contains(uniquesites,'p'));
    savepath=fullfile(plotparam.savepath,'summary_lfp_');
   fontsize=11;
   set(figsess,'position',[50,50,1700,600],'color',[1 1 1]);
end
if plotlfp && strcmp(subj,'patra') && constrainpairs
    %constrain lfp pair for each da site recorded loop through each session
    %to get optimal lfp pair for each da site, already called assignpairs
    %in 'plotmultiple.m' and suppleid here
    sites=assignedsites;
    uniquesites=unique({sites(1:end).site});
    cnsites=uniquesites(contains(uniquesites,'c'));
    psites=uniquesites(contains(uniquesites,'p'));
    savepath=fullfile(plotparam.savepath,'summary_lfp_');
   fontsize=11;
   set(figsess,'position',[50,50,1700,600],'color',[1 1 1]);
end
if uiponly
    %only uip sensors
    uips=find(contains({sites.type},'uip'));
    sites=sites(uips);
end
plotnums=1:length(daregions);
targlsites={};
lsites={};
%get corresponding lfp paired site for each da site for xinfo data
if xflag
    targlsites=plotparam.lfpchs;
    lsites=getlfpsites(sessnums,targlsites);
end
if grpreg
    %group regions
    plotnums=1;
end
condlabels={};
condfull=[];
numsforlab=3;
if ~isempty(trtable)
    numsforlab=5;
end
for ic=1:length(ttypes)
    for iic=1:length(ttypes{ic})
        if iic==1
            if length(ttypes{ic}{1})<numsforlab
                condlabels{ic}=ttypes{ic}{1};
            else
                condlabels{ic}=ttypes{ic}{1}(1:numsforlab);
            end
        else
            if length(ttypes{ic}{iic})<numsforlab
                condlabels{ic}=[condlabels{ic} '_' ttypes{ic}{iic}];
            else
                condlabels{ic}=[condlabels{ic} '_' ttypes{ic}{iic}(1:numsforlab)];
            end                
        end
    end
    condfull=[condfull condlabels{ic} '_vs_'];
end
for ic=1:length(ttypes)
    if ~grpreg
        %separate plots  for different conditions & regions
        for ir=1:length(daregions)
            plotvar{ir}=[dalabels{ir}];
          %  plotvar{ir+length(ttypes)*(ic-1)}=[dalabels{ir} '_' condlabels{ic}];
        end
    else  
        %separate plots for different conditions
        plotvar{ic}=[dalabels{:}];
    end
end
%plot da differences (max, mean) targ window big vs small for 
%CN Left, CN right, P left, P right
clf(figsess,'reset');
set(figsess,'color',[1 1 1]);
axpos={};
if ~ovrdays || ~ovrsess
for ip=plotnums
    axa{ip}=subplot(1,length(plotnums),ip);   hold(axa{ip},'on');
    set(axa{ip},'units','pixels');
    axpos{ip}=get(axa{ip},'position');
    if ~scattypes
   % set(axa{ip},'xtick',1:length(condlabels),'xticklabel',condlabels,'ticklabelinterpreter','none','xticklabelrotation',45);
        set(axa{ip},'xtick',1:length(condlabels),'xticklabel',condlabels,'ticklabelinterpreter','none');

    set(axa{ip},'xlim',[0 length(condlabels)+1]);
    end
    %set(axa{ip},'xTickLabelRotation',90)
    set(axa{ip},'position',[(ip-1)*axsiz(1)+ip*(mar/2)+off,125,axsiz(1),axsiz(2)])
    ax=text(axa{ip},10,axsiz(2),plotvar{ip},'units','pixels','fontweight','bold');
end
end
if ovrdays || ovrsess
    for ip=plotnums
    axa{ip}=subplot(1,length(plotnums),ip);   hold(axa{ip},'on');
    set(axa{ip},'units','pixels');
    axpos{ip}=get(axa{ip},'position');
    %set(axa{ip},'xTickLabelRotation',90)
    set(axa{ip},'position',[(ip-1)*axsiz(1)+ip*(mar/2)+off,125,axsiz(1),axsiz(2)])
    ax=text(axa{ip},10,axsiz(2),plotvar{ip},'units','pixels','fontweight','bold');
    end
end
if plotsens
    %cn vs put
for ip=1
    axa{ip}=subplot(1,1,ip);   hold(axa{ip},'on');
    set(axa{ip},'units','pixels');
    axpos{ip}=get(axa{ip},'position');
    if ~scattypes
   % set(axa{ip},'xtick',1:length(condlabels),'xticklabel',condlabels,'ticklabelinterpreter','none','xticklabelrotation',45);
        set(axa{ip},'xtick',1:2,'xticklabel',{'cn','put'},'ticklabelinterpreter','none');

    set(axa{ip},'xlim',[0 3]);
    end
    %set(axa{ip},'xTickLabelRotation',90)
    set(axa{ip},'position',[(ip-1)*axsiz(1)+ip*(mar/2)+off,125,axsiz(1),axsiz(2)])
end
end
titletext=[event ' | ' targval];  
set(axa{1},'units','pixels')
axpos=get(axa{1},'position');
text(axa{1},-100,axpos(4)+160,titletext,'units','pixels','fontweight','bold','fontsize',fontsize,'interpreter','none');
labelx={};
ranj=[];        %random jitter x axis for each site
%Session legend label
if ~avgsites
axsesslegend=axes;   hold(axsesslegend,'on');
%  set(axsesslegend,'units','pixels');
legpos=get(axsesslegend,'position');
set(axsesslegend,'position',[.1,.9,.85,.09])
ax=text(axa{ip},10,axsiz(2),plotvar{ip},'units','pixels','fontweight','bold');
set(axsesslegend,'box','off','visible','off')
text(axsesslegend,-.4,.5,'sess # :','interpreter','tex','fontweight','bold','color',[0 0 0],'fontsize',fontsize);
for is=1:length(sessnums)
   % text(axa{1},is*75,axpos(4)+150,[num2str(sessnums(is)) ' ' sessmarktxt{is}],'interpreter','tex','units','pixels','fontweight','bold','color',[0 0 0]);
    scatter(axsesslegend,is,.5,marksize(is),sessmark{is},'markeredgecolor',[0 0 0]);
    text(axsesslegend,is-.5,.5,num2str(sessnums(is)),'interpreter','tex','fontweight','bold','color',[0 0 0],'fontsize',fontsize);
    ranj(is)=rand(1,1)*.8-.5;
end
end
if simple && ~scattypes
    ranj=ranj.*.8+.10;
    sessmark=repmat({'o'},1,length(sessnums));
    marksize=repmat(100,1,length(sessnums));
end
if avgsites
    for is=1:length(uniquesites)
    ranj(is)=rand(1,1)*.8-.5;
    end
end
countavgs={};
pairlab={};
sitenum=0;
countc=0;
countp=0;
cursessids=[];
curda=[];
dasites=[];
labelsx={};         %xinfo legend labels for da/lfp pairs
sitecolors={};
empcount=1;
emptysites={};
%{
%get labels for valid site nums/pairs (ie if values exist for site)
%get lfp pair labels, may be more than just da sites since different
%pairs could exist per session
for ise=1:length(sessnums)
    targses=find(strcmp({trialgrps.sessid},num2str(sessnums(ise))));
    trialinfo=trialgrps(targses).trialinfo;
    cursessids=find([sites.sessnum]==sessnums(ise));
    curda={sites(cursessids).probeid};
    dasites={sites(cursessids).site};
    for ida=1:length(curda)
    lfpsites={};
    daregion=contains(dasites(ida),'c');     %1=='c'
    dasite=dasites(ida);
    lfpsite='';
       flagempty=0;
    %get ch # for current sess and probe id
    targsi=find([sites.sessnum]==sessnums(ise) & ...
        strcmp({sites.probeid},curda(ida)));
    %if found ch# continue, otherwise skip loop
    if ~isempty(targsi)
        ich=sites(targsi(1)).ch;        %ch # defined by fscv recording ch
        if plotlfp
            %only if betatm given instead of datm (not for xinfo data)
            %Channel # as defined in betatm (here labeled as datm variable)
            ich=[];
            if ~isempty(datm{targses})
            for iich=1:length(datm{targses}{1})
                if strcmp(datm{targses}{1}{iich}.site,curda{ida})
                    ich=iich;
                end
            end
            end
        end
        %get lfp pair for da site, for xinfo data
        daid=find(strcmp(dapair,curda(ida)));
        lfptarg='';
        if ~plotlfp && ~isempty(daid)
            %not getting betatm values             
            lfptarg=lfppair{daid(1)};           %JUST FIRST GOOD LFP PAIR
            if sessnums(ise)==127 && strcmp(curda(ida),'cl6')
                lfptarg='cl1-cl4';
            end
        end
        datarg={};
        if xflag
            %get xinfo data
            targrow=find((contains({xinfo.siteda},curda(ida)) & ...
                contains({xinfo.sitelfp},lfptarg) & ...
                strcmp({xinfo.event},event)) & ...
                strcmp({xinfo.sessionid},num2str(sessnums(ise))) &...
                contains({xinfo.sessiontype},'big')==1);
            if plotlfp && contains(targval,'xx')
                targrow=find((contains({xinfo.sitelfp},curda(ida)) & ...
                strcmp({xinfo.event},event)) & ...
                strcmp({xinfo.sessionid},num2str(sessnums(ise))) &...
                contains({xinfo.sessiontype},'big')==1);
                if ~isempty(targrow)
                    targrow=targrow(1);
                end
            end
            if isempty(targrow)
               % find another lfppair lfptarg
               lfpsitesinsess=unique({xinfo(find(contains({xinfo.sessionid},num2str(sessnums(ise))))).sitelfp});
               lfptarg=intersect({lfppair{daid}},lfpsitesinsess);
               if length(lfptarg)>1                  
                    lfptarg=lfptarg(1);                   
               end
               disp(['new lfp targ for sess # ' num2str(sessnums(ise)) ' : ' lfptarg{:}]); 
               targrow=find((contains({xinfo.siteda},curda(ida)) & ...
                contains({xinfo.sitelfp},lfptarg) & ...
                strcmp({xinfo.event},event)) & ...
                strcmp({xinfo.sessionid},num2str(sessnums(ise))) &...
                contains({xinfo.sessiontype},'big')==1);
                if isempty(targrow)
                    %still empty
                    display(['targrow still empty for sess # ' num2str(sessnums(ise)) ' da site : ' curda{ida}]);
                    flagempty=1;
                end
                 targl=find([lsites.sessnum]==sessnums(ise) & ...
                        strcmp({lsites.probeid},lfptarg));
                    if plotlfp && contains(targval,'xx')
                        targl=find([lsites.sessnum]==sessnums(ise) & ...
                        strcmp({lsites.probeid},curda(ida)));
                    end
                    if ~isempty(targl)
                    lfpsite=lsites(targl).site;
                    else
                        lfpsite='xx';
                    end
                    dasite=dasite{:};        
            else   
             datarg=getfield(xinfo(targrow),'daall');    %'dapos' or 'daneg' types get below
             if isempty(datarg)
                flagempty=1;
             else
                    %get lfp site pair name iiregardles of lfp type value, if
                    %exits in xinfo
                    targl=find([lsites.sessnum]==sessnums(ise) & ...
                        strcmp({lsites.probeid},lfptarg));
                    if plotlfp && contains(targval,'xx')
                        targl=find([lsites.sessnum]==sessnums(ise) & ...
                        strcmp({lsites.probeid},curda(ida)));
                    end
                    lfpsite=lsites(targl).site;
                    dasite=dasite{:};            
             end
            end
        else
            %datm or betatm variables
            if ~isempty(datm{targses}) && ~ isempty(ich)
                datarg=datm{targses}{1}{ich};           %get datm values
                if isempty(datarg)
                    flagempty=1;
                else
                   dasite=dasite{:};
                    lfpsite=''; 
                end
            else
                flagempty=1;
            end
        end
        if ~flagempty
            sitenum=countc+countp+1;
            if daregion
                %CN
                countc=countc+1;
                %{
                if avgsites
                %sitenum fixed by unique sites
                countc=find(strcmp(cnsites,dasite));
                end
                %}
                pid=1;
                sitecolor=markc(find(strcmp(cnsites,dasite)),:);
                if grpreg
                    pid=1;      %group regions together, just separate conditions to separate plots
                    sitecolor=markc(find(strcmp(cnsites,dasite)),:);       %separate color group for p
                else
                    sitenum=countc;
                end
            else
                %put
                countp=countp+1;
                %{
                if avgsites
                %sitenum fixed by unique sites
                countp=find(strcmp(psites,dasite));
                end
                %}
                pid=2;
                sitecolor=markp(find(strcmp(psites,dasite)),:);
                if grpreg
                    pid=1;      %group regions together, just separate conditions to separate plots
                    sitecolor=markp(find(strcmp(psites,dasite)),:);       %separate color group for p
                else
                    sitenum=countp;
                end
            end         
            sitecolors(pid).color(sitenum,:)=sitecolor;   
            if iscell(dasite)
                dasite=dasite{:};
            end
             labelsx{pid}(sitenum).da=dasite;
             labelsx{pid}(sitenum).lfp=lfpsite;
            if ~isempty(labelsx{pid}(sitenum).da) || ~isempty(labelsx{pid}(sitenum).lfp)
                sitecolors(pid).label{sitenum}=[labelsx{pid}(sitenum).da newline labelsx{pid}(sitenum).lfp];
            else
                sitecolors(pid).label{sitenum}='';
            end
        end
    end
    end
end
%end loop through all da/sess


pairlab={};
for iix=1:length(labelsx)
    if ~isempty(labelsx{iix})
    for iiy=1:length(labelsx{iix})
        if ~isempty(labelsx{iix}(iiy).da)
            if ~isempty(labelsx{iix}(iiy).lfp)
                pairlab{iix}{iiy}=[labelsx{iix}(iiy).da newline labelsx{iix}(iiy).lfp];
            else
                pairlab{iix}{iiy}=labelsx{iix}(iiy).da;
            end
        else
            pairlab{iix}{iiy}='';
        end
    end
    pairlab{iix}=unique(pairlab{iix});
    end
end
        
%make legend with lfp pairs below da sites
countc=0;
countp=0;
sitenum2=0;
axpos=get(axa{1},'position');
ax=text(axa{1},-100,axpos(4)+80,'site id : ','units','pixels','fontweight','bold','color',[0 0 0],'fontsize',fontsize,'interpreter','none');
for iix=1:length(pairlab)
    for iiy=1:length(pairlab{iix})
        is=iiy;
        if iix>1
            is=(iix-1)*length(pairlab{iix-1})+iiy;
        end
        labell{is}=pairlab{iix}{iiy};  
        axpos=get(axa{1},'position');
        daregion=contains(labell{is},'c');     %1=='c'
        pid=1;
        sitenum2=countc+countp+1;
        if ~isempty(labell{is})
        if daregion
            %CN
            countc=countc+1;          
            sitenum=find(contains(sitecolors(pid).label,pairlab{iix}{iiy}));
            if ~isempty(sitenum)
            if ~grpreg
                pid=1;      %group regions together, just separate conditions to separate plots
                ax=text(axa{1},iiy*90-50,axpos(4)+90,labell{is},'units','pixels','fontweight','bold','color',sitecolors(pid).color(sitenum(1),:));
            else
                ax=text(axa{1},iiy*90-50,axpos(4)+90,labell{is},'units','pixels','fontweight','bold','color',sitecolors(pid).color(sitenum(1),:));
            end
            end
        else
            %put
            countp=countp+1;
            if ~grpreg                
                pid=2;
                sitenum2=countp;
                sitenum=find(contains(sitecolors(pid).label,pairlab{iix}{iiy}));
                if ~isempty(sitenum)
                ax=text(axa{1},iiy*90-50,axpos(4)+40,labell{is},'units','pixels','fontweight','bold','color',sitecolors(pid).color(sitenum(1),:));
                end
            else
                sitenum=find(contains(sitecolors(pid).label,[pairlab{iix}{iiy} ]));
                if ~isempty(sitenum)
                ax=text(axa{1},(iiy-countc)*90-50,axpos(4)+40,labell{is},'units','pixels','fontweight','bold','color',sitecolors(pid).color(sitenum(1),:));
                end
            end
        end 
        end
    end
end
%}

avgda={};
yvals={};
svals={};       %session ids for corresponding yvals
sitevals={};        %organized by unique sites, ie sitevals(sessnum).vals{icond} , sitevals(sessnum).site = siteda , sitevals(sessnum).sessnum = sessnum
nn=0;
sitenum=0;
countc=0;
countp=0;
cursessids=[];
curda=[];
dasites=[];
countnum=1;
for ise=1:length(sessnums)
targses=find(strcmp({trialgrps.sessid},num2str(sessnums(ise))));
trialinfo=trialgrps(targses).trialinfo;
avgy=nan(1,2);
ciy=nan(1,2);
day='';
date='';
sensortype='';
cursessids=find([sites.sessnum]==sessnums(ise));
curda={sites(cursessids).probeid};
if sessnums(ise)==127 && plotlfp
    %remove cl1-cl5 for this session
    curda=curda(~contains(curda,'cl1-cl5'));
end
dasites={sites(cursessids).site};
if ((~xflag && ~isempty(datm{targses})) || xflag) && ~isempty(curda)
if ~plotlfp
    day=sites(cursessids(1)).day;
    date=sites(cursessids(1)).date{:};
end
for ida=1:length(curda)
%each da separate subplot  
lfpsites={};
daregion=contains(dasites(ida),'c');     %1=='c'
dasite=dasites(ida);
sitecolor=[0 0 1];
sitenum=countc+countp+1;
siteid=find(strcmp(uniquesites,dasite));
emptyflag=0;
    
    
%get ch # for current sess and probe id
targsi=find([sites.sessnum]==sessnums(ise) & ...
    strcmp({sites.probeid},curda(ida)));
if contains(targval,'xx')
targrow=find((contains({xinfo.sitelfp},curda(ida)) & ...
    strcmp({xinfo.event},event)) & ...
    strcmp({xinfo.sessionid},num2str(sessnums(ise))) &...
    contains({xinfo.sessiontype},'big')==1);   
if isempty(targrow)
    emptyflag=1;
end
end
%if found ch# continue, otherwise skip loop
if ~isempty(targsi) && ~emptyflag
    if ~plotlfp
    day=sites(targsi).day;
    sensortype=sites(targsi).type;

    end
    if daregion
    %CN
    countc=countc+1;
    pid=1;
    if grpreg
    else
        sitenum=countc;
    end
    else
        %put
        countp=countp+1;
        pid=2;
        if grpreg
            pid=1;      %group regions together, just separate conditions to separate plots
        else
            sitenum=countp;
        end
    end 
    ich=sites(targsi(1)).ch;        %ch # defined by fscv recording ch
    lfpsite='';
    if plotlfp
        ich=[];
        %only if betatm given instead of datm (not for xinfo data)
        %Channel # as defined in betatm (here labeled as datm variable)
        if ~isempty(datm{targses})
        for iich=1:length(datm{targses}{1})
            if strcmp(datm{targses}{1}{iich}.site,curda{ida})
                ich=iich;
            end
        end
        end
        yaxlabel='beta';
        if ~isempty(ich)
                lfpsite=labelsx{pid}(sitenum).lfp;      %from pre-label analysis above
        end
    else
        %get lfp pair for da site, for xinfo data
       % daid=find(strcmp(dapair,curda(ida)));  

        %lfptarg=lfppair{daid(1)};           %JUST FIRST GOOD LFP PAIR
        lfpsite=labelsx{pid}(sitenum).lfp;      %from pre-label analysis above
    end
    lfptarg='';
    if ~isempty(lfpsite)
        lfpidx=find(contains({lsites.site},lfpsite));
        if ~isempty(lfpidx)
            lfptarg=lsites(lfpidx(1)).probeid;          %get general label
        else
            lfpsitesinsess=unique({xinfo(find(contains({xinfo.sessionid},num2str(sessnums(ise))))).sitelfp});
            daid=find(strcmp(dapair,curda(ida)));
               lfptarg=intersect({lfppair{daid}},lfpsitesinsess);
               if length(lfptarg)>1
                   lfptarg=lfptarg(1);
               end
        end
    end
dapercond={};
datarg={};           %get all da targeted values for trial condition groups
for icond=1:length(ttypes)
dvals={};
%get trials corresponding to groups of trial conditions in current group ttypes{icond}
%get trial types (big/small/break) for current group
targtrialtypes=[any(strcmp('big',ttypes{icond}))  ...
    any(strcmp('small',ttypes{icond})) ...
    any(strcmp('targetbreak',ttypes{icond})) ...
    any(strcmp('fixbreak',ttypes{icond})) ]; %get logical array for if want big/small/targ/fix
targtrialtypes=find(targtrialtypes==1);
nottypes=find(~strcmp(ttypes{icond},'big') & ...
    ~strcmp(ttypes{icond},'small') &...
    ~strcmp(ttypes{icond},'targetbreak') &...
   ~strcmp(ttypes{icond},'fixbreak'));
if isempty(sessiontypes)
    sessiontypes={'big','small'};
end
%get trial #'s for current condition, ie. not trial type identifier
%get da values for current group condition
%initialize each sess, as different avg  even if same
yvals{pid}{sitenum}{icond}=[];

lfpsite=[];
for itype=targtrialtypes
    trialnums=[];
    davals=[];
    if xflag && ~isempty(lfptarg)
        %get xinfo data
        targrow=find((contains({xinfo.siteda},curda(ida)) & ...
            contains({xinfo.sitelfp},lfptarg) & ...
            strcmp({xinfo.event},event)) & ...
            strcmp({xinfo.sessionid},num2str(sessnums(ise))) &...
            contains({xinfo.sessiontype},sessiontypes(itype))==1);        
        if contains(targval,'xx')
        targrow=find((contains({xinfo.sitelfp},curda(ida)) & ...
            strcmp({xinfo.event},event)) & ...
            strcmp({xinfo.sessionid},num2str(sessnums(ise))) &...
            contains({xinfo.sessiontype},sessiontypes(itype))==1);   
            targrow=targrow(1);
        end
        datarg=getfield(xinfo(targrow),'daall');    %'dapos' or 'daneg' types get below
        yaxlabel=targval;
    elseif ~xflag
        if ~isempty(ich)
        datarg=datm{targses}{itype}{ich};           %get all da targeted values for current trial type big/small
        end
    end
for tt=nottypes
    if isempty(trtable)
        %get list from trialinfo that has coarse groups defined
        curcond=ttypes{icond}{tt};
        targt=find(contains(trialinfo(itype).names,curcond)==1);
        trialnums=[trialnums trialinfo(itype).nums{targt}];  
    else
        %first get good trials for current trial type big/small/targ
        goodtrialgrp=find(contains(trialinfo(itype).names,'all')==1);
        goodtrials=trialinfo(itype).nums{goodtrialgrp};  
        %get trial list
        listses=find(strcmp({trtable.sessid},num2str(sessnums(ise))));
        listtyp=find(contains({trtable(listses).trorg.type},sessiontypes(itype)));
        trlisttarg=find(contains({trtable(listses).trorg(listtyp).grp.label},ttypes{icond}{tt}));
        trialnumtemp=trtable(listses).trorg(listtyp).grp(trlisttarg).trials;
        trialnumtemp=intersect(trialnumtemp,goodtrials);
        trialnums=[trialnums trialnumtemp];
    end
end
trialnums=unique(trialnums);        %trial nums for current trial type big/small
if ~isempty(datarg)
    trialids=[];
    if xflag
        trialnums=intersect(trialnums,getfield(xinfo(targrow),datype));     %get dapos/daall/neg trials
        trialids=find(ismember(datarg.trials,trialnums)==1); %trial ids for targeted metric
    else
        trialids=find(ismember(datarg.trialnums,trialnums)==1); %trial ids for targeted metric
    end
    if ~xcflag
        %not xcov variable, just get parameters as stored in xinfo or
        %datm/betatm
        davals=[];
        if ~contains(targval,'xx')
            davals=getfield(datarg,targval);    %all da values for good trials for sess/type
        end
        if contains(targval,'xx')
            %recompute lfp win mean based on win provided in xcov parameter
            %rather than using supplied in datm variable
            alnidx=getfield(datarg,'mididx');  
            xwinids=xwin(1)*ratefscv+alnidx+offset(1)*ratefscv:xwin(2)*ratefscv+alnidx+offset(2)*ratefscv;        
            tempdata=[];
            if plotlfp
                tempdata=getfield(datarg,'lfptracesaln');
            else
                tempdata=getfield(datarg,'datracesaln');
            end            
            davals=nanmean(tempdata(:,xwinids),2)';
        end
        if contains(targval,'ts') || contains(targval,'delt')
            alnts=getfield(datarg,'mididx')/ratefscv;        %get alignment idx (all signals shifted relative to this event and this represent that idx)
            davals=davals./ratefscv;
               if contains(targval,'damax')
                %remove any below alnts
                    damaxts=getfield(datarg,'damaxts');
                    damaxdata=getfield(datarg,'damax');     %only positive values counted
                    dapostrls=find(damaxdata>0);
                    badtrls=find(damaxts<alnts.*ratefscv);
                    trialids=trialids(~ismember(trialids,badtrls));
                    trialids=trialids(ismember(trialids,dapostrls));%only positive values counted
               end
               if contains(targval,'lfpmin')
                %remove any below alnts
                    damaxts=getfield(datarg,'lfpmints');
                    badtrls=find(damaxts<alnts.*ratefscv);
                    trialids=trialids(~ismember(trialids,badtrls));
               end
            if ~contains(targval,'delt')
                davals=davals-alnts;
            end
        end    
    else
        davals=getxcovpa(xinfo(targrow).xcovda,xinfo(targrow).tsx,targval);       %get xcov parameters
    end
    yvals{pid}{sitenum}{icond}=[yvals{pid}{sitenum}{icond} davals(trialids)];         %da values for specific group of conditions (ie. trial nums) 
end
end
%finish getting group values
sitevals(countnum).vals{icond}=yvals{pid}{sitenum}{icond};
sitevals(countnum).site=dasite;
sitevals(countnum).sitelfp=lfpsite;
sitevals(countnum).sessnum=sessnums(ise);    
sitevals(countnum).date=date;
sitevals(countnum).day=day;
sitevals(countnum).type=sensortype;
end
%finish getting condition types
if isempty(yvals{pid}{sitenum}{1})
    sitevals(countnum).vals=[];
end
end
countnum=countnum+1;
end
%finish da channels
end
end
%finish sessions

%plotting yvals
empties=cellfun('isempty',{sitevals.vals});
sitevals(empties)=[];
if ovrdays
    %remove sites that were only measured one day
    sitesall=unique([sitevals.site]);
    for is=1:length(sitesall)
        sids=find(contains([sitevals.site],sitesall{is}));
        if length(sids)<=1
            sitevals(sids)=[];
        end
    end  
    if uiponly
        %only uip sensors
        uips=find(contains({sitevals.type},'uip'));
        sitevals=sitevals(uips);
    end
end
numplots=length(sitevals);
allsites=unique([sitevals.site]);
if avgsites
    numplots=length(allsites);
end
ise=1;
countp=1;
countc=1;
sitecount=1;
avgda={};
diffgrps={};
for sitenum=1:numplots
avgy=nan(1,2);
ciy=nan(1,2);
dasite='';
sessnum=0;
dapercond={};
plotvals=[];
pid=1;
day=0;
sensortype='';
if ~avgsites
    dasite=sitevals(sitenum).site;
    sessnum=sitevals(sitenum).sessnum;
    plotvals=sitevals(sitenum).vals;
    day=sitevals(sitenum).day;
    sensortype=sitevals(sitenum).type;
    removeflag=0;
    for icond=1:length(plotvals)
        if length(plotvals{icond})<mintrials
            %not enough samples
            disp(['site ' dasite{:} ' sess ' num2str(sessnum) ' samples < ' num2str(mintrials) ', removing ' ]);
            removeflag=1;
        end
    end
    if removeflag
    for icond=1:length(plotvals)
        plotvals{icond}=nan;
    end
    end
else
    siteids=find(strcmp(allsites,allsites{sitenum}));
    dasite=sitevals(siteids(1)).site;
    sessnum=sitevals(siteids(1)).sessnum;
    day=sitevals(siteids(1)).day;
        sensortype=sitevals(siteids(1)).type;

    for iis=1:length(siteids)
        for icond=1:length(sitevals(1).vals)
            if iis==1
                plotvals{icond}=sitevals(siteids(iis)).vals{icond};
            else
                plotvals{icond}=[plotvals{icond} sitevals(siteids(iis)).vals{icond}];
            end

        end
    end
end
ise=find(sessnums==sessnum);
if contains(dasite,'p') && ~grpreg
    pid=2;
    sitecount=countp;
    countp=countp+1;
elseif ~grpreg
    sitecount=countc;
    countc=countc+1;
end
if ~isempty(plotvals)
for icond=1:length(plotvals)
dapercond{icond}=plotvals{icond};
if plotv
    dapercond{icond}=var(dapercond{icond},'omitnan');
end
end
sitecolor=sitecolors(pid).color(sitecount,:);
sitecolorid=find(contains(sitecolors(pid).label,dasite));
sitecolor=sitecolors(pid).color(sitecolorid(1),:);
if simple && ~scattypes
    sitecolor=[0 0 0];
end
if plotsensors
    sitecolor=[0 0 0];
    if strcmp(sensortype,'uip')
        sitecolor=[1 0 0];
    end
end
%normalize values
danormval=1;
if normda
    danormval=abs(nanmedian([plotvals{:}]));        %median of all values for all plot conditions for current probe
    lfpmeanval=abs(nanmean([plotvals{:}]));
    lfpmaxval=max(abs([plotvals{:}]));
    lfpminval=min(abs([plotvals{:}]));
     lfpstdval=nanstd(abs([plotvals{:}]));
   
    for icond=1:length(ttypes)   
        if ~plotlfp && ~minmaxnorm
            plotvals{icond}=plotvals{icond}./danormval;
        end
        if plotlfp
           % plotvals{icond}=(plotvals{icond}-lfpmeanval)./lfpstdval;
            plotvals{icond}=(plotvals{icond}-lfpmeanval)./(lfpmaxval-lfpminval);
        end
        if minmaxnorm && ~plotlfp
            plotvals{icond}=(plotvals{icond}-lfpmeanval)./(lfpmaxval-lfpminval);
        end
    end
end

%plot probe values if not plotting mean over all or not scatter
for icond=1:length(ttypes)
    curplotval=plotvals{icond};
    ciy(icond)=nanstd(curplotval)./sqrt(length(curplotval(~isnan(curplotval))))*1.96;
    avgy(icond)=nanmean(curplotval);
    avgda{pid}(icond,sitecount)=avgy(icond);
    is=find(strcmp(uniquesites,sites(targsi).site)==1); %idx for unique site name
    if ~meanda && ~scattypes && ~ovrdays
        %plot ci line for each da value
    aci=plot(axa{pid},[icond icond]+ranj(ise),[avgy(icond)-ciy(icond)  avgy(icond)+ciy(icond)],'-','linewidth',1,'color',sitecolor);  
    aci.Color(4)=.5;
    end
end

if ovrdays || ovrsess
    %plot as data vs time
    alphac=0.75;
    lwid=.75;
    alphaf=0;
    edgecolor=sitecolor;

    for icond=1:length(ttypes)
        if ~isnan(avgy(icond))
        scatter(axa{pid},day,...
        avgy(icond),100,'o','markeredgecolor',edgecolor,...
            'MarkerEdgeAlpha',alphac,'linewidth',lwid,'markerfacecolor',[0 0 0],'markerfacealpha',alphaf); 
            aci=plot(axa{pid},[day day],[avgy(icond)-ciy(icond)  avgy(icond)+ciy(icond)],'-','linewidth',1,'color',sitecolor);  
    aci.Color(4)=.5;
        end
    end   
end

if ~ovrdays
    %if not plotting over time
diffy=avgy(1)-avgy(2);
if gain
    diffy=(avgy(1)-avgy(2))/(avgy(1)+avgy(2));
end
diffgrps{pid}(sitecount)=diffy;
%get signifiance for difference between first pair of coniditon groups
sig=zeros(1,floor(length(ttypes)/2));
if mod(length(ttypes),2)==0 && length(ttypes)>=2
    for itt=1:2:length(ttypes)
        [sig((itt+1)/2),pval]=ttest2(dapercond{itt},dapercond{itt+1});
        if isnan(sig((itt+1)/2))
            sig((itt+1)/2)=0;
        end
    end
    %sig for pairs
end
if ~scattypes && ~plotsens
    %plot as bar graph with types on x-axis
    alphac=0.3;
    lwid=.75;
    alphaf=0;
    for icond=1:length(ttypes)
        if icond+1<=length(ttypes) && mod(icond,2)==1
            itt=(icond+1)/2;
            if sig(itt)
                alphac=1;
                lwid=2;
                alphaf=0;
                if simple
                    alphac=.35;
                    alphaf=.25;
                end   
            else
                alphac=0.3;
                lwid=.75;
                alphaf=0;
            end
        end
        if ~isnan(avgy(icond))
        scatter(axa{pid},icond+ranj(ise),...
        avgy(icond),marksize(ise),sessmark{ise},'markeredgecolor',sitecolor,...
            'MarkerEdgeAlpha',alphac,'linewidth',lwid,'markerfacecolor',[0 0 0],'markerfacealpha',alphaf); 
        end
    end   
end
if scattypes && length(ttypes)==2
    %scatter plot only if length ttypes ==2
    %scatter plot with 1 type on y-axis and 2nd type on x-axis
    facecolor=sitecolor;
    alphaf=0.5;
    if sig
        alphac=1;
        lwid=1;     
    else
        alphac=0.7;
        alphaf=.1;
        lwid=1;
        facecolor=[1 1 1];
    end
    if ~any(isnan(avgy))
    if ~simple
        scatter(axa{pid},avgy(2),...
            avgy(1),marksize(ise),sessmark{ise},'markerfacecolor',sitecolor,...
            'MarkerEdgeAlpha',alphac,'MarkerFaceAlpha',alphac,'linewidth',lwid);
    else
        if contains(dasite,'c')
        scatter(axa{pid},avgy(2),...
            avgy(1),100,'o','markerfacecolor',facecolor,'markeredgecolor',sitecolor,...
            'MarkerEdgeAlpha',alphac,'MarkerFaceAlpha',alphaf,'linewidth',lwid);
        else
        scatter(axa{pid},avgy(2),...
            avgy(1),100,'sq','markerfacecolor',facecolor,'markeredgecolor',sitecolor,...
            'MarkerEdgeAlpha',alphac,'MarkerFaceAlpha',alphaf,'linewidth',lwid);     
        end
    end
    end
end
if plotsens && length(ttypes)==2
    %compare cn vs put sensitivity 
    %nomralized values differences cond 1 vs cond2
    plotax=1;
    if contains(dasite,'p')
        plotax=2;
    end
    if sig
        alphac=1;
        lwid=1.5;
        if simple
            alphac=.5;
        end   
    else
        alphac=0.3;
        lwid=.75;
    end
    diffy=avgy(1)-avgy(2);
    if gain
        diffy=(avgy(1)-avgy(2))/(avgy(1)+avgy(2));
    end
    if simple
        [cmap,num,typ]=brewermap(2,'RdBu');
        sitecolor=cmap(1,:);
        if diffy<0
            sitecolor=cmap(2,:);
        end
    end
    if ~any(isnan(avgy))
    if simple
        if sig
           scatter(axa{1},plotax+ranj(ise),...
        abs(diffy),marksize(ise),sessmark{ise},'markeredgecolor',sitecolor,...
            'MarkerEdgeAlpha',alphac,'linewidth',lwid,'markerfacecolor',sitecolor,'markerfacealpha',.25);
        else
           scatter(axa{1},plotax+ranj(ise),...
        abs(diffy),marksize(ise),sessmark{ise},'markeredgecolor',sitecolor,...
        'MarkerEdgeAlpha',alphac,'linewidth',lwid);            
        end
    else
        scatter(axa{1},plotax+ranj(ise),...
        diffy,marksize(ise),sessmark{ise},'markeredgecolor',sitecolor,...
        'MarkerEdgeAlpha',alphac,'linewidth',lwid);
    end  
    end
end
if ~meanda && ~scattypes && ~plotsens
%plot line between types for each da value pairs to show trend
aline=plot(axa{pid},(1:length(ttypes))+ranj(ise),avgy,'-','linewidth',1,'color',sitecolor);      
aline.Color(4)=.3;
end
if plotlines
%plot line between types for each da value pairs to show trend
if ~any(isnan(avgy))
    for itt=1:2:length(ttypes)      
        if sig((itt+1)/2)            
        aline=plot(axa{pid},(itt:itt+1)+ranj(ise),avgy(itt:itt+1),'-','linewidth',1,'color',sitecolor);      
        aline.Color(4)=.15;
        end
    end
end
end
end
end
sitecount=sitecount+1;
set(findall(axa{pid},'-property','FontSize'),'FontSize',fontsize)
end
%finish plotting

if ~ovrdays
if plotsens
    yaxlabel=[yaxlabel ' ' [ttypes{1}{:}]];
    for icond=2:length(ttypes)
        yaxlabel=[yaxlabel ' - ' [ttypes{icond}{:}]];
    end
end
%finish getting all values for yvals{sitenum}{icond} & avgda{pid}(icond,ida)
if normda
    yaxlabel= [yaxlabel ' norm'];
end
if meanda && ~scattypes && ~plotsens
for ip=plotnums
    %mean changes for all values from one type to another type
    meansda=nanmean(avgda{ip},2);
    stdsda=nanstd(avgda{ip},[],2);
    cisda=stdsda./sqrt(size(avgda{ip},2)).*1.96;
    cilines=meansda-cisda ;
    cilinespos=meansda+cisda;
    for iip=1:length(cilines)
        aci{iip}=plot(axa{ip},[iip iip],[cilines(iip) cilinespos(iip)],'-','linewidth',3,'color',[0 0 0]);  
        aci{iip}.Color(4)=.75;
    end
    %aci2=plot(axa{ip},[2 2],[cilines(2) cilinespos(2)],'-','linewidth',2,'color',[0 0 0]);  
    %aci.Color(4)=.75;     aci2.Color(4)=.75;
    if ~simple
        %trend line
        if length(ttypes)==2
        aline=plot(axa{ip},(1:length(ttypes)),meansda,'-','linewidth',1,'color',[0 0 0]);  
        aline.Color(4)=.7;
        end
    end
    if simple
        medianda=nanmedian(reshape(avgda{ip},1,size(avgda{ip},1)*size(avgda{ip},2)));       %plot median horizontal line
        amed=plot(axa{ip},[0 length(cilines)+.5],[medianda medianda],'-','linewidth',1,'color',[0 0 0]);  
        amed.Color(4)=.5;
    %plot bar graph
    %{
        ylimax=get(axa{ip},'ylim');
    ax2 = axes('Position',get(axa{ip},'Position'),'XAxisLocation','top',...
    'YAxisLocation','right','Color','none','XColor','k','YColor','k');
    hold(ax2, 'all');  %   <--------------------------------
    barsimp=bar(axa{ip},(1:length(ttypes)),meansda,'FaceAlpha',0.15,'FaceColor',[0 0 0],'EdgeColor','none');
            set(axa{ip},'ylim',ylimax)
    %}
    end

end
end
if plotsens
    %bar graph sensitiivties cn vs put, absolute values
 %mean changes for all values from one type to another type
 for ip=1:2
     %CN VS PUT ip
    meansda(ip)=nanmean(abs(diffgrps{ip}));
    stdsda(ip)=nanstd(abs(diffgrps{ip}));
    cisda(ip)=stdsda(ip)./sqrt(length(diffgrps{ip})).*1.96;
    cilines(ip)=meansda(ip)-cisda(ip) ;
    cilinespos(ip)=meansda(ip)+cisda(ip);
    aci=plot(axa{1},[ip ip],[cilines(ip) cilinespos(ip)],'-','linewidth',3,'color',[0 0 0]);  
    aci.Color(4)=.75;   
 end
    if simple
        %plot bar graph
        ylimax=get(axa{1},'ylim');
       % ylimax=[-max(abs(ylimax)) max(abs(ylimax))];
        ax2 = axes('Position',get(axa{1},'Position'),'XAxisLocation','top',...
        'YAxisLocation','right','Color','none','XColor','k','YColor','k');
        hold(ax2, 'all');  %   <--------------------------------
        barsimp=bar(axa{1},(1:length(meansda)),meansda,.5,'FaceAlpha',0.05,'FaceColor',[0 0 0],'EdgeColor',[0 0 0]);
                set(axa{1},'ylim',ylimax)
    end    
    
end
if ~scattypes
ylabel(axa{1},yaxlabel);
else
ylabel(axa{1},yaxlabel,'interpreter','none');
end
if scattypes && length(ttypes)==2
    maxlim=max(nanmean([avgda{:}],2))+3*max(nanstd([avgda{:}],[],2));
    minlim=min(nanmean([avgda{:}],2))-3*max(nanstd([avgda{:}],[],2));
    for ip=plotnums
        xlim(axa{ip},[minlim maxlim]);
        ylim(axa{ip},[minlim maxlim]);
        %plot unitary line
        aci=plot(axa{ip},[minlim maxlim],[minlim maxlim],'--','linewidth',1,'color',[0 0 0]);  
        aci.Color(4)=.25; 
        tlab1=ttypes{1}{1};
        if length(ttypes{1})>1
        for it=2:length(ttypes{1})
            tlab1=[tlab1 '/' ttypes{1}{it}];
        end
        end
        tlab2=ttypes{2}{1};
        if length(ttypes{2})>1
        for it=2:length(ttypes{2})
            tlab2=[tlab2 '/' ttypes{2}{it}];
        end
        end
        ylabel(axa{ip},[yaxlabel ' '  tlab1],'interpreter','none');
        xlabel(axa{ip},[yaxlabel ' '  tlab2],'interpreter','none');
    end
end
end
groupvals.sitevals=sitevals;
groupvals.cond=ttypes;
groupvals.metric=targval;
groupvals.event=event;
groupvals.signal='da';
if plotlfp
    groupvals.signal='lfp';
end
savename=[savepath event '_' targval '_' condfull];
if gain
    savename=[savename '_gain'];
end
if plotsens
    savename=[savename '_sens'];
end
if normda
    savename=[savename '_z'];
end
if meanda
    savename=[savename '_avg'];
end
if scattypes
    savename=[savename '_scat'];
end
if xflag
    savename=[savename '_xdat'];
end
if plotv
    savename=[savename '_v'];
end
if grpreg
    savename=[savename '_grp'];
end
if simple
    savename=[savename '_simp'];
end
if avgsites
    savename=[savename '_avgsites'];
end
if ovrdays
    savename=[savename '_overdays'];
end
if minmaxnorm
    savename=[savename '_minmaxnorm'];
end
if patraonly
    savename=[savename '_patra'];
end
if cleoonly
    savename=[savename '_cleo'];
end
if constrainpairs
        savename=[savename '_constrainedlfps'];
end
if uiponly
            savename=[savename '_uipsonly'];
end
savefig(figsess,savename);
saveas(figsess,savename,'jpg')
print(figsess,savename,'-painters','-depsc');

end
