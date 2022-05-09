function plotdasummary2(datm,plotparam,varargin)
%3/25/19 
%plot da targeted values for grouped conditions eg. {{'big','left'},{'small','left','aftersm'}}
%sites defined from excel sheet indicating different sites based on sessnum
%& site label from getsites function
%plot bar scatter plots of properties of da peaks or win during
%event/targval types, values from datm variables
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
mark=[0 0.1 .7; 0.8 0 .5; 0 .7 0;.7 0 0; 0.4 0 .7];
sessmark={'.','+','o','^','s','*','v','p','h','x'};
sessmarktxt={'.','+','o','\Delta','sq','*','\nabla','p','h','x'};
mark=[linspecer(10,'qualitative'); linspecer(10,'qualitative')];
sessmark={'sq','o','^','+','v','*','p','h','x','d','<','>','o','^','v','<','>','+','*','p','h','x','d','<','>'};
sessmarktxt={' sq','\o','\Delta','+','\nabla','*','p','h','\times','d','<','>','.','\Delta','\nabla','<','>','+','*','p','h','x','d','<','>'};
marksize=[200 200 150 200 200 200 200 200 200 200 200 200 50 50 50 50 50 50 50 50 50 50 50 200 200 200 200 200 200 200 200 200 200 200 200 200];
cmapp=[1 .5 0;1 .55 0; 1 .6 0; 1 .65 0; 1 .7 0; 1 .75 0; 1 .8 0; 1 .85 0; 1 .9 0];
cmapc=[1 .5 0;1 .55 0; 1 .6 0; 1 .65 0; 1 .7 0; 1 .75 0; 1 .8 0; 1 .85 0; 1 .9 0];
cmap2=cool;
cmap1=hot;
markc=cmap1(1:2:end,:);
markp=cmap2(1:2:end,:);

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
avgsites=0;
subj='patra';
while argnum<=length(varargin)
    switch varargin{argnum}
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
            
    end
    argnum=argnum+1;
end
targdasites=plotparam.dasites;
sites=getsites(sessnums,targdasites,subj);
uniquesites=unique({sites(1:end).site});
cnsites=uniquesites(contains(uniquesites,'c'));
psites=uniquesites(contains(uniquesites,'p'));
[dapair,lfppair]=getsitepairs(targdasites,subj);

if plotlfp
    targdasites=plotparam.lfpchs;
    sites=getlfpsites(sessnums,targdasites,subj);
    uniquesites=unique({sites(1:end).site});
    cnsites=uniquesites(contains(uniquesites,'c'));
    psites=uniquesites(contains(uniquesites,'p'));
    savepath=fullfile(plotparam.savepath,'summary_lfp_');
   fontsize=11;
   set(figsess,'position',[50,50,1700,600],'color',[1 1 1]);
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
if simple
    ranj=ranj.*.5;
    sessmark=repmat({'o'},1,length(plotparam.sessnums));
    marksize=repmat(100,1,length(plotparam.sessnums));
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
    sitenum=countc+countp+1;
    if avgsites
        %sitenum fixed by unique sites
        sitenum=find(strcmp(uniquesites,dasite));
    end
    if daregion
        %CN
        countc=countc+1;
        if avgsites
        %sitenum fixed by unique sites
        countc=find(strcmp(cnsites,dasite));
        end
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
        if avgsites
        %sitenum fixed by unique sites
        countp=find(strcmp(psites,dasite));
        end
        pid=2;
        sitecolor=markp(find(strcmp(psites,dasite)),:);
        if grpreg
            pid=1;      %group regions together, just separate conditions to separate plots
            sitecolor=markp(find(strcmp(psites,dasite)),:);       %separate color group for p
        else
            sitenum=countp;
        end
    end     
    %get ch # for current sess and probe id
    targsi=find([sites.sessnum]==sessnums(ise) & ...
        strcmp({sites.probeid},curda(ida)));
    %if found ch# continue, otherwise skip loop
    if isempty(targsi)
        emptysites(empcount).sessnum=sessnums(ise);
        emptysites(empcount).da=curda{ida};
        empcount=empcount+1;
    end
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
        end
        datarg={};
        if xflag
            %get xinfo data
            targrow=find((contains({xinfo.siteda},curda(ida)) & ...
                contains({xinfo.sitelfp},lfptarg) & ...
                strcmp({xinfo.event},event)) & ...
                strcmp({xinfo.sessionid},num2str(sessnums(ise))) &...
                contains({xinfo.sessiontype},'big')==1);
            if isempty(targrow)
               % find another lfppair lfptarg
               lfpsitesinsess=unique({xinfo(find(contains({xinfo.sessionid},num2str(sessnums(ise))))).sitelfp});
               lfptarg=intersect({lfppair{daid}},lfpsitesinsess);
               disp(['new lfp targ for sess # ' num2str(sessnums(ise)) ' : ' lfptarg{:}]); 
               targrow=find((contains({xinfo.siteda},curda(ida)) & ...
                contains({xinfo.sitelfp},lfptarg) & ...
                strcmp({xinfo.event},event)) & ...
                strcmp({xinfo.sessionid},num2str(sessnums(ise))) &...
                contains({xinfo.sessiontype},'big')==1);
                if isempty(targrow)
                    %still empty
                    error(['targrow still empty for sess # ' num2str(sessnums(ise)) ' da site : ' curda{ida}]);
                end
            end   
             datarg=getfield(xinfo(targrow),'daall');    %'dapos' or 'daneg' types get below
             if isempty(datarg)
                emptysites(empcount).sessnum=sessnums(ise);
                emptysites(empcount).da=curda{ida};
                empcount=empcount+1;
             else
                    %get lfp site pair name iiregardles of lfp type value, if
                    %exits in xinfo
                    targl=find([lsites.sessnum]==sessnums(ise) & ...
                        strcmp({lsites.probeid},lfptarg));
                    lfpsite=lsites(targl).site;
                    labelsx{pid}(sitenum).da=dasite{:};
                    labelsx{pid}(sitenum).lfp=lfpsite;

             end
        else
            if ~isempty(datm{targses}) && ~ isempty(ich)
            datarg=datm{targses}{1}{ich};           %get datm values
            if isempty(datarg)
                emptysites(empcount).sessnum=sessnums(ise);
                emptysites(empcount).da=curda{ida};
                empcount=empcount+1;
            else
               labelsx{pid}(sitenum).da=dasite{:};
                labelsx{pid}(sitenum).lfp=''; 
            end
            else
                 emptysites(empcount).sessnum=sessnums(ise);
                emptysites(empcount).da=curda{ida};
                empcount=empcount+1;
            end
        end

        if ~isempty(datarg)
        sitecolors(pid).color(sitenum,:)=sitecolor;     
        if ~isempty(labelsx{pid}(sitenum).da) || ~isempty(labelsx{pid}(sitenum).lfp)
        sitecolors(pid).label{sitenum}=[labelsx{pid}(sitenum).da newline labelsx{pid}(sitenum).lfp];
        else
                    sitecolors(pid).label{sitenum}='';
        end
        else
                                sitecolors(pid).label{sitenum}='';
        end
    end
    end
end
%end loop through all da/sess
pairlab={};
for iix=1:length(labelsx)
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
            if ~grpreg
                pid=1;      %group regions together, just separate conditions to separate plots
                ax=text(axa{1},iiy*90-100,axpos(4)+90,labell{is},'units','pixels','fontweight','bold','color',sitecolors(pid).color(sitenum(1),:));
            else
                ax=text(axa{1},iiy*90-100,axpos(4)+90,labell{is},'units','pixels','fontweight','bold','color',sitecolors(pid).color(sitenum(1),:));
            end
        else
            %put
            countp=countp+1;
            if ~grpreg                
                pid=2;
                sitenum2=countp;
                sitenum=find(contains(sitecolors(pid).label,pairlab{iix}{iiy}));
                ax=text(axa{1},iiy*90-100,axpos(4)+40,labell{is},'units','pixels','fontweight','bold','color',sitecolors(pid).color(sitenum(1),:));
            else
                sitenum=find(contains(sitecolors(pid).label,[pairlab{iix}{iiy} ]));
                ax=text(axa{1},(iiy-countc)*90-100,axpos(4)+40,labell{is},'units','pixels','fontweight','bold','color',sitecolors(pid).color(sitenum(1),:));
            end
        end 
        end
    end
end


avgda={};
yvals={};
nn=0;
sitenum=0;
countc=0;
countp=0;
cursessids=[];
curda=[];
dasites=[];
for ise=1:length(sessnums)
targses=find(strcmp({trialgrps.sessid},num2str(sessnums(ise))));
trialinfo=trialgrps(targses).trialinfo;
avgy=nan(1,2);
ciy=nan(1,2);
cursessids=find([sites.sessnum]==sessnums(ise));
curda={sites(cursessids).probeid};
if sessnums(ise)==127 && plotlfp
    %remove cl1-cl5 for this session
    curda=curda(~contains(curda,'cl1-cl5'));
end
dasites={sites(cursessids).site};
if (~xflag && ~isempty(datm{targses})) || xflag
for ida=1:length(curda)
%each da separate subplot  
lfpsites={};
daregion=contains(dasites(ida),'c');     %1=='c'
dasite=dasites(ida);
sitecolor=[0 0 1];
sitenum=countc+countp+1;
if avgsites
    sitenum=find(strcmp(uniquesites,dasite));
end
if daregion
    %CN
    countc=countc+1;
    if avgsites
        countc=find(strcmp(cnsites,dasite));
    end
    pid=1;
    if grpreg
    else
        sitenum=countc;
    end
else
    %put
    countp=countp+1;
    if avgsites
        countp=find(strcmp(psites,dasite));
    end
    pid=2;
    if grpreg
        pid=1;      %group regions together, just separate conditions to separate plots
    else
        sitenum=countp;
    end
end     
    
%get ch # for current sess and probe id
targsi=find([sites.sessnum]==sessnums(ise) & ...
    strcmp({sites.probeid},curda(ida)));
%if found ch# continue, otherwise skip loop
if ~isempty(targsi)
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
        lfptarg=lsites(lfpidx(1)).probeid;          %get general label
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
if ~avgsites
    %initialize each sess, as different avg  even if same
yvals{pid}{sitenum}{icond}=[];
end
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
        %{
        if exclusive && icond>1
            %make sure current trial nums do not include those from
            %previous conditions, make mutually exclusive
            trialnumtemp=intersect(trialnumtemp,goodtrials);
            trialnumtemp=trialnumtemp(find(~ismember(trialnumtemp,prevtrials)));
            prevtrials=[prevtrials trialnumtemp];            
        end
        if exclusive && icond==1
            prevtrials=trialnumtemp;
        end
        %}
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
        davals=getfield(datarg,targval);    %all da values for good trials for sess/type
        if contains(targval,'ts') || contains(targval,'delt')
            alnts=getfield(datarg,'mididx')/ratefscv;        %get alignment idx (all signals shifted relative to this event and this represent that idx)
            davals=davals./ratefscv;
               if contains(targval,'damax')
                %remove any below alnts
                    damaxts=getfield(datarg,'damaxts');
                    badtrls=find(damaxts<alnts.*ratefscv);
                    trialids=trialids(~ismember(trialids,badtrls));
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
    yvals{pid}{sitenum}{icond}=[yvals{pid}{sitenum}{icond} davals(trialids)];         %da values for specific group of conditions (ie. trial nums) or site (ie. if avgsites)
end
end
%finish getting group values
dapercond{icond}=yvals{pid}{sitenum}{icond};
if plotv
    dapercond{icond}=var(dapercond{icond},'omitnan');
end
%get different condition values
end
if ~isempty(dapercond{1})
 %sitenumc=find(contains(sitecolors(pid).label,[dasite{:} newline lfpsite]));
sitecolor=sitecolors(pid).color(sitenum,:);
if simple
    sitecolor=[0 0 0];
end
%normalize values
danormval=1;
if normda
    danormval=nanmedian([dapercond{:}]);        %median of all values for all plot conditions for current probe
    for icond=1:length(ttypes)        
        yvals{pid}{sitenum}{icond}=yvals{pid}{sitenum}{icond}./danormval;
    end
end
%plot probe values if not plotting mean over all or not scatter
for icond=1:length(ttypes)
    curplotval=yvals{pid}{sitenum}{icond};
    ciy(icond)=nanstd(curplotval)./sqrt(length(curplotval(~isnan(curplotval))))*1.96;
    avgy(icond)=nanmean(curplotval);
    avgda{pid}(icond,sitenum)=avgy(icond);
    is=find(strcmp(uniquesites,sites(targsi).site)==1); %idx for unique site name
    if ~meanda && ~scattypes
        %plot ci line for each da value
    aci=plot(axa{pid},[icond icond]+ranj(ise),[avgy(icond)-ciy(icond)  avgy(icond)+ciy(icond)],'-','linewidth',1,'color',sitecolor);  
    aci.Color(4)=.5;
    end
end
%get signifiance for difference between first pair of coniditon groups
sig=0;
if length(ttypes)==2
[sig,pval]=ttest2(dapercond{1},dapercond{2});
end
if isnan(sig)
    sig=0;
end
if ~scattypes
    %plot as bar graph with types on x-axis
    if sig
        alphac=1;
        lwid=2;
        if simple
            alphac=.35;
        end   
    else
        alphac=0.3;
        lwid=.75;
    end
    if simple
        if sig
           scatter(axa{pid},(1:length(ttypes))+ranj(ise),...
        avgy,marksize(ise),sessmark{ise},'markeredgecolor',sitecolor,...
            'MarkerEdgeAlpha',alphac,'linewidth',lwid,'markerfacecolor',[0 0 0],'markerfacealpha',.25);
        else
           scatter(axa{pid},(1:length(ttypes))+ranj(ise),...
        avgy,marksize(ise),sessmark{ise},'markeredgecolor',sitecolor,...
        'MarkerEdgeAlpha',alphac,'linewidth',lwid);            
        end
    else
        scatter(axa{pid},(1:length(ttypes))+ranj(ise),...
        avgy,marksize(ise),sessmark{ise},'markeredgecolor',sitecolor,...
        'MarkerEdgeAlpha',alphac,'linewidth',lwid);
    end   
end
if scattypes && length(ttypes)==2
    %scatter plot only if length ttypes ==2
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
%plot line between types for each da value pairs to show trend
aline=plot(axa{pid},(1:length(ttypes))+ranj(ise),avgy,'-','linewidth',1,'color',sitecolor);      
aline.Color(4)=.3;
end
if plotlines
%plot line between types for each da value pairs to show trend
    if sig
    aline=plot(axa{pid},(1:length(ttypes))+ranj(ise),avgy,'-','linewidth',1,'color',sitecolor);      
    aline.Color(4)=.15;
    end
end
end
set(findall(axa{pid},'-property','FontSize'),'FontSize',fontsize)
end
end
end
end
%finish getting all values for yvals{sitenum}{icond} & avgda{pid}(icond,ida)
if normda
    yaxlabel= [yaxlabel ' norm'];
end
if meanda && ~scattypes
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
    %plot bar graph
        ylimax=get(axa{ip},'ylim');
    ax2 = axes('Position',get(axa{ip},'Position'),'XAxisLocation','top',...
    'YAxisLocation','right','Color','none','XColor','k','YColor','k');
    hold(ax2, 'all');  %   <--------------------------------
    barsimp=bar(axa{ip},(1:length(ttypes)),meansda,'FaceAlpha',0.15,'FaceColor',[0 0 0],'EdgeColor','none');
            set(axa{ip},'ylim',ylimax)
    end

end
end
ylabel(axa{1},yaxlabel,'interpreter','none');
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

savename=[savepath event '_' targval '_' condfull];
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

savefig(figsess,savename);
saveas(figsess,savename,'jpg')
print(figsess,savename,'-painters','-depsc');

end
