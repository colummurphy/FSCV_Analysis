function groupvals=plotbehsummary(xinfo,binfos,plotparam,varargin)
%04/28/20, same as plotdasumary2 but for behavioral meas, only care about
%session #'s not sites, input xbinfos for xinfo argument
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
mar=140;
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
normall=0;
fontsize=14;
trialgrps=plotparam.trialgrps;
savepath=fullfile(plotparam.savepath,'summary_da_');
sessnums=plotparam.sessnums;
conditions={'left','right'};
condlabels={'L','R'};
behvars={'lick'};
ratefscv=10;
plotvar={};
event='targ';
sessiontypes={'big','small'};
countsig={};
yaxlabel='\Delta[DA] (nM) ';
normda=0;
mintrials=6;
meanda=0;
scattypes=0;
plotlfp=0;
datype='goodtrials';        %all good trials for xinfo data
xflag=0;            %get xinfo data
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
xwin=[0 4];
ovrdays=0;
ovrsess=0;
wins={};

while argnum<=length(varargin)
    switch varargin{argnum}
        case 'beh'
            argnum=argnum+1;
            behvars=varargin{argnum};
        case 'wins'
            argnum=argnum+1;
            wins=varargin{argnum};      %windows for each variable
        case 'minmaxnorm'
            minmaxnorm=1;
        case 'sens'
            plotsens=1;         %ie. cond1/cond2 for cn vs put
        case 'avgsites'
            %average per site across sessions
            avgsites=1;       
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
        case 'scattypes'
            scattypes=1;        %scatter plot of data
        case 'ttypes'
            argnum=argnum+1;
            ttypes=varargin{argnum};        %condition types grouped eg. {{'big','left'},{'small','left','aftersm'}}
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
        case 'gain'
            plotsens=1;
            gain=1;
        case 'cleopatra'
            cleopatra=1;        %both subjects;
        case 'patra'
            patraonly=1;
            sessnums=sessnums(sessnums>35);
        case 'cleo'
            sessnums=sessnums(sessnums<35);
            cleoonly=1;
        case 'normall'
            normall=1;
            
    end
    argnum=argnum+1;
end
targdasites=plotparam.dasites;
if length(sessnums)>23
    fontsize=11;
end
savepath=fullfile(plotparam.savepath,filesep,'beh',filesep);
if ~isdir(savepath)
    mkdir(savepath);
end
savepath=[savepath 'summary_beh_'];
set(figsess,'position',[50,50,1700,600],'color',[1 1 1]);

plotnums=1:length(behvars);
if length(plotnums)>=3
    set(figsess,'position',[50,50,1700,900],'color',[1 1 1]);
end
numsforlab=3;
if ~isempty(trtable)
    numsforlab=5;
end
condfull=[];
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
    if ic~=length(ttypes)
    condfull=[condfull condlabels{ic} '_vs_'];
    else
        condfull=[condfull condlabels{ic}];
    end
end

%plot da differences (max, mean) targ window big vs small for 
%CN Left, CN right, P left, P right
clf(figsess,'reset');
set(figsess,'color',[1 1 1]);
axpos={};
if ~ovrdays || ~ovrsess
for ip=plotnums
    if length(plotnums)<4
        axa{ip}=subplot(1,length(plotnums),ip);   hold(axa{ip},'on');
    else
        axa{ip}=subplot(2,3,ip);   hold(axa{ip},'on');
    end        
    set(axa{ip},'units','pixels');
    axpos{ip}=get(axa{ip},'position');
    if ~scattypes
   % set(axa{ip},'xtick',1:length(condlabels),'xticklabel',condlabels,'ticklabelinterpreter','none','xticklabelrotation',45);
        set(axa{ip},'xtick',1:length(condlabels),'xticklabel',condlabels,'ticklabelinterpreter','none');

    set(axa{ip},'xlim',[0 length(condlabels)+1]);
    end
    %set(axa{ip},'xTickLabelRotation',90)
    if length(plotnums)<4
        set(axa{ip},'position',[(ip-1)*axsiz(1)+ip*(mar/2),125,axsiz(1),axsiz(2)])
    else
        set(axa{ip},'position',[(mod(ip-1,3))*axsiz(1)+(mod(ip-1,3)+1)*(mar/2),475-ceil((ip-3)/3)*(axsiz(2)+90),axsiz(1),axsiz(2)])
    end
    ax=text(axa{ip},10,axsiz(2),behvars{ip},'units','pixels','fontweight','bold');
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
        set(axa{ip},'xtick',1:length(behvars),'xticklabel',behvars,'ticklabelinterpreter','none');

    set(axa{ip},'xlim',[0 length(behvars)]);
    end
    %set(axa{ip},'xTickLabelRotation',90)
    set(axa{ip},'position',[(ip-1)*axsiz(1)+ip*(mar/2)+off,125,axsiz(1),axsiz(2)])
end
end
savelabel=[];
titletext=[];
for ip=plotnums
    if ip==1
        titletext=[event ' | ' behvars{ip}];
        savelabel=[event '_' condfull '_' behvars{ip} '_' num2str(wins{ip}(1)) '-' num2str(wins{ip}(2)) 's'];
    else
        titletext=[titletext ' | ' behvars{ip}];
        savelabel=[savelabel '_' behvars{ip} '_' num2str(wins{ip}(1)) '-' num2str(wins{ip}(2)) 's'];
    end
end
savelabel(strfind(savelabel,'.'))='d';      %replace dots with 'd'
set(axa{1},'units','pixels')
axpos=get(axa{1},'position');
text(axa{1},0,axpos(4)+110,titletext,'units','pixels','fontweight','bold','fontsize',fontsize,'interpreter','none');
labelx={};
ranj=[];        %random jitter x axis for each site
%Session legend label
if ~avgsites
axsesslegend=axes;   hold(axsesslegend,'on');
%  set(axsesslegend,'units','pixels');
legpos=get(axsesslegend,'position');
set(axsesslegend,'position',[.1,.90,.85,.09])
set(axsesslegend,'box','off','visible','off')
text(axsesslegend,-.5,1,'sess # :','interpreter','tex','fontweight','bold','color',[0 0 0],'fontsize',fontsize);
randomnumbers=rand(length(sessnums),1);
for is=1:length(sessnums)
   % text(axa{1},is*75,axpos(4)+150,[num2str(sessnums(is)) ' ' sessmarktxt{is}],'interpreter','tex','units','pixels','fontweight','bold','color',[0 0 0]);
    scatter(axsesslegend,is+.5,1,marksize(is),sessmark{is},'markeredgecolor',[0 0 0]);
    text(axsesslegend,is-.2,1,num2str(sessnums(is)),'interpreter','tex','fontweight','bold','color',[0 0 0],'fontsize',fontsize);
    ranj(is)=randomnumbers(is)*.8-.5;
end
end
if simple && ~scattypes && ~plotsens
    ranj=ranj.*.8+.10;
    sessmark=repmat({'o'},1,length(sessnums));
    marksize=repmat(100,1,length(sessnums));
end
if simple && ~scattypes && plotsens
    ranj=ranj.*.6+.05;
    sessmark=repmat({'o'},1,length(sessnums));
    marksize=repmat(100,1,length(sessnums));
end
if avgsites
    for is=1:length(uniquesites)
    ranj(is)=rand(1,1)*.8-.5;
    end
end
yvals={};
sitevals={};        %organized by unique sites, ie sitevals(sessnum).vals{icond} , sitevals(sessnum).site = siteda , sitevals(sessnum).sessnum = sessnum
curda=[];
countnum=1;
for ise=1:length(sessnums)
targses=find(strcmp({trialgrps.sessid},num2str(sessnums(ise))));
trialinfo=trialgrps(targses).trialinfo;
avgy=nan(1,2);
ciy=nan(1,2);
 
for bi=1:length(behvars)
        yaxlabel=behvars{bi};
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

%get trial #'s for current condition, ie. not trial type identifier
%get da values for current group condition
%initialize each sess, as different avg  even if same
yvals{ise}{bi}{icond}=[];

lfpsite=[];
for itype=targtrialtypes
    trialnums=[];
    davals=[];
    targrow=[];
    %get xbinfo data
    if ~contains(behvars{bi},'rt')
        %not reaction time, use xbinfos
    targrow=find((strcmp({xinfo.sitelfp},behvars{bi}) & ...
        strcmp({xinfo.event},event) & ...
        strcmp({xinfo.sessionid},num2str(sessnums(ise))) &...
        contains({xinfo.sessiontype},sessiontypes(itype)))==1);    
    end
    if length(targrow)>1
        trialcells={xinfo(targrow).goodtrials};
        maxtrialnums=0;
        chotarg=1;
        for it=1:length(trialcells)
            trialnums=length(trialcells{it});
            if trialnums>maxtrialnums
                maxtrialnums=trialnums;
                chotarg=it;
            end
        end
        targrow=targrow(it);
    end
    datarg={};
    if ~contains(behvars{bi},'rt')
        datarg=getfield(xinfo(targrow),'daall');    %'dapos' or 'daneg' types get below
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
    trialids=[];
    davals=[];
    alnidx=300;
    xwinids=[];
    tempdata=[];
    btarg=[];
    if ~contains(behvars{bi},'rt')
        %xbinfos data
        tempdata=getfield(datarg,'lfptracesaln');
        alnidx=getfield(datarg,'mididx');  
        xwinids=wins{bi}(1)*ratefscv+alnidx:wins{bi}(2)*ratefscv+alnidx;        
        trialnums=intersect(trialnums,getfield(xinfo(targrow),datype));     %get dapos/daall/neg trials
        trialids=find(ismember(datarg.trials,trialnums)==1); %trial ids for targeted metric
    else
        %reaction time from binfos
        btarg=find(contains({binfos.sessiontype},sessiontypes{itype}) &...
            strcmp({binfos.sessionid},num2str(sessnums(ise))) &...
            strcmp({binfos.evt},event)); 
        trialsl=binfos(btarg).seltrialsl;
    end
    rate=ratefscv;
    switch behvars{bi}
        case 'lrt'
            %left reaction time
            trialsl=binfos(btarg).seltrialsl;
            trialnums=intersect(trialnums,trialsl);     %get dapos/daall/neg trials
            trialids=find(ismember(trialsl,trialnums)==1); %trial ids for targeted metric
            temp=binfos(btarg).target_lrt(trialids);
            davals=temp.*1000;
        case 'rrt'
            %right reaction time
            trialsr=binfos(btarg).seltrialsr;
            trialnums=intersect(trialnums,trialsr);     %get dapos/daall/neg trials
            trialids=find(ismember(trialsr,trialnums)==1); %trial ids for targeted metric
            temp=binfos(btarg).target_rrt(trialids);
            davals=temp.*1000;            
        case 'eye'
            btarg=find(contains({binfos.sessiontype},sessiontypes{itype}) &...
                strcmp({binfos.sessionid},num2str(sessnums(ise))) &...
                strcmp({binfos.evt},event)); 
            etrialsall=getfield(datarg,'trials');                 
            etrialsr=binfos(btarg).seltrialsr;      %right side only
            targtrials=intersect(etrialsr, trialnums);
            etrialsel=find(ismember(etrialsall,targtrials));
            temp=datarg.lfptracesaln(etrialsel,:);
            eyewin=-temp(:,xwinids); 
            eye=nanmean(eyewin,2);            
            %Z SCORE EYE 072019 from plottracesel
            glitchwidth=2;
            basewin=[alnidx-7*rate:alnidx+4*rate];          %fix baseline window for eye below
            %zscore eye wrt baseline as done in olivier et al 2014 frontier
            %remove blinks
            extwin=[alnidx-3*rate:alnidx+4*rate];
            relwin=intersect(xwinids,extwin)-extwin(1)+1;
            datawin=-temp(:,extwin); %flip so original sign
            datawin=deglitchinterpsimple(datawin,3.5,glitchwidth);
            basedata=-temp(:,basewin);   %normalizing to pre-fix and targ window allows normalization trial to trial to look at faster changes within 1 sec 06/04 without phase trend eg. 67/68/69
            basedata=-basedata;
            for ix=1:size(basedata,1)
                [xx, xb]=rmoutliers(basedata(ix,:));
                basedata(ix,xb)=nan;
            end
            baseline=nanstd(basedata,[],2);          %normalized to each trial so independent of phase changes
            dataz=(datawin-nanmean(basedata,2))./baseline;
            datawin=dataz;
            eyewin=dataz(:,relwin);     %user selected windo
            eye=nanmean(eyewin,2);
            eye(isoutlier(eye))=nan;  
           % if sessnums(ise)==114
            %    eye(1:end)=nan;
           % end
            davals=eye';
        case 'lick'
            licktrials=datarg.trials;
            targtrials=intersect(licktrials, trialnums);
            trialsel=find(ismember(licktrials,targtrials));
            licks=tempdata(trialsel,xwinids); 
            lick=nanmean(licks,2);
            lick(isoutlier(lick))=nan;  
            davals=lick';
        case 'pulse'
            %rrstd
            hrtrials=datarg.trials;
            targtrials=intersect(hrtrials, trialnums);
            trialsel=find(ismember(hrtrials,targtrials));
            hrwin=tempdata(trialsel,xwinids); 
            hr=nanmean(hrwin,2);
           rrtemp=1./tempdata(trialsel,:).*60; 
           rrwin=rrtemp(:,xwinids);
           rr=nanmean(rrwin,2);
            rrstd=nanstd(rrwin,[],2).*1000; %in milliseconds
            %outofrange=find(rrstd<.01);
           % rrstd(outofrange)=nan;
            rmssd=rms(diff(rrtemp(:,xwinids),1,2),2);  
            rrstd(isoutlier(rrstd))=nan;        %remove outliers using median criterai
            if nanmean(rrstd)>30
                %remoe all session too large, something wrong
                disp(['sess ' num2str(sessnums(ise)) ' unstable pulse']);
                rrstd(1:end)=nan;
            end
            davals=rrstd';    
    end
    if contains(behvars{bi},'ts') || contains(behvars{bi},'delt')
        alnts=getfield(datarg,'mididx')/ratefscv;        %get alignment idx (all signals shifted relative to this event and this represent that idx)
        davals=davals./ratefscv;
           if contains(behvars{bi},'damax')
            %remove any below alnts
                damaxts=getfield(datarg,'damaxts');
                damaxdata=getfield(datarg,'damax');     %only positive values counted
                dapostrls=find(damaxdata>0);
                badtrls=find(damaxts<alnts.*ratefscv);
                trialids=trialids(~ismember(trialids,badtrls));
                trialids=trialids(ismember(trialids,dapostrls));%only positive values counted
           end
           if contains(behvars{bi},'lfpmin')
            %remove any below alnts
                damaxts=getfield(datarg,'lfpmints');
                badtrls=find(damaxts<alnts.*ratefscv);
                trialids=trialids(~ismember(trialids,badtrls));
           end
        if ~contains(behvars{bi},'delt')
            davals=davals-alnts;
        end
    end    
    yvals{ise}{bi}{icond}=[yvals{ise}{bi}{icond} davals];         %da values for specific group of conditions (ie. trial nums) 
end
%finish getting group values
sitevals(countnum).vals{icond}=yvals{ise}{bi}{icond};
sitevals(countnum).beh=behvars{bi};
sitevals(countnum).sessnum=sessnums(ise);    
end
%finish getting condition types
if isempty(yvals{ise}{bi}{1})
    sitevals(countnum).vals=[];
end
countnum=countnum+1;
end
end
%finish da channels
%finish sessions

%plotting yvals
empties=cellfun('isempty',{sitevals.vals});
sitevals(empties)=[];

numplots=length(behvars);

ise=1;
countp=1;
countc=1;
sitecount=1;
avgda={};
diffgrps={};
for bi=1:numplots
    dasite=behvars{bi};
    for ise=1:length(sessnums)
        avgy=nan(1,2);
        ciy=nan(1,2);
        dapercond={};
        plotvals=[];
        pid=1;
        sessnum=sessnums(ise);
        plotvals=yvals{ise}{bi};
        removeflag=0;
        for icond=1:length(plotvals)
            if length(plotvals{icond})<mintrials
                %not enough samples
                disp(['site ' dasite ' sess ' num2str(sessnum) ' samples < ' num2str(mintrials) ', removing ' ]);
                removeflag=1;
            end
        end
        if removeflag
        for icond=1:length(plotvals)
            plotvals{icond}=nan;
        end
        end
       
if ~isempty(plotvals)
        for icond=1:length(plotvals)
            dapercond{icond}=plotvals{icond};
            if plotv
                dapercond{icond}=var(dapercond{icond},'omitnan');
            end
        end
        sitecolor=markc(ise,:);      
        if simple && ~scattypes
            sitecolor=[0 0 0];
        end
%normalize plotvals values zscore
danormval=1;
if normda && ~plotsens
    if (~strcmp(behvars{bi},'pulse') && ~contains(behvars{bi},'rt')) || normall
    meanval=nanmean([plotvals{:}]);        %mean of all values for all plot conditions for current prob   
    stdval=nanstd([plotvals{:}]);
    for icond=1:length(ttypes)   
        plotvals{icond}=(plotvals{icond}-meanval)./stdval;    
    end
    end
end

%plot probe values if not plotting mean over all or not scatter
for icond=1:length(ttypes)
    curplotval=plotvals{icond};
    ciy(icond)=nanstd(curplotval)./sqrt(length(curplotval(~isnan(curplotval))))*1.96;
    avgy(icond)=nanmean(curplotval);
    avgda{bi}(icond,ise)=avgy(icond);
    if ~meanda && ~scattypes && ~ovrdays
        %plot ci line for each da value
        aci=plot(axa{bi},[icond icond]+ranj(ise),[avgy(icond)-ciy(icond)  avgy(icond)+ciy(icond)],'-','linewidth',1,'color',sitecolor);  
        aci.Color(4)=.5;
    end
end

%if plotting sens between  2 conditions only
diffy=avgy(1)-avgy(2);
if gain
    diffy=(avgy(1)-avgy(2))/(avgy(1)+avgy(2));
end
diffgrps{bi}(ise)=diffy;
%get signifiance for difference between first pair of coniditon groups
sig=zeros(1,floor(length(ttypes)/2));
if mod(length(ttypes),2)==0 && length(ttypes)>=2
    for itt=1:2:length(ttypes)
        [sig((itt+1)/2),pval]=ttest2(dapercond{itt},dapercond{itt+1});
        if isnan(sig((itt+1)/2))
            sig((itt+1)/2)=0;
        end
         %store significance and non-sig results for each combo of conditions tested
            countsig{bi}{(itt+1)/2}(ise).ttype=[[ttypes{itt}{:} '.'] [ttypes{itt+1}{:}]];      %first condition vs 2nd cond group
            countsig{bi}{(itt+1)/2}(ise).dasite=behvars{bi};
            countsig{bi}{(itt+1)/2}(ise).sessnum=sessnum;
            countsig{bi}{(itt+1)/2}(ise).sig=sig((itt+1)/2);
            countsig{bi}{(itt+1)/2}(ise).p=pval;
            countsig{bi}{(itt+1)/2}(ise).pos=nanmean(dapercond{itt})<nanmean(dapercond{itt+1});     %positive slope from cond 1 to cond 2
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
        scatter(axa{bi},icond+ranj(ise),...
        avgy(icond),marksize(ise),sessmark{ise},'markeredgecolor',sitecolor,...
            'MarkerEdgeAlpha',alphac,'linewidth',lwid,'markerfacecolor',[0 0 0],'markerfacealpha',alphaf); 
        end
    end  
end
if ~meanda && ~scattypes && ~plotsens
%plot line between types for each da value pairs to show trend
aline=plot(axa{bi},(1:length(ttypes))+ranj(ise),avgy,'-','linewidth',1,'color',sitecolor);      
aline.Color(4)=.3;
end
if plotlines
%plot line between types for each da value pairs to show trend
if ~any(isnan(avgy))
    for itt=1:2:length(ttypes)      
        if sig((itt+1)/2)            
        aline=plot(axa{bi},(itt:itt+1)+ranj(ise),avgy(itt:itt+1),'-','linewidth',1,'color',sitecolor);      
        aline.Color(4)=.15;
        end
    end
end
end

end
end
sitecount=sitecount+1;
yaxlabel=dasite;
if plotsens
    yaxlabel=[yaxlabel ' ' [ttypes{1}{:}]];
    for icond=2:length(ttypes)
        yaxlabel=[yaxlabel ' - ' [ttypes{icond}{:}]];
    end
end
%finish getting all values for yvals{sitenum}{icond} & avgda{pid}(icond,ida)
if normda
    yaxlabel= [yaxlabel ' z'];
end
if strcmp(behvars{bi},'pulse') && ~normall
    yaxlabel='rrstd (ms)';
end
if strcmp(behvars{bi},'lrt') && ~normall
    yaxlabel='reaction time, left (ms)';
end
if strcmp(behvars{bi},'rrt') && ~normall
    yaxlabel='reaction time, right (ms)';
end
if ~scattypes
ylabel(axa{bi},yaxlabel);
else
ylabel(axa{bi},yaxlabel,'interpreter','none');
end
set(findall(axa{bi},'-property','FontSize'),'FontSize',fontsize)

end
%finish plotting individual session meas


if meanda && ~scattypes && ~plotsens
for ip=plotnums
    %mean changes for all values from one type to another type, across
    %behaviors
    meansda=nanmean(avgda{ip},2);
    stdsda=nanstd(avgda{ip},[],2);
    cisda=stdsda./sqrt(size(avgda{ip},2)).*1.96;
    cilines=meansda-cisda ;
    cilinespos=meansda+cisda;
    for iip=1:length(cilines)
        aci{iip}=plot(axa{ip},[iip iip],[cilines(iip) cilinespos(iip)],'-','linewidth',3,'color',[0 0 0]);  
        aci{iip}.Color(4)=.75;
        %label number sig sites / total for each tested pair
        if mod(iip,2)
            ylims=get(axa{ip},'ylim');
            totalsig=length(find([countsig{ip}{(iip+1)/2}.sig]>0));
            totalsites=length(countsig{ip}{(iip+1)/2});
            %ax=text(axa{iip},iip,ylims(2),[num2str(totalsig) '/' num2str(totalsites)],'fontweight','bold');
            totalsigpos=length(find([countsig{ip}{(iip+1)/2}.sig]>0 & [countsig{ip}{(iip+1)/2}.pos]==0));
            totalsigneg=length(find([countsig{ip}{(iip+1)/2}.sig]>0 & [countsig{ip}{(iip+1)/2}.pos]==1));
            ax=text(axa{ip},iip,ylims(2)+(ylims(2)-ylims(1))/10,['+' num2str(totalsigpos) '/' num2str(totalsites) sprintf('\n') '-' num2str(totalsigneg) '/' num2str(totalsites)]);
        end
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
end
end

if plotsens
    %bar graph sensitiivties cn vs put, absolute values
 %mean changes for all values from one type to another type
 %/04/2020 rather than total average, separate averae for positive/negative
    for ip=plotnums
    meansdapos(ip)=nanmean((diffgrps{ip}(diffgrps{ip}>0)));
    meansdaneg(ip)=nanmean((diffgrps{ip}(diffgrps{ip}<0)));
        stdsdapos(ip)=nanstd(diffgrps{ip}(diffgrps{ip}>0));
        stdsdaneg(ip)=nanstd(diffgrps{ip}(diffgrps{ip}<0));
    cisdapos(ip)=stdsdapos(ip)./sqrt(length(diffgrps{ip}(diffgrps{ip}>0))).*1.96;
    cisdaneg(ip)=stdsdaneg(ip)./sqrt(length(diffgrps{ip}(diffgrps{ip}<0))).*1.96;
    
    aci=plot(axa{1},[ip ip],[meansdapos(ip)-cisdapos(ip) meansdapos(ip)+cisdapos(ip)],'-','linewidth',3,'color',[0 0 0]);    
        aci.Color(4)=.75;  
        aci=plot(axa{1},[ip ip],[meansdaneg(ip)-cisdaneg(ip) meansdaneg(ip)+cisdaneg(ip)],'-','linewidth',3,'color',[0 0 0]);    
            aci.Color(4)=.75;  
              
        %label number sig sites / total for each tested pair
            ylims=get(axa{ip},'ylim');
            totalsig=length(find([countsig{ip}{1}.sig]>0));
            totalsites=length(countsig{ip}{1});
            %ax=text(axa{iip},iip,ylims(2),[num2str(totalsig) '/' num2str(totalsites)],'fontweight','bold');
            csites=find(contains({countsig{ip}{1}.dasite},'c'));
            psites=find(contains({countsig{ip}{1}.dasite},'p'));
            siteids=csites;
            if ip==2
                siteids=psites;
            end
            totalsigpos=length(find([countsig{ip}{1}(siteids).sig]>0 & [countsig{ip}{1}(siteids).pos]==0));
            totalsigneg=length(find([countsig{ip}{1}(siteids).sig]>0 & [countsig{ip}{1}(siteids).pos]==1));    
            ax=text(axa{ip},ip,ylims(2),['+' num2str(totalsigpos) '/' num2str(length(siteids)) sprintf('\n') '-' num2str(totalsigneg) '/' num2str(length(siteids))]);      
    end
    if simple
        %plot bar graph
        ylimax=get(axa{1},'ylim');
       % ylimax=[-max(abs(ylimax)) max(abs(ylimax))];
        ax2 = axes('Position',get(axa{1},'Position'),'XAxisLocation','top',...
        'YAxisLocation','right','Color','none','XColor','k','YColor','k');
        hold(ax2, 'all');  %   <--------------------------------
        %barsimp=bar(axa{1},(1:length(meansda)),meansda,.5,'FaceAlpha',0.05,'FaceColor',[0 0 0],'EdgeColor',[0 0 0]);
        %        set(axa{1},'ylim',ylimax)
                [cmap,num,typ]=brewermap(2,'RdBu');
        barsimp=bar(axa{1},(1:length(meansdapos)),meansdapos,.5,'FaceAlpha',0.05,'FaceColor',cmap(1,:),'EdgeColor',[0 0 0]);
        barsimp=bar(axa{1},(1:length(meansdaneg)),meansdaneg,.5,'FaceAlpha',0.05,'FaceColor',cmap(2,:),'EdgeColor',[0 0 0]);
        set(axa{1},'ylim',ylimax)
    end     
end


groupvals.sitevals=sitevals;
groupvals.cond=ttypes;
groupvals.metric=targval;
groupvals.event=event;
savename=[savepath savelabel];
if gain
    savename=[savename '_gain'];
end
if plotsens
    savename=[savename '_sens'];
end
if normda
    savename=[savename '_zSel'];
end
if meanda
    savename=[savename '_avg'];
end
if scattypes
    savename=[savename '_scat'];
end

if plotv
    savename=[savename '_v'];
end
if simple
    savename=[savename '_simp'];
end
if avgsites
    savename=[savename '_avgsites'];
end

if patraonly
    savename=[savename '_patra'];
end
if normall
    savename=[savename '_zAll'];
end

savefig(figsess,savename);
saveas(figsess,savename,'jpg')
print(figsess,savename,'-painters','-depsc');

end
