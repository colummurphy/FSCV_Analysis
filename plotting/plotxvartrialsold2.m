function xcovdata=plotxvartrials(data,eventmarks,eventnames,plotparam,varargin)
%plotxvar
%data.lfp is betalfp data each channel in {ilfp}
%data.da is da each channel in {ida}
%ax{1} is axes for trial by trial plot, ax{2} for mean
%eventmarks{1} corresponds to eventnames{1} label, defines times for
%windows for xvar in seconds
xcovdata={};
fscvchs=plotparam.chnums;
lfpchs=plotparam.lfpchs;
rates=plotparam.xrates;
shuffle=0;
cluster=0;
vartype='none';         %xvar type, 'none', 'coeff', 'unbiased'
seltrials=1:size(data.lfp{1},1);
argnum=1;
label='';
peak=0;
stdthres=1;     %# of std's above std
unsquared=0;
barlags=0;
poscluster=1;
sortxvar=0;
burstgrp=0;
burstalign={};      %in nlx samples
burstwin=0.25;          %+/- seconds
burstthres=2;     %#stds 
noburst=0;
noburstwin=.75;
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'type'
            argnum=argnum+1;
            vartype=varargin{argnum};
        case 'trials'
            argnum=argnum+1;
            seltrials=varargin{argnum};
        case 'label'
            argnum=argnum+1;
            label=varargin{argnum};
        case 'shuffle'
            shuffle=1;
            label=[label '_shuffle'];
        case 'cluster'
            cluster=1;
            argnum=argnum+1;            
            %next argumen should be pos or neg to define what type of
            %changes
            if strcmp(varargin{argnum},'neg')               
                poscluster=0;
                label=[label '_clustneg'];
            else
                label=[label '_clustpos'];
            end
        case 'peaks'
            peak=1;     %look at peaks, not averages in window, for pos/neg clustering
            label=[label '_peaks'];
        case 'barlags'
            barlags=1;      %plot bar graph of lags
            
        case 'stds'
            %user defined std thres for da peak or mean changes
            argnum=argnum+1;
            stdthres=varargin{argnum};
            label=[label '_' num2str(stdthres) 'peaks'];
        case 'unsquared'
            unsquared=1;
            label=[label '_unsq'];
        case 'sortxvar'
            sortxvar=1;
        case 'burstgroup'
            %select trials where beta burst occurs at specified time points
            %listed for each trial & specify # stds above trial mean for
            %identify burst amptliude
            argnum=argnum+1;
            burstgrp=1;
            burstalign=varargin{argnum};
            argnum=argnum+1;
            burstthres2=varargin{argnum};
            if ~isempty(burstthres2)
                burstthres=burstthres;
            end
            label=[label '_burstgrp'];      
        case 'noburst'
            noburst=1;
            argnum=argnum+1;
            burstalign=varargin{argnum};
            argnum=argnum+1;
            burstthres2=varargin{argnum};
            if ~isempty(burstthres2)
                burstthres=burstthres;
            end
            label=[label '_noburstgrp'];   
    end
    argnum=argnum+1;
end
origseltrials=seltrials;
lfpchs=plotparam.lfpchs;
%get lfp ch groups p & c
pgroup=find(contains(lfpchs,'p')==1);
cgroup=find(contains(lfpchs,'c')==1);
maxplots=max(length(pgroup),length(cgroup));

fcvlfpgroups={};
%group fscv/lfp channels based on region p/c
for ii=1:length(fscvchs)
    ich=fscvchs(ii);
    fsite=plotparam.sites(ich);
    isp=contains(fsite,'p');
    fcvlfpgroups{ii}.fsite=ich;
    if isp
        fcvlfpgroups{ii}.lsites=pgroup;
    else
        fcvlfpgroups{ii}.lsites=cgroup;
    end
end
    
%set up trial by trial xcov plots
figsize=[1200 800];
pad=30;
sizet=[round(figsize(1)/(maxplots)-pad*maxplots/2) figsize(2)-100];

if ~ispc
    figsize=[1000 700];
    sizet=[round(figsize(1)/(maxplots+2)) figsize(2)-100];
    pad=25;
end

figxvar=figure('visible','off');     %figure for each channel
if ispc
figxvar=figure('visible','on');     %figure for each channel
end
set(figxvar,'Color', [1 1 1],'Position',[20,50,figsize(1),figsize(2)]);
set(0,'CurrentFigure',figxvar);    %set figure handle to current figure
axxvar=axes(figxvar);
hold(axxvar,'on');
axt={};
numrows=1;
if maxplots>6
    numrows=2; sizet(2)=sizet(2)/2;
end
for ilfp=1:maxplots
    axt{ilfp}=subplot(numrows,min(maxplots,6),ilfp);
    set(axt{ilfp},'units','pixel');
   if ilfp<=6
   set(axt{ilfp},'position',[pad + sizet(1)*(ilfp-1)+pad*(ilfp), ...
        figsize(2)-sizet(2)-pad,sizet(1),sizet(2)]);
   else
  set(axt{ilfp},'position',[pad + sizet(1)*(ilfp-7)+pad*(ilfp-6), ...
        figsize(2)-sizet(2)*2-pad*2,sizet(1),sizet(2)]);       
    end
end

figxvarsort={};
axtsort={};
figxvarsortneg={};
axtsortneg={};
if sortxvar
    figxvarsort=figure('visible','off');     %figure for each channel
    if ispc
    figxvarsort=figure('visible','on');     %figure for each channel
    end
    set(figxvarsort,'Color', [1 1 1],'Position',[20,50,figsize(1),figsize(2)]);
    set(0,'CurrentFigure',figxvarsort);    %set figure handle to current figure
    axxvarsort=axes(figxvarsort);
    hold(axxvarsort,'on');
    for ilfp=1:maxplots
        axtsort{ilfp}=subplot(numrows,min(maxplots,6),ilfp);
        set(axtsort{ilfp},'units','pixel');
       if ilfp<=6
       set(axtsort{ilfp},'position',[pad + sizet(1)*(ilfp-1)+pad*(ilfp), ...
            figsize(2)-sizet(2)-pad,sizet(1),sizet(2)]);
       else
      set(axtsort{ilfp},'position',[pad + sizet(1)*(ilfp-7)+pad*(ilfp-6), ...
            figsize(2)-sizet(2)*2-pad*2,sizet(1),sizet(2)]);       
        end
    end
    
    figxvarsortneg=figure('visible','off');     %figure for each channel
    if ispc
    figxvarsortneg=figure('visible','on');     %figure for each channel
    end
    set(figxvarsortneg,'Color', [1 1 1],'Position',[20,50,figsize(1),figsize(2)]);
    set(0,'CurrentFigure',figxvarsortneg);    %set figure handle to current figure
    axxvarsortneg=axes(figxvarsortneg);
    hold(axxvarsortneg,'on');
    for ilfp=1:maxplots
        axtsortneg{ilfp}=subplot(numrows,min(maxplots,6),ilfp);
        set(axtsortneg{ilfp},'units','pixel');
       if ilfp<=6
       set(axtsortneg{ilfp},'position',[pad + sizet(1)*(ilfp-1)+pad*(ilfp), ...
            figsize(2)-sizet(2)-pad,sizet(1),sizet(2)]);
       else
      set(axtsortneg{ilfp},'position',[pad + sizet(1)*(ilfp-7)+pad*(ilfp-6), ...
            figsize(2)-sizet(2)*2-pad*2,sizet(1),sizet(2)]);       
        end
    end  
end
%set up average plots
sizep=[250 250];
sizep=[round(figsize(1)/(maxplots+1)-maxplots*pad/10) round(figsize(2)/(3))];
if ~ispc
    sizep=[round(figsize(1)/(maxplots+2)) round(figsize(2)/(4))];
end
%figxvarmean=figure; 
figxvarmean=figure('visible','off'); 
if ispc
    figxvarmean=figure('visible','on');
end
set(figxvarmean, 'Color', [1 1 1],'Position',[20,50,figsize(1),figsize(2)]);
set(0,'CurrentFigure',figxvarmean);    %set figure handle to current figure
axvarmean=axes(figxvarmean);
ax={};
hold(axvarmean,'on');
set(ax,'units','pixel');
for ichda=1:length(fscvchs)
    for ilfp=1:maxplots
        ax{ichda,ilfp}=subplot(length(fscvchs),maxplots,...
            ilfp+(ichda-1)*maxplots);
        set(ax{ichda,ilfp},'units','pixels');
       set(ax{ichda,ilfp},'position',[pad + sizep(1)*(ilfp-1)+pad*(ilfp), ...
            figsize(2)-sizep(2)*ichda-pad*2*ichda,sizep(1),sizep(2)]);
       % plot(ax{ichda,ilfp},1,1);
    end
end
%set up bar plots, don't plot amplitudes anymore 11/25
axvarbar={};
figxvarbar=[];
axb={};
if barlags==1
figxvarbar=figure('visible','off'); 
if ispc
    figxvarbar=figure('visible','on');
end
set(figxvarbar, 'Color', [1 1 1],'Position',[20,50,figsize(1),figsize(2)]);
set(0,'CurrentFigure',figxvarbar);    %set figure handle to current figure
axvarbar=axes(figxvarbar);
for ii=1:length(fscvchs)
    axb{ii}=subplot(1,length(fscvchs),ii);
end
end
markcolor=[0 0 0 ];
xcovshuf={};
xcovdata={};
tsx={};
bardata={};
bararray=[];
barshuf={};
for ievent=1:length(eventmarks)
pathxvar{ievent}=[plotparam.pathname 'xvar_' eventnames{ievent} '_' ...
    vartype '_' label '\'];
if ~isdir(pathxvar{ievent})
    mkdir(pathxvar{ievent});
end
 postrials={};
 negtrials={};
 xcovdataneg={};
 xcovdata={};
 
 bid=1;
for ichda=1:length(fscvchs)
    ida=fscvchs(ichda);
    if cluster==1 
        %cluster signals based on positive/negative change during defined
        %event window in eventmarks{ievent} for given da channel
        das=[];
        baseline=[];
        dasub=[];
        postmean=[];
        postpeak=[];
        for itrial=1:size(data.da{ida},1)
            das(itrial,:)=fillmissing(data.da{ida}(itrial,:),'pchip');
            trialwinids=round(eventmarks{ievent}(1,itrial)./rates(1)*rates(2)):round(eventmarks{ievent}(2,itrial)./rates(1)*rates(2));
            if isempty(trialwinids)
                %ie event makrs go reverse or someting else
                trialwinids(1)=-10;
            end
            baseids=trialwinids(1)-2:trialwinids(1)+2;
            baseids=baseids(baseids>0);
            if any(trialwinids<1)
                baseline(itrial)=nan;
                dasub(itrial,:)=das(itrial,:)-baseline(itrial);
                postmean(itrial)=nan;
                postpeak(itrial)=nan;
            else
                baseline(itrial)=mean(das(itrial,baseids));
                dasub(itrial,:)=das(itrial,:)-baseline(itrial);
                postmean(itrial)=mean(dasub(itrial,trialwinids));
                postpeak(itrial)=max(dasub(itrial,trialwinids));
            end
        end                
        stdda=nanstd(data.da{ida}(:,230:370),[],2); %get stds only +/-7s of rew
        stdda2=nanstd(stdda,[],1)*stdthres;
        postrials{ida}=find(postmean>stdda2);    %get trials where post mean is > 1 std
        negtrials{ida}=find(postmean<-stdda2); 
        if peak==1
            %look at pos peaks, not mean in window
            postrials{ida}=find(postpeak>stdda2);
            negtrials{ida}=[];
        end
        if poscluster==1
            %select trials wiht positive increase
            seltrials=postrials{ida};
        else
            seltrials=negtrials{ida};
        end
    %end cluster signals
    end
    
    %only for p/c groups 11/22/2018
    for ii=1:length(fcvlfpgroups{ichda}.lsites)            
        %cross var with defined lfp chs
        %only for p / c groups
        cla(axt{ii});
        ilfp=fcvlfpgroups{ichda}.lsites(ii);
        %itplot=length(fcvlfpgroups{ichda}.lsites)*(ichda-1)+ii;
        xdata={data.lfp{ilfp} data.da{ida}};
        if unsquared
            xdata={sqrt(data.lfp{ilfp}) data.da{ida}};
        end
        
        %if burst present at selected time point +/- burstwin s
        %(ie beta env > thres at
        %target), get only these trials  
        if burstgrp==1
            btrials=[];
            for itrace=1:size(data.lfp{ilfp},1)
                bwin=burstalign{ievent}(itrace)-burstwin*rates(1):...
                    burstalign{ievent}(itrace)+burstwin*rates(1);
                if any(bwin<1) || any(bwin>size(data.lfp{ilfp},2))
                    continue
                end
                %find "bursts" of beta
                sqenvdata=data.lfp{ilfp}(itrace,:);
                thres=std(sqenvdata)*burstthres;
                idsbeta=find(sqenvdata>=thres);
                betainalign=intersect(idsbeta,bwin);    %beta in search win

                %make sure bursts outside targeted align event
                %are smaller than thres with padding
                winevent=eventmarks{ievent}(1,itrace)+250:eventmarks{ievent}(2,itrace)-250;
                winout=setdiff(winevent,bwin);
                betaoutside=intersect(idsbeta,winout);
                %if bursts outside below selectional fraction than ok
                if length(betainalign)>40 && length(betaoutside)/length(winout)<0.05
                    %should be greater than min width for burst (40 ms)
                    btrials=[btrials itrace];
                end  
            end
            seltrials=intersect(origseltrials,btrials); %get beta burst trials from indicated trials group only
            if length(seltrials)<7
                warning('too few trials after burst selection of trials, skip');
                continue
            end
        end
        if noburst==1
            %control for urst broup, when no bursts at targ aln
            %using wider ewindow
            btrials=[];
            for itrace=1:size(data.lfp{ilfp},1)
                bwin=burstalign{ievent}(itrace)-noburstwin*rates(1):...
                    burstalign{ievent}(itrace)+noburstwin*rates(1);
                if any(bwin<1) || any(bwin>size(data.lfp{ilfp},2))
                    continue
                end
                %find "bursts" of beta
                sqenvdata=data.lfp{ilfp}(itrace,:);
                thres=std(sqenvdata)*burstthres;
                idsbeta=find(sqenvdata>=thres);
                betainalign=intersect(idsbeta,bwin);    %beta in search win

                %if bursts outside below selectional fraction than ok
                if length(betainalign)<2 
                    %should be greater than min width for burst (40 ms)
                    btrials=[btrials itrace];
                end  
            end
            seltrials=intersect(origseltrials,btrials); %get beta burst trials from indicated trials group only
            if length(seltrials)<7
                warning('too few trials after burst selection of trials, skip');
                continue
            end
        end
        %xdata={betalfp{ilfp} da{ida}};
        sitenames={plotparam.lfpchs{ilfp} plotparam.sites{ida}};
        xcovshuf{ii}=xvardata(xdata,plotparam.xrates,eventmarks{ievent},plotparam,axt{ii},seltrials,...
            'eventsnlxids','sitename',sitenames,'fbands',plotparam.fbands,'type',vartype,'shuffle','eventnames',eventnames{ievent});

        if ii==1
        xcovdata{ii}=xvardata(xdata,plotparam.xrates,eventmarks{ievent},plotparam,axt{ii},seltrials,...
            'eventsnlxids','sitename',sitenames,'fbands',plotparam.fbands,'type',vartype,'eventnames',eventnames{ievent},'firstplot');
        else
        xcovdata{ii}=xvardata(xdata,plotparam.xrates,eventmarks{ievent},plotparam,axt{ii},seltrials,...
            'eventsnlxids','sitename',sitenames,'fbands',plotparam.fbands,'type',vartype,'eventnames',eventnames{ievent});
        end
        if isempty(xcovdata{ii})
            %did not process because too few trials, skip
            continue
        end
        if sortxvar
                if ii==1
                    plotxtrials(axtsort{ii},xcovdata{ii},'leftlabel','sortpos');
                    plotxtrials(axtsortneg{ii},xcovdata{ii},'leftlabel','sortneg');
                else
                    plotxtrials(axtsort{ii},xcovdata{ii},'sortpos');
                    plotxtrials(axtsortneg{ii},xcovdata{ii},'sortneg');
                end
        end        
        
        if shuffle==1
            xcovdata{ii}=xcovshuf{ii};
        end
        xcovdata{ii}.sitename=plotparam.lfpchs{ilfp};
        xcovdata{ii}.sitenameda=plotparam.sites{ida};
        if barlags==1
            %get lags with peak greater than .5 std
            %targs=find(xcovdata{ii}.xcovlag(:,2)>.5*nanstd(xcovdata{ii}.xcovlag(:,2)));
           %targs=find(~isnan(xcovdata{ii}.xcovlag(:,2)));
            %lagpeaks=xcovdata{ii}.xcovlag(targs,1);
            %11/24/2018, get all, not just > stds..
            %refind > thres of 1.5 stds
            abovethrespos=[];
            abovethresneg=[];
            countpos=1; countneg=1;
            abovethrespshuf=[];
            abovethresnshuf=[];
            countpshuf=1; countnshuf=1;
            lagpeakspos=[];
            lagpeaksneg=[];
            for itrial=1:size(xcovdata{ii}.xcovda,1)
                stdxcov(itrial)=nanstd(xcovdata{ii}.xcovda(itrial,:));
                lagpeakspos(itrial)=xcovdata{ii}.xcovlag(itrial,2);
                if lagpeakspos(itrial)>stdxcov(itrial)*2.5
                    abovethrespos(countpos)=itrial;
                    countpos=countpos+1;
                end
                lagpeaksneg(itrial)=xcovdata{ii}.xcovlaganti(itrial,2);
                if abs(lagpeaksneg(itrial))>stdxcov(itrial)*2.5
                    abovethresneg(countneg)=itrial;
                    countneg=countneg+1;
                end
                stdxcovshuf(itrial)=nanstd(xcovshuf{ii}.xcovda(itrial,:));
                lagpeak=xcovshuf{ii}.xcovlag(itrial,2);
                if lagpeak>stdxcovshuf(itrial)*2.5
                    abovethrespshuf(countpshuf)=itrial;
                    countpshuf=countpshuf+1;
                end
                lagpeakneg=xcovshuf{ii}.xcovlaganti(itrial,2);
                if abs(lagpeakneg)>stdxcovshuf(itrial)*2.5
                    abovethresnshuf(countnshuf)=itrial;
                    countnshuf=countnshuf+1;
                end

            end
            lagpeaks=xcovdata{ii}.xcovlag(abovethrespos,1);
            %bardata{ichda,ii}.coeffs=xcovdata{ii}.xcovlag(targs,2);
            bardata{ichda,ii}.coeffs=xcovdata{ii}.xcovlag(abovethrespos,2);
            bardata{ichda,ii}.lags=lagpeaks;    %lags for site pair
            meand=nanmean(bardata{ichda,ii}.coeffs);
            ci=nanstd(bardata{ichda,ii}.coeffs)./sqrt(length(bardata{ichda,ii}.coeffs))*1.96;
            cidata=[meand-ci meand+ci];
            bardata{ichda,ii}.coeffci=cidata;
            meand=nanmean(bardata{ichda,ii}.lags);
            ci=nanstd(bardata{ichda,ii}.lags)./sqrt(length(bardata{ichda,ii}.lags))*1.96;
            cidata=[meand-ci meand+ci];
            bardata{ichda,ii}.lagci=cidata;
            
            %negative lags
            %THIS is getting lags for peaks that are negative not lags that
            %are negative
           % targs=find(xcovdata{ii}.xcovlag(:,2)<-.5*nanstd(xcovdata{ii}.xcovlag(:,2)));
           % bardata{ichda,ii}.coeffsneg=xcovdata{ii}.xcovlag(targs,2);
            bardata{ichda,ii}.coeffsneg=xcovdata{ii}.xcovlaganti(abovethresneg,2);
            bardata{ichda,ii}.lagsneg=xcovdata{ii}.xcovlaganti(abovethresneg,1);
            meand=nanmean(bardata{ichda,ii}.coeffsneg);
            ci=nanstd(bardata{ichda,ii}.coeffsneg)./sqrt(length(bardata{ichda,ii}.coeffsneg))*1.96;
            cidata=[meand-ci meand+ci];
            bardata{ichda,ii}.coeffcineg=cidata;
            meand=nanmean(bardata{ichda,ii}.lagsneg);
            ci=nanstd(bardata{ichda,ii}.lagsneg)./sqrt(length(bardata{ichda,ii}.lagsneg))*1.96;
            cidata=[meand-ci meand+ci];
            bardata{ichda,ii}.lagcineg=cidata;
            
            bardata{ichda,ii}.label=[plotparam.sites{ida} ', ' plotparam.lfpchs{ilfp}];
            bardata{ichda,ii}.event=['xvar_' eventnames{ievent} '_' ...
                vartype '_' label];
            bardata{ichda,ii}.bid=bid;
            
            %Shuffled get lags with peak greater than .5 std
            if isempty(abovethrespshuf)
                abovethrespshuf=1:size(xcovdata{ii}.xcovda,1);
            end
            lagpeaks=xcovshuf{ii}.xcovlag(abovethrespshuf,1);
            barshuf{ichda,ii}.coeffs=xcovshuf{ii}.xcovlag(abovethrespshuf,2);
            barshuf{ichda,ii}.lags=lagpeaks;    %lags for site pair
            meand=nanmean(barshuf{ichda,ii}.coeffs);
            ci=nanstd(barshuf{ichda,ii}.coeffs)./sqrt(length(barshuf{ichda,ii}.coeffs))*1.96;
            cidata=[meand-ci meand+ci];
            barshuf{ichda,ii}.coeffci=cidata;
            meand=nanmean(barshuf{ichda,ii}.lags);
            ci=nanstd(barshuf{ichda,ii}.lags)./sqrt(length(barshuf{ichda,ii}.lags))*1.96;
            cidata=[meand-ci meand+ci];
            barshuf{ichda,ii}.lagci=cidata;
            if isempty(abovethresnshuf)
                abovethresnshuf=1:size(xcovdata{ii}.xcovda,1);
            end            
            %negative lags
            barshuf{ichda,ii}.coeffsneg=xcovshuf{ii}.xcovlaganti(abovethresnshuf,2);
            barshuf{ichda,ii}.lagsneg=xcovshuf{ii}.xcovlaganti(abovethresnshuf,1);
            meand=nanmean(barshuf{ichda,ii}.coeffsneg);
            ci=nanstd(barshuf{ichda,ii}.coeffsneg)./sqrt(length(barshuf{ichda,ii}.coeffsneg))*1.96;
            cidata=[meand-ci meand+ci];
            barshuf{ichda,ii}.coeffcineg=cidata;
            meand=nanmean(barshuf{ichda,ii}.lagsneg);
            ci=nanstd(barshuf{ichda,ii}.lagsneg)./sqrt(length(barshuf{ichda,ii}.lagsneg))*1.96;
            cidata=[meand-ci meand+ci];
            barshuf{ichda,ii}.lagcineg=cidata;
            
            barshuf{ichda,ii}.label=[bardata{ichda,ii}.label '-s'];
            barshuf{ichda,ii}.event=bardata{ichda,ii}.event;
            currsize=size(bararray,1);
            idsfill=currsize+(1:length(bardata{ichda,ii}.lags));  
            bararray(idsfill,1)=bid;        %fill with id of sitepair
            bararray(idsfill,2)=bardata{ichda,ii}.lags;
            [H,pp] = ttest2(bardata{ichda,ii}.lags, barshuf{ichda,ii}.lags,'tail','both');
            bardata{ichda,ii}.plags=pp;
            [H,pp] = ttest2(bardata{ichda,ii}.lagsneg, barshuf{ichda,ii}.lagsneg,'tail','both');
            bardata{ichda,ii}.plagsneg=pp;
                   %     chivalues(ii)=dg_chi2test2([NB ND],1,'binwidths',[TB TD]);

            bid=bid+1;
        end
        xmod=xcovdata{ii}.xcovda;
        tsx=xcovdata{ii}.tsx;
        xmodmean=nanmean(xmod,1);
        nums=length(~isnan(mean(xmod,2)));
        xmodshuf=xcovshuf{ii}.xcovda; tsxshuf=xcovshuf{ii}.tsx;
        xmodshufmean=nanmean(xmodshuf,1);
        numshuf=length(~isnan(mean(xmodshuf,2)));

        xmodci=nanstd(xmod,[])./sqrt(nums)*1.96;
        xmodshufci=nanstd(xmodshuf,[])./sqrt(numshuf)*1.96;
        maxxmod=find(abs(xmodmean)==max(abs(xmodmean)));
        peakxmod=xmodmean(maxxmod);
        lagxmod=tsx(maxxmod);
        cixmod=xmodci(maxxmod);
        peakxmodshuf=xmodshufmean(maxxmod);
        cixmodshuf=xmodshufci(maxxmod);
        if ~isempty(peakxmodshuf)
            if peakxmod+cixmod<peakxmodshuf-cixmodshuf || ...
                    peakxmod-cixmod>peakxmodshuf+cixmodshuf
                %if sig different
                markcolor=[1 0 0];
            else
                markcolor=[0 0 0];
            end
        cla(ax{ichda,ii});
        plot(ax{ichda,ii},tsx,xmodmean,'marker','none','linewidth',2,'color',markcolor)
        xlabel(ax{ichda,ii},'lag (s)')
        if ii==1
            ylabel(ax{ichda,ii},'trial averaged x-cov')
        end
        title(ax{ichda,ii},[eventnames{ievent} ' | xcov(lfp ' lfpchs{ilfp} ' , DA ' plotparam.sites{ida} ')'])
        hold(ax{ichda,ii},'on'); 
        plot(ax{ichda,ii},tsx,xmodmean+xmodci,'--','linewidth',1,'color',markcolor)
        plot(ax{ichda,ii},tsx,xmodmean-xmodci,'--','linewidth',1,'color',markcolor)
        xlim(ax{ichda,ii},[min(tsx) max(tsx)]);
        
        end       
    end 
    save([pathxvar{ievent} 'xvar_' plotparam.sites{ida}],'xcovdata');
   savename=['lfp_x_da_' plotparam.sites{ida}];
    saveas(figxvar,[pathxvar{ievent} savename],'jpg')
    if sortxvar
        saveas(figxvarsort,[pathxvar{ievent} savename 'sorted'],'jpg')
        saveas(figxvarsortneg,[pathxvar{ievent} savename 'sortedneg'],'jpg')
    end
    %if cluster==1
    %    save([pathxvar{ievent} 'xvarneg_' plotparam.sites{ida}],'xcovdataneg', 'postrials','negtrials');
    %    saveas(figxvarn,[pathxvar{ievent} savename 'n'],'jpg')
   % end
    
end
%%
%bar graph mean lags for site pairs for current event window
if barlags==1
    cla(axb{1}); cla(axb{2});
    hold(axb{1},'on'); hold(axb{2},'on')
    
    %plot(axvarbar,bararray(:,1),bararray(:,2),'o')

    for ichda=1:length(fscvchs)
            bid=1;
    labels=[];
       % plotid=(ichda-1)*2+1;
        plotid=ichda;
                    title(axb{plotid},bardata{1,1}.event)
        cmark=[0 0 0];
        cmarksig=[1 0 0];
            numsitepairs=length(fscvchs)*length(fcvlfpgroups{ichda}.lsites);
        for ii=1:length(fcvlfpgroups{ichda}.lsites)
            %plot lags both positive and negative peaks
            if bardata{ichda,ii}.plags<0.05
                cmark=[1 0 0];
            else
                cmark=[ 0 0 0];
            end
            if bardata{ichda,ii}.plagsneg<0.05
                cmarkneg=[1 0 0];
            else
                cmarkneg=[ 0 0 0];
            end
            xdata=repmat(bid,1,length(bardata{ichda,ii}.lags));
            ranjitter=rand(1,length(bardata{ichda,ii}.lags))*.15+.25;
            scatter(axb{plotid},xdata-ranjitter,bardata{ichda,ii}.lags,65,'o','markeredgecolor',cmark,'MarkerEdgeAlpha',.3,'linewidth',1.5);
            line(axb{plotid},'xdata',[bid bid]-.1,'ydata',bardata{ichda,ii}.lagci);                     
            xdata=repmat(bid,1,length(bardata{ichda,ii}.lagsneg));   
                        ranjitter=rand(1,length(bardata{ichda,ii}.lagsneg))*.15+.25;
            scatter(axb{plotid},xdata+ranjitter,bardata{ichda,ii}.lagsneg,65,'sq','markeredgecolor',cmarkneg,'MarkerEdgeAlpha',.3,'linewidth',1.5);
            line(axb{plotid},'xdata',[bid bid]+.1,'ydata',bardata{ichda,ii}.lagcineg);
            ylim(axb{plotid},[-3 3]);
            
            %plots coeffs
            %{
            xdata=repmat(bid,1,length(bardata{ichda,ii}.coeffs));
            plot(axb{plotid+1},xdata,bardata{ichda,ii}.coeffs,'o','color',cmark);
            line(axb{plotid+1},'xdata',[bid bid],'ydata',bardata{ichda,ii}.coeffci);
            xdata=repmat(bid,1,length(bardata{ichda,ii}.coeffsneg));
            plot(axb{plotid+1},xdata,bardata{ichda,ii}.coeffsneg,'sq','color',cmarkneg);
            line(axb{plotid+1},'xdata',[bid bid],'ydata',bardata{ichda,ii}.coeffcineg);
            ylim(axb{plotid+1},[-3 3]);
                        %}

            labels{bid}=bardata{ichda,ii}.label;            
            bid=bid+1;
            
            %plot shuffeld same
            %{
            xdata=repmat(bid,1,length(barshuf{ichda,ii}.lags));
            plot(axb{plotid},xdata-.25,barshuf{ichda,ii}.lags,'o','color',cmark);
            line(axb{plotid},'xdata',[bid bid]-.1,'ydata',barshuf{ichda,ii}.lagci);

            xdata=repmat(bid,1,length(barshuf{ichda,ii}.lagsneg));                        
            plot(axb{plotid},xdata+.25,barshuf{ichda,ii}.lagsneg,'sq','color',cmarkneg);
            line(axb{plotid},'xdata',[bid bid]+.1,'ydata',barshuf{ichda,ii}.lagcineg);
                        %{
            xdata=repmat(bid,1,length(barshuf{ichda,ii}.coeffs));
            plot(axb{plotid+1},xdata,barshuf{ichda,ii}.coeffs,'o','color',cmark);
            line(axb{plotid+1},'xdata',[bid bid],'ydata',barshuf{ichda,ii}.coeffci);
                        
            xdata=repmat(bid,1,length(barshuf{ichda,ii}.coeffsneg));
            plot(axb{plotid+1},xdata,barshuf{ichda,ii}.coeffsneg,'sq','color',cmarkneg);
            line(axb{plotid+1},'xdata',[bid bid],'ydata',barshuf{ichda,ii}.coeffcineg);
                        %}
            labels{bid}=barshuf{ichda,ii}.label;
            bid=bid+1;    
       %}
        end
        set(axb{plotid},'xtick',1:numsitepairs,'xticklabel',labels);
        set(axb{plotid},'xTickLabelRotation',90)
        %{
        set(axb{plotid+1},'xtick',1:numsitepairs*2,'xticklabel',labels);
        set(axb{plotid+1},'xTickLabelRotation',90)    
        ylabel(axb{plotid+1},'peak xcov')
        %}
        ylabel(axb{plotid},'lag at peak xcov (s)')
 
    end
            save([pathxvar{ievent} 'xvar_bardata_da_' plotparam.sites{ida}], 'barshuf', 'bardata');
        saveas(figxvarbar,[pathxvar{ievent} 'lfp_x_da_bar_da_' plotparam.sites{ida}],'jpg')
                saveas(figxvarbar,[pathxvar{ievent} savename 'bar' '.fig'])

   % names = {'CRHS'; 'ELLY'; 'LGWD'; 'ECFS'; 'THMS'};
%set(gca,'xtick',[1:5],'xticklabel',names)
end

if ispc
set(findall(figxvarmean,'-property','FontSize'),'FontSize',12)
else
 set(findall(figxvarmean,'-property','FontSize'),'FontSize',8)
end
   
%plot mean xvar
savename=['lfp_x_da_' plotparam.sites{ida}];
saveas(figxvarmean,[pathxvar{ievent} savename 'mean'],'bmp')
    
end


close all;
end

