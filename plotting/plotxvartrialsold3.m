function xcovdata=plotxvartrials(data,eventmarks,eventnames,plotparam,varargin)
%plotxvar
%data.lfp is betalfp data each channel in {ilfp}
%data.da is da each channel in {ida}
%ax{1} is axes for trial by trial plot, ax{2} for mean
%eventmarks{1} corresponds to eventnames{1} label, defines times for
%windows for xvar in seconds
%12/06/2018 took out all shuffle functions, saved old ver as
%plotxvartrialsold2.m for reference
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
numtypes=1;
numlfptypes=1;
lfplabels={};
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
            label=[label '_'];
        case 'dasignaltypes'
            cluster=1;
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
            numlfptypes=2;
            lfplabels={'burstgrp','noburstgrp'};
            burstalign=varargin{argnum};
            argnum=argnum+1;
            burstthres2=varargin{argnum};
            if ~isempty(burstthres2)
                burstthres=burstthres;
            end
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
    
origlabel=label;
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
    vartype '_' label filesep];
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
    %Middle of window defined as subtraction point/baseline
    das=[];
    sigtypes={};
    stdall=[];      %mean std deviation baseline points averaged for all trials as global baseline
    basedastds=[];
    for itrial=1:size(data.da{ida},1)
        das(itrial,:)=fillmissing(data.da{ida}(itrial,:),'pchip');
        trialwinids=round(eventmarks{ievent}(1,itrial)./rates(1)*rates(2)):round(eventmarks{ievent}(2,itrial)./rates(1)*rates(2));
         if ~(any(trialwinids<1) || length(trialwinids)<7)
             basedastds(itrial)=nanstd(das(itrial,trialwinids(1:median(1:length(trialwinids)))));
         else
             basedastds(itrial)=nan;
         end
    end
    stdall=nanmean(basedastds);
    signames={};
   % das=[];
    baseline=[];
    dasub=[];
    postmean=[];
    postpeak=[];        %max change after middle of window
    targda=[];          %change at middle of window
    postrials{ida}=[];
    negtrials{ida}=[];
    pausetrials{ida}=[];        %trials with transient suppression of da at middle of window & rebound
    nopausetrials{ida}=[];      %trials without pause/rebound da, but increases
    sustainedtrials{ida}=[];        %trials without pause, sustained
    for itrial=1:size(data.da{ida},1)
       % das(itrial,:)=fillmissing(data.da{ida}(itrial,:),'pchip');
        trialwinids=round(eventmarks{ievent}(1,itrial)./rates(1)*rates(2)):round(eventmarks{ievent}(2,itrial)./rates(1)*rates(2));
        if isempty(trialwinids)
            %ie event makrs go reverse or someting else
            trialwinids(1)=-10;
        end
        if any(trialwinids<1) || length(trialwinids)<7
            baseline(itrial)=nan;
            dasub(itrial,:)=das(itrial,:)-baseline(itrial);
            postpeak(itrial)=nan;
        else
            mididx=median(1:length(trialwinids));
            baseids=trialwinids(mididx)-2:trialwinids(mididx);
            baseids=baseids(baseids>0);
            postids=trialwinids(mididx)+1:trialwinids(end);
            baseline(itrial)=nanmean(das(itrial,baseids));
            dasub(itrial,:)=das(itrial,:)-baseline(itrial);
            %look at max abs to see if neg/pos change
            postda=max(abs(dasub(itrial,postids)));       
            maxid=find(abs(dasub(itrial,postids))==postda)+trialwinids(mididx);
            if length(maxid)>1
                maxid=maxid(1);
            end
            if isnan(postda) || isempty(maxid)
                postpeak(itrial)=nan;
            else
                maxda=dasub(itrial,maxid);
                postpeak(itrial)=maxda;
                stdda=nanstd(das(itrial,trialwinids(1):trialwinids(mididx))); %get stds only +/-7s of rew
                if maxda>=stdall*2
                    postrials{ida}=[postrials{ida} itrial];
                end
                if maxda<0 && abs(maxda)>stdall*2
                    negtrials{ida}=[negtrials{ida} itrial];
                end
            end

            %look for transient da pauses & no pauses
            quartidx=median(1:mididx);
            baseids=trialwinids(1):trialwinids(quartidx);
            baseids=baseids(baseids>0);
            targids=trialwinids(mididx-2:mididx+2);     %look at midx for puase
            baseline(itrial)=nanmean(das(itrial,baseids));
            dasub(itrial,:)=das(itrial,:)-baseline(itrial);
            %look at change at target point to see if transient
            %deflection
            targda=dasub(itrial,targids);       
            if any(isnan(targda)) || isempty(targda)
                %
            else
                %targda(itrial)=targda;
                stdda=nanstd(das(itrial,trialwinids(1):trialwinids(quartidx))); %get stds only +/-7s of rew
                deflects=any(targda<-stdall*.25) && nanmean(dasub(itrial,trialwinids(quartidx:mididx-2)))>0;       %if decreases below baseline
                maxinwin=max(dasub(itrial,trialwinids(mididx:end)),[],'omitnan');
                mininwin=min(dasub(itrial,trialwinids(mididx:end)),[],'omitnan');
                rebounds=maxinwin>stdall*.25;
                nopauses=all(targda>-stdall*.25);
                sustains=all(targda>=-stdall) && maxinwin<stdall*2 && mininwin>-stdall*2;
                if isnan(maxinwin) || isempty(maxinwin) || stdda==0
                    %not a pause/rebound or sustained response do nothing
                elseif deflects && rebounds
                    pausetrials{ida}=[pausetrials{ida} itrial];
                elseif nopauses && rebounds
                    %no deflection, but might increase
                    nopausetrials{ida}=[nopausetrials{ida} itrial];
                elseif sustains
                    %no deflection, sustained
                    sustainedtrials{ida}=[sustainedtrials{ida} itrial];
                end
            end

        end
    end                
    if poscluster==1
        %select trials with positive increase
        seltrials=postrials{ida};
    else
        seltrials=negtrials{ida};
    end
    %end cluster signals
    sigtypes{1}=postrials{ida};
    sigtypes{2}=negtrials{ida};
    sigtypes{3}=pausetrials{ida};
    sigtypes{4}=nopausetrials{ida};
    sigtypes{5}=sustainedtrials{ida};
    signames={'dapos','daneg','dapause','danopause','dasust'};
    numtypes=length(signames);
    end
    
    for sigtype=1:numtypes
    %scroll different types of targeted clustered signals
        pathxvar{ievent}=[plotparam.pathname 'xvar_' eventnames{ievent} '_' ...
        vartype filesep];
        if ~isdir(pathxvar{ievent})
            mkdir(pathxvar{ievent});
        end
        if numtypes>1
            %use signames labels for da changes groups
            label=[origlabel signames{sigtype}];
            seltrials=sigtypes{sigtype};
        end
        for lfptype=1:numlfptypes
               %scroll lfp groups
                %REDUNDANT BECUASE FIND BURSTS EVERY DA CHANNEL....
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
            if burstgrp==1                
            %if burst present at selected time point +/- burstwin s
            %(ie beta env > thres at
            %target), get only these trials  
                btrials=[];
                nbtrials=[];    %no burst trials
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
                     %if bursts outside below selectional fraction than ok
                    if length(betainalign)<2 
                        %should be greater than min width for burst (40 ms)
                        nbtrials=[nbtrials itrace];
                    end  
                end
                lfptypes{1}=intersect(origseltrials,btrials); %get beta burst trials from indicated trials group only
                lfptypes{2}=intersect(origseltrials,nbtrials); %get no beta burst trials from indicated trials group only
            end

               if ~isempty(lfplabels)
                   %indicates that we should cluster lfps by burst
                   label=[origlabel lfplabels{lfptype}];
                   seltrials=lfptypes{lfptype};
               end
                if length(seltrials)<7
                    warning('too few trials after burst selection of trials, skip');
                    continue
                end
                %xdata={betalfp{ilfp} da{ida}};
                sitenames={plotparam.lfpchs{ilfp} plotparam.sites{ida}};
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
                %xcovdata{ii}.sitename=plotparam.lfpchs{ilfp};
                %xcovdata{ii}.sitenameda=plotparam.sites{ida};
               
                xmod=xcovdata{ii}.xcovda;
                tsx=xcovdata{ii}.tsx;
                xmodmean=nanmean(xmod,1);
                nums=length(~isnan(mean(xmod,2)));

                xmodci=nanstd(xmod,[])./sqrt(nums)*1.96;
                maxxmod=find(abs(xmodmean)==max(abs(xmodmean)));
                peakxmod=xmodmean(maxxmod);
                lagxmod=tsx(maxxmod);
                cixmod=xmodci(maxxmod);

                cla(ax{ichda,ii});
                plot(ax{ichda,ii},tsx,xmodmean,'marker','none','linewidth',2,'color',markcolor)
                xlabel(ax{ichda,ii},'lag (s)')
                if ii==1
                    ylabel(ax{ichda,ii},'trial averaged x-cov')
                end
                title(ax{ichda,ii},[eventnames{ievent} ' | ' label ' | xcov(lfp ' lfpchs{ilfp} ' , DA ' plotparam.sites{ida} ')'])
                hold(ax{ichda,ii},'on'); 
                plot(ax{ichda,ii},tsx,xmodmean+xmodci,'--','linewidth',1,'color',markcolor)
                plot(ax{ichda,ii},tsx,xmodmean-xmodci,'--','linewidth',1,'color',markcolor)
                xlim(ax{ichda,ii},[min(tsx) max(tsx)]);
            end 
            save([pathxvar{ievent} 'xvar_' label '_' plotparam.sites{ida}],'xcovdata');
            savename=['lfp_x_da_' label '_' plotparam.sites{ida}];
            saveas(figxvar,[pathxvar{ievent} savename ],'jpeg')
            %saveas(fig2,[pathlfp{3} trialname '.tif'],'tif')
            if sortxvar
                saveas(figxvarsort,[pathxvar{ievent} savename 'sorted'],'jpeg')
                saveas(figxvarsortneg,[pathxvar{ievent} savename 'sortedneg'],'jpeg')
            end  
        end
    end
end

if ispc
    set(findall(figxvarmean,'-property','FontSize'),'FontSize',12)
else
    set(findall(figxvarmean,'-property','FontSize'),'FontSize',8)
end
   
%plot mean xvar
savename=['lfp_x_da_' label '_' plotparam.sites{ida}];
saveas(figxvarmean,[pathxvar{ievent} savename 'mean'],'bmp')
    
end


close all;
end

