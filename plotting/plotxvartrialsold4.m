function xcovdata=plotxvartrials(data,eventmarks,eventnames,plotparam,varargin)
%USED PRIOR TO 1/26/2019
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
        pade=2;             %window samples (10 hz rate)
        basesamples=3;      %window samples (10 hz rate)
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
            %label=[label '_'];
        case 'cluster'
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
sizep=[round(figsize(1)/(maxplots+1)-maxplots*pad/10) round(figsize(2)/3)];
if ~ispc
    sizep=[round(figsize(1)/(maxplots+2)) round(figsize(2)/(4))];
end
%figxvarmean=figure; 
figxvarmean=figure('visible','off'); 
if ispc
    figxvarmean=figure('visible','on');
end
set(figxvarmean, 'Color', [1 1 1],'Position',[20,50,figsize(1),round(figsize(2)/2)]);
set(0,'CurrentFigure',figxvarmean);    %set figure handle to current figure
axvarmean=axes(figxvarmean);
ax={};
hold(axvarmean,'on');
set(ax,'units','pixel');
for ilfp=1:maxplots
    ax{ilfp}=subplot(1,maxplots,ilfp);
    set(ax{ilfp},'units','pixels');
    set(ax{ilfp},'position',[pad + sizep(1)*(ilfp-1)+pad*(ilfp), ...
            pad*2,sizep(1),sizep(2)]);
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
    %interfix, intertarg, different event windows..
pathxvar{ievent}=[plotparam.pathname 'xvar_' eventnames{ievent} '_' ...
    vartype label filesep];
if ~isdir(pathxvar{ievent})
    mkdir(pathxvar{ievent});
end

for ichda=1:length(fscvchs)
    xcovdata={};        %saved in separate file for each da ch
    ida=fscvchs(ichda);
    
    %selected trials based on da characteristics in window
    targalninc=[];      %transient da increase aligned to targ (center wind)
    targalndec=[];      %transient pause aligned to targ (center of window)
    sustinc=[];         %sustained increase without local valley
    sustdec=[];     %sustained decrease without local rise
    sust=[];
    trialsinc=[];       %alninc and sustinc
    trialsdec=[];
    tracealninc=[];
    tracealndec=[];
    tracesustinc=[];
    tracesustdec=[];
    tracesust=[];
    traceinc=[];    %traces inc aln & sust
    tracedec=[];
    tsinc=[];   %delta timestamp of maxima for alninc and sustinc
    tsdec=[];
    betatarginc=[];
    betatargdec=[];
    tracebetainc=[];
    tracebetadec=[];
    tsbetainc=[];   %delta timestamp of maxima for inc
    tsbetadec=[];
    tsbetaincdapeak=[]; %timestamp of maxima da when beta burst inc
    for iax=1:length(ax)
        cla(ax{iax});
        cla(axt{iax});
    end
for ii=1:length(fcvlfpgroups{ichda}.lsites)            
    %cross var with defined lfp chs for p / c groups
    cla(axt{ii});
    ilfp=fcvlfpgroups{ichda}.lsites(ii);
    %itplot=length(fcvlfpgroups{ichda}.lsites)*(ichda-1)+ii;
    xdata={data.lfp{ilfp} data.da{ida}};
    if unsquared
        xdata={sqrt(data.lfp{ilfp}) data.da{ida}};
    end
    sitenames={plotparam.lfpchs{ilfp} plotparam.sites{ida}};
    if ii==1
        xcovdata{ii}=xvardata(xdata,plotparam.xrates,eventmarks{ievent},plotparam,axt{ii},seltrials,...
            'sitename',sitenames,'fbands',plotparam.fbands,'type',vartype,'eventnames',eventnames{ievent},'firstplot');
    else
        xcovdata{ii}=xvardata(xdata,plotparam.xrates,eventmarks{ievent},plotparam,axt{ii},seltrials,...
            'sitename',sitenames,'fbands',plotparam.fbands,'type',vartype,'eventnames',eventnames{ievent});
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
    %data for plotting means
    xmod=xcovdata{ii}.xcovda;
    tsx=xcovdata{ii}.tsx;
    xmodmean=nanmean(xmod,1);
    nums=length(~isnan(mean(xmod,2)));
    xmodci=nanstd(xmod,[])./sqrt(nums)*1.96;    
    
    plot(ax{ii},tsx,xmodmean,'marker','none','linewidth',2,'color',markcolor)
    xlabel(ax{ii},'lag (s)')
    if ii==1
        ylabel(ax{ii},'trial averaged x-cov')
        title(ax{ii},[eventnames{ievent} ' | ' lfpchs{ilfp} ' , ' plotparam.sites{ida}])
    else
    title(ax{ii},[lfpchs{ilfp} ' , ' plotparam.sites{ida}])
    end
    hold(ax{ii},'on'); 
    plot(ax{ii},tsx,xmodmean+xmodci,'--','linewidth',1,'color',markcolor)
    plot(ax{ii},tsx,xmodmean-xmodci,'--','linewidth',1,'color',markcolor)
    xlim(ax{ii},[min(tsx) max(tsx)]);
    
    if cluster
        %group data based on different signal characteristics
        %cluster signals based on positive/negative change during defined
        %event window in eventmarks{ievent} for given da channel
        %Middle of window defined as subtraction point/baseline
        %NO NOW CHANGED SO EVENTMARKS only identifies align idx
        %/1/10/19, 
        %FIXED ALSO WINIDS referring to tid rather than itrial
        %MADE SO WINIDS refers to centered range around cue iddx
        %since contains only seltrials
        %cluster based on da properties only if first lfp channel
        %cluster based on lfp properties only if first da channel 
        if ii==1
        das=xcovdata{ii}.dadata;            %all trials
        winids=xcovdata{ii}.winids;     %only seltrials
        seltrials=xcovdata{ii}.seltrials;      
        sigtypes={};
        signames={};
        stdall=[];      %mean std deviation baseline points averaged for all trials as global baseline
        basedastds=[];
        for itrial=1:length(seltrials)
            %get std of signal baseline in first half of each signal win
            curda=das(seltrials(itrial),winids(itrial,:));
            basedastds(itrial)=nanstd(curda(1:median(1:length(curda))));    
        end
        stdall=nanmean(basedastds);     %mean of all stds
        baseline=[];
        dasub=[];
        posttargda=nan(1,length(seltrials));      %mean da post mididx for all sel trials
        maxdachange=nan(1,length(seltrials));     %max inflection of da post mididx for all sel trials
        postmean=[];
        postpeak=[];        %max change after middle of window
        targda=[];          %change at middle of window
        traceda=[];     %all baseline sub da for cur win
        xcovdata{ii}.daavgstds=stdall;%da baseline
        for itrial=1:length(seltrials)
            curda=das(seltrials(itrial),:);    %da trace
            trialwinids=winids(itrial,:);
            mididx=median(1:length(trialwinids));         %middle of window. ie targ period
            %baseids=winids(mididx)-2:winids(mididx);
            baseids=trialwinids(mididx)-pade-basesamples:trialwinids(mididx)-pade;
            baseids=baseids(baseids>0);
            targids=trialwinids(mididx-pade:mididx+pade);     %look at midx for puase
            baseline=nanmean(curda(baseids));
            curda=curda-baseline;
            %look at change at target point to see if transient
            %deflection
            targda=curda(targids);       
            if any(isnan(targda)) || isempty(targda)
                %
            else
                stdda=nanstd(curda(baseids)); %get stds curr trial baseline
                %local targ dec if below baseline thres & pos for remainder
                avgtargda=nanmean(targda);
                targdec=any(targda<=-stdall) ;
                %&& ...       nanmean(curda(trialwinids(mididx+pade:end)))>avgtargda+stdall*.25;       
                targinc=any(targda>=stdall*.75);
                %&& ... nanmean(curda(trialwinids(mididx+pade:end)))<avgtargda;
                if sum(curda(trialwinids))~=0
                    posttargda(itrial)=nanmean(curda(trialwinids(mididx+pade:end)));
                end
                %find max change in da transient
                maxdachange(itrial)=nan;                
                maxdaid=find(curda(trialwinids(mididx-pade:end))==max(curda(trialwinids(mididx-pade:end)),[],'omitnan'))+mididx-pade-1;
                if ~isempty(maxdaid) && sum(curda(trialwinids))~=0
                    mindaid=find(curda(trialwinids(pade+1:maxdaid))==min(curda(trialwinids(pade+1:maxdaid)),[],'omitnan'))+pade;
                    if ~isempty(mindaid)
                        maxdachange(itrial)=curda(trialwinids(maxdaid))-curda(trialwinids(mindaid));                
                    end    
                end
                
                traceda(itrial,:)=curda(trialwinids);       %store all sub da win traces
                sustpos=nanmean(curda(trialwinids(mididx+pade:end)))>stdall*.75;
                sustneg=nanmean(curda(trialwinids(mididx+pade:end)))<-stdall;
                nopause=any(targda<=stdall*.25) && any(targda>=-stdall*.25) &&...
                    nanmean(curda(trialwinids(mididx+pade:end)))<stdall*.75 &&...
                    nanmean(curda(trialwinids(mididx+pade:end)))>-stdall*.75;
                if stdda==0
                    %not a pause/rebound or sustained response do nothing
                elseif targdec
                    targalndec=[targalndec seltrials(itrial)];
                    tracealndec=[tracealndec; curda(trialwinids)];
                    tstarg=find(curda(trialwinids(mididx-pade:end))==min(targda,[],'omitnan'))+mididx-pade-1;
                    tsdec=[tsdec tstarg];  
                    trialsdec=[trialsdec; seltrials(itrial)];
                    tracedec=[tracedec; curda(trialwinids)];
                elseif targinc
                    targalninc=[targalninc seltrials(itrial)];
                     tracealninc=[tracealninc; curda(trialwinids)];
                     %get max delta ts from center of window
                     tstarg=find(curda(trialwinids(mididx-pade:end))==max(targda,[],'omitnan'))+mididx-pade-1;
                     tsinc=[tsinc tstarg];
                     trialsinc=[trialsinc; seltrials(itrial)];
                     traceinc=[traceinc; curda(trialwinids)];
               elseif sustpos 
                    sustinc=[sustinc seltrials(itrial)];
                    tracesustinc=[tracesustinc; curda(trialwinids)];
                    %get max delta ts from center of window
                    if ~targinc
                        %if not already stored maxima
                        tstarg=find(curda(trialwinids(mididx-pade:end))==max(curda(trialwinids(mididx-pade:end)),[],'omitnan'))+mididx-pade-1;
                        tsinc=[tsinc tstarg]; 
                        trialsinc=[trialsinc; seltrials(itrial)];
                        traceinc=[traceinc; curda(trialwinids)];
                    end
                elseif sustneg 
                    sustdec=[sustdec seltrials(itrial)];
                    tracesustdec=[tracesustdec; curda(trialwinids)];
                    if ~targdec
                        %if not already stored
                        tstarg=find(curda(trialwinids(mididx-pade:end))==min(curda(trialwinids(mididx-pade:end)),[],'omitnan'))+mididx-pade-1;
                        tsdec=[tsdec tstarg]; 
                        trialsdec=[trialsdec; seltrials(itrial)];
                        tracedec=[tracedec; curda(trialwinids)];
                    end
                elseif ~(targdec && targinc && sustpos && sustneg) && nopause
                    %no deflection, sustained
                    sust=[sust seltrials(itrial)];
                    tracesust=[tracesust; curda(trialwinids)];
                end
            end
            dasub(itrial,:)=curda;      %store subtracted da signal
        end
        end
        
        %if ichda==1
            %only do all lfp groups when for all da channels even if
            %repeat?
            blfp=xcovdata{ii}.lfpdata;
            winids=xcovdata{ii}.winids;
            seltrials=xcovdata{ii}.seltrials;      
            sigtypes={};
            signames={};
            basestds=[];
            basemean=[];
            betatarginc=[];
            betatargdec=[];
            tracebetainc=[];
            tracebetadec=[];
            tsbetainc=[];   %delta timestamp of maxima for inc
            tsbetadec=[];
            tsbetaincdapeak=[];
            tsbetaincdadiffpeak=[];
            for itrial=1:length(seltrials)
                %get std of signal baseline in first half of each signal win
                curlfp=blfp(seltrials(itrial),winids(itrial,:));
                basemean(itrial)=nanmean(curlfp(1:median(1:length(curlfp))));    
                basestds(itrial)=nanstd(curlfp(1:median(1:length(curlfp))));
            end
            meanbaseline=nanmean(basemean);     %mean of all baseline periods
            meanstds=nanmean(basestds);     %mean of all baseline periods
            xcovdata{ii}.betabaseline=meanbaseline;
            xcovdata{ii}.betastds=meanstds;
            for itrial=1:length(seltrials)
                curlfp=blfp(seltrials(itrial),:);    %da trace
                trialwinids=winids(itrial,:);
                mididx=median(1:length(trialwinids));         %middle of window. ie targ period
                targids=trialwinids(mididx-pade:mididx+pade);     %look at midx for puase
                %look at change at target point to see if transient
                targlfp=curlfp(targids);       
                targburst=nanmean(targlfp)>meanbaseline+meanstds ;
                %&&...nanmean(curlfp(trialwinids(mididx+pade+1:end)))<=meanbaseline+meanstds*2;
                targsup=nanmean(targlfp)<meanbaseline-meanstds ;
                %&&... nanmean(curlfp(trialwinids(mididx+pade+1:end)))>=meanbaseline-meanstds*2;
                if targburst
                    betatarginc=[betatarginc seltrials(itrial)];
                    tracebetainc=[tracebetainc; curlfp(trialwinids)];
                    tstarg=find(curlfp(trialwinids(mididx-pade:end))==max(targlfp,[],'omitnan'))+mididx-pade-1;
                    tsbetainc=[tsbetainc tstarg];
                    %just looking for maxima usually outputs first or last
                    %idx, need to look for local peak
                    %following beta burst ts tstarg 
                    tstargda=find(dasub(itrial,trialwinids(tstarg:end))==max(dasub(itrial,trialwinids(mididx:end)),[],'omitnan'))+tstarg-1;
                    diffda=find(diff(dasub(itrial,trialwinids(tstarg:end)),1,2)>0)+tstarg-1;
                    tstargdadiff=nan;
                    if ~isempty(diffda)
                    tstargdadiff=diffda(1);
                    end
                    tsbetaincdapeak=[tsbetaincdapeak tstargda];    
                    tsbetaincdadiffpeak=[tsbetaincdadiffpeak tstargdadiff];  
                elseif targsup
                    betatargdec=[betatargdec seltrials(itrial)];
                    tracebetadec=[tracebetadec; curlfp(trialwinids)]; 
                    tstarg=find(curlfp(trialwinids(mididx-pade:end))==min(targlfp,[],'omitnan'))+mididx-pade-1;
                    tsbetadec=[tsbetadec tstarg];
                end
            end
        %end        
    %store da traces groups in all lfp chs
    %get sorted targ da's for separate group sortedda
    [xx, sortid]=sort(posttargda);
    sortedtrials=seltrials(sortid);
    [xx, sortidmax]=sort(maxdachange);
    sortedtrialsmax=seltrials(sortidmax);
    xcovdata{ii}.posttargda=posttargda;       %actual values
    xcovdata{ii}.maxdachange=maxdachange;
   
    xcovdata{ii}.grouptrials.sortedda=sortedtrials;     %trial nums
    xcovdata{ii}.grouptrials.sorteddamax=sortedtrialsmax;
    xcovdata{ii}.grouptrials.targalndec=targalndec;
    xcovdata{ii}.grouptrials.targalninc=targalninc;
    xcovdata{ii}.grouptrials.sustinc=sustinc;
    xcovdata{ii}.grouptrials.sustdec=sustdec;
    xcovdata{ii}.grouptrials.sust=sust;
    xcovdata{ii}.grouptrials.inc=trialsinc;
    xcovdata{ii}.grouptrials.dec=trialsdec;

    xcovdata{ii}.grouptrialsts.inc=tsinc;
    xcovdata{ii}.grouptrialsts.dec=tsdec;
    
    xcovdata{ii}.grouptracesda.all=traceda;
    xcovdata{ii}.grouptracesda.targalndec=tracealndec;
    xcovdata{ii}.grouptracesda.targalninc=tracealninc;
    xcovdata{ii}.grouptracesda.sustinc=tracesustinc;
    xcovdata{ii}.grouptracesda.sustdec=tracesustdec;
    xcovdata{ii}.grouptracesda.sust=tracesust;
    xcovdata{ii}.grouptracesda.inc=traceinc;
    xcovdata{ii}.grouptracesda.dec=tracedec;
        
    %xcov grouped by da properties
    xcovdata{ii}.groupxcov.sortedda=...
        xcovdata{ii}.xcovda(sortid,:);
    xcovdata{ii}.groupxcov.sorteddamax=...
        xcovdata{ii}.xcovda(sortidmax,:);
    xcovdata{ii}.groupxcov.targalndec=...
        xcovdata{ii}.xcovda(ismember(seltrials,targalndec)==1,:);
    xcovdata{ii}.groupxcov.targalninc=...
        xcovdata{ii}.xcovda(ismember(seltrials,targalninc)==1,:);
    xcovdata{ii}.groupxcov.sustinc=...
        xcovdata{ii}.xcovda(ismember(seltrials,sustinc)==1,:);
    xcovdata{ii}.groupxcov.sustdec=...
        xcovdata{ii}.xcovda(ismember(seltrials,sustdec)==1,:);
    xcovdata{ii}.groupxcov.sust=...
        xcovdata{ii}.xcovda(ismember(seltrials,sust)==1,:);
    xcovdata{ii}.groupxcov.inc=...
        xcovdata{ii}.xcovda(ismember(seltrials,trialsinc)==1,:);
    xcovdata{ii}.groupxcov.dec=...
        xcovdata{ii}.xcovda(ismember(seltrials,trialsdec)==1,:);
    
    %xcoeff near 0 lag for different groups
    numpts=size(xcovdata{ii}.xcovda,2);
    midids=median(1:numpts)-2:median(1:numpts)+2;
    coefmean=nanmean(xcovdata{ii}.xcovda(:,midids),2);
    xcovdata{ii}.groupxcoef.sortedda=nanmean(xcovdata{ii}.xcovda(sortid,midids),2);
    xcovdata{ii}.groupxcoef.sorteddamax=nanmean(xcovdata{ii}.xcovda(sortidmax,midids),2);    
    xcovdata{ii}.groupxcoef.targalndec=nanmean(xcovdata{ii}.groupxcov.targalndec(:,midids),2);
    xcovdata{ii}.groupxcoef.targalninc=nanmean(xcovdata{ii}.groupxcov.targalninc(:,midids),2);
    xcovdata{ii}.groupxcoef.sustinc=nanmean(xcovdata{ii}.groupxcov.sustinc(:,midids),2);
    xcovdata{ii}.groupxcoef.sustdec=nanmean(xcovdata{ii}.groupxcov.sustdec(:,midids),2);
    xcovdata{ii}.groupxcoef.sust=nanmean(xcovdata{ii}.groupxcov.sust(:,midids),2);
    xcovdata{ii}.groupxcoef.inc=nanmean(xcovdata{ii}.groupxcov.inc(:,midids),2);
    xcovdata{ii}.groupxcoef.dec=nanmean(xcovdata{ii}.groupxcov.dec(:,midids),2);

    %peak max/min lags for different gorups
    %only those lags w/ peak above thres?
    %find(xcovdata{ii}.xcovlag(:,2)>nanmean(xcovdata{ii}.xcovlag(:,2)))
    peaklags=xcovdata{ii}.xcovlag(:,1);     %lag ts for positive max
    antilags=xcovdata{ii}.xcovlaganti(:,1); %lag ts for negative max
    xcovdata{ii}.groupxlags.sortedda=peaklags(sortid);
    xcovdata{ii}.groupxlags.sorteddamax=peaklags(sortidmax);
    xcovdata{ii}.groupxlags.targalndec=peaklags(ismember(seltrials,targalndec)==1);
    xcovdata{ii}.groupxlags.targalninc=peaklags(ismember(seltrials,targalninc)==1);
    xcovdata{ii}.groupxlags.sustinc=peaklags(ismember(seltrials,sustinc)==1);
    xcovdata{ii}.groupxlags.sustdec=peaklags(ismember(seltrials,sustdec)==1);
    xcovdata{ii}.groupxlags.sust=peaklags(ismember(seltrials,sust)==1);
    xcovdata{ii}.groupxlags.inc=peaklags(ismember(seltrials,trialsinc)==1);
    xcovdata{ii}.groupxlags.dec=peaklags(ismember(seltrials,trialsdec)==1);
    
    xcovdata{ii}.groupxlagsanti.sortedda=antilags(sortid);
    xcovdata{ii}.groupxlagsanti.sorteddamax=antilags(sortidmax);
    xcovdata{ii}.groupxlagsanti.targalndec=antilags(ismember(seltrials,targalndec)==1);
    xcovdata{ii}.groupxlagsanti.targalninc=antilags(ismember(seltrials,targalninc)==1);
    xcovdata{ii}.groupxlagsanti.sustinc=antilags(ismember(seltrials,sustinc)==1);
    xcovdata{ii}.groupxlagsanti.sustdec=antilags(ismember(seltrials,sustdec)==1);
    xcovdata{ii}.groupxlagsanti.sust=antilags(ismember(seltrials,sust)==1);
    xcovdata{ii}.groupxlagsanti.inc=antilags(ismember(seltrials,trialsinc)==1);
    xcovdata{ii}.groupxlagsanti.dec=antilags(ismember(seltrials,trialsdec)==1);
    
    %store lfp traces groups in all da chs
    xcovdata{ii}.grouptrials.betatarginc=betatarginc;
    xcovdata{ii}.grouptrials.betatargdec=betatargdec;  
    
    xcovdata{ii}.grouptrialsts.incbeta=tsbetainc;
    xcovdata{ii}.grouptrialsts.decbeta=tsbetadec;
    
    %find subsequent DA peak or initial DA incline point 
    %following beta burst in trial
    xcovdata{ii}.grouptrialsts.incbetada=tsbetaincdapeak;
    xcovdata{ii}.grouptrialsts.incbetadadiff=tsbetaincdadiffpeak;
    
    xcovdata{ii}.grouptraceslfp.betatarginc=tracebetainc;
    xcovdata{ii}.grouptraceslfp.betatargdec=tracebetadec;   
    
    xcovdata{ii}.groupxcov.betatarginc=...
        xcovdata{ii}.xcovda(ismember(seltrials,betatarginc)==1,:);
    xcovdata{ii}.groupxcov.betatargdec=...
        xcovdata{ii}.xcovda(ismember(seltrials,betatargdec)==1,:);   
    
    %xcoeff near 0 lag for different lfp groups
    xcovdata{ii}.groupxcoef.betatarginc=nanmean(xcovdata{ii}.groupxcov.betatarginc(:,midids),2);
    xcovdata{ii}.groupxcoef.betatargdec=nanmean(xcovdata{ii}.groupxcov.betatargdec(:,midids),2);
    
    %peak max/min lags for different lfp groups
    xcovdata{ii}.groupxlags.betatarginc=peaklags(ismember(seltrials,betatarginc)==1);
    xcovdata{ii}.groupxlags.betatargdec=peaklags(ismember(seltrials,betatargdec)==1);
    
    xcovdata{ii}.groupxlagsanti.betatarginc=antilags(ismember(seltrials,betatarginc)==1);
    xcovdata{ii}.groupxlagsanti.betatargdec=antilags(ismember(seltrials,betatargdec)==1);

    end               
    %end cluster signals    
            
end
%end ii lfp types

%save all xvardata for this da channel
if ispc
    set(findall(figxvarmean,'-property','FontSize'),'FontSize',10)
else
    set(findall(figxvarmean,'-property','FontSize'),'FontSize',8)
end
%save all plots fir this da channel
savename=['lfp_x_da_' plotparam.sites{ida}];
saveas(figxvar,[pathxvar{ievent} savename ],'jpeg')
if sortxvar
    saveas(figxvarsort,[pathxvar{ievent} savename 'sorted'],'jpeg')
    saveas(figxvarsortneg,[pathxvar{ievent} savename 'sortedneg'],'jpeg')
end 
saveas(figxvarmean,[pathxvar{ievent} savename 'mean'],'jpeg')

%plot trialbytrial xcov & means for clustered signals
if isfield(xcovdata{ii},'groupxcov')
groupnames=fieldnames(xcovdata{ii}.groupxcov);
numgroups=length(fieldnames(xcovdata{ii}.groupxcov));
for ifi=1:numgroups
    for iax=1:length(ax)
        cla(ax{iax});
        cla(axt{iax});
    end
    curgroup=groupnames{ifi};
    for ilf=1:size(xcovdata,2)
        xmod=getfield(xcovdata{ilf}.groupxcov,curgroup);
        tsx=xcovdata{ilf}.tsx;
        xmodmean=nanmean(xmod,1);
        nums=length(~isnan(mean(xmod,2)));
        xmodci=nanstd(xmod,[])./sqrt(nums)*1.96;    
        if nums>=5
            %only if more than 5 valid signals
        if ilf==1
            plotxtrials(axt{ilf},xcovdata{ilf},'leftlabel','groupxcov',curgroup);
        else
            plotxtrials(axt{ilf},xcovdata{ilf},'groupxcov',curgroup);
        end
        plot(ax{ilf},tsx,xmodmean,'marker','none','linewidth',2,'color',markcolor)
        xlabel(ax{ilf},'lag (s)')
        if ilf==1
            ylabel(ax{ilf},'trial averaged x-cov')
            title(ax{ilf},[eventnames{ievent} ' | ' curgroup ' | ' xcovdata{ilf}.sitename{1} ' , ' xcovdata{ilf}.sitename{2}])
        else
            title(ax{ilf},[xcovdata{ilf}.sitename{1} ' , ' xcovdata{ilf}.sitename{2}])
        end
        hold(ax{ilf},'on'); 
        plot(ax{ilf},tsx,xmodmean+xmodci,'--','linewidth',1,'color',markcolor)
        plot(ax{ilf},tsx,xmodmean-xmodci,'--','linewidth',1,'color',markcolor)
        xlim(ax{ilf},[min(tsx) max(tsx)]);
        end
        if ~isfield(xcovdata{ilf}.grouptraceslfp,curgroup)
            %if current group does not have lfp traces already stored,
            %retrieve
            %add lfp sel da traces
            blfp=xcovdata{ilf}.lfpdata;
            winids=xcovdata{ilf}.winids;
            seltrials=xcovdata{ilf}.seltrials; 
            sel=getfield(xcovdata{ilf}.grouptrials,curgroup);
            sellfptraces=[];
            tsincdabeta=[];     %only for da increases (inc) group
            for it=1:length(sel)
                curtrace=blfp(sel(it),winids(ismember(seltrials,sel(it)),:));
                sellfptraces=[sellfptraces; curtrace];
                if strcmp(curgroup,'inc')
                    %find max beta burst in trial for da inc group
                    trialwin=winids(ismember(seltrials,sel(it)),:);
                    mididx=median(1:length(trialwin));
                    tstarg=find(curtrace==max(curtrace(pade+1:end-pade),[],'omitnan'));
                    if ~isempty(tstarg)
                        if max(curtrace(pade+1:end-pade),[],'omitnan')<=xcovdata{ilf}.betabaseline+xcovdata{ilf}.betastds
                            %if not above threshold, does not count as
                            %burst
                            tstarg=nan;
                        end
                    end
                    tsincdabeta=[tsincdabeta tstarg];
                end
            end
            xcovdata{ilf}.grouptraceslfp=...
                setfield(xcovdata{ilf}.grouptraceslfp,curgroup,...
                sellfptraces);
            if strcmp(curgroup,'inc')
                xcovdata{ilf}.grouptrialsts.incdabeta=tsincdabeta;     %ts of beta burst for da inc trials
            end
        end
        if ~isfield(xcovdata{ilf}.grouptracesda,curgroup)
            %add da sel lfp traces
            das=xcovdata{ilf}.dadata;
            winids=xcovdata{ilf}.winids;
            seltrials=xcovdata{ilf}.seltrials;
            sel=getfield(xcovdata{ilf}.grouptrials,curgroup);
            seldatraces=[];
            for it=1:length(sel)                
                curda=das(sel(it),winids(ismember(seltrials,sel(it)),:));
                mididx=median(1:length(curda));         %middle of window. ie targ period
                baseids=mididx-pade-basesamples:mididx-pade;
                baseids=baseids(baseids>0);
                targids=mididx-pade:mididx+pade;     %look at midx for puase
                baseline=nanmean(curda(baseids));
                curda=curda-baseline;
                seldatraces=[seldatraces; curda];
            end
            xcovdata{ilf}.grouptracesda=...
                setfield(xcovdata{ilf}.grouptracesda,curgroup,...
                seldatraces);

        end
    end
    %save plots for each cluster group for this da channel
    savename=['lfp_x_da_' curgroup '_' plotparam.sites{ida}];
    saveas(figxvar,[pathxvar{ievent} savename ],'jpeg')
    saveas(figxvarmean,[pathxvar{ievent} savename 'mean'],'jpeg')
    
    %save all xvar data
    save([pathxvar{ievent} 'xvar_' plotparam.sites{ida}],'xcovdata');

end    
end
%end plotting clusters

end
%end ichda
end
%end eventtypes

close all;

end

