function xcovdata=plotxvartrials(data,eventmarks,plotparam,varargin)
%1/26/2019
%UPDATED
%plotxvar
%data.lfp is betalfp data each channel in {ilfp}
%data.da is da each channel in {ida}
%ax{1} is axes for trial by trial plot, ax{2} for mean
%eventmarks{1} corresponds to eventnames{1} label, defines times for
%windows for xvar in seconds
%12/06/2018 took out all shuffle functions, saved old ver as
%plotxvartrialsold2.m for reference
%get sessnum
sessidfind=strfind(plotparam.pathname,'analysis');
sessfind=strfind(plotparam.pathname,filesep);
sessid=plotparam.pathname(sessfind(end-1)-3:sessfind(end-1)-1);
if ~isnumeric(str2num(sessid))
    %delete one character from front (assuming maximum 3 digits)
    sessid=sessid(2:end);
    if ~isnumeric(str2num(sessid))
        sessid=sessid(2:end);
    end
end
sessnum=str2num(sessid);
alnts=30;           %align to 30s for all data
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
        xwin=[-1 2];
        isbeh=0;
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'xwin'
            argnum=argnum+1;
            xwin=varargin{argnum};
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
        case 'beh'
            isbeh=1;
    end
    argnum=argnum+1;
end
xwintemp=xwin;      %store not for outcome
origseltrials=seltrials;
lfpchs=plotparam.lfpchs;
%get lfp ch groups p & c
fcvlfpgroups={};
maxplots=5;
if ~isbeh
pgroup=find(contains(lfpchs,'p')==1);
cgroup=find(contains(lfpchs,'c')==1);
maxplots=max(length(pgroup),length(cgroup));
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
else
    lfpchs=data.sitesbeh;
    maxplots=length(lfpchs);
    for ii=1:length(fscvchs)
    ich=fscvchs(ii);
    fcvlfpgroups{ii}.fsite=ich;
    fcvlfpgroups{ii}.lsites=1:maxplots;
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
if maxplots>6
    sizep=[round(figsize(1)/(6)-maxplots*pad/10) round(figsize(2)/3)];
end
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
        if maxplots>6
            if ilfp==1
               set(figxvarmean,'Position',[20,50,figsize(1),figsize(2)]);
            end
                ax{ilfp}=subplot(2,6,ilfp);
    set(ax{ilfp},'units','pixels');    
            
            %{
               if ilfp<=6
       set(ax{ilfp},'position',[pad + sizep(1)*(ilfp-1)+pad*(ilfp), ...
            figsize(2)-sizep(2)-pad*2,sizep(1),sizep(2)]);
               else
      set(ax{ilfp},'position',[pad + sizep(1)*(ilfp-7)+pad*(ilfp-6), ...
            figsize(2)-sizep(2)*2-pad*2,sizep(1),sizep(2)]);       
               end
            %}
        else
       ax{ilfp}=subplot(1,maxplots,ilfp);
    set(ax{ilfp},'units','pixels');    
            set(ax{ilfp},'position',[pad + sizep(1)*(ilfp-1)+pad*(ilfp), ...
            pad*2,sizep(1),sizep(2)]);
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
rewardidx=plotparam.alignidx;



for ievent=1:length(eventmarks)
    xwin=xwintemp;
    %interfix, intertarg, different event windows..
    evname=eventmarks{ievent}.name;
    evts=eventmarks{ievent}.tcsc;
    evidsda=eventmarks{ievent}.da;
pathxvar{ievent}=[plotparam.pathname 'xvar_' evname '_' ...
    vartype label filesep];
if isbeh
    pathxvar{ievent}=[plotparam.pathname 'xvar_beh_' evname '_' ...
    vartype label filesep];
end
if ~isdir(pathxvar{ievent})
    mkdir(pathxvar{ievent});
end


for ichda=1:length(fscvchs)
    xcovdata={};        %saved in separate file for each da ch
    ida=fscvchs(ichda);
        seltrials=data.trials{ida};
    da=data.da{ida};
    fsda=data.rates(2);           %rates 2 is the da rate
%SHIFT SIGNALS FIRST PRIOR TO SUBMITTING TO XVARDATA BASED ON INDICATED ALN
pade=round(.2*fsda);       %padding before subseq event
basesamples=round(.3*fsda);        %.3 s for baseline average;
basepad=1;              %+/- .1 s around aling idx for baseline signal
rewext=50;       %extend baseline for reward baseline if too nan%changed from 10 to 50 04/27/2019
alnidxda=round(alnts*fsda);       %rewardidx, ie middle of data recording file align new idx here
idsaln=evidsda;     %fscv sample ids for events for da signal
%if not aligning to outcome, shift signals to indicated event mark
if ~contains(evname,'outcome')
    meanaln=nanmean(idsaln);
    stdaln=nanstd(idsaln);
    badtrls=(idsaln<meanaln-stdaln*20 | idsaln>meanaln+stdaln*20 | isnan(idsaln));
    da(badtrls,:)=nan;          %set signal to nan for bad trials
    daaln=[];
    for itrial=1:size(da,1)
        if ~isnan(idsaln(itrial))
        shiftpos=round((alnidxda-idsaln(itrial))); %numshift to alnidx
        daaln(itrial,:)=circshift(da(itrial,:),shiftpos,2);
        else
            daaln(itrial,:)=nan(1,size(da,2));
        end
        if isnan(nanmean(da(itrial,:)))            
            badtrls(itrial)=1;
        end
    end        
    da=daaln; %replace original data with shifted data
end
%baseline subtract data
baseids=alnidxda-basepad:alnidxda+basepad;   
if contains(evname,'outcome') 
    %use pretarg as baseline 04/27/2019
    %da signal around rew, nan extends a lot so  extend base window
    baseids=alnidxda-basepad-rewext:alnidxda+basepad;
    %extend xwin;
    xwin=[-1 5];
end
if contains(evname,'targ')
    %make window longer according to sess #
    xwin=[-1 3.5];
    if sessnum<66
        xwin=[-1 2];
    end
end
baseline=nanmean(da(:,baseids),2);
dasub=da-baseline; %baseline subtracted DA
da=dasub;           %replace original da data with subtracted data

%selected trials based on da characteristics in window
sustinc=[];         %sustained increase without local valley
sustdec=[];     %sustained decrease without local rise
sust=[];
meansda=[];
maxdas=[];
trialsinc=[];       %alninc and sustinc
trialsdec=[];
tsinc=[];   %delta timestamp of maxima for alninc and sustinc
tsdec=[];    
sortedtrials=[];
sortedtrialsmax=[];
stdall=nan;
for iax=1:length(ax)
    cla(ax{iax});
    cla(axt{iax});
end

for ii=1:length(fcvlfpgroups{ichda}.lsites)            
    %cross var with defined lfp chs for p / c groups
    cla(axt{ii});
    ilfp=fcvlfpgroups{ichda}.lsites(ii);
    lfp=data.lfp{ilfp};
    fs=data.rates(1);           %rates 1 is the lfp rate
    if isbeh
        fs=data.ratesbeh(ilfp);
    end
%SHIFT SIGNALS FIRST PRIOR TO SUBMITTING TO XVARDATA BASED ON INDICATED ALN
pade=round(.2*fs);       %padding before subseq event
basesamples=round(.3*fs);        %.3 s for baseline average;
alnidx=round(alnts*fs);       %rewardidx, ie middle of data recording file align new idx here
idsaln=round(evts.*fs);     %fscv sample ids for events for da signal
basepad=round(.1*fs);
%if not aligning to outcome, shift signals to indicated event mark
if ~contains(evname,'outcome')
    meanaln=nanmean(idsaln);
    stdaln=nanstd(idsaln);
    badtrls=(idsaln<meanaln-stdaln*20 | idsaln>meanaln+stdaln*20 | isnan(idsaln));
    lfp(badtrls,:)=nan;          %set signal to nan for bad trials
    daaln=[];
    for itrial=1:size(lfp,1)
        if ~isnan(idsaln(itrial))
            shiftpos=round((alnidx-idsaln(itrial)));
            daaln(itrial,:)=circshift(lfp(itrial,:),shiftpos,2);
        else
            daaln(itrial,:)=nan(1,size(lfp,2));
        end
        if isnan(nanmean(lfp(itrial,:)))
            badtrls(itrial)=1;
        end
    end
    if strcmp(lfpchs{ilfp},'eyedist')
        %background sub
        baseids=alnidx-pade-basepad:alnidx-basepad;   
        baseline=nanmean(daaln(:,baseids),2);
        dasub=daaln-baseline; %baseline subtracted DA
        daaln=dasub;           %replace original da data with subtracted data
    end
    lfp=daaln;    %replace original data with shifted data
end
    xdata={lfp da};
    sitenames={lfpchs{ilfp} plotparam.sites{ida}};
    plotparam.alnts=alnts;
    xrates=data.rates;
    if isbeh
        xrates=[fs data.rates(2)];      %beh rate & da rate
    end
    %get xvar's of data provided in xdata
    if ii==1
        xcovdata{ii}=xvardata(xdata,xrates,plotparam,axt{ii},seltrials,...
            'sitename',sitenames,'fbands',plotparam.fbands,'type',vartype,'eventnames',evname,'firstplot','xwin',xwin);
    else
        xcovdata{ii}=xvardata(xdata,xrates,plotparam,axt{ii},seltrials,...
            'sitename',sitenames,'fbands',plotparam.fbands,'type',vartype,'eventnames',evname,'xwin',xwin);
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
        title(ax{ii},[evname ' | ' lfpchs{ilfp} ' , ' plotparam.sites{ida}])
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
        das=xcovdata{ii}.dadata;            %all trials shifted data, downconverted
        winids=xcovdata{ii}.winids;     %only seltrials
        seltrials=xcovdata{ii}.seltrials;      
        basepad=1;          %.1 s +/- for baseline sub
        rewext=4;
        %get std of signal baseline in first half of each signal win
        curda=das(seltrials,winids);
       % mididx=round(median(1:length(winids)));
        mididx=alnidxda-winids(1)+1;
        basedastds=nanstd(curda(:,1:round(median(1:size(curda,2)))),[],2);
        stdall=nanmean(basedastds);     %mean of all stds
        %baseline subtract data
        baseids=mididx-basepad:mididx+basepad;   
        if contains(evname,'outcome') 
            baseids=mididx-basepad-rewext:mididx+basepad;
        end
        baseline=nanmean(curda(:,baseids),2);
        dasub=curda-baseline; %baseline subtracted DA
        curda=dasub;           %replace original da data with subtracted data
        [maxdas,maxdasts]=max(curda(:,mididx:end),[],2,'omitnan');     %max post aln
        abovethres=find(maxdas>stdall);
        trialsinc=seltrials(abovethres);
        tsinc=maxdasts+mididx-1;
        [mindas,mindasts]=min(curda(:,mididx:end),[],2,'omitnan');
        belowthres=find(mindas<-stdall);
        trialsdec=seltrials(belowthres);
        tsdec=mindasts+mididx-1;
        rebounds=intersect(abovethres,belowthres);  %trials that show min/max are rebound usu
        %REBOUNDS NOT USED YET
        decsonly=setdiff(belowthres,rebounds);  %dec da trials without rebound only
        meansda=nanmean(curda(:,mididx:end),2);   
        sustpos=find(meansda>stdall*.75);
        sustinc=seltrials(sustpos);
        sustneg=find(meansda<-stdall);
        sustdec=seltrials(sustneg);             
        sust=find(meansda<stdall*.75 & meansda>-stdall*.75);
            %get sorted targ da's for separate group sortedda
        [xx, sortid]=sort(meansda);
        sortedtrials=seltrials(sortid);
        [xx, sortidmax]=sort(maxdas);
        sortedtrialsmax=seltrials(sortidmax);
        
    xcovdata{ii}.daavgstds=stdall;%da baseline        
    seltrials=xcovdata{ii}.seltrials;      
    xcovdata{ii}.posttargda=meansda;       %actual values
    xcovdata{ii}.maxdachange=maxdas;   
    xcovdata{ii}.grouptrials.sortedda=sortedtrials;     %trial nums
    xcovdata{ii}.grouptrials.sorteddamax=sortedtrialsmax;
    xcovdata{ii}.grouptrials.sustinc=sustinc;
    xcovdata{ii}.grouptrials.sustdec=sustdec;
    xcovdata{ii}.grouptrials.sust=sust;
    xcovdata{ii}.grouptrials.inc=trialsinc;
    xcovdata{ii}.grouptrials.dec=trialsdec;
    xcovdata{ii}.grouptrialsts.inc=tsinc;
    xcovdata{ii}.grouptrialsts.dec=tsdec;
        
    %xcov grouped by da properties
    xcovdata{ii}.groupxcov.sortedda=...
        xcovdata{ii}.xcovda(sortid,:);
    xcovdata{ii}.groupxcov.sorteddamax=...
        xcovdata{ii}.xcovda(sortidmax,:);
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
    xcovdata{ii}.groupxcoef.sortedda=nanmean(xcovdata{ii}.groupxcov.sortedda(:,midids),2);
    xcovdata{ii}.groupxcoef.sorteddamax=nanmean(xcovdata{ii}.groupxcov.sorteddamax(:,midids),2);        
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
    xcovdata{ii}.groupxlags.sustinc=peaklags(ismember(seltrials,sustinc)==1);
    xcovdata{ii}.groupxlags.sustdec=peaklags(ismember(seltrials,sustdec)==1);
    xcovdata{ii}.groupxlags.sust=peaklags(ismember(seltrials,sust)==1);
    xcovdata{ii}.groupxlags.inc=peaklags(ismember(seltrials,trialsinc)==1);
    xcovdata{ii}.groupxlags.dec=peaklags(ismember(seltrials,trialsdec)==1);
    
    xcovdata{ii}.groupxlagsanti.sortedda=antilags(sortid);
    xcovdata{ii}.groupxlagsanti.sorteddamax=antilags(sortidmax);
    xcovdata{ii}.groupxlagsanti.sustinc=antilags(ismember(seltrials,sustinc)==1);
    xcovdata{ii}.groupxlagsanti.sustdec=antilags(ismember(seltrials,sustdec)==1);
    xcovdata{ii}.groupxlagsanti.sust=antilags(ismember(seltrials,sust)==1);
    xcovdata{ii}.groupxlagsanti.inc=antilags(ismember(seltrials,trialsinc)==1);
    xcovdata{ii}.groupxlagsanti.dec=antilags(ismember(seltrials,trialsdec)==1);        
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
            title(ax{ilf},[evname ' | ' curgroup ' | ' xcovdata{ilf}.sitename{1} ' , ' xcovdata{ilf}.sitename{2}])
        else
            title(ax{ilf},[xcovdata{ilf}.sitename{1} ' , ' xcovdata{ilf}.sitename{2}])
        end
        hold(ax{ilf},'on'); 
        plot(ax{ilf},tsx,xmodmean+xmodci,'--','linewidth',1,'color',markcolor)
        plot(ax{ilf},tsx,xmodmean-xmodci,'--','linewidth',1,'color',markcolor)
        xlim(ax{ilf},[min(tsx) max(tsx)]);
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

