function [plotdata,da,info]=setTrialAxAvg(data,plotparam,hax,seltrials,varargin)
%01/2019 updates as with setTrialAx
%now must provide data generated from setTrialAx
%04/17/2018
%plot trialbytrial data on hax axis for given window (plotparam.win)
%data is generated in compileTrials..m script 
%is a structure with PCA components and behavior variables
%assume plotting data.da variable for signal of interest
%varargin arguments:
%markers for events (eg fix or target) overlaid on color
%output subtracted da signals
da=[];
da=data;
cla(hax);
plotnum=0;
numtrials=size(da,1);
if isempty(seltrials)
    seltrials=1:numtrials;
end
numtrials=size(da(seltrials,:),1);
argnum=1;
markers={};
%markers{1} = fix cue appearance
%markers{2} = target cue apperance
%markers{3} = eye on fix cue
iscsc=0;        %not csc data, ie fscv time stamps
cscname=[];
lfp2=0;
logscale=0;
avg=0;
zz=0;
phys=0;
pulse=0;
eye=0;
blink=0;
noplot=0;
sitename={};
rclfp=[];
basetrial=0;
spect=0;        %plot avg spectrum
filtband='';       %default filter band is beta low
filtlick=[0 10];        %freq bands of lick filter
filtbetal=[12 18];        %low beta%    changed from 13 to 12 07/07/18
filtbetah=[18 33];        %high beta
filttheta=[6 10];           %theta
filtdelta=[0.5 4];           %delta
filtgammal=[40 58]; 
filtbeta=[13 30];
    bandfilt=[0 0];
ulabel=[];
ylimit=[];
triallabel=[];
    envwidth=0;
plotnum=0;
islfp=0;
iscsc=0;
samplerate=[];
scales=[];
means=[];
ts=[];
if isfield(plotparam,'triallabel')
    triallabel=plotparam.triallabel;
end
while argnum <= length(varargin)
    switch varargin{argnum}        
        case 'cscale'
            argnum=argnum+1;
            plotparam.cscale=varargin{argnum};
        case 'plotnum'
            argnum=argnum+1;
            plotnum=varargin{argnum};
        case 'filt'
            %get filt band (filter here 08/04/2018)
            argnum=argnum+1;        %fscv ch site name
            filtband=varargin{argnum};     %flag to plot lfp2
       case 'bandfilt'
            %provide actual freq limits
            argnum=argnum+1;
            bandfilt=varargin{argnum};
        case 'fbands'
            %provide pass bands for each site
            %sitename must be defined by this point
            argnum=argnum+1;
            fbands=varargin{argnum};
            bandid=find(strcmp(sitename,fbands)==1);
            bandfilt=fbands{bandid(1)+1};
        case 'lfp2'
            %plot 2nd filtered band of lfps
            lfp2=1;     %flag to plot lfp2
        case 'lfp3'
            lfp2=2;
        case 'log'
            logscale=1;
        case 'avg'
            avg=1;      %plot mean
        case 'zscore'
            zz=1;      %plot zscore
        case 'spectrum'
            spect=1;
        case 'label'
            argnum=argnum+1;
            ulabel=varargin{argnum};
       case 'sitename'
            argnum=argnum+1;        %fscv ch site name
            sitename=varargin{argnum};
       case 'smoothwin'
            %provide smooth width for enveloped signal
            argnum=argnum+1;
            envwidth=varargin{argnum};  %in seconds
        case 'ylimit'
            argnum=argnum+1;
            ylimit=varargin{argnum};
        case 'fixedbasetrial'
            %uniform baseline for all trials that is all trials
            basetrial=1;
       case 'rclfp'
            %lfp already reconverted and provided here filtered
            argnum=argnum+1;
            rclfp=varargin{argnum};
        case 'plotnum'
            argnum=argnum+1;
            plotnum=varargin{argnum};
        case 'info'
            %user provides sampler ate & ts
            argnum=argnum+1;
            infos=varargin{argnum};
            samplerate=infos.samplespersec;     %for extracted phys signals
            ts=infos.ts;
        case 'scales'
            argnum=argnum+1;
            scalesmeans=varargin{argnum};            
            if ~isempty(scalesmeans)
            scales=scalesmeans(1);
            if length(scalesmeans)>1
            means=scalesmeans(2);
            end
            end
        otherwise           
            %name provided, csc signal
            iscsc=1;
            cscname=varargin{argnum};
            cscid=find(ismember(plotparam.cscNames,cscname));
            if ~(contains(cscname,'eye') || contains(cscname,'lick') || ...
                    contains(cscname,'pulse') || contains(cscname,'blink'))
                islfp=1;
            end
            %samplerate=plotparam.ratelfp;
            sitename=cscname;
    end
    argnum = argnum + 1;
end
if basetrial==1
    baseline='alltrial';
end
win=plotparam.win;      %in fscv samples
rewardidx=plotparam.alignidx;
idsfixeye=plotparam.samplesfixeye{1};
idsfix=plotparam.samplesfix{1};
idstarg=plotparam.samplestarg{1};
idstargeye=plotparam.samplestargeye{1};

vertplotwidth=plotparam.vertplotwidth;
vertplotwidth2=plotparam.vertplotwidth2;
colortrial=plotparam.colormap;
colorm=plotparam.markercolors;
samplespersec=plotparam.samplespersec;
chnum=plotparam.chnum;
trialtype=[plotparam.trialtype ' | ' triallabel];

if ~isempty(samplerate)
    samplespersec=ceil(samplerate);
    if isfield(data,'relts')
    tslfp=data.relts;           %original relative ts's for csc
    dcrate=ceil(length(tslfp)./length(rates));      %downconvert rate for ts 
    ts=tslfp(1:dcrate:end);
    end
    winneg=round((plotparam.alignidx-win(1))./plotparam.samplespersec*samplespersec);
    winpos=round((win(end)-plotparam.alignidx)./plotparam.samplespersec*samplespersec);
    rewardidx=round(plotparam.alignidx.*samplespersec./plotparam.samplespersec);
    winnew=rewardidx-winneg:rewardidx+winpos; 
    win=winnew;             %replace fscv domain window for plotting ids
    %convert marker id's to new sample rate
    idsfixeye=round(plotparam.tfixeye.*samplespersec);
    idsfix=round(plotparam.tfix.*samplespersec);
    idstarg=round(plotparam.ttarg.*samplespersec);
    idstargeye=round(plotparam.ttargeye.*samplespersec);
end


idsrew=repmat(rewardidx,1,size(da,1));
alnevt=plotparam.alnevt;
idsaln=[];
switch alnevt
    case 'fix'
        idsaln=idsfix;
    case 'fixeye'
        idsaln=idsfixeye;
    case 'targ'
        idsaln=idstarg;
    case 'targeye'
        idsaln=idstargeye;
    case 'outcome'
        %do not shift
        idsaln=idsrew;
end


info.alnevt=alnevt;
info.idsaln=idsaln;
alnidx=rewardidx;
pade=round(.2*samplespersec);       %padding before subseq event
basesamples=round(.3*samplespersec);        %.3 s for baseline average;
basepad=1;              %+/- .1 s around aling idx for baseline signal
outdata=da;         %output signal, downconverted LFP, before shifting
rewext=5;       %extend baseline for reward baseline if too nan
info.alnevt=alnevt;
info.idsaln=idsaln;
alnidx=rewardidx;
%shift signals according to alnevt idsaln
badtrls=ones(1,size(da,1));
idsfixeye2=nan(1,length(idsaln));
idsfix2=nan(1,length(idsaln));
idstarg2=nan(1,length(idsaln));
idstargeye2=nan(1,length(idsaln));
idsrew2=nan(1,length(idsaln));

if ~isempty(idsaln) && ~strcmp(alnevt,'outcome')
    meanaln=nanmean(idsaln);
    stdaln=nanstd(idsaln);
    badtrls=(idsaln<meanaln-stdaln*20 | idsaln>meanaln+stdaln*20 | isnan(idsaln));
    da(badtrls,:)=nan;          %set signal to nan for bad trials
   % alnidx=nanmean(idsaln(~badtrls));       %new alignment idx for all trials
    daaln=nan(size(da,1),size(da,2));
    for itrial=1:size(da,1)
        if badtrls(itrial)==0 && ~isempty(idsaln(itrial)) && ~isnan(idsaln(itrial))
        %numshift to align to meanalnwin
        shiftpos=round((alnidx-idsaln(itrial)));
        %shift da signal to indicated alignment idx
        daaln(itrial,:)=circshift(da(itrial,:),shiftpos,2);
        %shift all markers also based on alignment idx specified
        idsfixeye2(itrial)=idsfixeye(itrial)+shiftpos;
        idsfix2(itrial)=idsfix(itrial)+shiftpos;
        idstarg2(itrial)=idstarg(itrial)+shiftpos;
        idstargeye2(itrial)=idstargeye(itrial)+shiftpos;
        idsrew2(itrial)=rewardidx+shiftpos;
        if isnan(nanmean(da(itrial,:)))
            %if all signal is nan, designate bad trial so don't process
            badtrls(itrial)=1;
        end
        end
    end
    %replace original markers with newly shifted events
    idsfixeye=idsfixeye2;
    idsfix=idsfix2;
    idstarg=idstarg2;
    idstargeye=idstargeye2;
    idsrew=idsrew2;
    %replace original data with shifted data
    da=daaln;
    %replace seltrials without badtrls
    %{
    goodtrials=find(badtrls==0);
    seltrials2=seltrials(ismember(seltrials,goodtrials));
    seltrials=seltrials2;
    numtrials=length(seltrials);
    %}
end
%baseline subtraction for DA fscv data only
if ~iscsc || strcmp(cscname,'eyedist')
    %since concistent minima occurs at eye on targ/fix, use aln event +/-   
    %rather than prior
    %eye distance also needs to be baselien subtracted since cumulative of
    %signal across initial windows when exported signals
    %baseline immediately prior to aln idx w/ pad
    baseids=alnidx-basepad:alnidx+basepad;   
    if strcmp(alnevt,'outcome') && ~strcmp(cscname,'eyedist')
        %da signal around rew, nan extends a lot so we have to extend base
        %window
        baseids=alnidx-basepad-rewext:alnidx+basepad+rewext;
    end
    baseline=nanmean(da(:,baseids),2);
    %baseline subtracted DA
    dasub=da-baseline;
    da=dasub;
end

interval=plotparam.interval;        %interval xtick marks def 2.5s
cscale=plotparam.cscale;
cminshow=plotparam.cminshow;

xticks=1:interval*samplespersec:size(da(:,win),2);
mintime=round((win(1)-alnidx))/samplespersec;
maxtime=round((win(end)-alnidx))/samplespersec;
xticklabels=mintime:interval:maxtime;
xticklabels=round(xticklabels.*10)./10;
xticklabels=num2str(xticklabels');

    
        
pdata=[];
pdata=da(seltrials,win(1):win(end));
plotdata.data=pdata;
numshift=zeros(size(da,1),1);
if logscale==1
   pdata=10*log10(pdata);
end
%remove any trials whose mean is above threshold
%with significant outlier signals
pstds=nanstd(pdata,0,1);
numdatapertimepoint=[];     %# of valid data points for each time point (ie not nan)
pci=[];     %confidence intervals calculated based on non nan samples
validtrials=[];
for ii=1:size(pdata,2)
    numdatapertimepoint(ii)=size(find(~isnan(pdata(:,ii))==1),1);
    if numdatapertimepoint(ii)<=.1*size(pdata,1)
        %if less than 5% of all trials then exclude this time point in
        %average, signal not meaningfull
        pdata(:,ii)=nan;
    end
    pci(ii)=pstds(ii)/sqrt(numdatapertimepoint(ii))*1.96;        %95%ci
end
for it=1:size(pdata,1)
    if length(~isnan(pdata(it,:)))/size(pdata,2)>.75
    %if 75% of data or great is not nan, then valid
    validtrials=[validtrials ...
        seltrials(it)];
    end
end
means1=nanmean(pdata,1);
means2=nanmean(means1);
thresoutlier=means2*3;
meanstrials=nanmean(pdata,2);
outtrial=find(abs(meanstrials)>abs(thresoutlier));
if iscsc    
pdata(outtrial,:)=nan;
end
%remove glitches based on max mean or std
maxpeak=means2*20;
%maxpeak=nanmean(stddata)*10;
pcis=[];
stddata=nanstd(pdata,[],1);
if iscsc
    pdata(pdata>abs(maxpeak))=nan;
else
    maxpeak=nanmean(stddata)*10;
    pdata(pdata>abs(maxpeak))=nan;
end
if avg==0 && zz==0
   % pdata=pdata;
elseif avg==1 && zz==0
    pdata=nanmean(pdata,1);
    pcis(1,:)=pdata+pci;
    pcis(2,:)=pdata-pci;
end
if zz==1
    %pdata1=zscore(pdata,0,2);       %generates all nan if any nan values
    pdata=zscorenan(pdata,2);        
end
if zz==1 && avg==1
    pdata=nanmean(pdata,1);
        %pcis(1,:)=pdata+pci;
   % pcis(2,:)=pdata-pci;
end

if isempty(scales) || isempty(means) || strcmp(cscname,'eyedist')
    stds=nanstd(pdata,[],2);
means=nanmean(pdata,2);
    scales=nanmean(stds);
    means=nanmean(means);
    if strcmp(sitename,'eye')
        %eye d, must remove blinks from getting scales
        stds=nanstd(eyedata,[],2);
        means=nanmean(eyedata,2);
        means=nanmean(means);
        scales=nanmean(stds);
    end
end
info.scales=[scales means];

plot(hax, pdata');
hold(hax,'on')
%plot(hax,pcis','color',[0.7 0.7 0.7],'linewidth',0.5,'linestyle','--')
%plot error bars
if avg==1 && zz==0
    plot(hax, pdata','color',[0 0 0],'linewidth',1);
    plot(hax,pcis','color',[0.7 0.7 0.7],'linewidth',0.5,'linestyle','-')
end
origpos=getpixelposition(hax);      %get original position 
plotdata.mean=pdata;
plotdata.ci=pcis;
plotdata.stds=stddata;
plotdata.trialsvalid=validtrials;
%plot event markers
if ~isempty(ylimit)
cscale=ylimit;
end
if ~isnan(cscale(1)) 
    ylim(hax,cscale); 
else
    cscale(1)=-scales*3;
    cscale(2)=scales*4;
    if avg
    cscale(2)=means+scales;
    cscale(1)=means-scales;
    end
    if iscsc
            cscale(2)=means+scales*10;
            cscale(1)=means-scales*2;
        if means<0 && strcmp(cscname,'eyed')
            cscale(1)=means-scales*3;
            cscale(2)=means+scales*4;
        end
        if avg
            cscale(2)=means+scales*1.5;
            cscale(1)=means-scales*1;
            if means<0 && (strcmp(cscname,'eyed') || strcmp(cscname,'eyex'))
            cscale(1)=means-scales*2;
            cscale(2)=means+scales*3;
            end
            if ~isempty(bandfilt)
                if bandfilt(1)>30
                    %gamma
                    cscale(2)=means+scales*2.5;
                    cscale(1)=means-scales*1.5;
                end
            end
        end
    end    
    if ~isempty(cscale) && cscale(2)>cscale(1) && ~any(isnan(cscale))
    ylim(hax,cscale)
    end
end
info.cscale=cscale;
 %plot event markers
line(hax,[nanmean(idsfix(seltrials)-win(1)+mean(numshift(seltrials))) nanmean(idsfix(seltrials)-win(1)+mean(numshift(seltrials)))],[cscale(1) cscale(2)],'Color',colorm(1,:)*.35,'LineWidth',.5,'LineStyle','--')
line(hax,[nanmean(idsfixeye(seltrials)-win(1)+mean(numshift(seltrials))) nanmean(idsfixeye(seltrials)-win(1)+mean(numshift(seltrials)))],[cscale(1) cscale(2)],'Color',colorm(2,:)*.45,'LineWidth',0.5,'LineStyle','--')
if ~strcmp(plotparam.trialtype,'fixbreak')
    if ~strcmp(plotparam.alnevt,'targ')
        line(hax,[nanmean(idstarg(seltrials)-win(1)+mean(numshift(seltrials))) nanmean(idstarg(seltrials)-win(1)+mean(numshift(seltrials)))],[cscale(1) cscale(2)],'Color',colorm(3,:)*.35,'LineWidth',.5,'LineStyle','--')
    end
    if ~strcmp(plotparam.alnevt,'targeye')
        line(hax,[nanmean(idstargeye(seltrials)-win(1)+mean(numshift(seltrials))) nanmean(idstargeye(seltrials)-win(1)+mean(numshift(seltrials)))],[cscale(1) cscale(2)],'Color',colorm(4,:)*.35,'LineWidth',.5,'LineStyle','--')
    end
end
if ~strcmp(plotparam.alnevt,'outcome')
   %not aligned to outcome, mark
line(hax,[nanmean(idsrew(seltrials)-win(1)+mean(numshift(seltrials))) nanmean(idsrew(seltrials)-win(1)+mean(numshift(seltrials)))],[cscale(1) cscale(2)],'Color',colorm(4,:)*.35,'LineWidth',.5,'LineStyle','--')
end

%organize plot
set(hax,'XTick',xticks)
set(hax,'xticklabel',xticklabels)
set(hax,'tickdir','out','box','off')
xlabel(hax,'time (s)')
xlim(hax,[0 win(end)-win(1)])
clabel='[DA] (nM)';
if zz==1
    clabel='z score da';
end

titlename=[];
if ~iscsc
    titlename=[trialtype ' DA, ch ' num2str(chnum)];    
       if ~isempty(sitename)
            %titlename=[trialtype ' DA, ' sitename{:} ' | ' alnevt];
            titlename=[trialtype ' | ' alnevt ' | ' sitename{:}];
        end
    %title(hax,[trialtype ' trials, DA, ch ' num2str(chnum)])
else
    %title(hax,[trialtype ' trials, lfp beta, ' cscname])
    %band=bandfilt;
    %title(hax,['lfp ' num2str(band(1)) '-' num2str(band(2)) ' Hz, ' cscname])
    titlename=[trialtype ' lfp ' num2str(bandfilt(1)) '-' num2str(bandfilt(2)) ' Hz, ' cscname ' | ' alnevt];
    clabel='\beta-lfp (\muV^2)';
end
if ~isempty(ulabel)
    titlename=[titlename ', ' ulabel];
end
%if plotnum<=2
if avg
    title(hax,titlename)
end

hold(hax,'off')
%plotdata=pdata;

end




