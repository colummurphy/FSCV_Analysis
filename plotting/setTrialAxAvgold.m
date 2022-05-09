function [plotdata,da]=setTrialAxAvg(data,plotparam,hax,seltrials,varargin)
%04/17/2018
%plot trialbytrial data on hax axis for given window (plotparam.win)
%data is generated in compileTrials..m script 
%is a structure with PCA components and behavior variables
%assume plotting data.da variable for signal of interest
%varargin arguments:
%markers for events (eg fix or target) overlaid on color
%output subtracted da signals
da=[];
if isfield(data,'da')
da=data.da;
else
    da=data.lfp{1};
end
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
if isfield(plotparam,'triallabel')
    triallabel=plotparam.triallabel;
end
while argnum <= length(varargin)
    if isnumeric(varargin{argnum})
        markers{argnum}=varargin{argnum}; %store marker events
    else
        %if string, then specifies what field (ie. da, ph, or
        %lfp cl1 p1, eyedetc. to use in data
        switch varargin{argnum}
            case 'da'
                da=getfield(data,varargin{argnum});
            case 'ph'
                da=getfield(data,varargin{argnum});
            case 'm'
                da=getfield(data,varargin{argnum});
             case 'bg'
                da=getfield(data,varargin{argnum});
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
            otherwise
                cscid=find(ismember(data.cscNames,varargin{argnum}));
               % da=10*log10(data.lfp{cscid});
                if isequal(varargin{argnum},'eyed') || isequal(varargin{argnum},'blink')
                   cscid=find(ismember(data.cscNames,'eyed'));
                    da=-data.lfp{cscid};
                    phys=1;
                    eye=1;
                    if isequal(varargin{argnum},'blink')
                        blink=1;
                    end
                else
                   if isempty(cscid)
                       return;
                   end
                   da=(data.lfp{cscid}).*1e9;
               end
               if isequal(varargin{argnum},'pulse')
                   pulse=1;
               end

                plotparam.cscale(1)=nan;
                iscsc=1;
                cscname=varargin{argnum};
                sitename=cscname;
        end
    end
    argnum = argnum + 1;
end
win=plotparam.win;      %in fscv samples


rewardidx=plotparam.alignidx;
shiftdata=[];
if isfield(plotparam,'shift')
    %shift data to targ or fix marker
    shiftdata=plotparam.shift;
end
vertplotwidth=plotparam.vertplotwidth;
vertplotwidth2=plotparam.vertplotwidth2;
colortrial=plotparam.colormap;
colorm=plotparam.markercolors;
win=plotparam.win;      %in fscv samples
samplespersec=plotparam.samplespersec;
chnum=plotparam.chnum;
trialtype=[plotparam.trialtype ' | ' triallabel];

if isempty(da)
    return;
end

if iscsc
    samplespersec=ceil(data.ratelfp);
    tslfp=data.relts;
    ts1=win(1)./plotparam.samplespersec;
    ts2=win(end)./plotparam.samplespersec;
    tslfprounded=round(tslfp.*1000)/1000;
    win1id=find(tslfprounded==ts1);
    win2id=find(tslfprounded==ts2);
    win=win1id:win2id;          %get window of samples in nlx ts id's
    %convert marker id's to new sample rate
    markers{1}=markers{1}.*samplespersec./plotparam.samplespersec;
    markers{2}=markers{2}.*samplespersec./plotparam.samplespersec;
    markers{3}=markers{3}.*samplespersec./plotparam.samplespersec;
    rewardidx=rewardidx.*samplespersec./plotparam.samplespersec;
if isempty(rclfp)
        %reconverted lfp not provided, then reconvert here   
    switch filtband
        case 'betal'
            bandfilt=filtbetal;            
        case 'betah'
            bandfilt=filtbetah;
        case 'theta'
            bandfilt=filttheta;
        case 'delta'
            bandfilt=filtdelta;
        case 'gamma'
            bandfilt=filtgammal;
        case 'beta'
            bandfilt=filtbeta;
    end
                info.filtfreq=bandfilt;
        
    if sum(bandfilt)>0
        %if another filter band specified 05/03/2018
        %filter at band selected
        for ii=1:size(da,1)
            datemp=filterLFP(da(ii,:),samplespersec,bandfilt);
            if envwidth==0
                winlength=round(samplespersec*.5/mean(bandfilt));
            else
                winlength=envwidth.*samplespersec; %user provided envelop width
            end
            datemp=datemp.^2;   %get power V^2
            da(ii,:)=smoothwin(datemp,winlength);   %smoothing 
        end
    end
    else
        da=rclfp;
    end
    
end


origeye=[];
if eye==1
    eyedata=[];
    origeye=da;
    for itrial=1:size(da,1)
        eyedata(itrial,:)=deglitchnanamp(da(itrial,:),2.5e-3,30);
    end

    
    da=eyedata;
    %take out infiinite slopes (ie blinks that can amplify avg)
end

if blink==1
    da=origeye;
        %calculate rate of blinks
    eyez=zscore(da,0,2);
    movingwindow=2;     %5 second long moving window
    stepsize=.1;       %0.5 second step size
   % thres=mean(eyez); thres=mean(thres);
   thres=3;
    movwin=round(movingwindow*samplespersec);
    stepwin=round(samplespersec*stepsize);
    findmaxwin=round(samplespersec*.1);    
    numovbins=movwin/stepwin;       %number overlapping bins
    bins=[];
    histbins=[];
    rates=[];
    pulsets=[];
    for ibin=1:numovbins
        firstid=stepwin*(ibin-1)+1;
        bins(ibin,:)=firstid:movwin:size(eyez,2);
    end
        %find peaks for blinks
    for itrial=1:size(eyez,1)
        [pks,locs,w,p] = findpeaks(eyez(itrial,:),'minpeakwidth',...
            round(samplespersec*.05),'minpeakdistance',round(samplespersec*.05)...
            ,'MinPeakprominence',thres);
        for ibin=1:numovbins
            histbins(ibin,:)=histc(locs,bins(ibin,:));
        end
        bb=reshape(bins+movwin/2-1,1,numel(bins));  %middle of each bin and single vector
        [aa,a1]=sort(bb);
        countpeaksperbin=reshape(histbins,1,numel(bins));
        countpeaksperbin=countpeaksperbin(a1);
        rates(itrial,:)=countpeaksperbin./(movwin/samplespersec);
        %rates(itrial,:)=countpeaksperbin./(movwin/samplespersec/60);
        pulsets(itrial,:)=round(bb./samplespersec*10)/10;
    end
    samplespersec=1/stepsize;
     tslfp=pulsets(1,:);
    ts1=plotparam.win(1)./plotparam.samplespersec;
    ts2=plotparam.win(end)./plotparam.samplespersec;
    tslfprounded=round(tslfp.*1000)/1000;
    win1id=find(tslfprounded==ts1);
    win2id=find(tslfprounded==ts2);
    win=win1id:win2id;          %get window of samples in nlx ts id's
    %convert marker id's to new sample rate
    markers{1}=markers{1}./ceil(data.ratelfp)*samplespersec;
    markers{2}=markers{2}./ceil(data.ratelfp)*samplespersec;
    markers{3}=markers{3}./ceil(data.ratelfp)*samplespersec;   
        rewardidx=rewardidx./ceil(data.ratelfp)*samplespersec;

    da=rates;
end

if pulse==1
    %calculate HR
    pulsed=zscore(da,0,2);
    movingwindow=3;     %5 second long moving window
    stepsize=.1;       %0.5 second step size
    thres=mean(pulsed); thres=mean(thres);
    movwin=round(movingwindow*samplespersec);
    stepwin=round(samplespersec*stepsize);
    findmaxwin=round(samplespersec*.1);    
    numovbins=movwin/stepwin;       %number overlapping bins
    bins=[];
    histbins=[];
    rates=[];
    pulsets=[];
    for ibin=1:numovbins
        firstid=stepwin*(ibin-1)+1;
        bins(ibin,:)=firstid:movwin:size(pulsed,2);
    end
    for itrial=1:size(pulsed,1)
    [pks,locs,w,p] = findpeaks(pulsed(itrial,:),'minpeakdistance',findmaxwin,'MinPeakHeight',1);
        for ibin=1:numovbins
            histbins(ibin,:)=histc(locs,bins(ibin,:));
        end
        bb=reshape(bins+movwin/2-1,1,numel(bins));  %middle of each bin and single vector
        [aa,a1]=sort(bb);
        countpeaksperbin=reshape(histbins,1,numel(bins));
        countpeaksperbin=countpeaksperbin(a1);
        rates(itrial,:)=countpeaksperbin./(movwin/samplespersec/60);
        %rates(itrial,:)=countpeaksperbin./(movwin/samplespersec/60);
        if mean(rates(itrial,:))<40
            %not heart rate
            rates(itrial,:)=nan(1,length(countpeaksperbin));
        end
        pulsets(itrial,:)=round(bb./samplespersec*10)/10;
    end
    samplespersec=1/stepsize;
     tslfp=pulsets(1,:);
    ts1=plotparam.win(1)./plotparam.samplespersec;
    ts2=plotparam.win(end)./plotparam.samplespersec;
    tslfprounded=round(tslfp.*1000)/1000;
    win1id=find(tslfprounded==ts1);
    win2id=find(tslfprounded==ts2);
    win=win1id:win2id;          %get window of samples in nlx ts id's
    %convert marker id's to new sample rate
    markers{1}=markers{1}./ceil(data.ratelfp)*samplespersec;
    markers{2}=markers{2}./ceil(data.ratelfp)*samplespersec;
    markers{3}=markers{3}./ceil(data.ratelfp)*samplespersec;   
        rewardidx=rewardidx./ceil(data.ratelfp)*samplespersec;

    da=rates;
end

baseline=plotparam.baseline;
pade=round(.2*samplespersec);       %padding before subseq event
basesamples=round(.5*samplespersec);        %.3 s for baseline average;
if iscsc
basesamples=round(.3*samplespersec);        %.3 s for baseline average lfp;
end
postrewdur=round(samplespersec*6.5); %5 s post reward interval to average
postrewdurshort=round(samplespersec*3); %3s post rew
postrewdurim=round(samplespersec*1.5); %1.5 post rew
avgpostrewpad=round(samplespersec*1);       %1s avg for postrew baseline
targimdur=round(samplespersec*1);      %1s post targ immediate
rewpredur=round(samplespersec*1.5);     %pre reward perio 1.5s

%re normalize data to targeted baseline period
baselinedata=[];
if ~iscsc
if isequal(baseline,'fix')    
    for xx=1:numtrials
        baselinemean=nanmean(da(seltrials(xx),markers{1}(seltrials(xx)):markers{2}(seltrials(xx))));
        if length(markers)>=3
            %use eye on fix cue to fix cue apperance before hand as
            %baseline
            baselinemean=nanmean(da(seltrials(xx),markers{1}(seltrials(xx)):markers{3}(seltrials(xx))));
            if isnan(baselinemean)
                %seaerch around points to find non-nan value
                baselinemean=nanmean(da(seltrials(xx),markers{1}(seltrials(xx))-5:markers{3}(seltrials(xx))+5));
            end
            if markers{3}(seltrials(xx))<markers{1}(seltrials(xx))
                %if for some reason no first marker found or at next trial
                baselinemean=nanmean(da(seltrials(xx),markers{3}(seltrials(xx))-25:markers{3}(seltrials(xx))));
            end
        end
        baselinedata(xx)=baselinemean;
        baselinerep=repmat(baselinemean,1,length(da(seltrials(xx),:)));
       % if ~iscsc
            %only fscv subtract baseline
            da(seltrials(xx),:)=da(seltrials(xx),:)-baselinerep;
       % end
    end
end
if isequal(baseline,'fixeye')    
    %1s post fix eye
    for xx=1:numtrials
        baselinemean=nanmean(da(seltrials(xx),markers{3}(seltrials(xx)):markers{3}(seltrials(xx))+10));
            if isnan(baselinemean)
                %seaerch around points to find non-nan value
                baselinemean=nanmean(da(seltrials(xx),markers{3}(seltrials(xx))-10:markers{3}(seltrials(xx))+10));
            end
            if markers{3}(seltrials(xx))<markers{1}(seltrials(xx))
                %if for some reason no first marker found or at next trial
                baselinemean=nanmean(da(seltrials(xx),markers{3}(seltrials(xx))-25:markers{3}(seltrials(xx))));
            end
        
        baselinedata(xx)=baselinemean;
        baselinerep=repmat(baselinemean,1,length(da(seltrials(xx),:)));
       % if ~iscsc
            %only fscv subtract baseline
            da(seltrials(xx),:)=da(seltrials(xx),:)-baselinerep;
       % end
    end
end
if isequal(baseline,'postcue')    
    %1s post fix first cue
    for xx=1:numtrials
        baselinemean=nanmean(da(seltrials(xx),markers{1}(seltrials(xx)):markers{1}(seltrials(xx))+10));
            if isnan(baselinemean)
                %seaerch around points to find non-nan value
                baselinemean=nanmean(da(seltrials(xx),markers{1}(seltrials(xx))-10:markers{1}(seltrials(xx))+10));
            end
            if markers{3}(seltrials(xx))<markers{1}(seltrials(xx))
                %if for some reason no first marker found or at next trial
                baselinemean=nanmean(da(seltrials(xx),markers{1}(seltrials(xx))-25:markers{1}(seltrials(xx))));
            end
        
        baselinedata(xx)=baselinemean;
        baselinerep=repmat(baselinemean,1,length(da(seltrials(xx),:)));
       % if ~iscsc
            %only fscv subtract baseline
            da(seltrials(xx),:)=da(seltrials(xx),:)-baselinerep;
       % end
    end
end
if isequal(baseline,'pretarg')    
    %.5s pretarg
    for xx=1:numtrials
        tmark=markers{2}(seltrials(xx))-5:markers{2}(seltrials(xx));
        if any(tmark<1)
            baselinemean=nan;
        else
        baselinemean=nanmean(da(seltrials(xx),tmark));
        end
        baselinedata(xx)=baselinemean;
        baselinerep=repmat(baselinemean,1,length(da(seltrials(xx),:)));
      %  if ~iscsc
            %only fscv subtract baseline
            da(seltrials(xx),:)=da(seltrials(xx),:)-baselinerep;
     %   end
    end
end
if isequal(baseline,'precue')    
    %1s pre fix first cue
    for xx=1:numtrials
        basetids=markers{1}(seltrials(xx))-10:markers{1}(seltrials(xx));
        if any(basetids<1)
            basetids=1:10;
        end
        baselinemean=nanmean(da(seltrials(xx),basetids));
            if isnan(baselinemean)
                %seaerch around points to find non-nan value
                baselinemean=nanmean(da(seltrials(xx),markers{1}(seltrials(xx))-20:markers{1}(seltrials(xx))));
            end
            if markers{3}(seltrials(xx))<markers{1}(seltrials(xx))
                %if for some reason no first marker found or at next trial
                baselinemean=nanmean(da(seltrials(xx),markers{1}(seltrials(xx))-25:markers{1}(seltrials(xx))));
            end
        
        baselinedata(xx)=baselinemean;
        baselinerep=repmat(baselinemean,1,length(da(seltrials(xx),:)));
        %if ~iscsc
            %only fscv subtract baseline
            da(seltrials(xx),:)=da(seltrials(xx),:)-baselinerep;
       % end
    end
end
if isequal(baseline,'prereward')
   % rewardidx=median(win);       %assume 0 seconds is reward period;
   %2 s before reward to 0.5 s before
        for xx=1:numtrials
        baselinemean=nanmean(da(seltrials(xx),rewardidx-pade-basesamples:rewardidx-pade));
        if length(markers)<3
            %if for some reason no first marker found or at next trial
            baselinemean=nanmean(da(seltrials(xx),rewardidx-25:rewardidx));
        end
        baselinedata(xx)=baselinemean;
        baselinerep=repmat(baselinemean,1,length(da(seltrials(xx),:)));
       % if ~iscsc
            da(seltrials(xx),:)=da(seltrials(xx),:)-baselinerep;
      %  end
        end
end
if isequal(baseline,'postreward')
   % rewardidx=median(win);       %assume 0 seconds is reward period;
   %2 s before reward to 0.5 s before
        for xx=1:numtrials
        baselinemean=nanmean(da(seltrials(xx),rewardidx+pade:rewardidx+pade+avgpostrewpad));
        
        baselinedata(xx)=baselinemean;
        baselinerep=repmat(baselinemean,1,length(da(seltrials(xx),:)));
        da(seltrials(xx),:)=da(seltrials(xx),:)-baselinerep;
        end
end
if isequal(baseline,'alltrial')
   %entire trial is baseline
        for xx=1:numtrials
        baselinemean=nanmean(da(seltrials(xx),markers{1}(seltrials(xx)):rewardidx+postrewdurshort));
        baselinedata(xx)=baselinemean;
        baselinerep=repmat(baselinemean,1,length(da(seltrials(xx),:)));
        da(seltrials(xx),:)=da(seltrials(xx),:)-baselinerep;
        end
end
end



interval=plotparam.interval;        %interval xtick marks def 2.5s
cscale=plotparam.cscale;
cminshow=plotparam.cminshow;

xticks=1:interval*samplespersec:size(da(:,win),2);
%mintime=round((min(win)-median(win))./samplespersec);
mintime=round((win(1)-rewardidx))/samplespersec;
maxtime=round((win(end)-rewardidx))/samplespersec;
%maxtime=round((max(win)-median(win))./samplespersec);
xticklabels=mintime:interval:maxtime;
xticklabels=round(xticklabels.*10)./10;
xticklabels=num2str(xticklabels');

    
        
pdata=[];
pdata=da(seltrials,win(1):win(end));
plotdata.data=pdata;
numshift=zeros(size(da,1),1);
if ~isempty(shiftdata)
%shift to targ/fix instead of default reward/outcome time
    if strcmp(shiftdata,'targ')
        %shift to targ cue appear, markers{2}
        avgtarg=round(mean(markers{2}));
        datemp=da;
        tempshift=[];
        for itrial=1:size(datemp,1)
            numshift(itrial)=avgtarg-markers{2}(itrial);
            tempshift(itrial,:)=circshift(datemp(itrial,:),[0 numshift(itrial)]);
        end
        pdata=tempshift(seltrials,win(1):win(end));
    end
        if strcmp(shiftdata,'fix')
        %shift to fix cue appear, markers{1}
        avgtarg=round(mean(markers{1}));
        datemp=da;
        tempshift=[];
        for itrial=1:size(datemp,1)
            numshift(itrial)=avgtarg-markers{1}(itrial);
            tempshift(itrial,:)=circshift(datemp(itrial,:),[0 numshift(itrial)]);
        end
        pdata=tempshift(seltrials,win(1):win(end));
    end
end
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
means=nanmean(pdata,1);
means2=nanmean(means);
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
means=nanmean(pdata,2);
stds=nanstd(pdata,0,2);
means=nanmean(means);
stds=nanmean(stds);
if ~isempty(ylimit)
cscale=ylimit;
end
if ~isnan(cscale(1)) 
    ylim(hax,cscale); 
else
    cscale(1)=means-2*stds;
    cscale(2)=means+3*stds;
    if iscsc
        cscale(1)=means-2*stds;
        cscale(2)=means+20*stds;
        if avg
            cscale(1)=means-3*stds;
            cscale(2)=means+8*stds;
        end
    end
    if phys
        cscale(1)=means-2*stds;
        cscale(2)=means+5*stds;
    end
    if pulse
        cscale(1)=means-5*stds;
        cscale(2)=means+5*stds;
    end
    if ~isempty(cscale) && cscale(2)>cscale(1) && ~any(isnan(cscale))
    ylim(hax,cscale)
    end
end

line(hax,[mean(markers{1}(seltrials)-win(1)+mean(numshift(seltrials))) mean(markers{1}(seltrials)-win(1)+mean(numshift(seltrials)))],[cscale(1) cscale(2)],'Color',colorm(1,:)*.35,'LineWidth',.5,'LineStyle','--')
line(hax,[mean(markers{3}(seltrials)-win(1)+mean(numshift(seltrials))) mean(markers{3}(seltrials)-win(1)+mean(numshift(seltrials)))],[cscale(1) cscale(2)],'Color',colorm(3,:)*.45,'LineWidth',0.5,'LineStyle','--')

if ~strcmp(plotparam.trialtype,'fixbreak')
    line(hax,[mean(markers{2}(seltrials)-win(1)+mean(numshift(seltrials))) mean(markers{2}(seltrials)-win(1)+mean(numshift(seltrials)))],[cscale(1) cscale(2)],'Color',colorm(2,:)*.45,'LineWidth',0.5,'LineStyle','--')
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
            titlename=[trialtype ' DA, ' sitename{:}];
        end
    %title(hax,[trialtype ' trials, DA, ch ' num2str(chnum)])
else
    %title(hax,[trialtype ' trials, lfp beta, ' cscname])
    %band=bandfilt;
    %title(hax,['lfp ' num2str(band(1)) '-' num2str(band(2)) ' Hz, ' cscname])
    titlename=[trialtype ' lfp ' num2str(bandfilt(1)) '-' num2str(bandfilt(2)) ' Hz, ' cscname];
    clabel='\beta-lfp \times 10^9 (V^2)';
end
if ~isempty(ulabel)
    titlename=[titlename ', ' ulabel];
end
%if plotnum<=2
    title(hax,titlename)
%end

hold(hax,'off')
%plotdata=pdata;

end




