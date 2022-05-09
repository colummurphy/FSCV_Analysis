function [plotdata,da,info]=setTrialAx(data,plotparam,hax,seltrials,varargin)
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

plotnum=0;
numtrials=size(da,1);
if isempty(seltrials)
    seltrials=1:numtrials;
end
info={};
numtrials=size(da(seltrials,:),1);
argnum=1;
abssig=0;
markers={};
%markers{1} = fix cue appearance
%markers{2} = target cue apperance
%markers{3} = eye on fix cue
iscsc=0;        %not csc data, ie fscv time stamps
cscname=[];
pulse=0;
lick=0;
lfp2=0;
islfp=0;
eyex=0;
getfft=0;
logscale=0;
eye=0;
blink=0;
noplot=0;
sitename={};
rclfp=[];
basetrial=0;
filtband='';       %default filter band is beta low
filtlick=[0 10];        %freq bands of lick filter
filtbetal=[12 18];        %low beta%    changed from 13 to 12 07/07/18
filtbetah=[18 33];        %high beta
filttheta=[6 10];           %theta
filtdelta=[0.5 4];           %delta
filtgammal=[40 58]; 
filtbeta=[13 30];
    bandfilt=[0 0];
    fftfreqlim=[5 50];
    ffttapers=[1.8 1];      %default tapers unless other wise spec
    fftwin=[1 0.25];       %1 s window use .25 s step as averaging across targeted window for fft to keep same f length
    rmdc=0;
    norm=0;
    envwidth=0;
    if isfield(plotparam,'fftwin')
        fftwin=plotparam.fftwin;
    end
    colortrial=plotparam.colormap;
triallabel=[];
if isfield(plotparam,'triallabel')
    triallabel=plotparam.triallabel;
end
while argnum <= length(varargin)
    if isnumeric(varargin{argnum})
        markers{argnum}=varargin{argnum}; %store marker events
        %markers{1} is fix cue appear
        %markers{2} is target cue appear
        %markers{3} is eye fix
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
            case 'sitename'
                argnum=argnum+1;        %fscv ch site name
                sitename=varargin{argnum};
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
            case 'smoothwin'
                %provide smooth width for enveloped signal
                argnum=argnum+1;
                envwidth=varargin{argnum};  %in seconds
            case 'getfft'
                getfft=1;
            case 'tapers'
                argnum=argnum+1;
                ffttapers=varargin{argnum};
            case 'log'
                logscale=1;
            case 'noplot'
                noplot=1;   %don't plot
            case 'fftwin'
                argnum=argnum+1;
                fftwin=varargin{argnum};
            case 'colormap'
                argnum=argnum+1;
                colortrial=varargin{argnum};    %user specified color map
            case 'fixedbasetrial'
                %uniform baseline for all trials that is all trials
                basetrial=1;
            case 'triallabel'
                argnum=argnum+1;
                triallabel=varargin{argnum};
            case 'rclfp'
                %lfp already reconverted and provided here filtered
                argnum=argnum+1;
                rclfp=varargin{argnum};
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
                   da=(data.lfp{cscid});
                   islfp=1;
               end
               if isequal(varargin{argnum},'pulse')
                   pulse=1;
               end
               if isequal(varargin{argnum},'eyex')
                   eyex=1;
               end
               if isequal(varargin{argnum},'lick')
                   lick=1;
               end
               %if logscale==1
                %   da=10*log10(da);
              % end
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
baseline=plotparam.baseline;
shiftdata=[];
if isfield(plotparam,'shift')
    %shift data to targ or fix marker
    shiftdata=plotparam.shift;
end
if basetrial==1
    baseline='alltrial';
end
colorsize=[500 500];
vertplotwidth=[500 500];
colorm=[0 0 0];
chnum=0;
if ~noplot
colorsize=plotparam.colorsize;
vertplotwidth=plotparam.vertplotwidth;
vertplotwidth2=plotparam.vertplotwidth2;
colorm=plotparam.markercolors;
chnum=plotparam.chnum;

end
samplespersec=ceil(plotparam.samplespersec);
trialtype=[plotparam.trialtype ' | ' triallabel];

%get baseline averages for selected window periods
%markers{1}=fix cue appearance
%markes{2}= target appearance
%markers{3}=fix eye




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
            %square & envelope signal based on default win length
            if envwidth==0
                winlength=round(samplespersec*.5/mean(bandfilt));
            else
                winlength=envwidth.*samplespersec; %user provided envelop width
            end
            datemp=datemp.^2;   %get power V^2
            datemp=smoothwin(datemp,winlength);   %smoothing 
            da(ii,:)=datemp.*1e9;
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
  %filter data
        sitename='eye';
   % da=eyedata;
    %take out infiinite slopes (ie blinks that can amplify avg)
end
if eyex==1
   datatemp=[];
  for itrial=1:size(da,1)
      datatemp=deglitchnanamp(da(itrial,:),2.5e-3,30);
      da(itrial,:)=datatemp;
   end
  %filter data
        sitename='eyex';
   % da=eyedata;
    %take out infiinite slopes (ie blinks that can amplify avg)
end
if lick==1
    %filter
    for ii=1:size(da,1)
    datemp=filterLFP(da(ii,:),samplespersec,[0 100]);            
    da(ii,:)=datemp;
    end
end
if blink==1
    sitename='blink';
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
        sitename='pulse';

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



pade=round(.2*samplespersec);       %padding before subseq event
basesamples=round(.5*samplespersec);        %.3 s for baseline average;

%baseline subtraction for DA fscv data onluy
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
     %   end
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
      %  if ~iscsc
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
      %  if ~iscsc
            %only fscv subtract baseline
            da(seltrials(xx),:)=da(seltrials(xx),:)-baselinerep;
     %   end
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
       % if ~iscsc
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
   basetrial=1;
        for xx=1:numtrials
        baselinemean=nanmean(da(seltrials(xx),markers{1}(seltrials(xx)):rewardidx+postrewdurshort));
        baselinedata(xx)=baselinemean;
        baselinerep=repmat(baselinemean,1,length(da(seltrials(xx),:)));
        da(seltrials(xx),:)=da(seltrials(xx),:)-baselinerep;
        end
end
end
%CALCULAte averages for given intervals for each trial
%if fft flag, get avg spectra
%ie post-reward 0.2 - 2.5 s, post-target to reward makresr{3} to -0.2 s
%fix to target (markers{2} to marks{3} relative to baseline which is 0.5 s
%average before previous event
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
baselinedata=[];
postfixwin=[];          %average window after fix to target appear
postfixpeak=[];         %find peak after fix to target appear
fixpeakidx=[];
posttargwin=[];      %average window after target appear to reward or 0s
posttargpeak=[];      %find peak after target appear until reward or 0s
posttargpeakdiff=[];      %LFP subtract baseline (same for DA)
targpeakidx=[];
targimwin=[];
targimpeak=[];
targimpeakdiff=[];%LFP subtract baseline (same for DA)
targimpeakidx=[];       %immediatly post targ 1s
postrewwin=[];      %average window after reward 0 to 5s
postrewpeak=[];      %find peak after reward 0 to 5 s
postrewpeakdiff=[];      %LFP subtract baseline (same for DA)

rewpeakidx=[];
rewshortwin=[];
rewshortpeak=[];     %peak 0 to 3 s post rew 
rewshortpeakdiff=[];     %LFP subtract baseline (same for DA)
rewshortpeakidx=[];
rewimwin=[];        %immediate 1.5s after rew
rewimpeak=[];
rewimpeakdiff=[];%LFP subtract baseline (same for DA)
targpeakvfix=[];
rewimpeakidx=[];
padfft=samplespersec*.5;        %500 ms padding +/- for peak fft
fftwinids=fftwin(1)*samplespersec;
fftstepids=fftwin(2)*samplespersec;
%peaks are averaged around 0.2 s +/- window _NO
%windows are offset to 0.2 - to subsequent event
%reward changes seem to occur  up to 6 s after reward onset

%initialize storage variable plotdata with event modulated signals
%store site name
if isempty(sitename)
    sitename=['ch' num2str(chnum)];
end
plotdata.site=sitename;
plotdata.trialnums=seltrials;
%datawin=data.lfp{cscid}(1,1:padfft*2+1);
%[S,f]=chronux_mtspectrumc(datawin',ffttapers,0,samplespersec,fftfreqlim,0,1);  
%defaultpeakfftsize=length(f);
if iscsc
datawin=data.lfp{cscid}(1,1:samplespersec*5).*1e6;     %convert to microvolts first            
[S,winmid,f]=lfp_mtspecgram2(datawin',fftwin,ffttapers,0,samplespersec,fftfreqlim,0,1,rmdc,norm);
%[S,f]=chronux_mtspectrumc(datawin',ffttapers,0,samplespersec,fftfreqlim,0,1);  
plotdata.fftf=f;        %default freq variable
end

for xx=1:numtrials
    %markers{3} eye fixed, markers{2} target appears, markers{1} fix cue
    tid=seltrials(xx);
    baseidx=markers{1}(tid);
    if markers{1}(tid)>rewardidx
        %something wrong with fix cue appearance ts since higher than
        %reward point
        if markers{3}(tid)<rewardidx
            baseidx=markers{3}(tid);
        else
            baseidx=markers{2}(tid)-30;
        end
    end
    baselineids=baseidx-pade-basesamples:baseidx-pade;
    if basetrial==1
        %uniform baseline for all events that is entire trial
        baselineids=baseidx:rewardidx+postrewdurshort;
    end
    if any(baselineids<1)
        %any outside window, something wrong with marker timestamps
        postfixwin(xx)=nan;
        postfixpeak(xx)=nan;
        fixpeakidx(xx)=nan;
        plotdata.basewin(xx)=nan;
        if islfp && getfft
        plotdata.fixwinfft{xx}=nan;
        plotdata.fixwinfftavg(xx,:)=nan(1,length(plotdata.fftf));
        plotdata.fixpostfftavg(xx,:)=nan(1,length(plotdata.fftf));
        plotdata.targprefftavg(xx,:)=nan(1,length(plotdata.fftf));
            plotdata.prefixwinfft{xx}=nan;
            plotdata.prefixwinfftavg(xx,:)=nan(1,length(plotdata.fftf));
        end
    else
        fixbaseline=nanmean(da(tid,baselineids));
        plotdata.basewin(xx)=fixbaseline;
        if iscsc
            fixbaseline=0;     %10/05/2018 change so LFPs also baseline
            %sub otherwise cannot account for decreases...for peak changes
            
        end
        if markers{2}(tid)-pade*2<=baseidx
            markers{2}(tid)=baseidx+1*samplespersec;
        end
        postfixids=baseidx:markers{2}(tid)-pade;        %12/19/2018 - pade
        postfixdwellids=baseidx+pade:markers{2}(tid)-pade;        %12/19/2018 add pade
        
        postfixwin(xx)=nanmean(da(tid,postfixdwellids))-fixbaseline;
        %abssig=abs(da(tid,postfixids)-fixbaseline);
        sigda=da(tid,postfixids)-fixbaseline;
        %postfixmax=max(abssig);
        postfixmax=max(sigda);
        nanratio=length(find(isnan(da(tid,postfixids))==1))/length(da(tid,postfixids));
        if isnan(postfixmax) || nanratio>.5
            postfixpeak(xx)=nan;
            fixpeakidx(xx)=nan;
        else
            postfixpeakidx=find(sigda==postfixmax)+baseidx-1;        %absolute index
            if length(postfixpeakidx>1)
                postfixpeakidx=postfixpeakidx(1);
            end
            postfixpeak(xx)=da(tid,postfixpeakidx)-fixbaseline; %relative max change
            fixpeakidx(xx)=postfixpeakidx-baseidx;     %relative idx from preceding event
        end
        if islfp && getfft
            %get fft spectra task modulated
            datawin=data.lfp{cscid}(tid,postfixids).*1e6;     %convert to microvolts first                        
            if fftwin(1)*samplespersec<size(datawin,2)
                [S,winmid,f]=lfp_mtspecgram2(datawin',fftwin,ffttapers,0,samplespersec,fftfreqlim,0,1,rmdc,norm);
                plotdata.fixwinfft{xx}=S;
                plotdata.fixwinfftavg(xx,:)=mean(S,1);
                plotdata.fixpostfftavg(xx,:)=S(1,:);     %1 s window immediately after event
            else
                plotdata.fixwinfft{xx}=nan;
                plotdata.fixwinfftavg(xx,:)=nan(1,length(plotdata.fftf));
                plotdata.fixpostfftavg(xx,:)=nan(1,length(plotdata.fftf));
            end
            %[S,f]=chronux_mtspectrumc(datawin',ffttapers,0,samplespersec,fftfreqlim,0,1);  
            pretargetids=markers{2}(tid)-fftwinids:markers{2}(tid);
            datawin=data.lfp{cscid}(tid,pretargetids).*1e6;     %convert to microvolts first            
            if fftwin(1)*samplespersec<size(datawin,2)
                [S,winmid,f]=lfp_mtspecgram2(datawin',fftwin,ffttapers,0,samplespersec,fftfreqlim,0,1,rmdc,norm);
                plotdata.targprefftavg(xx,:)=S(end,:);        %immediately prior target
            else
                plotdata.targprefftavg(xx,:)=nan(1,length(plotdata.fftf));
            end
            baselinefftids=baseidx-pade-fftwinids:baseidx-pade;
            if any(baselinefftids<1)
                plotdata.prefixwinfft{xx}=nan;
                plotdata.prefixwinfftavg(xx,:)=nan(1,length(plotdata.fftf));
            else
                datawin=data.lfp{cscid}(tid,baselinefftids).*1e6;     %convert to microvolts first            
                [S,winmid,f]=lfp_mtspecgram2(datawin',fftwin,ffttapers,0,samplespersec,fftfreqlim,0,1,rmdc,norm);
                plotdata.prefixwinfft{xx}=S;
                plotdata.prefixwinfftavg(xx,:)=mean(S,1);
            end
        end
    end
    
    targbaselineids=markers{2}(tid)-pade-basesamples:markers{2}(tid)-pade;
    if basetrial==1
        %uniform baseline for all events that is entire trial
        targbaselineids=baseidx:rewardidx+postrewdurshort;
    end
    if any(targbaselineids<1) || strcmp(plotparam.trialtype,'fixbreak') ||...
            isnan(fixbaseline)
        %if fix break trials then no target/reward, skip this
        %if any outside window, something wrong wiht marker timestamps
        %skip both target/reward values
        posttargwin(xx)=nan;
        posttargpeak(xx)=nan;
                posttargpeakdiff(xx)=nan;
        targpeakidx(xx)=nan;
        targimwin(xx)=nan;
        targimpeak(xx)=nan;
                targimpeakdiff(xx)=nan;
        targimpeakidx(xx)=nan;
        postrewwin(xx)=nan;
        postrewpeak(xx)=nan;
                postrewpeakdiff(xx)=nan;
        rewpeakidx(xx)=nan;
        rewshortpeakidx(xx)=nan;
        rewshortpeak(xx)=nan;
                rewshortpeakdiff(xx)=nan;
        rewshortwin(xx)=nan;
                    rewimpeak(xx)=nan;
         rewimpeakdiff(xx)=nan;
targpeakvfix(xx)=nan;
            rewimpeakidx(xx)=nan;
            rewimwin(xx)=nan;
            rewprewin(xx)=nan;
            rewprepeak(xx)=nan;
             rewprepeakdiff(xx)=nan;
            rewprepeakidx(xx)=nan;
            rewpeakvfix(xx)=nan;
            rewimpeakvfix(xx)=nan;
            rewshortpeakvfix(xx)=nan;
        if islfp && getfft
        plotdata.targwinfft{xx}=nan;
        plotdata.targwinfftavg(xx,:)=nan(1,length(plotdata.fftf));
        plotdata.targpostfftavg(xx,:)=nan(1,length(plotdata.fftf));    
        plotdata.targimwinfft{xx}=nan;
        plotdata.targimwinfftavg(xx,:)=nan(1,length(plotdata.fftf));
        plotdata.targimpostfftavg(xx,:)=nan(1,length(plotdata.fftf));  
        plotdata.rewprefftavg(xx,:)=nan(1,length(plotdata.fftf));
        plotdata.rewwinfft{xx}=nan;
        plotdata.rewwinfftavg(xx,:)=nan(1,length(plotdata.fftf));
        plotdata.rewpostfftavg(xx,:)=nan(1,length(plotdata.fftf));
        end
    else        
        targbaseline=nanmean(da(tid,targbaselineids));
        pbaseline=targbaseline;     %for peaks for lfps need to subtract to find min
        if iscsc
            targbaseline=0;
        end
        targids=markers{2}(tid)+pade:rewardidx-pade;     %'dwell' target win w/ padding 
        targimids=markers{2}(tid):markers{2}+targimdur;     %immediate post
        rewpreids=rewardidx-rewpredur:rewardidx;
        posttargwin(xx)=nanmean(da(tid,targids))-targbaseline;
        targimwin(xx)=nanmean(da(tid,targimids))-targbaseline;
        rewprewin(xx)=nanmean(da(tid,rewpreids))-targbaseline;
        %
        targdasub=da(tid,targids)-pbaseline;
        targdasubvfix=da(tid,targids)-fixbaseline;
        targdasubim=da(tid,targimids)-pbaseline;
        prerewdasub=da(tid,rewpreids)-pbaseline;   
        if abssig
            %look at absolute changes
            targdasub=abs(targdasub);
            targdasubvfix=abs(targdasubvfix);
            targdasubim=abs(targdasubim);
            prerewdasub=abs(prerewdasub);       
        end
            
        posttargmax=max(targdasub);
        targimmax=max(targdasubim);
        %targimmaxpos=max(absimsig);
        targmaxvfix=max(targdasubvfix);
        rewpremax=max(prerewdasub);
        %posttargmax=max(sigda);
       nanratio=length(find(isnan(da(tid,targids))==1))/length(da(tid,targids));
       nantrial=0;
       if length(da(tid,targids))==0
           nantrial=1;
       end
       if nanratio>0.5
           nantrial=1;
       end
       if isnan(posttargmax)
           nantrial=1;
       end
        if nantrial
            targpeakvfix(xx)=nan;
            posttargpeak(xx)=nan;
            posttargpeakdiff(xx)=nan;
            targpeakidx(xx)=nan;
            targimpeakidx(xx)=nan;
            targimpeak(xx)=nan;
            targimpeakdiff(xx)=nan;
            rewprepeak(xx)=nan;
             rewprepeakdiff(xx)=nan;           
            rewprepeakidx(xx)=nan;
        else
            targpeakvfixidx=find(targdasubvfix==targmaxvfix)+markers{2}(tid)-1;        %sig vs fix baseline
            posttargpeakidx=find(targdasub==posttargmax)+markers{2}(tid)-1;        %absolute index
            posttargimpeakidx=find(targdasubim==targimmax)+markers{2}(tid)-1;        %absolute index
            prewprepeakidx=find(prerewdasub==rewpremax)+rewardidx-rewpredur-1; 
            if length(posttargpeakidx)>1
                posttargpeakidx=posttargpeakidx(1);
            end
            if length(posttargimpeakidx)>1
                posttargimpeakidx=posttargimpeakidx(1);
            end
            if length(prewprepeakidx)>1
                prewprepeakidx=prewprepeakidx(1);
            end
            if length(targpeakvfixidx)>1
                targpeakvfixidx=targpeakvfixidx(1);
            end
            %for actual stored values use 0 as basleine for LFPs, but keep
            %difference for DA
            targpeakvfix(xx)=da(tid,targpeakvfixidx)-fixbaseline;
            posttargpeak(xx)=da(tid,posttargpeakidx)-targbaseline; %relative max change
            posttargpeakdiff(xx)=da(tid,posttargpeakidx)-pbaseline;
            targpeakidx(xx)=posttargpeakidx-markers{2}(tid);     %relative idx from preceding event

            if isempty(posttargimpeakidx)
                targimpeak(xx)=nan;
                targpeakvfix(xx)=nan;
                targimpeakdiff(xx)=nan;
                targimpeakidx(xx)=nan;    
            else
                targimpeak(xx)=da(tid,posttargimpeakidx)-targbaseline;
                targpeakvfix(xx)=nan;
                targimpeakdiff(xx)=da(tid,posttargimpeakidx)-pbaseline;
                targimpeakidx(xx)=posttargimpeakidx-markers{2}(tid);
            end
            if isempty(prewprepeakidx)
                rewprepeak(xx)=nan;
                rewprepeakdiff(xx)=nan;
                rewprepeakidx(xx)=nan;
            else
                rewprepeak(xx)=da(tid,prewprepeakidx)-targbaseline;
                rewprepeakdiff(xx)=da(tid,prewprepeakidx)-pbaseline;
                rewprepeakidx(xx)=prewprepeakidx-rewardidx+rewpredur;
            end
        end
        if islfp && getfft
            %get fft spectra task modulated
            datawin=data.lfp{cscid}(tid,targids).*1e6;     %convert to microvolts first     
            if fftwin(1)*samplespersec>length(datawin)
                 plotdata.targwinfft{xx}=nan;
                  plotdata.targwinfftavg(xx,:)=nan(1,length(plotdata.fftf));
                plotdata.targpostfftavg(xx,:)=nan(1,length(plotdata.fftf));
                display(['targ fft| trial ' num2str(seltrials(xx)) ' fftwin > data length; skipping, nans set']);
              else
                [S,winmid,f]=lfp_mtspecgram2(datawin',fftwin,ffttapers,0,samplespersec,fftfreqlim,0,1,rmdc,norm);
                plotdata.targwinfft{xx}=S;
                plotdata.targwinfftavg(xx,:)=mean(S,1);
                plotdata.targpostfftavg(xx,:)=S(1,:);
            end

            prerewids=rewardidx-fftwinids:rewardidx;
            datawin=data.lfp{cscid}(tid,prerewids).*1e6;     %convert to microvolts first 
            [S,winmid,f]=lfp_mtspecgram2(datawin',fftwin,ffttapers,0,samplespersec,fftfreqlim,0,1,rmdc,norm);
            plotdata.rewprefftavg(xx,:)=S(end,:);       %last point is reward pre
            
            datawin=data.lfp{cscid}(tid,targimids).*1e6;     %convert to microvolts first     
            if fftwin(1)*samplespersec>length(datawin)
                 plotdata.targimwinfft{xx}=nan;
                  plotdata.targimwinfftavg(xx,:)=nan(1,length(plotdata.fftf));
                plotdata.targimpostfftavg(xx,:)=nan(1,length(plotdata.fftf));
                display(['targ fft| trial ' num2str(seltrials(xx)) ' fftwin > data length; skipping, nans set']);
              else
                [S,winmid,f]=lfp_mtspecgram2(datawin',fftwin,ffttapers,0,samplespersec,fftfreqlim,0,1,rmdc,norm);
                plotdata.targimwinfft{xx}=S;
                plotdata.targimwinfftavg(xx,:)=mean(S,1);
                plotdata.targimpostfftavg(xx,:)=S(1,:);
            end
        end
        
        %reward peak relative to pre reward baseline
        rewbaseids=rewardidx-pade-basesamples:rewardidx-pade;
        if basetrial==1
            %uniform baseline for all events that is entire trial
            rewbaseids=baseidx:rewardidx+postrewdurshort;
        end
        rewbaseline=nanmean(da(tid,rewbaseids));
        pbaseline=rewbaseline;
        rewbaselineminidx=find(da(tid,rewbaseids)==min(da(tid,rewbaseids)));
        rewbaselinemin=min(da(tid,rewbaseids));         
        if iscsc
            rewbaseline=0;
            rewbaselinemin=0;
        end
        if isnan(rewbaselinemin)
            disp(['nan trial ' num2str(tid)]);
            rewbaselinemin=0;
        end
        rewids=rewardidx:rewardidx+postrewdur;
        postrewwin(xx)=nanmean(da(tid,rewids))-rewbaseline;
        rewdasub=da(tid,rewids)-pbaseline;
        if abssig
            rewdasub=abs(rewdasub);
        end
        %change to non-nan peak max
        sigda=da(tid,rewids)-pbaseline;
        postrewmax=max(rewdasub);     %change to abs 10/05/2018
        %postrewmax=max(sigda);
        nanratio=length(find(isnan(da(tid,rewids))==1))/length(da(tid,rewids));
        if isnan(postrewmax) || nanratio>0.5
            postrewpeak(xx)=nan;
            postrewpeakdiff(xx)=nan;
            rewpeakidx(xx)=nan;
            rewpeakvfix(xx)=nan;            
        else
            postrewpeakidx=find(rewdasub==postrewmax)+rewardidx-1;        %absolute index
            if length(postrewpeakidx)>1
                postrewpeakidx=postrewpeakidx(1);
            end
            postrewpeak(xx)=da(tid,postrewpeakidx)-rewbaseline; %relative max change
            postrewpeakdiff(xx)=da(tid,postrewpeakidx)-pbaseline; %relative max change
            %postrewpeakmax(xx)=da(tid,postrewpeakidx)-rewbaselinemin; %relative max change minus min
            rewpeakidx(xx)=postrewpeakidx-rewardidx;     %relative idx from preceding event
            rewpeakvfix(xx)=da(tid,postrewpeakidx)-fixbaseline;
        end
        
        if islfp && getfft
            %get fft spectra task modulated
            datawin=data.lfp{cscid}(tid,rewids).*1e6;     %convert to microvolts first                 
            [S,winmid,f]=lfp_mtspecgram2(datawin',fftwin,ffttapers,0,samplespersec,fftfreqlim,0,1,rmdc,norm);
            plotdata.rewwinfft{xx}=S;
            plotdata.rewwinfftavg(xx,:)=mean(S,1);
            plotdata.rewpostfftavg(xx,:)=S(1,:);
        end
        
        %reward peak short window (3s)
        rewids=rewardidx:rewardidx+postrewdurshort;
        rewshortwin(xx)=nanmean(da(tid,rewids))-rewbaseline;
        rewdasub=da(tid,rewids)-rewbaseline;
        if abssig            
            rewdasub=abs(rewdasub);
        end
        postrewmax=max(rewdasub);
        if isnan(postrewmax)
            rewshortpeak(xx)=nan;
            rewshortpeakdiff(xx)=nan;
            rewshortpeakidx(xx)=nan;
            rewshortpeakvfix(xx)=nan;
        else
            peakvfixidx=find(rewdasub==postrewmax)+rewardidx-1;        %absolute index
            if length(peakvfixidx)>1
                peakvfixidx=peakvfixidx(1);
            end
            rewshortpeak(xx)=da(tid,peakvfixidx)-rewbaseline; %relative max change
            rewshortpeakdiff(xx)=da(tid,peakvfixidx)-pbaseline; %relative max change
            rewshortpeakidx(xx)=peakvfixidx-rewardidx;     %relative idx from preceding event
            rewshortpeakvfix(xx)=da(tid,peakvfixidx)-fixbaseline;
        end
        
        %reward peak relative to pre baseline immediate     
        rewids=rewardidx:rewardidx+postrewdurim;
        rewimwin(xx)=nanmean(da(tid,rewids))-rewbaseline;
        rewdasub=abs(da(tid,rewids)-rewbaseline);
        postrewmax=max(rewdasub);
        if isnan(postrewmax)
            rewimpeak(xx)=nan;
           rewimpeakdiff(xx)=nan;
            rewimpeakidx(xx)=nan;
            rewimpeakvfix(xx)=nan;
        else
            peakvfixidx=find(rewdasub==postrewmax)+rewardidx-1;        %absolute index
            if length(peakvfixidx)>1
                peakvfixidx=peakvfixidx(1);
            end
            rewimpeak(xx)=da(tid,peakvfixidx)-rewbaseline; %relative max change
            rewimpeakdiff(xx)=da(tid,peakvfixidx)-pbaseline; %relative max change            
            rewimpeakidx(xx)=peakvfixidx-rewardidx;     %relative idx from preceding event
            rewimpeakvfix(xx)=da(tid,peakvfixidx)-fixbaseline;
        end
        
    end
    
end

plotdata.fixwin=postfixwin;
plotdata.fixwinci=nanstd(plotdata.fixwin)/sqrt(length(find(~isnan(postfixwin))))*1.96;  %95% confidence interval
plotdata.fixwinavg=nanmean(plotdata.fixwin);
plotdata.fixpeak=postfixpeak;
plotdata.fixpeakci=nanstd(plotdata.fixpeak)/sqrt(length(find(~isnan(postfixpeak))))*1.96;  %95% confidence interval
plotdata.fixpeakavg=nanmean(plotdata.fixpeak);
plotdata.fixpeakts=fixpeakidx./samplespersec;
plotdata.fixpeaktsci=nanstd(plotdata.fixpeakts)/sqrt(length(find(~isnan(fixpeakidx))))*1.96;  %95% confidence interval
plotdata.fixpeaktsavg=nanmean(plotdata.fixpeakts);

plotdata.targwin=posttargwin;
plotdata.targwinci=nanstd(plotdata.targwin)/sqrt(length(find(~isnan(posttargwin))))*1.96;  %95% confidence interval
plotdata.targwinavg=nanmean(plotdata.targwin);
plotdata.targpeakvprefix=targpeakvfix;
plotdata.targpeakvprefixci=nanstd(plotdata.targpeakvprefix)/sqrt(length(find(~isnan(targpeakvfix))))*1.96;  %95% confidence interval
plotdata.targpeakvprefixavg=nanmean(plotdata.targpeakvprefix);
plotdata.targpeak=posttargpeak;
plotdata.targpeakci=nanstd(plotdata.targpeak)/sqrt(length(find(~isnan(posttargpeak))))*1.96;  %95% confidence interval
plotdata.targpeakavg=nanmean(plotdata.targpeak);
plotdata.targpeakdiff=posttargpeakdiff;
plotdata.targpeakdiffci=nanstd(plotdata.targpeakdiff)/sqrt(length(find(~isnan(posttargpeakdiff))))*1.96;  %95% confidence interval
plotdata.targpeakdiffavg=nanmean(plotdata.targpeakdiff);
plotdata.targpeakts=targpeakidx./samplespersec;
plotdata.targpeaktsci=nanstd(plotdata.targpeakts)/sqrt(length(find(~isnan(targpeakidx))))*1.96;  %95% confidence interval
plotdata.targpeaktsavg=nanmean(plotdata.targpeakts);


plotdata.targimwin=targimwin;
plotdata.targimwinci=nanstd(plotdata.targimwin)/sqrt(length(find(~isnan(targimwin))))*1.96;  %95% confidence interval
plotdata.targimwinavg=nanmean(plotdata.targimwin);
plotdata.targimpeak=targimpeak;
plotdata.targimpeakci=nanstd(plotdata.targimpeak)/sqrt(length(find(~isnan(targimpeak))))*1.96;  %95% confidence interval
plotdata.targimpeakavg=nanmean(plotdata.targimpeak);
plotdata.targimpeakdiff=targimpeakdiff;
plotdata.targimpeakdiffci=nanstd(plotdata.targimpeakdiff)/sqrt(length(find(~isnan(targimpeakdiff))))*1.96;  %95% confidence interval
plotdata.targimpeakdiffavg=nanmean(plotdata.targimpeakdiff);
plotdata.targimpeakts=targimpeakidx./samplespersec;
plotdata.targimpeaktsci=nanstd(plotdata.targimpeakts)/sqrt(length(find(~isnan(targimpeakidx))))*1.96;  %95% confidence interval
plotdata.targimpeaktsavg=nanmean(plotdata.targimpeakts);

plotdata.rewprewin=rewprewin;
plotdata.rewprewinci=nanstd(plotdata.rewprewin)/sqrt(length(find(~isnan(rewprewin))))*1.96;  %95% confidence interval
plotdata.rewprewinavg=nanmean(plotdata.rewprewin);
plotdata.rewprepeak=rewprepeak;
plotdata.rewprepeakci=nanstd(plotdata.rewprepeak)/sqrt(length(find(~isnan(rewprepeak))))*1.96;  %95% confidence interval
plotdata.rewprepeakavg=nanmean(plotdata.rewprepeak);
plotdata.rewprepeakdiff=rewprepeakdiff;
plotdata.rewprepeakdiffci=nanstd(plotdata.rewprepeakdiff)/sqrt(length(find(~isnan(rewprepeakdiff))))*1.96;  %95% confidence interval
plotdata.rewprepeakdiffavg=nanmean(plotdata.rewprepeakdiff);
plotdata.rewprepeakts=rewprepeakidx./samplespersec;
plotdata.rewprepeaktsci=nanstd(plotdata.rewprepeakts)/sqrt(length(find(~isnan(rewprepeakidx))))*1.96;  %95% confidence interval
plotdata.rewprepeaktsavg=nanmean(plotdata.rewprepeakts);

        
plotdata.rewwin=postrewwin;
plotdata.rewwinci=nanstd(plotdata.rewwin)/sqrt(length(find(~isnan(postrewwin))))*1.96;  %95% confidence interval
plotdata.rewwinavg=nanmean(plotdata.rewwin);
plotdata.rewpeak=postrewpeak;
plotdata.rewpeakci=nanstd(plotdata.rewpeak)/sqrt(length(find(~isnan(postrewpeak))))*1.96;  %95% confidence interval
plotdata.rewpeakavg=nanmean(plotdata.rewpeak);
plotdata.rewpeakts=rewpeakidx./samplespersec;
plotdata.rewpeaktsci=nanstd(plotdata.rewpeakts)/sqrt(length(find(~isnan(rewpeakidx))))*1.96;  %95% confidence interval
plotdata.rewpeaktsavg=nanmean(plotdata.rewpeakts);
plotdata.rewpeakdiff=postrewpeakdiff;
plotdata.rewpeakdiffci=nanstd(plotdata.rewpeakdiff)/sqrt(length(find(~isnan(postrewpeakdiff))))*1.96;  %95% confidence interval
plotdata.rewpeakdiffavg=nanmean(plotdata.rewpeakdiff);

plotdata.rewpeakvprefix=rewpeakvfix;
plotdata.rewpeakvprefixci=nanstd(plotdata.rewpeakvprefix)/sqrt(length(find(~isnan(rewpeakvfix))))*1.96;  %95% confidence interval
plotdata.rewpeakvprefixavg=nanmean(plotdata.rewpeakvprefix);
plotdata.rewimpeakvprefix=rewimpeakvfix;
plotdata.rewimpeakvprefixci=nanstd(plotdata.rewimpeakvprefix)/sqrt(length(find(~isnan(rewimpeakvfix))))*1.96;  %95% confidence interval
plotdata.rewimpeakvprefixavg=nanmean(plotdata.rewimpeakvprefix);
plotdata.rewshortpeakvprefix=rewshortpeakvfix;
plotdata.rewshortpeakvprefixci=nanstd(plotdata.rewshortpeakvprefix)/sqrt(length(find(~isnan(rewshortpeakvfix))))*1.96;  %95% confidence interval
plotdata.rewshortpeakvprefixavg=nanmean(plotdata.rewshortpeakvprefix);

            
plotdata.rewshortpeak=rewshortpeak;
plotdata.rewshortpeakci=nanstd(plotdata.rewshortpeak)/sqrt(length(find(~isnan(rewshortpeak))))*1.96;  %95% confidence interval
plotdata.rewshortpeakavg=nanmean(plotdata.rewshortpeak);
plotdata.rewshortpeakdiff=rewshortpeakdiff;
plotdata.rewshortpeakdiffci=nanstd(plotdata.rewshortpeakdiff)/sqrt(length(find(~isnan(rewshortpeakdiff))))*1.96;  %95% confidence interval
plotdata.rewshortpeakdiffavg=nanmean(plotdata.rewshortpeakdiff);
plotdata.rewshortpeakts=rewshortpeakidx./samplespersec;
plotdata.rewshortpeaktsci=nanstd(plotdata.rewshortpeakts)/sqrt(length(find(~isnan(rewshortpeakidx))))*1.96;  %95% confidence interval
plotdata.rewshortpeaktsavg=nanmean(plotdata.rewshortpeakts);
plotdata.rewshortwin=rewshortwin;
plotdata.rewshortwinci=nanstd(plotdata.rewshortwin)/sqrt(length(find(~isnan(rewshortwin))))*1.96;  %95% confidence interval
plotdata.rewshortwinavg=nanmean(plotdata.rewshortwin);

plotdata.rewimpeak=rewimpeak;
plotdata.rewimpeakci=nanstd(plotdata.rewimpeak)/sqrt(length(find(~isnan(rewimpeak))))*1.96;  %95% confidence interval
plotdata.rewimpeakavg=nanmean(plotdata.rewimpeak);
plotdata.rewimpeakdiff=rewimpeakdiff;
plotdata.rewimpeakdiffci=nanstd(plotdata.rewimpeakdiff)/sqrt(length(find(~isnan(rewimpeakdiff))))*1.96;  %95% confidence interval
plotdata.rewimpeakdiffavg=nanmean(plotdata.rewimpeakdiff);
plotdata.rewimpeakts=rewimpeakidx./samplespersec;
plotdata.rewimpeaktsci=nanstd(plotdata.rewimpeakts)/sqrt(length(find(~isnan(rewimpeakidx))))*1.96;  %95% confidence interval
plotdata.rewimpeaktsavg=nanmean(plotdata.rewimpeakts);
plotdata.rewimwin=rewimwin;
plotdata.rewimwinci=nanstd(plotdata.rewimwin)/sqrt(length(find(~isnan(rewimwin))))*1.96;  %95% confidence interval
plotdata.rewimwinavg=nanmean(plotdata.rewimwin);

plotdata.baseline=baseline;
plotdata.baselinedata=baselinedata;

if noplot==0
    cla(hax);
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
    if logscale==1
        da=10.*log10(da);
    end
    pdata=da(seltrials,win(1):win(end));
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
    %imagetrials=image(hax, pdata,'cdatamapping','scaled');
    imagetrials=image(pdata,'parent',hax,'cdatamapping','scaled');
    hold(hax,'on')

    origpos=getpixelposition(hax);      %get original position 
    %set color map to specified 
    colormap(hax,colortrial)
    %set artifact (on nan values) mask ontop of image using alpha data properties
    artTime=isnan(pdata);   %find artifact points (nan periods)
    artTime=abs(artTime-1);         %make alpha data mask by inverting 1/0's
    artTime2=artTime;
    maskGray=artTime2==0;             %find Zero indices representing artifact mask
    maskGray=maskGray*.15;            %make gray rather than white default by making non-zero
    artTime=artTime+maskGray;
    set(imagetrials, 'AlphaData', artTime);

    %plot event markers
    for jj=1:5:numtrials
        line(hax,[markers{1}(seltrials(jj))-win(1)+numshift(seltrials(jj))],jj,'Color',colorm(1,:),'LineWidth',0.5,'LineStyle',':','marker','s','markersize',3,'markerfacecolor','k')
        if ~strcmp(plotparam.trialtype,'fixbreak')
            line(hax,[markers{2}(seltrials(jj))-win(1)+numshift(seltrials(jj))],jj,'Color',colorm(2,:),'LineWidth',0.5,'LineStyle',':','marker','.','markersize',3)
        end
        line(hax,[markers{3}(seltrials(jj))-win(1)+numshift(seltrials(jj))],jj,'Color',colorm(3,:),'LineWidth',0.5,'LineStyle',':','marker','.','markersize',2)

    end
    
    %organize plot
    set(hax,'XTick',xticks)
    set(hax,'xticklabel',xticklabels)
    set(hax,'tickdir','out','box','off')
    xlabel(hax,'time (s)')
    clabel='[DA] (nM)';
    if ~iscsc
        title(hax,[trialtype ' DA ch ' num2str(chnum)])
        if ~isempty(sitename)
            title(hax,[trialtype ' DA ' sitename])
        end
    else
        title(hax,[trialtype ' lfp beta ' cscname])
        band=[];
      band=bandfilt;
        title(hax,{[trialtype] ['lfp ' num2str(band(1)) '-' num2str(band(2)) ' Hz '] [cscname]})
        clabel='\beta-lfp \times 10^9 (V^2)';
    end
    matlabver=version('-release');
    matlabver=str2num(matlabver(regexp(matlabver,'\d')));
    if matlabver>2013
        %colorbar function introduced in matlab 2014
        h1=colorbar(hax,'southoutside');
        cpos = getpixelposition(h1);
        ylabelbar=ylabel(h1,clabel); ypos = getpixelposition(ylabelbar);
        cpos(4) = 15; 
        set(h1,'Units','Pixels','Position', [cpos(1) origpos(2)-75 cpos(3)*2/3 cpos(4)]);
    end
    if plotnum==1
        %if first plot plot labels
        ylabel(hax,'sorted trial #')
       % origpos(1)=origpos(1)+10;
        set(hax, 'Units','Pixels','Position',  [origpos(1) origpos(2) colorsize(1) colorsize(2)]);
    else
        set(hax,'Ytick',[])
        set(hax, 'Units','Pixels','Position',  [origpos(1) origpos(2) colorsize(1) colorsize(2)]);
    end
    if ~isempty(cscale)
        if ~isnan(cscale(1)) 
            caxis(hax,cscale); 
        else
            stds=nanstd(da,0,2);
            means=nanmean(da,2);
            stds=nanmean(stds);
            means=nanmean(means);
            cscale(1)=means-0.5*stds;
            cscale(2)=means+1.75*stds;
            if iscsc
                cscale(1)=means-.25*stds;
                cscale(2)=means+3*stds;
            end
            if ~any(isnan(cscale))
            caxis(hax,cscale)
            end
        end
        if ~iscsc
            if matlabver>2013
            set(h1, 'ylim', [cminshow cscale(2)])
            end
        end
    end
        set(hax,'YDir','reverse')        %flip y-axis values so first trial on top

    hold(hax,'off')
    
%end plotting functions
end

end




