function [plotdata,outdata,info]=setTrialAx(data,plotparam,hax,seltrials,varargin)
%1/26/2019 STILL NEED TO UPDATE OT USE ttarg, tfix vars for nlx signals
%NO SMOOTHING, CONVERT TO ENVELOPE
%Update massive 1/24/2019
%use new trialbytrial(1).ttarg, tfix, etc. for nlx signals markers rather
%than samplesfix, etc.
%baseline sub around aln event +/-1 rather than before with padding
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
smoothdata=0;
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
ploteyedist=0;
ploteyev=0;
getfft=0;
logscale=0;
eye=0;
blink=0;
noplot=0;
sitename={};
winphys=0.1;            %smoothing window for lick/etc.phys sign
rclfp=[];
basetrial=0;
info=[];
alnevt='targeye';
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
means=[];
scales=[];
if isfield(plotparam,'triallabel')
    triallabel=plotparam.triallabel;
end
shiftdata=[];
if isfield(plotparam,'shift')
    %shift data to targ or fix marker
    shiftdata=plotparam.shift;
end
colorsize=[500 500];
vertplotwidth=[500 500];
colorm=[0 0 0];
chnum=0;
samplespersec=ceil(plotparam.samplespersec);
tslfp=[];
%fscv default event markers in da samples
win=plotparam.win;      %in fscv samples
rewardidx=plotparam.alignidx;
idsfixeye=plotparam.samplesfixeye{1};
idsfix=plotparam.samplesfix{1};
idstarg=plotparam.samplestarg{1};
idstargeye=plotparam.samplestargeye{1};

while argnum <= length(varargin)        
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
            smoothdata=1;
        case 'smooth'
            smoothdata=1;            
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
            getfft=0;
        case 'info'
            argnum=argnum+1;
            info=varargin{argnum};
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
            cscid=find(ismember(data.cscNames,varargin{argnum}));
            if strcmp(varargin{argnum},'eyedist') || strcmp(varargin{argnum},'eyev')
                cscid=find(ismember(data.cscNames,'eyex'));
            end
           % da=10*log10(data.lfp{cscid});
            if isequal(varargin{argnum},'eyed') || isequal(varargin{argnum},'blink')
               cscid=find(ismember(data.cscNames,'eyed'));
               if isfield(data,'lfp')
                da=-data.lfp{cscid};
               end
                phys=1;
                eye=1;
                if isequal(varargin{argnum},'blink')
                    blink=1;
                end
            elseif strcmp(varargin{argnum},'eyedist') || ...
                    strcmp(varargin{argnum},'eyev') || ...
                    strcmp(varargin{argnum},'eyex')
                %eye x data
                if isfield(data,'lfp')
                da=(data.lfp{cscid});
                end
            else
                %lfp
                if isfield(data,'lfp')
               da=(data.lfp{cscid}).*1e6;           %convert to microvolts
                end
               islfp=1;
           end
           if isequal(varargin{argnum},'pulse')
               pulse=1;
           end
           if isequal(varargin{argnum},'eyex')
               eyex=1;
           end
           if isequal(varargin{argnum},'eyedist')
               eyex=1;
               ploteyedist=1;
           end
           if isequal(varargin{argnum},'eyev')
               eyex=1;
               ploteyev=1;
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
    argnum = argnum + 1;
end
trialtype=[plotparam.trialtype ' | ' triallabel];
if ~noplot
colorsize=plotparam.colorsize;
colorm=plotparam.markercolors;
chnum=plotparam.chnum;
end
if ~isempty(rclfp)
    %csc data provided from previous export, already reconverted, get
    %new sample rates
    samplespersec=ceil(data.ratelfp);       %default dc  csc sample rate
    if ~isempty(info)
        samplespersec=info.samplespersec;
    end
    tslfp=data.relts;           %original relative ts's for csc
    dcrate=ceil(length(tslfp)./length(rclfp));      %downconvert rate for ts 
    tslfpd=tslfp(1:dcrate:end);
    winneg=round(plotparam.alignidx-plotparam.win(1))/plotparam.samplespersec;
    winpos=round(plotparam.win(end)-plotparam.alignidx)/plotparam.samplespersec;
    rewardidx=round(plotparam.alignidx.*samplespersec./plotparam.samplespersec);
    winnew=rewardidx-round(winneg*samplespersec):rewardidx+round(winpos*samplespersec); 
    win=winnew;             %replace fscv domain window for plotting ids
    %convert marker id's to new sample rate
    idsfixeye=round(plotparam.tfixeye.*samplespersec);
    idsfix=round(plotparam.tfix.*samplespersec);
    idstarg=round(plotparam.ttarg.*samplespersec);
    idstargeye=round(plotparam.ttargeye.*samplespersec);
    da=rclfp;
end
    
if iscsc && isempty(rclfp)
    samplespersec=ceil(data.ratelfp);
    tslfp=data.relts;           %original relative ts's for csc
    winneg=round(plotparam.alignidx-plotparam.win(1))/plotparam.samplespersec;
    winpos=round(plotparam.win(end)-plotparam.alignidx)/plotparam.samplespersec;
    rewardidx=round(plotparam.alignidx.*samplespersec./plotparam.samplespersec);
    winnew=rewardidx-round(winneg*samplespersec):rewardidx+round(winpos*samplespersec); 
    win=winnew;             %replace fscv domain window for plotting ids
    %convert marker id's to new sample rate
    idsfixeye=round(plotparam.tfixeye.*samplespersec);
    idsfix=round(plotparam.tfix.*samplespersec);
    idstarg=round(plotparam.ttarg.*samplespersec);
    idstargeye=round(plotparam.ttargeye.*samplespersec);
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
            dafilt=filterLFP(da(ii,:),samplespersec,bandfilt);
            %square & envelope signal based on default win length
            
            if envwidth==0
                winlength=round(samplespersec*.5/mean(bandfilt));
            else
                winlength=envwidth.*samplespersec; %user provided envelop width
            end
            
            dasquared=dafilt.^2;   %get power V^2 (uv^2)
            if smoothdata
                datemp=envwave(dasquared);
                %Returns nan automatically if no extrama found warning
                if isnan(datemp)
                    datemp=dasquared;       %use un-env data
                end
                daenv=smoothwin(datemp,winlength);   %smoothing 
            else
            %REMOVE SMOOTHING 1/30/2019, replace with envelope
                daenv=envwave(datemp);
            end
            da(ii,:)=daenv;      %convert to microvolts squared
        end
    end    
end
info.envwidth=envwidth;
origeye=[];
if isequal(cscname,'eyed')
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
eyed=[];
eyeusacs=[];
eyev=[];
if eyex==1 && isempty(rclfp)
   datatemp=[];
  for itrial=1:size(da,1)
      datatemp=deglitchnanamp(da(itrial,:),2.5e-3,30);
      %datatemp=da(itrial,:);
      da(itrial,:)=datatemp;
      smootheye=smoothwin(datatemp,20);
      eyevel=diff(smootheye);
      eyevel(isnan(eyevel))=0;
      abseyevel=abs(eyevel);
      eyedist=cumtrapz(abseyevel);
      thres=nanmean(abseyevel);
      maxlim=nanstd(abseyevel);
      [pks,locs,w,p] = findpeaks(abseyevel,'minpeakwidth',...
            round(samplespersec*.01),'minpeakdistance',round(samplespersec*.05)...
            ,'MinPeakprominence',thres);
        sacslogic=zeros(1,length(datatemp));
        locsbelowlim=locs(pks<=maxlim);
        sacslogic(locsbelowlim)=1;
      eyeusacs(itrial,:)=sacslogic;
      eyed(itrial,:)=eyedist;
      eyev(itrial,:)=eyevel;
   end
  %filter data
    sitename='eyex';
    if ploteyedist
        sitename='eyedist';
        da=eyed;
    end
    if ploteyev
        sitename='eyev';
        da=abs(eyev);
    end
end
if lick==1 && isempty(rclfp)
    %filter
    for ii=1:size(da,1)
    datemp=filterLFP(da(ii,:),samplespersec,[0 100]);       
    datemp=smoothwin(datemp,winphys);   %smoothing
    da(ii,:)=datemp;
    end
end
if blink==1  && isempty(rclfp)
    sitename='blink';
    if isempty(rclfp)
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
    tslfp=data.relts;           %original relative ts's for csc
    dcrate=ceil(length(tslfp)./length(rates));      %downconvert rate for ts 
    tslfp=tslfp(1:dcrate:end);
    winneg=round(plotparam.alignidx-plotparam.win(1))/plotparam.samplespersec;
    winpos=round(plotparam.win(end)-plotparam.alignidx)/plotparam.samplespersec;
    rewardidx=round(plotparam.alignidx.*samplespersec./plotparam.samplespersec);
    winnew=rewardidx-round(winneg*samplespersec):rewardidx+round(winpos*samplespersec); 
    win=winnew;             %replace fscv domain window for plotting ids
    %convert marker id's to new sample rate
    idsfix=round(idsfix./ceil(data.ratelfp)*samplespersec);
    idstarg=round(idstarg./ceil(data.ratelfp)*samplespersec);
    idsfixeye=round(idsfixeye./ceil(data.ratelfp)*samplespersec);  
    idstargeye=round(idstargeye./ceil(data.ratelfp)*samplespersec);  
    da=rates;
    end
end

if pulse==1 && isempty(rclfp)
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
    samplespersec=1/stepsize;       %DEFINE NEW SAMPLE RATE
    tslfp=pulsets(1,:);
    tslfp=data.relts;           %original relative ts's for csc
    dcrate=ceil(length(tslfp)./length(rates));      %downconvert rate for ts 
    tslfp=tslfp(1:dcrate:end);
    winneg=round(plotparam.alignidx-plotparam.win(1))/plotparam.samplespersec;
    winpos=round(plotparam.win(end)-plotparam.alignidx)/plotparam.samplespersec;
    rewardidx=round(plotparam.alignidx.*samplespersec./plotparam.samplespersec);
    winnew=rewardidx-round(winneg*samplespersec):rewardidx+round(winpos*samplespersec); 
    win=winnew;             %replace fscv domain window for plotting ids
    %convert marker id's to new sample rate
    idsfix=round(idsfix./ceil(data.ratelfp)*samplespersec);
    idstarg=round(idstarg./ceil(data.ratelfp)*samplespersec);
    idsfixeye=round(idsfixeye./ceil(data.ratelfp)*samplespersec);   
    idstargeye=round(idstargeye./ceil(data.ratelfp)*samplespersec);  
    da=rates;
end
info.ts=tslfp;
info.samplespersec=samplespersec;
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


pade=round(.2*samplespersec);       %padding before subseq event
basesamples=round(.3*samplespersec);        %.3 s for baseline average;
basepad=1;              %+/- .1 s around aling idx for baseline signal
outdata=da;         %output signal, downconverted LFP, before shifting, all trials
rewext=3;       %extend baseline for reward baseline if too nan
rewext2=5;      %max allowed extension
%info.alnevt=alnevt;
%info.idsaln=idsaln;
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
if ~iscsc
    %since concistent minima occurs at eye on targ/fix, use aln event +/-   
    %rather than prior
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
    if basetrial
        %mean entire trial for baseline sub
        baseline=nanmean(da(:,win(1):win(end)),2);
        baselineall=repmat(baseline,1,size(da,2));
        dasub=da-baselineall;
    end
    da=dasub;
end
%CALCULAte averages for given intervals for each trial
%if fft flag, get avg spectra
%ie post-reward 0.2 - 2.5 s, post-target to reward makresr{3} to -0.2 s
%fix to target (idstarg to marks{3} relative to baseline which is 0.5 s
%average before previous event
postrewdur=round(samplespersec*5); %5 s post reward interval to average
postrewdurshort=round(samplespersec*3); %3s post rew
postrewdurim=round(samplespersec*1.5); %1.5 post rew
targimdur=round(samplespersec*1);      %1s post targ immediate
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
if iscsc && isempty(rclfp) && getfft
datawin=data.lfp{cscid}(1,1:samplespersec*5).*1e6;     %convert to microvolts first            
[S,winmid,f]=lfp_mtspecgram2(datawin',fftwin,ffttapers,0,samplespersec,fftfreqlim,0,1,rmdc,norm);
%[S,f]=chronux_mtspectrumc(datawin',ffttapers,0,samplespersec,fftfreqlim,0,1);  
plotdata.fftf=f;        %default freq variable
end

%initialize all variables first QQQQQQQQQQQQQQQQQQQQQQ
postfixwin=nan(1,numtrials);
plotdata.basewin=nan(1,numtrials);     
postfixpeak=nan(1,numtrials);
fixpeakidx=nan(1,numtrials);
posttargwin=nan(1,numtrials);
posttargpeak=nan(1,numtrials);
posttargpeakabs=nan(1,numtrials);

targpeakidx=nan(1,numtrials);
targpeakabsidx=nan(1,numtrials);

targimwin=nan(1,numtrials);
targimpeak=nan(1,numtrials);
targimpeakidx=nan(1,numtrials);
postrewwin=nan(1,numtrials);
postrewpeak=nan(1,numtrials);
rewpeakidx=nan(1,numtrials);
rewshortpeakidx=nan(1,numtrials);
rewshortpeak=nan(1,numtrials);
rewshortwin=nan(1,numtrials);
rewimpeak=nan(1,numtrials);
rewimpeakidx=nan(1,numtrials);
rewimwin=nan(1,numtrials);  

if islfp && getfft
for xx=1:numtrials
    plotdata.fixwinfft{xx}=nan;
    plotdata.prefixwinfft{xx}=nan;
    plotdata.targwinfft{xx}=nan;
    plotdata.targimwinfft{xx}=nan;
    plotdata.rewwinfft{xx}=nan;
end
plotdata.fixwinfftavg=nan(numtrials,length(plotdata.fftf));
plotdata.fixpostfftavg=nan(numtrials,length(plotdata.fftf));
plotdata.targprefftavg=nan(numtrials,length(plotdata.fftf));
plotdata.prefixwinfftavg=nan(numtrials,length(plotdata.fftf));
plotdata.targwinfftavg=nan(numtrials,length(plotdata.fftf));
plotdata.targpostfftavg=nan(numtrials,length(plotdata.fftf));  
plotdata.targimwinfftavg=nan(numtrials,length(plotdata.fftf));
plotdata.targimpostfftavg=nan(numtrials,length(plotdata.fftf));
plotdata.rewwinfftavg=nan(numtrials,length(plotdata.fftf));
plotdata.rewpostfftavg=nan(numtrials,length(plotdata.fftf));
end

for xx=1:numtrials
%process only selected trials
%idsfixeye eye fixed, idstarg target appears, markers{1} fix cue
tid=seltrials(xx);
if badtrls(tid)==0 && ~isnan(idsfix(tid))
    %not bad trial as already designated above
    %%%%%
    %FIX CHARACTERISTICS
    %baseidx=idsfixeye(tid);     %TAKE OUT EYE ALN 2/2/2019
    baseidx=idsfix(tid);
    if eyex
        %eye x signal need to take before cue appears
        baseidx=idsfix(tid)-pade;
    end
    baselineids=baseidx-basepad:baseidx+basepad;
    if basetrial==1
        %uniform baseline for all events that is entire trial
        baselineids=idsfix(tid):idsrew(tid)+postrewdurshort;
    end
    fixbaseline=nanmean(da(tid,baselineids));
    if iscsc && ~ploteyedist
        %need baseline sub for eye distance, not for other eye x vars
        fixbaseline=0;     %10/05/2018 change so LFPs also baseline
    end
    plotdata.basewin(xx)=fixbaseline;
    %postfixids=baseidx:idstarg(tid)-pade;        %12/19/2018 - pade
    postfixids=baseidx+pade:idstarg(tid)-pade;        %12/19/2018 add pade
    postfixwin(xx)=nanmean(da(tid,postfixids))-fixbaseline;
    sigda=da(tid,postfixids)-fixbaseline;
    [postfixmax,postfixpeakidxrel]=max(sigda);
    nanratio=length(find(isnan(da(tid,postfixids))==1))/length(da(tid,postfixids));
    if ~any(isnan(postfixmax)) && nanratio<=.5
        %not too much nan
        postfixpeakidx=postfixpeakidxrel+baseidx+pade-1;        %absolute index
        if length(postfixpeakidx)>1
            postfixpeakidx=postfixpeakidx(1);
        end
        postfixpeak(xx)=da(tid,postfixpeakidx)-fixbaseline; %relative max change
        fixpeakidx(xx)=postfixpeakidx-baseidx;     %relative offset from base idx
    end
    if islfp && getfft
        %get fft spectra task modulated
        datawin=data.lfp{cscid}(tid,postfixids).*1e6;     %convert to microvolts first                        
        if fftwin(1)*samplespersec<size(datawin,2)
            [S,winmid,f]=lfp_mtspecgram2(datawin',fftwin,ffttapers,0,samplespersec,fftfreqlim,0,1,rmdc,norm);
            plotdata.fixwinfft{xx}=S;
            plotdata.fixwinfftavg(xx,:)=mean(S,1);
            plotdata.fixpostfftavg(xx,:)=S(1,:);     %1 s window immediately after event
        end
        %[S,f]=chronux_mtspectrumc(datawin',ffttapers,0,samplespersec,fftfreqlim,0,1);  
        pretargetids=idstarg(tid)-fftwinids:idstarg(tid);
        datawin=data.lfp{cscid}(tid,pretargetids).*1e6;     %convert to microvolts first            
        if fftwin(1)*samplespersec<=size(datawin,2)
            [S,winmid,f]=lfp_mtspecgram2(datawin',fftwin,ffttapers,0,samplespersec,fftfreqlim,0,1,rmdc,norm);
            plotdata.targprefftavg(xx,:)=S(end,:);        %immediately prior target
        end
        baselinefftids=baseidx-pade-fftwinids:baseidx-pade;
        if ~any(baselinefftids<1)
            datawin=data.lfp{cscid}(tid,baselinefftids).*1e6;     %convert to microvolts first            
            [S,winmid,f]=lfp_mtspecgram2(datawin',fftwin,ffttapers,0,samplespersec,fftfreqlim,0,1,rmdc,norm);
            plotdata.prefixwinfft{xx}=S;
            plotdata.prefixwinfftavg(xx,:)=mean(S,1);
        end
    end
    
    %%%%
    %TARGET CHARACTERISTICS
    %targidx=idstargeye(tid);        %TAKE OUT EYE ALN 2/2/2109
    targidx=idstargeye(tid);        %change back from targ to targeye 2/4/2019
    if eyex
        %eye x or eye dist signal need to take before cue appears
        targidx=idstarg(tid)-pade;
    end
    targbaselineids=targidx-basepad:targidx+basepad;
    if basetrial==1
        %uniform baseline for all events that is entire trial
        targbaselineids=baselineids;
    end
    if ~any(targbaselineids<1) && ~strcmp(plotparam.trialtype,'fixbreak') &&...
            ~isnan(fixbaseline) && ~any(isnan(targbaselineids))  
        targbaseline=nanmean(da(tid,targbaselineids));
        if iscsc && ~ploteyedist
            targbaseline=0;
        end
        targids=targidx+pade:idsrew(tid)-pade;     %'dwell' target win w/ padding 
        targimids=targidx:targidx+targimdur;     %immediate post
        posttargwin(xx)=nanmean(da(tid,targids))-targbaseline;
        targimwin(xx)=nanmean(da(tid,targimids))-targbaseline;
        targdasub=da(tid,targids)-targbaseline;
        targdasubim=da(tid,targimids)-targbaseline;
        [posttargmax, targmaxidxrel]=max(targdasub);
        [posttargmaxabs,targmaxabsidxrel]=max(abs(targdasub));
        posttargmaxabs=targdasub(targmaxabsidxrel);     %signed abs max
        [targimmax, targimmaxidxrel]=max(targdasubim);
       nanratio=length(find(isnan(da(tid,targids))==1))/length(da(tid,targids));
       if ~any(isnan(posttargmax)) && nanratio<=.5
            posttargpeakidx=targmaxidxrel+targidx+pade-1;        %absolute index
            posttargpeakabsidx=targmaxabsidxrel+targidx+pade-1;        %absolute index
            posttargimpeakidx=targimmaxidxrel+targidx-1;        %absolute index
            if length(posttargpeakidx)>1
                posttargpeakidx=posttargpeakidx(1);
            end
            if length(posttargpeakabsidx)>1
                posttargpeakabsidx=posttargpeakabsidx(1);
            end
            if length(posttargimpeakidx)>1
                posttargimpeakidx=posttargimpeakidx(1);
            end
            posttargpeakabs(xx)=da(tid,posttargpeakabsidx)-targbaseline;
            posttargpeak(xx)=da(tid,posttargpeakidx)-targbaseline; %relative max change
            targpeakidx(xx)=posttargpeakidx-targidx;     %relative idx from preceding event
            targpeakabsidx(xx)=posttargpeakabsidx-targidx;
            if ~isempty(posttargimpeakidx)
                targimpeak(xx)=da(tid,posttargimpeakidx)-targbaseline;
                targimpeakidx(xx)=posttargimpeakidx-targidx;
            end            
        end
        if islfp && getfft
            %get fft spectra task modulated
            datawin=data.lfp{cscid}(tid,targids).*1e6;     %convert to microvolts first     
            if fftwin(1)*samplespersec<=length(datawin)
                [S,winmid,f]=lfp_mtspecgram2(datawin',fftwin,ffttapers,0,samplespersec,fftfreqlim,0,1,rmdc,norm);
                plotdata.targwinfft{xx}=S;
                plotdata.targwinfftavg(xx,:)=mean(S,1);
                plotdata.targpostfftavg(xx,:)=S(1,:);
            end            
            datawin=data.lfp{cscid}(tid,targimids).*1e6;     %convert to microvolts first     
            if fftwin(1)*samplespersec<=length(datawin)
                [S,winmid,f]=lfp_mtspecgram2(datawin',fftwin,ffttapers,0,samplespersec,fftfreqlim,0,1,rmdc,norm);
                plotdata.targimwinfft{xx}=S;
                plotdata.targimwinfftavg(xx,:)=mean(S,1);
                plotdata.targimpostfftavg(xx,:)=S(1,:);
            end
        end
    end
    
    %%%%
    %OUTCOME CHARACTERISTICS
    %reward peak relative to pre reward baseline
    reidx=idsrew(tid);
     if eyex
        %eye x signal need to take before cue appears
        reidx=idsrew(tid)-pade;
    end
    rewbaseids=reidx-basepad:reidx+basepad;
    if basetrial==1
        %uniform baseline for all events that is entire trial
        rewbaseids=baselineids;
    end
    rewbaseline=nanmean(da(tid,rewbaseids));
    if isnan(rewbaseline)
        %try extending baseline for rew for DA signals
        if ~iscsc
            rewbaseids=reidx-basepad-rewext:reidx+basepad+rewext;
            rewbaseline=nanmean(da(tid,rewbaseids));
            if isnan(rewbaseline)
                rewbaseids=reidx-basepad-rewext2:reidx+basepad+rewext2;
                rewbaseline=nanmean(da(tid,rewbaseids));
            end
        end
    end        
    if iscsc && ~ploteyedist
        rewbaseline=0;
    end
    if ~isnan(rewbaseline)
        rewids=reidx+pade:reidx+postrewdur;
        postrewwin(xx)=nanmean(da(tid,rewids))-rewbaseline;
        rewdasub=da(tid,rewids)-rewbaseline;
        [postrewmax,postrewmaxidxrel]=max(rewdasub);     %change to abs 10/05/2018
        nanratio=length(find(isnan(da(tid,rewids))==1))/length(da(tid,rewids));
        if ~any(isnan(postrewmax)) && nanratio<=0.5
            postrewpeakidx=postrewmaxidxrel+reidx+pade-1;        %absolute index
            if length(postrewpeakidx)>1
                postrewpeakidx=postrewpeakidx(1);
            end
            postrewpeak(xx)=da(tid,postrewpeakidx)-rewbaseline; %relative max change
            rewpeakidx(xx)=postrewpeakidx-reidx;     %relative idx from preceding event
        end        
        if islfp && getfft
            %get fft spectra task modulated
            datawin=data.lfp{cscid}(tid,rewids).*1e6;     %convert to microvolts first  
            if fftwin(1)*samplespersec<=length(datawin)
                [S,winmid,f]=lfp_mtspecgram2(datawin',fftwin,ffttapers,0,samplespersec,fftfreqlim,0,1,rmdc,norm);
                plotdata.rewwinfft{xx}=S;
                plotdata.rewwinfftavg(xx,:)=mean(S,1);
                plotdata.rewpostfftavg(xx,:)=S(1,:);
            end
        end
        
        %reward peak short window (3s)
        rewids=reidx+pade:reidx+postrewdurshort;
        rewshortwin(xx)=nanmean(da(tid,rewids))-rewbaseline;
        rewdasub=da(tid,rewids)-rewbaseline;
        [postrewmax,postrewmaxidxrel]=max(rewdasub);
        if ~isnan(postrewmax)
            maxidx=postrewmaxidxrel+reidx+pade-1;        %absolute index
            if length(maxidx)>1
                maxidx=maxidx(1);
            end
            rewshortpeak(xx)=da(tid,maxidx)-rewbaseline; %relative max change
            rewshortpeakidx(xx)=maxidx-reidx;     %relative idx from preceding event
        end
        
        %reward peak relative to pre baseline immediate , no pad for im    
        rewids=reidx:reidx+postrewdurim;
        rewimwin(xx)=nanmean(da(tid,rewids))-rewbaseline;
        rewdasub=da(tid,rewids)-rewbaseline;
        [postrewmax,postrewmaxidxrel]=max(rewdasub);
        if ~isnan(postrewmax)
            maxidx=find(rewdasub==postrewmax)+reidx-1;        %absolute index
            if length(maxidx)>1
                maxidx=maxidx(1);
            end
            rewimpeak(xx)=da(tid,maxidx)-rewbaseline; %relative max change
            rewimpeakidx(xx)=maxidx-reidx;     %relative idx from preceding event
        end
        
    %end outcome    
    end
    
%only not badtrls    
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
plotdata.targpeak=posttargpeak;
plotdata.targpeakci=nanstd(plotdata.targpeak)/sqrt(length(find(~isnan(posttargpeak))))*1.96;  %95% confidence interval
plotdata.targpeakavg=nanmean(plotdata.targpeak);
plotdata.targpeakts=targpeakidx./samplespersec;
plotdata.targpeaktsci=nanstd(plotdata.targpeakts)/sqrt(length(find(~isnan(targpeakidx))))*1.96;  %95% confidence interval
plotdata.targpeaktsavg=nanmean(plotdata.targpeakts);

plotdata.targpeakabs=posttargpeakabs;
plotdata.targpeakci=nanstd(plotdata.targpeakabs)/sqrt(length(find(~isnan(posttargpeakabs))))*1.96;  %95% confidence interval
plotdata.targpeakavg=nanmean(plotdata.targpeakabs);
plotdata.targpeakabsts=targpeakabsidx./samplespersec;
plotdata.targpeakabstsci=nanstd(plotdata.targpeakabsts)/sqrt(length(find(~isnan(targpeakabsidx))))*1.96;  %95% confidence interval
plotdata.targpeakabstsavg=nanmean(plotdata.targpeakabsts);

plotdata.targimwin=targimwin;
plotdata.targimwinci=nanstd(plotdata.targimwin)/sqrt(length(find(~isnan(targimwin))))*1.96;  %95% confidence interval
plotdata.targimwinavg=nanmean(plotdata.targimwin);
plotdata.targimpeak=targimpeak;
plotdata.targimpeakci=nanstd(plotdata.targimpeak)/sqrt(length(find(~isnan(targimpeak))))*1.96;  %95% confidence interval
plotdata.targimpeakavg=nanmean(plotdata.targimpeak);
plotdata.targimpeakts=targimpeakidx./samplespersec;
plotdata.targimpeaktsci=nanstd(plotdata.targimpeakts)/sqrt(length(find(~isnan(targimpeakidx))))*1.96;  %95% confidence interval
plotdata.targimpeaktsavg=nanmean(plotdata.targimpeakts);
        
plotdata.rewwin=postrewwin;
plotdata.rewwinci=nanstd(plotdata.rewwin)/sqrt(length(find(~isnan(postrewwin))))*1.96;  %95% confidence interval
plotdata.rewwinavg=nanmean(plotdata.rewwin);
plotdata.rewpeak=postrewpeak;
plotdata.rewpeakci=nanstd(plotdata.rewpeak)/sqrt(length(find(~isnan(postrewpeak))))*1.96;  %95% confidence interval
plotdata.rewpeakavg=nanmean(plotdata.rewpeak);
plotdata.rewpeakts=rewpeakidx./samplespersec;
plotdata.rewpeaktsci=nanstd(plotdata.rewpeakts)/sqrt(length(find(~isnan(rewpeakidx))))*1.96;  %95% confidence interval
plotdata.rewpeaktsavg=nanmean(plotdata.rewpeakts);
            
plotdata.rewshortpeak=rewshortpeak;
plotdata.rewshortpeakci=nanstd(plotdata.rewshortpeak)/sqrt(length(find(~isnan(rewshortpeak))))*1.96;  %95% confidence interval
plotdata.rewshortpeakavg=nanmean(plotdata.rewshortpeak);
plotdata.rewshortpeakts=rewshortpeakidx./samplespersec;
plotdata.rewshortpeaktsci=nanstd(plotdata.rewshortpeakts)/sqrt(length(find(~isnan(rewshortpeakidx))))*1.96;  %95% confidence interval
plotdata.rewshortpeaktsavg=nanmean(plotdata.rewshortpeakts);
plotdata.rewshortwin=rewshortwin;
plotdata.rewshortwinci=nanstd(plotdata.rewshortwin)/sqrt(length(find(~isnan(rewshortwin))))*1.96;  %95% confidence interval
plotdata.rewshortwinavg=nanmean(plotdata.rewshortwin);

plotdata.rewimpeak=rewimpeak;
plotdata.rewimpeakci=nanstd(plotdata.rewimpeak)/sqrt(length(find(~isnan(rewimpeak))))*1.96;  %95% confidence interval
plotdata.rewimpeakavg=nanmean(plotdata.rewimpeak);
plotdata.rewimpeakts=rewimpeakidx./samplespersec;
plotdata.rewimpeaktsci=nanstd(plotdata.rewimpeakts)/sqrt(length(find(~isnan(rewimpeakidx))))*1.96;  %95% confidence interval
plotdata.rewimpeaktsavg=nanmean(plotdata.rewimpeakts);
plotdata.rewimwin=rewimwin;
plotdata.rewimwinci=nanstd(plotdata.rewimwin)/sqrt(length(find(~isnan(rewimwin))))*1.96;  %95% confidence interval
plotdata.rewimwinavg=nanmean(plotdata.rewimwin);

%plotdata.alnevt=alnevt;

%prepare for plotting
interval=plotparam.interval;        %interval xtick marks def 2.5s
cscale=plotparam.cscale;
cminshow=plotparam.cminshow;
xticks=1:interval*samplespersec:size(da(:,win),2);
%mintime=round((min(win)-median(win))./samplespersec);
mintime=round((win(1)-alnidx))/samplespersec;
maxtime=round((win(end)-alnidx))/samplespersec;
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
%separate from alnevt that sets default
if strcmp(shiftdata,'targ')
    %shift to targ cue appear, idstarg
    avgtarg=round(mean(idstarg));
    datemp=da;
    tempshift=[];
    for itrial=1:size(datemp,1)
        numshift(itrial)=avgtarg-idstarg(itrial);
        tempshift(itrial,:)=circshift(datemp(itrial,:),[0 numshift(itrial)]);
    end
    pdata=tempshift(seltrials,win(1):win(end));
end
if strcmp(shiftdata,'fix')
    %shift to fix cue appear, idsfix
    avgtarg=round(mean(idsfix));
    datemp=da;
    tempshift=[];
    for itrial=1:size(datemp,1)
        numshift(itrial)=avgtarg-idsfix(itrial);
        tempshift(itrial,:)=circshift(datemp(itrial,:),[0 numshift(itrial)]);
    end
    pdata=tempshift(seltrials,win(1):win(end));
end
end

if isempty(scales)
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

if noplot==0
    cla(hax);
   
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
        line(hax,[idsfix(seltrials(jj))-win(1)+numshift(seltrials(jj))],jj,'Color',colorm(1,:),'LineWidth',0.5,'LineStyle',':','marker','s','markersize',3)
        line(hax,[idsfixeye(seltrials(jj))-win(1)+numshift(seltrials(jj))],jj,'Color',colorm(2,:),'LineWidth',0.5,'LineStyle',':','marker','.','markersize',3)
        if ~strcmp(plotparam.trialtype,'fixbreak')
            if ~strcmp(plotparam.alnevt,'targ')
            line(hax,[idstarg(seltrials(jj))-win(1)+numshift(seltrials(jj))],jj,'Color',colorm(3,:),'LineWidth',0.5,'LineStyle',':','marker','o','markersize',3)
            end
            if ~strcmp(plotparam.alnevt,'targeye')
            line(hax,[idstargeye(seltrials(jj))-win(1)+numshift(seltrials(jj))],jj,'Color',colorm(4,:),'LineWidth',0.5,'LineStyle',':','marker','.','markersize',3)
            end
        end
       if ~strcmp(plotparam.alnevt,'outcome')
           %not aligned to outcome, mark
           line(hax,[idsrew(seltrials(jj))-win(1)+numshift(seltrials(jj))],jj,'Color',colorm(5,:),'LineWidth',0.5,'LineStyle',':','marker','sq','markersize',3)
       end

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
            title(hax,[trialtype ' | ' alnevt ' | ' sitename{:}])
        end
    else
        title(hax,[trialtype ' lfp beta ' cscname])
        band=[];
      band=bandfilt;
        title(hax,{[trialtype] ['lfp ' num2str(band(1)) '-' num2str(band(2)) ' Hz | ' alnevt] [cscname]})
        clabel='\beta-lfp (\muV^2)';
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
            cscale(1)=-scales*.5;
            cscale(2)=scales*4;
            if iscsc
                cscale(2)=scales*5;
                cscale(1)=scales*1;
                if means<0 && contains(sitename,'eye') 
                    cscale(1)=means-scales*2;
                    cscale(2)=means+scales*3;
                end
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




