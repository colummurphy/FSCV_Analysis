function [xcovdata]=xvardata(data,rates,plotparam,hax,seltrials,varargin)
%1/28/2019 UPDATES
%Data provided should already shifted to alignment event
%08/04/2018 cross covariance trialbytrial data
%supply pair of data{1} lfp and data{2} da to look at cross-variance between 2
%supply task events for window periods to look at trial by trial (in secon
%if only one list of events, then 2nd is assumed to be alignidx
%data is in the form of rows (trials) x columns (time)
        %markers{1} is fix cue appear
        %markers{2} is target cue appear
        %markers{3} is eye fix
%11/24/2018 look for max lags and max anticor lags (neg peaks)
%1/10/2019 SERIOUS PROBLEMS
%Used itrial id instead of tid to select winids
%middle of win ids no longer aligned since selecting min of trial id win

fontsize=12;
plotnum=0;
numtrials=size(data{1},1);
if isempty(data{1})
    numtrials=size(data{2},1);
end
if isempty(seltrials)
    seltrials=1:numtrials;
end
argnum=1;
noplot=0;
sitename={};
freqband=[0 0];     %freq band of filtered signal
xcovtype='coeff';       %default xcov type argument
eventnames={};
shuffle=0;
excludethres=.30;       %default percentage of data nan less than this ok
xcovdata={};
splitwinflag=0;     %split window in half to look at cov before/after peak separately
firstplot=0;
win=[-1 1];     %+/- 1 s from alignidx provided 1/10/19
single=0;
tcovflag=0;
%coeff - normalizes the sequence so that the covariances at
%zero lag are identically 1.0.
%unbiased - scales the raw covariance by 1/(M-abs(k)), where k
%is the index into the result.    
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'firstplot'
            firstplot=1;
        case 'cscale'
            argnum=argnum+1;
            plotparam.cscale=varargin{argnum};
        case 'sitename'
            argnum=argnum+1;        %fscv ch site name
            sitename=varargin{argnum};
        case 'fbands'
            %provide pass bands for each site
            %sitename must be defined by this point
            argnum=argnum+1;
            fbands=varargin{argnum};
            %assume first data is lfp
            bandid=find(strcmp(sitename{1},fbands)==1);
            if ~isempty(bandid)
                freqband=fbands{bandid+1};         
            end
        case 'noplot'
            noplot=1;   %don't plot
        case 'type'
            argnum=argnum+1;
            xcovtype=varargin{argnum};           
        case 'eventnames'
            argnum=argnum+1;
            eventnames=varargin{argnum};
        case 'shuffle'
            %shuffle trials of data 2, keep data 1 same
            shuffle=1;
            noplot=1;
        case 'excludethres'
            %percent of data nan, need to exluce if above this
            argnum=argnum+1;
            excludethres=varargin{argnum};
        case 'xwin'
            argnum=argnum+1;
            win=varargin{argnum};
        case 'splitwin'
            splitwinflag=1;     %split window in half to look at cov before/after peak separately  
        case 'single'
            single=1;           %single sample type provided
        case 'tcov'
            %trial covariance of pairs of time points
            tcovflag=1;
    end
    argnum = argnum + 1;
end
%resample data to higher sampling rate so equal sample rates
numsamples(1)=size(data{1},2);
numsamples(2)=size(data{2},2);
resampled={};
resampled{2}=data{2};
resampled{1}=data{1};
ts=[];
if single && ~isempty(resampled{1})
ts=[0:1/10:numsamples(1)/rates(1)-1/rates(1)];
else
ts=[0:1/rates(2):numsamples(2)/rates(2)-1/rates(2)];
end
if ~single
if rates(1)>rates(2)
    %downsample data 1 (lfp) to da (data 2 10 hz)
    resampled{1}=[];
    for itrial=1:numtrials
        downsamp = downsampleSingleBatch(data{1}(itrial,:),round(rates(1)/rates(2)));
        downsamp=[downsamp 0];
        resampled{1}(itrial,:)=downsamp;        
    end
end
elseif numsamples(1)>0
    %lfp provided, no da, downsample to 10 Hz default
        resampled{1}=[];
    for itrial=1:numtrials
        downsamp = downsampleSingleBatch(data{1}(itrial,:),round(rates(1)/10));
        downsamp=[downsamp 0];
        resampled{1}(itrial,:)=downsamp;        
    end
end
fs=rates(2);
xcovdata.rate=fs;
alnidx=plotparam.alnts.*fs;
%make sure all same length windows across trials;
winids=alnidx+win(1)*fs+1:alnidx+win(2)*fs+1;
badtrials=[];
meanslfp1=nanmean(resampled{1},1);
meanslfp=nanmean(meanslfp1);
thresoutlier=meanslfp*3.5;            %Max mean trial signal acceptable in LFP
maxpeak=abs(meanslfp*20);                      %max peak acceptable in lfp
meanstrialslfp=nanmean(resampled{1},2);
outtrial=find(abs(meanstrialslfp)>abs(thresoutlier));
resampled{1}(outtrial,:)=nan;           %MAKE NAN ANY LFP OUTLIERS
for itrial=1:numtrials
    %smooth out any nan or imag data if less than exclude thres
    dadata=[];
    lfpdata=[];
    if ~isempty(resampled{2})
        dadata=resampled{2}(itrial,:);
    end
    if ~isempty(resampled{1})
        lfpdata=resampled{1}(itrial,:);
    end
    if ~isempty(dadata)
   if any(imag(dadata(winids)))
       %if any imaginary component means data had been cut off from
       %start or end of recording, cannot smooth, skip trial, fill all nan
       dadata=nan(1,length(dadata));
       badtrials=[badtrials itrial];
   end
   if any(isnan(dadata(winids)))
        if length(find(isnan(dadata(winids))==1))...
                /length(dadata(winids))<=excludethres
            [dadata,TF] = fillmissing(dadata,'pchip');
        else
            %if > 25% of trial data is nan, make all nan 
            dadata=nan(1,length(dadata));
            badtrials=[badtrials itrial];
        end            
   end
    if any(abs(dadata(winids))>12*nanstd(resampled{2}(itrial,winids)))
        %sometimes filling data with pchip creates artificial
        %slopes and boundaries, should remove.
        dadata=nan(1,length(dadata));
        badtrials=[badtrials itrial];
    end
    %{
    if any(abs(dadata(winids))>abs(10*nanmean(resampled{2}(itrial,winids))))
        %sometimes filling data with pchip creates artificial
        %slopes and boundaries, should remove.
        dadata=nan(1,length(dadata));
        badtrials=[badtrials itrial];
    end
    %}
    if nanstd(dadata(winids))>10
        %sometimes filling data with pchip creates artificial
        %slopes and boundaries, should remove.
        %stdbad=[stdbad itrial];
        dadata=nan(1,length(dadata));
        badtrials=[badtrials itrial];
    end
    
    resampled{2}(itrial,:)=dadata;          
    end
    
    if ~isempty(lfpdata)
    if sum(lfpdata(winids))==0 
        resampled{1}(itrial,:)=nan(1,length(lfpdata));
        badtrials=[badtrials itrial];
    end
    if any(lfpdata(winids)>maxpeak)
        resampled{1}(itrial,:)=nan(1,length(lfpdata));
        badtrials=[badtrials itrial];
    end
    end
end
xcovdata.lfpdata=resampled{1};
xcovdata.dadata=resampled{2};
xcovdata.alnidx=alnidx;
xcovdata.winids=winids;
tsplot=ts(xcovdata.winids);
tsplot=tsplot-tsplot(1);
tsx=[-fliplr(tsplot(2:end)) tsplot];
if tcovflag
    tsx=tsplot;
end
xcovdata.tsx=tsx;
xcovdata.freqband=freqband;
xcovdata.eventnames=eventnames;
xcovdata.sitename=sitename;
origseltrials=seltrials;
if ~isempty(badtrials)
    %remove bad trials from seltrials
    seltrials(ismember(seltrials,badtrials))=[];
end
xcovdata.seltrials=seltrials;

xcovdata.xcovda=nan(length(seltrials),length(winids)*2-1);
xcovdata.xcovpre=[];
xcovdata.xcovpost=[];


if ~single
if length(seltrials)<3
    xcovdata=[];
    return;
end
xcovda=nan(length(seltrials),length(winids)*2-1);
pdata=[];
if tcovflag
    xcovda=nan(length(winids),length(winids));
    pdata=nan(length(winids),length(winids));
end
xcovpre=[];
xcovpost=[];
[a,b]=sort(rand(1,length(seltrials)));
[a,b1]=sort(rand(1,length(seltrials)));
shufseltrials=seltrials(b);     %shuffled selectedtrials
ipre=1;
%ONLY PROCESS SELTRIALS AND STORE ONLY SELTRIALS
if ~tcovflag
    %normal cross-covariance
for itrial=1:length(seltrials)
    tid=seltrials(itrial);       
    data1=resampled{1}(tid,winids);
    data2=resampled{2}(tid,winids);
    if shuffle==1
        data2=resampled{2}(shufseltrials(itrial),winids);
    end    
   [xcovtemp,~]=xcov(data1,data2,xcovtype);    
   xcovda(itrial,:)=xcovtemp;
   if splitwinflag
       %split window in halves to look at cov beofre/after peak separately
       idh=round(length(data1)/2);
       ipre=idh;
       ipost=idh+1;
       if length(1:idh)<length(idh+1:length(data1))
           %uneven windows, make even
           ipre=ipre+1;
       elseif length(1:idh)>length(idh+1:length(data1))
           ipre=ipre-1;
       end
       [prexcov,~]=xcov(data1(1:ipre),data2(1:ipre),xcovtype);    
       [postxcov,~]=xcov(data1(ipost:end),data2(ipost:end),xcovtype); 
       xcovpre(itrial,:)=prexcov;
       xcovpost(itrial,:)=postxcov;
   end
end
end

if tcovflag
    %trial series covariance of pairs of time points across signals
    avgwin=1;
    avgwin=0;
    for it=1:length(winids)
        sampids=winids(it)-avgwin:winids(it)+avgwin;
        sampids=sampids(sampids>0 & sampids<=winids(end));
        resamp1=nanmean(resampled{1}(seltrials,sampids),2);
        data1=resamp1-nanmean(resamp1);   %beta
        for itt=1:length(winids)
            sampids2=winids(itt)-avgwin:winids(itt)+avgwin;
            sampids2=sampids2(sampids2>0 & sampids2<=winids(end));
            resamp2=nanmean(resampled{2}(seltrials,sampids2),2);
            data2=resamp2-nanmean(resamp2); %da
            [xcovtemp,ptemp]=corr(data1,data2);    
            xcovda(it,itt)=xcovtemp;
            pdata(it,itt)=ptemp;
        end
    end
end


xcovdata.xcovda=xcovda;
xcovdata.pdata=pdata;
xcovdata.xcovpre=xcovpre;
xcovdata.xcovpost=xcovpost;

if ~tcovflag
    
tsplot=ts(xcovdata.winids(1,1:ipre));
tsplot=tsplot-tsplot(1);
tshalf=[-fliplr(tsplot(2:end)) tsplot];
xcovdata.tshalf=tshalf;

%get significant peak lag time points for each trial and their value at max
xcovlag=[];
xcovlaganti=[];     %neg peaks
pad=round(.05*size(xcovda,2));
actids=pad+1:size(xcovda,2)-pad;    %ids to search before padding
for itrial=1:length(seltrials)
    if any(isnan(xcovda(itrial,:)))
        %if signal nanned out because too many nans , store nan
        xcovlag(itrial,1)=nan;            %store lag at max in seconds
        xcovlag(itrial,2)=nan;
        xcovlaganti(itrial,1)=nan;
        xcovlaganti(itrial,2)=nan;
    else        
        %11/22/2018, find pos/neg peaks within each tiral
        [maxx,maxt]=max(xcovda(itrial,actids));     %max pos signal 11/22/2018
        maxt=maxt+actids-1;      %lag of max
        [minn,mint]=min(xcovda(itrial,actids));
        mint=mint+actids-1;        %lag of min
        if length(maxt)>1   
            %get closest max to zero lag
            cc=find(abs(maxt-median(1:size(xcovda,2)))==min(abs(maxt-median(1:size(xcovda,2)))));
            xcovlag(itrial,1)=tsx(maxt(cc));
            xcovlag(itrial,2)=xcovda(itrial,maxt(cc));
        else
            xcovlag(itrial,1)=tsx(maxt);            %store lag at max in seconds
            xcovlag(itrial,2)=xcovda(itrial,maxt);  %store value at max in xcov
        end
        if length(mint)>1 
            cc=find(abs(mint-median(1:size(xcovda,2)))==min(abs(mint-median(1:size(xcovda,2)))));
            xcovlaganti(itrial,1)=tsx(mint(cc));
            xcovlaganti(itrial,2)=xcovda(itrial,mint(cc));
        else
            xcovlaganti(itrial,1)=tsx(mint);            %store lag at max in seconds
            xcovlaganti(itrial,2)=xcovda(itrial,mint);  %store value at max in xcov
        end
    end
end
xcovdata.xcovlag=xcovlag;
xcovdata.xcovlaganti=xcovlaganti;

if ~isempty(xcovpre)
%get significant peak lag time points for each trial and their value at
%peak 
xcovlagpre=[];
for itrial=1:length(seltrials)
    if any(isnan(xcovpre(itrial,:))) 
        %if signal nanned out because too many nans , store nan
        xcovlagpre(itrial,1)=nan;            %store lag at max in seconds
        xcovlagpre(itrial,2)=nan;
    else
        maxx=max(abs(xcovpre(itrial,:)));        %max absolute of signal
        maxt=find(abs(xcovpre(itrial,:))==maxx);      %lag of max
        if length(maxt)>1
            xcovlagpre(itrial,1)=nan;
            xcovlagpre(itrial,2)=nan;
        else
        xcovlagpre(itrial,1)=tshalf(maxt);            %store lag at max in seconds
        xcovlagpre(itrial,2)=xcovpre(itrial,maxt);  %store value at max in xcov
        end
    end
end
%get lags for positive/negative coeff's
postargs=find(xcovlagpre(:,2)>0);
negtargs=find(xcovlagpre(:,2)<0);
xcovdata.lagsprepos=xcovlagpre(postargs,1);
xcovdata.lagspreneg=xcovlagpre(negtargs,1);
xcovdata.xcovlagpre=xcovlagpre;
end

if ~isempty(xcovpost)
%get significant peak lag time points for each trial and their value at
%peak 
xcovlagpost=[];
for itrial=1:length(seltrials)
    if any(isnan(xcovpost(itrial,:))) 
        %if signal nanned out because too many nans , store nan
        xcovlagpost(itrial,1)=nan;            %store lag at max in seconds
        xcovlagpost(itrial,2)=nan;
    else
        maxx=max(abs(xcovpost(itrial,:)));        %max absolute of signal
        maxt=find(abs(xcovpost(itrial,:))==maxx);      %lag of max
        if length(maxt)>1
            xcovlagpost(itrial,1)=nan;
            xcovlagpost(itrial,2)=nan;
        else
        xcovlagpost(itrial,1)=tshalf(maxt);            %store lag at max in seconds
        xcovlagpost(itrial,2)=xcovpost(itrial,maxt);  %store value at max in xcov
        end
    end
end
%get lags for positive/negative coeff's
postargs=find(xcovlagpost(:,2)>0);
negtargs=find(xcovlagpost(:,2)<0);
xcovdata.lagspostpos=xcovlagpost(postargs,1);
xcovdata.lagspostneg=xcovlagpost(negtargs,1);
xcovdata.xcovlagpost=xcovlagpost;
end

if shuffle==1
    xcovdata.shuftrials=shufseltrials;
end
end

if noplot==0    
    if firstplot
        plotxtrials(hax,xcovdata,'leftlabel');
    else
        plotxtrials(hax,xcovdata);
    end
end
end
end




