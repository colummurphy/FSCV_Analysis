function [xcovdata]=xvardata(data,rates,events,plotparam,hax,seltrials,varargin)
%PRE 1/28/2019
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
if isempty(seltrials)
    seltrials=1:numtrials;
end
%not supplied use all data, if same length
%alignidx{1}=repmat(1,1,numtrials);
%alignidx{2}=repmat(size(data{2},2),1,numtrials);    %default data{2} fscv domain
%if ~isempty(events)
   % alignidx{1}=events(1,:);        %in nlx sample rate 1000 hz
   % alignidx{2}=events(2,:);   
    alignidx=events;        %cue single align for all trials diff
%end


plotwin=0;
pade=5;     %
argnum=1;
iscsc=0;        %not csc data, ie fscv time stamps
cscname=[];
pulse=0;
lfp2=0;
logscale=0;
eye=0;
blink=0;
noplot=0;
sitename={};
freqband=[0 0];     %freq band of filtered signal
plotxcor=0;
xcovtype='coeff';       %default xcov type argument
eventnames={};
shuffle=0;
excludethres=.25;       %default percentage of data nan less than this ok
eventsampsflag=0;
xcovdata={};
splitwinflag=0;     %split window in half to look at cov before/after peak separately
firstplot=0;
win=[-1 1];     %+/- 1 s from alignidx provided 1/10/19
%coeff - normalizes the sequence so that the covariances at
%zero lag are identically 1.0.
%unbiased - scales the raw covariance by 1/(M-abs(k)), where k
%is the index into the result.    
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'firstplot'
            firstplot=1;
        case 'cor'
            plotxcor=1;
        case 'pade'
            argnum=argnum+1;
            pade=varargin{argnum};
        case 'cscale'
            argnum=argnum+1;
            plotparam.cscale=varargin{argnum};
        case 'plotnum'
            argnum=argnum+1;
            plotnum=varargin{argnum};
        case 'sitename'
            argnum=argnum+1;        %fscv ch site name
            sitename=varargin{argnum};
        case 'filtband'
            argnum=argnum+1;
            freqband=varargin{argnum};
        case 'filt'
            %get filt band (filter here 08/04/2018)
            argnum=argnum+1;        %fscv ch site name
            filtband=varargin{argnum};     %flag to plot lfp2    
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
        case 'plotwin'
            argnum=argnum+1;
            plotwin=varargin{argnum};   %how much time from xcov want to retain
        case 'log'
            logscale=1;
        case 'noplot'
            noplot=1;   %don't plot
        case 'type'
            argnum=argnum+1;
            xcovtype=varargin{argnum};           
        case 'eventnames'
            argnum=argnum+1;
            eventnames=varargin{argnum};
        case 'eventsnlxids'
            %user provided events variable in nlx samples domain ie. not fscv time
            %scale default so need to convert for both datas if uncommon
            %rates
            eventsampsflag=1;
        case 'shuffle'
            %shuffle trials of data 2, keep data 1 same
            shuffle=1;
            noplot=1;
        case 'excludethres'
            %percent of data nan, need to exluce if above this
            argnum=argnum+1;
            excludethres=varargin{argnum};
        case 'splitwin'
            splitwinflag=1;     %split window in half to look at cov before/after peak separately
            
    end
    argnum = argnum + 1;
end

if eventsampsflag==1
    %convert nlx samples domain to fscv time scales for common domain
    alignidx=round(alignidx./rates(1).*rates(2));
end

%resample data to higher sampling rate so equal sample rates
numsamples(1)=size(data{1},2);
numsamples(2)=size(data{2},2);
resampled={};
ts=[];
resampled{2}=data{2};
resampled{1}=data{1};
ts=[0:1/rates(2):numsamples(2)/rates(2)-1/rates(2)];
if rates(1)>rates(2)
    %downsample data 1 (lfp) to da (data 2 10 hz)
    resampled{1}=[];
    ts1=[0:1/rates(1):numsamples(1)/rates(1)-1/rates(1)];
    ts=[0:1/rates(2):numsamples(2)/rates(2)-1/rates(2)];
    for itrial=1:numtrials
        downsamp = downsampleSingleBatch(data{1}(itrial,:),round(rates(1)/rates(2)));
        downsamp=[downsamp 0];
        resampled{1}(itrial,:)=downsamp;        
    end
end
xcovdata.rate=rates(2);
%make sure all same length windows across trials;
%otherwise cut off end for minimum
triallength=[];
trialwinids=[];
for itrial=1:numtrials
 %   trialwinids=events(1,tid):events(2,tid);
    trialwinids(itrial,:)=alignidx(itrial)+win(1)*rates(2):alignidx(itrial)+win(2)*rates(2);
end
badtrials=[];
for itrial=1:numtrials
    %smooth out any nan or imag data if less than exclude thres
    dadata=resampled{2}(itrial,:);
    lfpdata=resampled{1}(itrial,:);
    curwin=trialwinids(itrial,:);
   % trialwinids=alignidx{1}(1,itrial):alignidx{2}(1,itrial);
    %trialwinids=alignidx{1}(1,itrial):alignidx{1}(1,itrial)+minnums-1;
    if any(curwin<1)  || isempty(curwin) || ...
        any(curwin>length(resampled{2}(itrial,:))) || ...
        any(curwin>length(resampled{1}(itrial,:)))
       dadata=nan(1,length(dadata));
       resampled{2}(itrial,:)=dadata; 
       resampled{1}(itrial,:)=nan(1,length(lfpdata));
       %REMOVE from seltrials
               %REMOVE FROM SELTRIALS
        if any(ismember(seltrials,itrial))
            seltrials(find(seltrials==itrial))=[];
            badtrials=[badtrials itrial];
        end
        disp([num2str(itrial) ' trial out of bounds']);
       continue;
    else
       if any(imag(dadata(curwin)))
           %if any imaginary component means data had been cut off from
           %start or end of recording, cannot smooth, skip trial, fill all
           %with nan
           dadata=nan(1,length(dadata));
       end
       if any(isnan(dadata(curwin)))
            if length(find(isnan(dadata(curwin))==1))...
                    /length(dadata(curwin))<=excludethres
       %if > 25% of trial data is nan, make all nan so not part of
      %cumulative calculation, exclusion criteria
                [dadata,TF] = fillmissing(dadata,'pchip');
                %itrial
                
            else
                dadata=nan(1,length(dadata));
            end            
       end
    end
    if any(abs(dadata(curwin))>12*nanstd(resampled{2}(itrial,curwin)))
        %sometimes filling data with pchip creates artificial
        %slopes and boundaries, should remove.
        dadata=nan(1,length(dadata));
    end
    if any(abs(dadata(curwin))>abs(10*nanmean(resampled{2}(itrial,curwin))))
        %sometimes filling data with pchip creates artificial
        %slopes and boundaries, should remove.
        dadata=nan(1,length(dadata));
    end
    if std(dadata(curwin))>10
        %sometimes filling data with pchip creates artificial
        %slopes and boundaries, should remove.
        dadata=nan(1,length(dadata));
    end
    if sum(lfpdata(curwin))==0
        resampled{1}(itrial,:)=nan(1,length(lfpdata));
    end
   resampled{2}(itrial,:)=dadata;          
end
xcovdata.lfpdata=resampled{1};
xcovdata.dadata=resampled{2};
if ~isempty(badtrials)
for ib=1:length(badtrials)
    if ismember(badtrials(ib),seltrials)
        disp(['removing trial # ' num2str(badtrials(ib))]);
        seltrials(find(seltrials==badtrials(ib)))=[];
    end
end
end

xcovdata.seltrials=seltrials;

%plotting test
%{
%figure; plot(ts, zscore(resampled{1}(1,:)))
%yyaxis right;
%plot(ts, zscore(resampled{2}(1,:)))
figure; plot(ts(trialwinids),resampled{1}(1,trialwinids))
hold on; yyaxis right; plot(ts(trialwinids),resampled{2}(1,trialwinids))
%}
xcovda=[];
xcovpre=[];
xcovpost=[];
xcovlag=[];
xcovlagpre=[];
xcovlagpost=[];
[a,b]=sort(rand(1,length(seltrials)));
[a,b1]=sort(rand(1,length(seltrials)));

shufseltrials=seltrials(b);     %shuffled selectedtrials
shufseltrials1=seltrials(b1);     %shuffled selectedtrials
xcovdata.winids=[];


ipre=1;
%ONLY PROCESS SELTRIALS AND STORE ONLY SELTRIALS
for itrial=1:length(seltrials)
    tid=seltrials(itrial);
    curwin=trialwinids(tid,:);
    %PROBLEM NOTED 1/10/19
    %SELECTED alignidx{1}(1,itrial) rather than tid    
    xcovdata.winids(itrial,:)=curwin;
    %figure; plot(resampled{2}(itrial,trialwinids))
    %hold on; yyaxis right; plot(resampled{1}(itrial,trialwinids))
    if any(curwin<1) || any(curwin>length(resampled{1}(tid,:)))
        data1=nan(1,length(curwin));
        data2=nan(1,length(curwin));
        
        warning([num2str(tid) ' trial out of bounds AGAN']);
    else
    data1=resampled{1}(tid,curwin);
    data2=resampled{2}(tid,curwin);
    end
    if shuffle==1
        shufwin=trialwinids(shufseltrials(itrial),:);
        if any(shufwin<1)
        data2=nan(1,length(shufwin));
        else
        data2=resampled{2}(shufseltrials(itrial),shufwin);
        end
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
   %xcovlag(itrial,:)=lagt;
end
tsplot=ts(xcovdata.winids(1,:));
tsplot=tsplot-tsplot(1);
tsx=[-fliplr(tsplot(2:end)) tsplot];
xcovdata.seltrials=seltrials;

xcovdata.xcovda=xcovda;
xcovdata.tsx=tsx;
xcovdata.xcovpre=xcovpre;
xcovdata.xcovpost=xcovpost;

tsplot=ts(xcovdata.winids(1,1:ipre));
tsplot=tsplot-tsplot(1);
tshalf=[-fliplr(tsplot(2:end)) tsplot];
xcovdata.tshalf=tshalf;



%get significant peak lag time points for each trial and their value at
%peak 
xcovlag=[];
xcovlaganti=[];     %neg peaks
pad=round(.05*size(xcovda,2));
actids=pad+1:size(xcovda,2)-pad;    %ids to search before padding
for itrial=1:length(seltrials)
    %tid=seltrials(itrial);
    if any(isnan(xcovda(itrial,:)))
        %if signal nanned out because too many nans , store nan
        xcovlag(itrial,1)=nan;            %store lag at max in seconds
        xcovlag(itrial,2)=nan;
        xcovlaganti(itrial,1)=nan;
        xcovlaganti(itrial,2)=nan;
    else
        %maxx=max(abs(xcovda(itrial,:)));        %max absolute of signal
        %maxt=find(abs(xcovda(itrial,:))==maxx);      %lag of max
        %11/22/2018, find pos/neg peaks within each tiral
        maxx=max(xcovda(itrial,actids));     %max pos signal 11/22/2018
        maxt=find(xcovda(itrial,:)==maxx);      %lag of max
        minn=min(xcovda(itrial,actids));
        mint=find(xcovda(itrial,:)==minn);      %lag of max
        if length(maxt)>1 || maxt<=pad || maxt>size(xcovda,2)-pad
            xcovlag(itrial,1)=nan;
            xcovlag(itrial,2)=nan;
            %if more than 1, get closest to 0
            if length(maxt)>1
                cc=find(abs(maxt-median(1:size(xcovda,2)))==min(abs(maxt-median(1:size(xcovda,2)))));
                xcovlag(itrial,1)=tsx(maxt(cc));
                xcovlag(itrial,2)=xcovda(itrial,maxt(cc));
            end
        else
            xcovlag(itrial,1)=tsx(maxt);            %store lag at max in seconds
            xcovlag(itrial,2)=xcovda(itrial,maxt);  %store value at max in xcov
        end
        if length(mint)>1 || mint<=pad || mint>size(xcovda,2)-pad
            xcovlaganti(itrial,1)=nan;
            xcovlaganti(itrial,2)=nan;
            %if more than 1, get closest to 0
            if length(mint)>1
                cc=find(abs(mint-median(1:size(xcovda,2)))==min(abs(mint-median(1:size(xcovda,2)))));
                xcovlaganti(itrial,1)=tsx(mint(cc));
                xcovlaganti(itrial,2)=xcovda(itrial,mint(cc));
            end
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
%{
%bin lags with negative or positive peaks
thres=nanstd(xcovlag(:,2),[])*.5;       %hreshold is 50% std
posids=find(xcovlag(:,2)>thres);
tbins=min(tsx):.25:max(tsx);
[npos, binspos]= histc(xcovlag(posids,1),tbins);
bins.countspos=npos;
bins.binidspos=binspos;
negids=find(xcovlag(:,2)<-thres);
[nneg, binsneg]= histc(xcovlag(negids,1),tbins);
bins.countsneg=nneg;
bins.binidsneg=binsneg;
bins.tbins=tbins;
%}

xcorda=[];
if plotxcor==1
    for itrial=1:length(seltrials)
        tid=seltrials(itrial);
        trialwinids=alignidx{1}(1,tid):alignidx{2}(1,tid);
        xcortemp=xcorr(resampled{1}(tid,trialwinids),resampled{2}(tid,trialwinids),xcovtype);
       xcorda(itrial,:)=xcortemp;
    end
    xcovdata.xcorda=xcorda;

end
%test plot
%{
tsplot=ts(trialwinids);
tsplot=tsplot-tsplot(1);
tsx=[-fliplr(tsplot(2:end)) tsplot];
plot(tsx,xcovda(3,:))
%}
%figure; plot(xcovda(1:9,:)')
%legend(nlx.cscNames{1:9})        
%xts(dacount)=winfirstts; 
%xcovdata.bins=bins;

if shuffle==1
    xcovdata.shuftrials=shufseltrials;
end
xcovdata.freqband=freqband;
xcovdata.eventnames=eventnames;
xcovdata.sitename=sitename;
if noplot==0    
    if firstplot
        plotxtrials(hax,xcovdata,'leftlabel');
    else
    plotxtrials(hax,xcovdata);
    end
end

end




