function detected=detectdauncorr(alldata,parameters,glitchids,varargin)
%prior to 2/9/2019, after this made so also correlated since now 3 pcs 
%11/06/2018 since peak can plateau and have variable peaks within it, thus
%prolonging the peak and making the center of peak less definable
%use instead max in differential pre peak (ie max slope) as peak, and make
%sure there is a negative differential subsequent as defined by width of
%peak. Store both max slop and max ts (abs peak)
%make sure stable baseline 2 s long (ie not rebound da signal where local valley)
%10/31/2018 need to make so it detects first maximum not after another hump
%on top of hump..
%similar to detectdatransientschunks.m
%detect dopamine peaks in data assuming that glitches removed & da
%projections from pca are clean & bad periods are nann'ed out
%do not care about correlation of CV at specific ts to DA
%weights used on each K CV template to find optimal sub point for bg sub
%parameters.maxDurationPeak defines +/-maximum length of "transient"
%ie what defines a transient
%parameters used from parameters variable:
%BGWin (scan offset from tParam.timeWin) only if tParam not empty
%peakthreshold = # of std's (calculated from entire given signal)
%to use for amplitude threshold
%RThres, correlation coefficient for da detection
%RThresOut, correlation coefficient threshold for movement
%detection/removal, fed into getpctdirect.m
%QaDAPH, Q threshold for removing signal w/ too much extraneous signals
%CV - Ats Matrix with all CV templates for DA/pH/BG/movement
%Cts - Concentration matrix (Cts) for scaling component signals
%meanBGThres - threshold for std of background current whether too unstable
%deflectionThres - higher threshold if large movement, allow time for
%dischargin current to stabilize
%K - k matrices for correlating detected CVs to DA/movement
%peakWindow - window +/- around detected da correlated signal to search for
%max
%maxDurationPeak = window +/- around detected peak to search for valleys
numpcs=3;       %da,ph,bg
if isfield(parameters,'numpcs')
    %how many components to use, just ph & da, or inlcude bg/move?
   numpcs=parameters.numpcs; 
end

samplerate=parameters.samplerate;
argnum=1;
relts=0;
absts=[];
winbaseline=10;
%default sliding window sizes & overlap
winb=120;   
winb=60;        %in seconds changed 09/28/2018, 2 min for bg sub point
winbov=round(0.75*winb);%overlap win, higher better to remove false neg's
winp=30;        %window size for sliding window for local peak from pcda  <10/31/2018
winp=10;
winpstep=5;       %step size for sliding window for local peak from pcda
winpstep=1;       %09/28/2018step size for sliding window for local peak from pcda

percentout=0.3;     %if >60% of window is bad then skip to next bg win sub
maxdurpeak=parameters.maxDurationPeak;
maxreturnlength=parameters.maxDurationPeak; %how long it takes for peak to decrease
BGavgi=mean(alldata,2); 
ImaxBG=max(BGavgi);
detected.maxBG=ImaxBG;
peakthreshold=parameters.peakthreshold; %# of std's for threshold of peak
Rthres=parameters.RThres;
tracewin=parameters.tracepad;       %window +/- around max to save trace
KDA=parameters.K(:,1);
KPH=parameters.K(:,2);
KBG=parameters.K(:,3);
KM=parameters.K(:,4);
BGavg=2;
cvavg=1;
searchwinbroad=100;     %window to search for half life
ats=parameters.CV;
cts=parameters.Cts;
wda=1;
wm=-1;       %weights for how much these bad signals matter in detecting optimal bg sub
wb=-0.6;
wp=-0.3;
badzone=[];
%11/02/2018 since already find artifacts in find sigartifacts this is redundant and false negatives
rthresout=parameters.RThresOut+0.1;    

while argnum <= length(varargin)
    switch varargin{argnum}       
        case 'timestamp'
            %output detected timestamps relative to indicated absolute time
            %stamp here (ie. neuralynx saved TS)
            warning('use ts option instead');
            argnum=argnum+1;
            relts=varargin{argnum};
        case 'ts'
            %better than using just start nlx timestamp above because
            %fscv.events(:,4) contains ts that were synchronized around
            %triggers from nlx
            argnum=argnum+1;
            absts=varargin{argnum};          
        case 'windowbg'
            %window iwin size (sliding bg window) in seconds
            argnum=argnum+1;
            winbg=varargin{argnum};
            winbov=winbg(2);       %overlap size
            winb=winbg(1);         %win size
        case 'windowpeak'
            %window sliding for peak in seconds
            argnum=argnum+1;
            winpeak=varargin{argnum};
            winpstep=winpeak(2);
            winp=winpeak(1);
        case 'badwin'
            %percentage of signal in bg window if correlated to K's not of DA
            %skip this window and move to next
            argnum=argnum+1;
            percentout=varargin{argnum};
        case 'badzone'
            %signals in these time periods correspond to sig movement
            %as detected by findsigartifacts
            %ignore these periods (given in samples)
            argnum=argnum+1;
            badzone=varargin{argnum};
        otherwise
            error('badoption');
    end
    argnum = argnum + 1;
end


%Get info on loaded FSCV data
%initialize holding variable detected for detected parameters
detected={}; 
detected.maxTSpeak=[];
detected.maxTS=[];
detected.prevalleyTS=[];
detected.postvalleyTS=[];
detected.datrace=[];
detected.tstrace=[];
detected.tids=[];
deflectionFlag=0;       %movement occured, wait until stable
count=0;
overlap=20;      %within (2s) overlap # of samples of previously stored, replace or skip
badoverlap=[10 25];      %after these samples of movement artifact do not search
            %check more after moremvent artifact than before
if isfield(parameters,'badoverlap')
    badoverlap=parameters.badoverlap;
end
%generate iwin intervals
fullsize=size(alldata);
%check that absts provided is equal in length to samples
if ~isempty(absts)
    if length(absts)~=fullsize(2)
        error('ts provided is not same length as signals provided');
    end
end
winbs=winb*samplerate;              %windows in samples
winbovs=winbov*samplerate;
winps=winp*samplerate;
winpsteps=winpstep*samplerate;

dacount=0;
for iwin=1:winbs-1:fullsize(2)-10
    %find optimal bg sub id in bg window
    %create window boundaries first
    iwinstartid=iwin-winbovs;
    iwinendid=iwinstartid+winbs-1+2*winbovs;
    idsiwin=iwinstartid:iwinendid;
    idswinbg=idsiwin(idsiwin>0 & idsiwin<fullsize(2));
    data=alldata(:,idswinbg);
    sizedata=size(data); 
    relid=(idswinbg(1)-1);      %id offset for each data in window being analyzed
    relativets=relts+(idswinbg(1)-1)/samplerate;
    %tstart=5; tend=sizedata(2)-10;  %default start/end times
   % windur=abs(tstart-tend);    %time period for finding correlating signals from BG point

    %find optimal bg sub id by:
    %maximize correlated da signals, minimize signals correlated to m/bg/ph
    %in this order
    count=1;
    candidateids=[];
    for tpoint=5:10:length(idswinbg)-5
        if ismember(tpoint-BGavg:tpoint+BGavg,glitchids)
            %if within glitch period skip
            continue
        end        
        bg=data(:,tpoint-BGavg:tpoint+BGavg);
        bgstd=std(bg,0,2);  %standard deviation of BG current versus time
        if max(bgstd)>parameters.meanBGThres 
            %check if mean BG is unstable, then go to next BG
            continue
        end
        if max(bgstd)>parameters.deflectionThres
            %if large charging current, maybe movement
            %related, causing time to recover
            deflectionFlag=1;
            continue
        end
        if deflectionFlag==1
            %need to wait for current to stabilize, skip until
            %then & check when BG is stable to below the normal
            %threshold
            if max(BGRs)<parameters.stableLevel
                deflectionFlag=0;
            end
            continue
        end  
        bgmean=mean(bg,2);
        %tile BGavg vector for entire meas matrix span
        refmatrix=repmat(bgmean,1,size(data,2));   
        subdata=data-refmatrix;
        
        %find # correlated points to templates of da, m, bg,ph
        rr=corr(subdata,KDA);     %detect positive corr's to DA, not abs
        rda=find(rr>Rthres);        %idx's of correlated DA in window
        if isempty(rda)
            %no correlated signals this sub idx, go to next indx
            continue
        end
        rm=corr(subdata,KM);     %detect positive corr's to m, not abs
        rmout=find(abs(rm)>parameters.RThresOut);
        rbg=corr(subdata,KBG);     %detect positive corr's to m, not abs
        rbout=find(abs(rbg)>parameters.RThresOut);
        rph=corr(subdata,KPH);     %detect positive corr's to m, not abs
        rpout=find(abs(rph)>parameters.RThresOut);
        %if most of window is ph/bg/movement passing percent threshold,skip
        if length(rpout)/size(data,2)>percentout || ...
                length(rbout)/size(data,2)>percentout || ...
                     length(rmout)/size(data,2)>percentout
                 continue
        end
        
        candidateids(count,1)=tpoint;        %store sample idx
        candidateids(count,2)=length(rda);      %store # of corr da's
        candidateids(count,3)=length(rmout);    %store # of corr bad
        candidateids(count,4)=length(rbout);
        candidateids(count,5)=length(rpout);
    
        count=count+1;
        
    end
    
    if isempty(candidateids)
        %no good bg sub points found
        continue
    end
    %find optimal subidx by calculating magnitude of idx performance from
    %scan
    mag=candidateids;
    mag(:,2)=mag(:,2).*wda;         %multiply weights by each frequency
    mag(:,3)=mag(:,3).*wm;
    mag(:,4)=mag(:,4).*wb;
    mag(:,5)=mag(:,5).*wp;
    mag(:,6)=mag(:,2)+mag(:,3)+mag(:,4)+mag(:,5);
    optimalidx=find(mag(:,6)==max(mag(:,6)));    
    optidx=mag(optimalidx,1);    

    %sub data by optimal idx
    bg=data(:,optidx-BGavg:optidx+BGavg);
    bgmean=mean(bg,2);
    refmatrix=repmat(bgmean,1,size(data,2));   
    subdata=data-refmatrix;
    
    %compute pct da with removed M/Q thres & glitches
    glitchidsrel=glitchids(glitchids>=idswinbg(1) & glitchids<=idswinbg(end));
    glitchidsrel=glitchidsrel-idswinbg(1)+1;
    ipcr = getpctdirect(subdata,'ats',ats,'cts',cts,...
        'qthres',parameters.QaDAPH,'rthresm',rthresout,...
        'rthresb',rthresout,'rthresp',rthresout,...
        'imaxbg',max(abs(bgmean)),'glitchids',glitchidsrel,'numpcs',numpcs);
    daproj=ipcr.DAiso;
    dastd=nanstd(daproj,0,2);       %09/28/2018
    %for current bg window & projected da, slide window to find maxima
    for iiwin=1:winpsteps:length(idswinbg)
        %make windows for sliding for finding local maxima
        iiwinstart=iiwin;
    	iiwinendid=iiwin+winps-1;        
        idspwin=iiwinstart:iiwinendid;
        idsinpwin=idspwin(idspwin>0 & idspwin<size(data,2));
        daprojwin=daproj(idsinpwin);
        %make longer window for see when returns to baseline >10s
        idspwinlong=iiwinstart:iiwinendid+maxreturnlength;
        idspwinlong=idspwinlong(idspwinlong>0 & idspwinlong<size(data,2));
        if length(idspwinlong)<50
            %if length of signal < 5 s skip
            continue
        end
        daprojwinlong=daproj(idspwinlong);
        winids=idsinpwin;
        %wints=relts+(winids(1)-1)/samplerate;
        %danorm=zscore(daprojwin,0,'omitnan');       %normalize to window
        %danormfull=zscore(daproj,0,'omitnan');      %normalize full sig
        %danormwin=danormfull(idsinpwin);
        nonnanids=find(~isnan(daprojwin));
        if length(nonnanids)<=30
            continue
        end
        dasub=nanmean(daprojwin(nonnanids(1:30)));       %subtracted to beginning of window
        %skip if std of baseline window too high, above dastd? so to
        %prevent detection of transient rebounds
        %%%%%
        dastdbase=std(daprojwin(nonnanids(1:30)));
        dastdbase2=nanstd(daprojwin);
        dasubwin=daprojwin-dasub;
        dasubwinlong=daprojwinlong-dasub;
        %dastd=nanstd(daprojwin,0,2);           %make consistent through
        %entire bg rather than change each cycle to find peak
        %look at subtracted signal both using subtraction at start and at
        %end (for downwards trends) 09/28/2018
        %Nevermind - just decrease step size to account for small
       %inflections during downard
       % dasubend=nanmean(daprojwin(nonnanids(end-4:end))); 
        %dasubwinrev=daprojwin-dasubend;
        mindiff=peakthreshold*dastd;            %peak threshold based on stds
        if mindiff>parameters.peakDiff
            %if std thres too high, use default peak conc change threshold
            mindiff=parameters.peakDiff;
        end
        
        %get ids of threshold crossing
        crossesthres=find(dasubwin>mindiff);     
        %crossesthres2=find(dasubwinrev>mindiff);  
        %crossesthres=[crossesthres crossesthres2];
        if isempty(crossesthres)
            %no threshold crossing, slide window
            continue
        end
        %get first time when crosses & find when returns to value
        %returnsthres=find(dasubwin(crossesthres(1)+1:end)...
         %   <=dasubwin(crossesthres(1)))+crossesthres(1);
        %if isempty(returnsthres)
            %does not return to thres value, find max-mindiff point
            %instead, if got too high
        searchwin=crossesthres(1)+1:crossesthres(1)+1+maxdurpeak;
        searchwin=searchwin(searchwin>0 & searchwin<length(dasubwin));
        if isempty(searchwin)
            continue
        end
        %find max within sliding window 3 s (know that rise time usually
        %faster than 3 s - 5 s at least even for ramp relative to peakthres
        localmaxfirst=find(dasubwin==max(dasubwin(searchwin)));
        %find when returns back to baseline below peakthres
        %this can be much longer if wide peak and want to ignore small
        %fluctuations on peak..10/31/2018 change to dasubwinlong
        returnsthres=find(dasubwinlong(localmaxfirst+1:end)...
            <=dasubwin(localmaxfirst)-mindiff)+localmaxfirst;
        %end
        if isempty(returnsthres)
            %does not return to any minimal value, slide window
            continue
        end
        %find maximum within these 2 points
        %change to find max locally wthin subwin rather than to
        %returnsthres(1), find out which one shorter
        maxend=length(dasubwin);
        if returnsthres(1)<length(dasubwin)
            maxend=returnsthres(1);
        end
        localmaxfirst=find(dasubwin==max(dasubwin(crossesthres(1):maxend)));
        damax=dasubwin(localmaxfirst);
        %find declineing point in long window 10/31/2018
        nextvalley=find(dasubwinlong(localmaxfirst+1:end)<=damax-mindiff)+localmaxfirst;
        if isempty(nextvalley) || length(nextvalley)<3
            %does not cross threshold in other direction, slide window
            %or recession too short
            continue
        end
        %get max until next valley
        %10/31/2018 change to get first max rather than center of peak if
        %small fluctuations at peak
        localmaxid=find(dasubwinlong==max(dasubwinlong(crossesthres(1):nextvalley(1))));
        peakda=damax;
        maxts=winids(localmaxfirst);
        if localmaxid~=localmaxfirst
            %if not same as first detected peak, make sure difference between
            %first detected peak is less than mindiff otherwise use this
            %later peak
            if dasubwinlong(localmaxid)>damax+mindiff
               peakda=dasubwinlong(localmaxid);
               maxts=idspwinlong(localmaxid);
            else
                localmaxid=localmaxfirst;
            end
        end
        
        %if in bad zone skip               
        if any(ismember(maxts-badoverlap(2)+relid:maxts+badoverlap(1)+relid,...
                badzone))
            continue
        end
        
        %11/06/2018
        %find actual first peak in signal, ie when slope is maximal
        %first smooth signal with 3 s moving window
        %look for peak diff nearby preceding, thnen look for peak neg diff after
        %then look for max DA signal in between these 2 diff peaks
        smoothda=smooth(dasubwinlong,30);
        diffda=diff(smoothda);
        maxdiffpre=find(diffda==max(diffda(1:localmaxfirst)));
        if isempty(maxdiffpre)
            continue
        end
        if diffda(maxdiffpre)<=0
            %non positive peak, something worng...
            %why would it be defined as a maxima in the first place then?
            continue
        end
        %when peak declines should occur within 2 s after peak point
        peakdeclineids=find(diffda(maxdiffpre:end)<=0)+maxdiffpre-1;
        if isempty(peakdeclineids)
            %signal does not plateau or decrease
            continue
        end
        firstdecline=find(diff(peakdeclineids)>1);
        peakdeclineidx=length(diffda);      %default
        if ~isempty(firstdecline)
            firstdecline=firstdecline(1);
            peakdeclineidx=peakdeclineids(firstdecline);
        end  
        %search local win between 2 differential peaks
        localwin=maxdiffpre:peakdeclineidx;
        %get peak da signal between these 2 diff peaks
        realfirstpeak=find(dasubwinlong(localwin)==...
            max(dasubwinlong(localwin)))+localwin(1)-1;
        if isempty(realfirstpeak)
                %occurs only when values are nan in localwin
                continue
        end
        localmaxid=realfirstpeak(1);
        %use max positive slop as peak' 11/06/2018
        localmaxidslope=maxdiffpre+1;
        peakda=dasubwinlong(localmaxid);
        maxts=idspwinlong(localmaxid);
        maxslopets=idspwinlong(localmaxidslope);
        
        %not changed from before 11/06
        %if already stored sample id, skip
        if isempty(absts)
            %if stored
            if any(ismember(maxslopets-overlap:maxslopets+overlap,...
                    round((detected.maxTS-relativets).*samplerate)+1))
                continue
            end
            if any(ismember(maxts-overlap:maxts+overlap,...
                    round((detected.maxTSpeak-relativets).*samplerate)+1))
                continue
            end
        else
            checkids=maxslopets-overlap+relid:maxslopets+overlap+relid;
            checkids=checkids(checkids>0 & checkids<=length(absts));
            if any(ismember(absts(checkids),...
                    detected.maxTS))
                continue
            end
            checkids=maxts-overlap+relid:maxts+overlap+relid;
            checkids=checkids(checkids>0 & checkids<=length(absts));
            if any(ismember(absts(checkids),...
                    detected.maxTSpeak))
                continue
            end
        end
        valleys=find(dasubwinlong<=peakda-mindiff);
        prepeakvalley=valleys(valleys<localmaxid);
        if isempty(prepeakvalley)
            %no rise interval, skip
            continue
        end
        prepeakid=prepeakvalley(end);
        prepeakts=idspwinlong(prepeakid);
        postpeakvalley=valleys(valleys>localmaxid);
        if isempty(postpeakvalley)
            %never reaches threshold on return down, skip
            continue
        end
        postpeakid=postpeakvalley(1);
        postpeakts=idspwinlong(postpeakid);
        
        %get da magnitude based on some baseline level
        %get minima da's in window
       % dasorted=sort(dasubwin,'ascend');
        %define baseline DA as median 10% of smallest signals in window
        %baselineda=median(dasorted(1:round(.2*searchwinbroad)),'omitnan');
        baselineda=0;
        %this baseline appraoch could significantly disrot if there is
        %big da decrease within window & inaccurately amplify actual signal
        %maxda=peakda-baselineda;
        maxda=peakda;
               
        %get half fall/rise times, look at broader window outside dasubwin
        %NEED TO REDEFINE ACCORDING TO SLOPES 
        broaderwin=maxts-searchwinbroad:maxts+searchwinbroad;
        broaderwin=broaderwin(broaderwin>0 & broaderwin<size(data,2));
        tsshalfs=find(daproj(broaderwin)-dasub-baselineda<=maxda/2)+broaderwin(1)-1;
        %initialize half lifes as nan assuming not defined within window
        halfrise=nan;     
        halffall=nan;
        %if empty, could mean interference in window preventing detection
        %of start/end half times, so do not disclude as peak (ie no skip)
        if ~isempty(tsshalfs)
            %signal increases/decreases to its half within broad win
            %rising interval
            halfrises=tsshalfs(tsshalfs<maxts);    %only those before maxts
            if ~isempty(halfrises)
                halfrise=maxts-halfrises(end);     %half rise time in samples
            end
            %falling interval
            halffalls=tsshalfs(tsshalfs>maxts);    %only those before maxts
            if ~isempty(halffalls)
                halffall=halffalls(1)-maxts;    %half fall time in samples   
            end
        end   
        
        %save da trace in defined tracewin period normalized to
        %pre-valley da of which defined peak da trheshold
        traceids=maxslopets+tracewin(1)*samplerate:maxslopets+tracewin(2)*samplerate;
        traceidsneg=traceids(traceids<=0);
        traceidsoverend=traceids(traceids>=size(data,2));
        if ~isempty(traceidsoverend)
            %skip if window exceeds data length
            continue
        end
        if ~isempty(traceidsneg)
            continue
        end
        traceidsinwin=traceids(traceids>0 & traceids<=size(data,2));
        %store as new detected signal
        dacount=dacount+1;      %increment counter of detected signals
        detected.stdDA(dacount)=mindiff;       %value used for threshold
        detected.maxDA(dacount)=maxda;       %da max (after baseline sub)
        detected.risehalftime(dacount)=halfrise./samplerate; 
        detected.fallhalftime(dacount)=halffall./samplerate; 

        %fill negative time values with imaginary value to
        %distinguish from nan artifact values
        dafillnegids=repmat(1j,length(traceidsneg),1);
        dafilloverids=repmat(1j,length(traceidsoverend),1);
        %use isreal() function to get only real values later
        detected.datrace(dacount,:) = [dafillnegids; ...
            daproj(traceidsinwin)'-dasub-baselineda;...
            dafilloverids];
        detected.tids(dacount,:)=traceids+relid;       %store fscv sample ids
        
        %absolute timestamps stored
        if isempty(absts)
            %use provided 'timestamp' as reference ts, referenced each bg
            %sub point            
            detected.prevalleyTS(dacount)=(prepeakts-1)./samplerate+relativets;
            detected.postvalleyTS(dacount)=(postpeakts-1)./samplerate+relativets;
            detected.maxTS(dacount)=(maxslopets-1)./samplerate+relativets;    
            detected.maxTSpeak(dacount)=(maxts-1)./samplerate+relativets;           
            detected.tstrace(dacount,:) = (traceids-1)./samplerate+relativets;  %store absolute time stamp
        else
            %more accurate NLX ts's are stored in fscv.events(:,4) where they
            %were incorporated in relation to shared triggers on both systems.
            %use absolute "ts" timestamps for fscv signal to store ts's
            detected.prevalleyTS(dacount)=absts(prepeakts+relid);
            detected.postvalleyTS(dacount)=absts(postpeakts+relid);
            detected.maxTS(dacount)=absts(maxslopets+relid);        
            detected.maxTSpeak(dacount)=absts(maxts+relid);        
            detected.tstrace(dacount,:) =absts(traceids+relid);       
        end
        
       %finish loop for current sliding window
     end
%finish bg sub loop
end

end