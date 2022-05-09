function detected=detectdatransientschunks(alldata,parameters,glitchids,tParam,varargin)
%chunk long data file into periods of overlapping windows and process
%checked on 05/21/2018 that bg subtraction has no effect on computed 
%da concentration changes but need to iteratively subtract because this
%affects magnitude of correlation coefficients, ie detected correlated DA
%use unsubtracted fscv data to determine whether transients are correlated
%to dopamine template based on its cv
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

%check matlab version if 2013 or earlier, does not support std with omitnan
%argument, need to use nanstd instead
matlabver=version('-release');
nums=regexp(matlabver,'\d');
matlabver=str2num(matlabver(nums));
old=0;
if matlabver<2014
    old=1;
end

samplerate=parameters.samplerate;
argnum=1;
relts=0;
winbaseline=10;
while argnum <= length(varargin)
    switch varargin{argnum}       
        case 'timestamp'
            %output detected timestamps relative to indicated absolute time
            %stamp here (ie. neuralynx saved TS)
            argnum=argnum+1;
            relts=varargin{argnum};
        otherwise
            error('badoption');
    end
    argnum = argnum + 1;
end


%Get info on loaded FSCV data
%initialize holding variable detected for detected parameters
percentout=0.3;     %if >60% of window is bad then skip to next bg win sub
detected={}; 
detected.maxTS=[];
detected.prevalleyTS=[];
detected.postvalleyTS=[];
detected.datrace=[];
detected.tstrace=[];
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
deflectionFlag=0;       %movement occured, wait until stable
count=0;
overlap=10;      %within overlap # of samples of previously stored, replace or skip


%generate chunk intervals
sizealldata=size(alldata);
chunksize=parameters.chunksize*samplerate;         %samples per chunk, default = 30s
chunkoverlap=round(.2*chunksize);         %5% overlap

for chunk=1:chunksize-1:sizealldata(2)-10
    chunkstartid=chunk-chunkoverlap;
    chunkendid=chunkstartid+chunksize-1+2*chunkoverlap;
    idschunk=chunkstartid:chunkendid;
    idschunkinwin=idschunk(idschunk>0 & idschunk<sizealldata(2));
    data=alldata(:,idschunkinwin);
    sizedata=size(data); 
    relativets=relts+(idschunkinwin(1)-1)/samplerate;

tstart=5; tend=sizedata(2)-10;  %default start/end times
windur=abs(tstart-tend);    %time period for finding correlating signals from BG point

for tpoint=tstart:10:tend 
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
    rr=corr(subdata,KDA);     %detect positive corr's to DA, not abs
    rm=corr(subdata,KM);     %detect positive corr's to m, not abs
    rmout=find(abs(rm)>parameters.RThresOut);
    rbg=corr(subdata,KBG);     %detect positive corr's to m, not abs
    rbout=find(abs(rbg)>parameters.RThresOut);
    rph=corr(subdata,KPH);     %detect positive corr's to m, not abs
    rpout=find(abs(rph)>parameters.RThresOut);
    rda=find(rr>Rthres);        %find idx's of correlated DA in window
    %if most of window is ph/bg/movement skip & try next sub point
    if length(rpout)/size(data,2)>percentout || ...
            length(rbout)/size(data,2)>percentout || ...
                 length(rmout)/size(data,2)>percentout
             continue
    end

    %only those local (within windur) to bg sub point
    rda=rda(rda>tpoint-windur & rda<tpoint+windur);
    if isempty(rda)
        %no correlated signals this sub idx, go to next indx
        continue
    end
    %get pct dopamine signal with removed M/Q thres & glitches
    glitchidsrel=glitchids(glitchids>=idschunkinwin(1) & glitchids<=idschunkinwin(end));
    glitchidsrel=glitchidsrel-(idschunkinwin(1)-1);
    ipcr = getpctdirect(subdata,'ats',ats,'cts',cts,...
        'qthres',parameters.QaDAPH,'rthresm',parameters.RThresOut,...
        'rthresb',parameters.RThresOut,'rthresp',parameters.RThresOut,...
        'imaxbg',max(abs(bgmean)),'glitchids',glitchidsrel);
    daproj=ipcr.DAiso;
    dastd=[];
    if old==0
        dastd=std(daproj,0,2,'omitnan');
    else
        dastd=nanstd(daproj,0,2);
    end
    stdwin=tpoint-winbaseline:tpoint+winbaseline;
    stdwin=stdwin(stdwin>0);
    if old==0
        dastd=std(daproj(stdwin),0,2,'omitnan');
    else
        dastd=nanstd(daproj(stdwin),0,2);
    end
    mindiff=peakthreshold*dastd;
    %if std too high, use default peak threshold instead
    if mindiff>parameters.peakDiff
        mindiff=parameters.peakDiff;
    end
    %for each detected correlated signal, determine whether it is a local
    %maximum, ie transient
    for idxcorr=1:length(rda)
        if idxcorr>1
        if rda(idxcorr)-1==rda(idxcorr-1) && rr(rda(idxcorr))<rr(rda(idxcorr-1))
            %if analyzed neighboring cv in last loop & previous cv more
            %correlated to DA, then skip this loop
            continue
        end
        end
        if isnan(daproj(rda(idxcorr)))
            %if nulled out portion of signal, ie glitch/movement/Qthres
            %skip
            continue
        end
        winmean=rda(idxcorr)-cvavg:rda(idxcorr)+cvavg;
        winmean=winmean(winmean>0 & winmean<sizedata(2));
        cvmean=mean(subdata(:,winmean),2);
        if corr(cvmean,KDA)<Rthres
            %if mean CV not correlated,use single sample cvmean
            %at known correlation idx
            cvmean=subdata(:,rda(idxcorr));
        end
        %determine local peak around correlated signal, should be nearby
        wins=(rda(idxcorr)-parameters.peakWindow-1:rda(idxcorr)+parameters.peakWindow);  
        win=wins(wins>0 & wins<sizedata(2));    %confined to recording period
        maxts=find(daproj==max(daproj(win)));  % max ts near correlated
        maxda=daproj(maxts);
        if isempty(maxda)
            continue
        end
        %get window for defining whether local peak exists
        maxwins=(rda(idxcorr)-parameters.maxDurationPeak-1:rda(idxcorr)+parameters.maxDurationPeak);  
        maxwin=maxwins(maxwins>0 & maxwins<sizedata(2));    %confined to recording period
        %get valley ts's around max based on threshold < max - std (absolute ts's)
        tssvalleys=find(daproj(maxwin)<maxda-mindiff)+maxwin(1)-1;
        if isempty(tssvalleys)
            %signal does not increase/decrease to threshold level within 
            %defined maxDurationPeak windows, skip
            continue
        end        
        %rising interval
        prepeaktss=tssvalleys(tssvalleys<maxts);    %only those before maxts
        if isempty(prepeaktss)
            %signal does not increase to max within defined interval, skip
            continue
        end
        prepeakts=prepeaktss(end);      %closest to max used to define width  
        stdprepeakwin=prepeakts-winbaseline:prepeakts;    %win for std basline
        stdprepeakwin=stdprepeakwin(stdprepeakwin>0);
        %calculate peak
        %get half life looking at larger window
        broadwins=maxts-searchwinbroad-1:maxts+searchwinbroad;  
        broadwin=broadwins(broadwins>0 & broadwins<sizedata(2));
        %determine relative amplitude of DA peak?
        %get minima da's in window
        dasorted=sort(daproj(broadwin),'ascend');
        %define baseline DA as mean 10% of smallest signals in window
        baselineda=0;
        if old==0
            baselineda=mean(dasorted(1:.2*searchwinbroad),'omitnan');
        else
            baselineda=nanmean(dasorted(1:.2*searchwinbroad));
        end
        %previous baseline biased & could inflate maxima, since decreases
        %in window could be possible
      %  baselineda=mean(dasorted,'omitnan');    %take average of entire signal
        %baseline from average of entire signal in period can be misleading
        %since many large fluctuations that are state dependent
      %  baselineids=maxts-searchwinbroad-1:maxts-1;
        %%baselineids=baselineids(baselineids>0);
       % baselineda=mean(dasorted(baselineids),'omitnan');
        %original baseline from minima in window is probably best to get
        %local da maxima        
        peakda=maxda-baselineda;

        %check that peak is above stds of baseline in pre valley
      %  peakthreshold=1;
      davalleythres=[];
        if old==0
            davalleythres=std(daproj(stdprepeakwin),0,2,'omitnan')*peakthreshold;
        else
            davalleythres=nanstd(daproj(stdprepeakwin),0,2)*peakthreshold;
        end
        if davalleythres>mindiff
            davalleythres=mindiff;
        end
        if maxda<davalleythres+daproj(prepeakts)
            continue
        end
            
%%%Maybe take out since we already check for glitches
        %if maxts-prepeakts<parameters.minPeakWidth
        %    continue
        %end     
        %get semi-valley ts's post max based on threshold < max - std*.5 (absolute ts's)
        tssvalleys2=find(daproj(maxwin)<maxda-mindiff*.5)+maxwin(1)-1;
        %falling interval
        postpeaktss=tssvalleys2(tssvalleys2>maxts);    %only those before maxts
        if isempty(postpeaktss)
            %does not decrease beyond threshold within defined interval, skip
            continue
        end
        postpeakts=postpeaktss(1);      %closest to max used to define width

        
        
        %shift daproj by baseline conc & find TS's below half peak conc
        tsshalfs=find(daproj(broadwin)-baselineda<=peakda/2)+broadwin(1)-1;
        %initialize half lifes as nan assuming not defined within window
        halfrisets=nan;     
        halffallts=nan;
        %if empty, could mean interference in window preventing detection
        %of start/end half times, so do not disclude as peak (ie no skip)
        if ~isempty(tsshalfs)
            %signal increases/decreases to its half within broad win
            %rising interval
            halfrise=tsshalfs(tsshalfs<maxts);    %only those before maxts
            if ~isempty(halfrise)
                halfrisets=maxts-halfrise(end);     %half rise time in samples
            end
            %falling interval
            halffall=tsshalfs(tsshalfs>maxts);    %only those before maxts
            if ~isempty(halffall)
                halffallts=halffall(1)-maxts;    %half fall time in samples   
            end
        end         

        %get da correlation info on max signal 
        winmax=maxts-cvavg:maxts+cvavg;
        winmax=winmax(winmax>0 & winmax<sizedata(2));
        cvmax=mean(subdata(:,winmax),2);    %cv at maxts
        rmax=corr(cvmax,KDA);  %da correlation at maxts
        rmaxm=corr(cvmax,KM);       %correlation to movement
        rmaxph=corr(cvmax,KPH);
        rmaxbg=corr(cvmax,KBG);
        %get signal at max if more correlated than curr idx & not
        %correlated to movement already be checked during 
        %getpctdirect run, but double check here
        if abs(rmaxph)>parameters.RThresOut || abs(rmaxbg)>parameters.RThresOut
            continue
        end
        if rmax<rr(rda(idxcorr)) || rmaxm>parameters.RThresOut
            %otherwise make CV to store, cvmean, originally correlated
            %non-max CV
            cvmax=cvmean;
            %R to store also default rda
            rmax=rr(rda(idxcorr));
         %   maxts=rda(idxcorr);
        end
        %check if this max ts already stored 
        %or within 1 s of already stored
        overlap=parameters.maxDurationPeak; 
        %set overlap to equal max duration of da otherwise overlapping
        %signals & false positives easily detected.
        if any(ismember(maxts-overlap:maxts+overlap,...
                round((detected.maxTS-relativets).*samplerate)+1))
            %if already stored, check if less correlated than current
            storedidx=find(ismember(round((detected.maxTS-relativets).*samplerate)+1,...
                maxts-overlap:maxts+overlap)==1);
            if length(storedidx)>1
                %if more than one overlap, get one closest in time
                differencet=abs(round((detected.maxTS(storedidx)-relativets)...
                    .*samplerate)+1-maxts);
                idclosestt=find(differencet==min(differencet));
                storedidx=storedidx(idclosestt(end));
            end
            if detected.R(storedidx)<rmax
                %if currently detected signal more highly correlated
                %than previously detected signal at same TS
                %Store and Replace
                %absolute timestamp
                detected.prevalleyTS(storedidx)=(prepeakts-1)./samplerate+relativets;
                detected.postvalleyTS(storedidx)=(postpeakts-1)./samplerate+relativets;
                detected.maxTS(storedidx)=(maxts-1)./samplerate+relativets;
                detected.stdDA(storedidx)=dastd;       %da std used for threshold
                detected.risetimestds(storedidx)=(maxts-prepeakts)./samplerate; 
                detected.falltimestds(storedidx)=(postpeakts-maxts)./samplerate;
                detected.maxDA(storedidx)=peakda;       %da max (after baseline sub)
                detected.risehalftime(storedidx)=halfrisets./samplerate; 
                detected.fallhalftime(storedidx)=halffallts./samplerate; 
                detected.R(storedidx)=rmax;
                detected.cv(:,storedidx)=cvmax;    
                %save da trace in defined tracewin period normalized to
                %pre-valley da of which defined peak da trheshold
                traceids=maxts+tracewin(1)*samplerate:maxts+tracewin(2)*samplerate;
                traceidsneg=traceids(traceids<=0);
                traceidsoverend=traceids(traceids>sizedata(2));
                traceidsinwin=traceids(traceids>0 & traceids<=sizedata(2));
                %fill negative time values with imaginary value to
                %distinguish from nan artifact values
                dafillnegids=repmat(1j,length(traceidsneg),1);
                dafilloverids=repmat(1j,length(traceidsoverend),1);
                %use isreal() function to get only real values later
                detected.datrace(storedidx,:) = [dafillnegids; ...
                    daproj(traceidsinwin)'-baselineda;...
                    dafilloverids];
                detected.tstrace(storedidx,:) = (traceids-1)./samplerate+relativets;
            end
        else   
            %check if valley ts's overlap with any previously stored valleys 
            %plateua-like dopamine signals could potentially share 2 dif
            %peaks, want to remove these types of repeats
            if any(ismember(postpeakts-overlap:postpeakts+overlap,...
                    round((detected.postvalleyTS-relativets).*samplerate)+1))
                continue
            end
            if any(ismember(prepeakts-overlap:prepeakts+overlap,...
                    round((detected.prevalleyTS-relativets).*samplerate)+1))
                continue
            end            
            
            %store as new detected signal
            count=count+1;      %increment counter of detected signals
            %absolute timestamp
            detected.prevalleyTS(count)=(prepeakts-1)./samplerate+relativets;
            detected.postvalleyTS(count)=(postpeakts-1)./samplerate+relativets;
            detected.maxTS(count)=(maxts-1)./samplerate+relativets;       
            detected.stdDA(count)=dastd;       %da std used for threshold
            detected.risetimestds(count)=(maxts-prepeakts)./samplerate; %relative ts ins amples
            detected.falltimestds(count)=(postpeakts-maxts)./samplerate;
            detected.maxDA(count)=peakda;       %da max (after baseline sub)
            detected.risehalftime(count)=halfrisets./samplerate; 
            detected.fallhalftime(count)=halffallts./samplerate; 
            detected.R(count)=rmax;
            detected.cv(:,count)=cvmax;   
            %save da trace in defined tracewin period normalized to
            %pre-valley da of which defined peak da trheshold
            traceids=maxts+tracewin(1)*samplerate:maxts+tracewin(2)*samplerate;
            traceidsneg=traceids(traceids<=0);
            traceidsoverend=traceids(traceids>sizedata(2));
            traceidsinwin=traceids(traceids>0 & traceids<=sizedata(2));
            %fill negative time values with imaginary value to
            %distinguish from nan artifact values
            dafillnegids=repmat(1j,length(traceidsneg),1);
            dafilloverids=repmat(1j,length(traceidsoverend),1);
            %use isreal() function to get only real values later
            detected.datrace(count,:) = [dafillnegids; ...
                daproj(traceidsinwin)'-baselineda;...
                dafilloverids];
            detected.tstrace(count,:) = (traceids-1)./samplerate+relativets;
        end
    %finish looping correlated signals for given bg sub point    
    end
%finish bg sub point
end
end

end