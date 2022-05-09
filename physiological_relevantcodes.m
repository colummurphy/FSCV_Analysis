%note that all physiological data measured at neurlaynx smapling frequency of 30 kHz, and downsampled here or elsewhere after to 10 Hz when used for FSCV analysis 
%Please make sure sampling windows are appropriate for your data

%Licking is just the magnitude of the x, y, and z inputs from
%accelerometer. I did not put that here since it's straightforward

%BEat to beat interval standard deviation (RRSTD)
 hrwin=hr(:,winslfp);       %get heart rate over a window (e.g 2 s? it really depends on analysis) for all trials (heart rate is just hte inverse of the time between pulses)
hrtemp=nanmean(hrwin,2);    %average heart rate in this window for all trials
rrtemp=1./hr.*60.*1e3;  %convert to beat to beat interval in s
rrdata=rrtemp(:,winslfp);   %get this over the window of interest for all trials

%Compute pulse timestamps, note that i had to transform time variable sampel period to
%match with FSCV, you do not need to do downsampling operatiosn done here
    sitename='pulse';
    %calculate HR
    pulsed=zscore(da,0,2);
    stepsize=.1;       %0.5 second step size
    newrate=1/stepsize;
    thres=mean(pulsed); thres=mean(thres);
    stepwin=round(samplespersec*stepsize);
    findmaxwin=round(samplespersec*.1);    
    bins=[];
    histbins=[];
    rates=[];
    pulsets=[];
    tslfp=data.relts;           %original relative ts's for csc
    rrint=[];
    hr=[];
    tsdown=tslfp(1:stepsize*samplespersec:end);
    hrdown=nan(size(pulsed,1),length(tsdown));
    rrdown=nan(size(pulsed,1),length(tsdown));
    for itrial=1:size(pulsed,1)
        [pks,locs,w,p] = findpeaks(pulsed(itrial,:),'minpeakdistance',findmaxwin,'MinPeakHeight',1);
        if ~isempty(locs)
        %get rr intervals
        intervals=diff(locs);
        interv=intervals./samplespersec;        %s between peaks
        hrtemp=[];
        rrinttemp=[];
        for iint=1:length(intervals)
            tsint=[];
            if iint==length(intervals)
                tsint=find(round(tslfp*samplespersec)>=locs(iint));
            else
                tsint=find(round(tslfp*samplespersec)<locs(iint+1) & round(tslfp*samplespersec)>=locs(iint));
            end
            rrinttemp(tsint)=interv(iint);
            hrtemp(tsint)=1/interv(iint)*60;
            if iint==1
                rrinttemp(1:tsint(1))=interv(iint);
                hrtemp(1:tsint(1))=1/interv(iint)*60;
            end
        end
        rrdown(itrial,:) = downsampleSingleBatch(rrinttemp,round(samplespersec/newrate));
        hrdown(itrial,:) = downsampleSingleBatch(hrtemp,round(samplespersec/newrate));
        hrdown(itrial,1)=nan;
        rrdown(itrial,1)=nan;
        end
    end


%deligtching of pupil diameter and eye position data if needed when
%calculating evoked pupillary response, without the blinks
      for itrial=1:size(da,1)
       eyedata(itrial,:)=deglitchnanamp(da(itrial,:),2.5e-3,30);
      end

