function samplesTemp=detectDA(dataFSCV,dParam,pcr,tParam)
%detectDA

%Get info on loaded FSCV data
sizeData=size(dataFSCV); 
samplesperscan=sizeData(1);
%initialize temp storage variables
samples_DA=zeros(1,sizeData(2)); 
samples_Q=nan(1,sizeData(2)); samples_M=zeros(1,sizeData(2));
%initialize holding variable samplesTemp for detected parameters
samplesTemp={}; samplesTemp.cv=nan(samplesperscan,sizeData(2));
samplesTemp.maxTS=zeros(1,sizeData(2));     samplesTemp.minTS=zeros(1,sizeData(2));
samplesTemp.duration=zeros(1,sizeData(2));  samplesTemp.maxI=zeros(1,sizeData(2));
samplesTemp.maxDA=zeros(1,sizeData(2));     samplesTemp.R=zeros(1,sizeData(2));
samplesTemp.maxCV=nan(samplesperscan,sizeData(2));   samplesTemp.noiseI=zeros(1,sizeData(2));
samplesTemp.noiseDA=zeros(1,sizeData(2));
        
numDA=0;        %counter of dopamine signals detected
numPeak=0;      %counter of peaks detected
artifact=0;     %counter of artifacts
deflectionFlag=0;       %deflection flag toggle for large BG current deviations
BGavgi=mean(dataFSCV,2); 
ImaxBG=max(BGavgi);
samplesTemp.maxBG=ImaxBG;
flag=[0 0 0 0];
CVavg=2;                %# samples to average +/-
BGavg=5;                %# samples to average +/-

t_start=5; t_end=sizeData(2)-10;
%if isempty(offset); offset=0; end
%use supplied time window for BG cursor, with dParam.BGwin param offsets
if ~isempty(tParam.timeWin)
    t_start=tParam.timeWin(1)*tParam.samplingRate-dParam.BGWin(1);
    t_end=tParam.timeWin(2)*tParam.samplingRate+dParam.BGWin(2);
    if t_start<1
        t_start=5;
    end
    if t_end>sizeData(2)
        t_end=sizeData(2)-10;
    end
end

for t_point=t_start:10:t_end 
        RR=0;
        BGR=dataFSCV(:,t_point-2:t_point+2);
        
        %standard deviation of BG current versus time
        BGRs=std(BGR,0,2);
        
        if max(BGRs)>dParam.meanBGThres 
            %check if mean BG is unstable, then go to next BG
            continue
        end
        if max(BGRs)>dParam.deflectionThres
            %if large depolarization current, maybe movement
            %related, causing time to recover
            deflectionFlag=1;
            continue
        end
        if deflectionFlag==1
            %need to wait for current to stabilize, skip until
            %then & check when BG is stable to below the normal
            %threshold
            if max(BGRs)<dParam.stableLevel
                deflectionFlag=0;
            end
            continue
        end


        BGR=mean(BGR,2);
        BGR=repmat(BGR,1,sizeData(2));
        CVIsubR=dataFSCV-BGR;
        diffI=diff(CVIsubR(dParam.Voxid,:),3);        %differential to find glitches along time and ignore these
        glitchTS=find(diffI>dParam.glitchThres);
       % CVIsubR = filterDisplay(samplesperscan,[],CVIsubR);
        %CV scan window extends - dParam.BGwin offset
        rescanIDinit=t_point-dParam.BGWin(1);       
        if rescanIDinit<5
            rescanIDinit=5;
        end
        %CV scan window extends + dParam.BGwin offset 
        rescanIDend=t_point+dParam.BGWin(2);    
        if rescanIDend>(sizeData(2)-2)
            rescanIDend=sizeData(2)-2;
        end
        %check if CV scan time window matches supplied scan window
        if ~isempty(tParam.timeWin)
            %if cv scan inititial time point < window we supplied use window
            if rescanIDinit<tParam.timeWin(1)*tParam.samplingRate
                rescanIDinit=tParam.timeWin(1)*tParam.samplingRate;
            end
            %if cv scan end time point > window we supplied use window
            if rescanIDend>tParam.timeWin(2)*tParam.samplingRate
                rescanIDend=tParam.timeWin(2)*tParam.samplingRate;
            end
        end
        
        %scan CV time window for matching CV's
        for t_point_rescan=rescanIDinit:5:rescanIDend    %scan BG-sub CVs
            if abs(mean(CVIsubR(:,t_point_rescan)))>0 && ~any(intersect(t_point_rescan,glitchTS))
                 CVIsubRscan=CVIsubR(:,t_point_rescan-1:t_point_rescan+1);
                 CVIsubRscan=mean(CVIsubRscan,2);   %check if current time shows correlated CV
                 Dproj=pcr.Vc'*CVIsubRscan;         %projections of CVmeas onto relevant PCs of Ats
                 E=CVIsubRscan-(pcr.Vc*Dproj);      %calculate residuals
                 Q=diag(E'*E)';                 %calculate sum of squares
                  %check if max CV current is along expected anodal scan
                  peakindex=find(CVIsubRscan==max(CVIsubRscan));
                  peakindex=peakindex(1);
                 if CVIsubRscan(dParam.Voxid)>dParam.noiseThres && CVIsubRscan(dParam.Voxid)<dParam.saturationLevel && peakindex>samplesperscan*0.4 && peakindex<samplesperscan*0.85
                     if Q < pcr.QaDAPH
                         %check if CV residual less than threshold
                         %correlate to template for DA
                            KscaledDA=pcr.K(:,1)*Dproj(1,:);       %ideal CV template to correlate for current projected conc
                            Rtemp=corr2(KscaledDA,CVIsubRscan);
                            RR = abs(Rtemp);
                     else
                         RR=nan;
                     end
                 else
                     RR=nan;
                 end
            else 
                RR=0;
                Q=pcr.Qa;
            end
            xRR=RR>=dParam.RThres;
            %store sample & TS if R >= Rthres & passes further checks
            if max(RR)>=dParam.RThres 
                if Q<samples_Q(t_point_rescan) || isnan(samples_Q(t_point_rescan))
                    %if Q not stored for this TS yet or
                    %assigned previously a larger value, store
                    %it with current value
                        samples_Q(t_point_rescan)=Q;
                end
                %check if current indexed sample has smaller R than
                %this scan so we should replace  or add it
                numDA=numDA+1;      %increment count of found DA CV's
                if samples_DA(t_point_rescan)~=0
                    %check & remove count if already counted
                    %this TS
                    numDA=numDA-1;
                end
                if samples_DA(t_point_rescan)<max(RR)
                    samples_DA(t_point_rescan)=max(RR);
                    %replace with higher RR
                    samplesTemp.cv(:,numDA)=CVIsubRscan;
                    samplesTemp.ts(:,numDA)=t_point_rescan;
                    %store CV
                end

           %find local peak dopamine changes
           %find peak time along I vs T +/- 10s
           %Use DA projected conc instead
            Dproj=pcr.Vc'*CVIsubR;               
            Cu=pcr.F*Dproj;         %concentrartions predicted of each L component (ie DA/pH/background)
            IDAproj=Cu(1,:);            %current da projection
            DprojDA=IDAproj./(0.0282*ImaxBG*1.24)*1000;    %convert current to concentration
           %look in both directions
           maxTSPos=0;
           windowPos=0;
           windowNeg=0;
           maxTSNeg=0;
           minTSPos=0;
           minTSNeg=0;
           maxTS=0;
           if t_point_rescan-dParam.peakWindow>1 && t_point_rescan+dParam.maxDurationPeak<sizeData(2)
               %determine window time frames for local peak
               %make sure scan window for peak is within recorded range
                 windowPos=(t_point_rescan+1:t_point_rescan+dParam.peakWindow);    %forward window
                 windowNeg=(t_point_rescan-dParam.peakWindow:t_point_rescan);      %back window
           else
                if t_point_rescan+dParam.maxDurationPeak>=sizeData(2)
                   windowNeg=t_point_rescan-dParam.peakWindow:t_point_rescan;
                   windowPos=t_point_rescan+1:sizeData(2);
                end
                if t_point_rescan-dParam.peakWindow<1
                    windowPos=t_point_rescan+1:t_point_rescan+dParam.peakWindow;
                    windowNeg=1:t_point_rescan;
                end
           end
           maxTSPos=find(DprojDA==max(DprojDA(windowPos)));  %find max at increasing Time points from current scan
          %focus on only forward direction from tscan, and
          %around max now
           if ismember(maxTSPos,samplesTemp.maxTS)
               %if already identified, skip, go to next loop
               continue
           end
          if maxTSPos>=sizeData(2)
               %if last TS, pad inwards
              maxTSPos=maxTSPos-2;
          end
          winMaxDepartPos=maxTSPos+1:maxTSPos+dParam.maxDurationPeak;     %window to find local decrease at increasing TS
          if max(winMaxDepartPos)>sizeData(2)
              winMaxDepartPos=maxTSPos+1:sizeData(2);
          end
           peakDepartPosTS=find(DprojDA(winMaxDepartPos)<=DprojDA(maxTSPos)-dParam.peakDiff);
           %find TS of [DA] less than max by predefined dParam.peakDiff
           %value which defines peak departure
           if isempty(peakDepartPosTS)
               %if no local decrease, skip, go to next loop
                peakDepartPosTS=0;
                continue
           end
           peakDepartPosTS=peakDepartPosTS-1+winMaxDepartPos(1);    %convert relative TS to absolute TS
           peakDepartPosTSPos=peakDepartPosTS(peakDepartPosTS>(maxTSPos+dParam.minPeakWidth)); 
           %make sure the time that it decreases is not due to
           %transient glitch
           if ~isempty(peakDepartPosTSPos)
                peakDepartPosTSPos=peakDepartPosTSPos(1);
           else
               %if no TS exists for local decrease > minPeakWidth around max
               continue
           end
           windowNeghalf=maxTSPos-dParam.maxDurationPeak:maxTSPos-1;   %check for local decrease from peak in reverse direction
           if min(windowNeghalf)<1
               windowNeghalf=1:maxTSPos-1;
           end
           halfTSNeg=find(DprojDA(windowNeghalf)<=DprojDA(maxTSPos)-dParam.peakDiff);
           if isempty(halfTSNeg)
                halfTSNeg=0;        %if zero, we skip this loop below
                continue
           end
           halfTSNeg=halfTSNeg-1+windowNeghalf(1);
           if ~isempty(halfTSNeg)
                halfTSNeg=halfTSNeg(end);
                else
                   halfTSNeg=0;
                   continue
           end
            minTSPos=peakDepartPosTSPos;        %set "valley" in forward direction to that value above minPeakWidth
            minTSNeg=halfTSNeg;             %set "valley" in reverse direction 
            if minTSPos==sizeData(2)
                minTSPos=minTSPos-1;
            end
           if maxTSPos<minTSPos && maxTSPos>minTSNeg && minTSNeg>=2
                maxTS=maxTSPos;
               if ~DprojDA(maxTS)>mean(DprojDA(minTSPos-1:minTSPos+1))+dParam.peakDiff
                   %If peak DA identified is not above peakdiff to
                   %identified local decrease [DA] valley
                   maxTS=0;
                   continue
               end
               if ~DprojDA(maxTS)>mean(DprojDA(minTSNeg-1:minTSNeg+1))+dParam.peakDiff
                   maxTS=0;
               end
           else
               maxTS=0;
           end
           if maxTS==0
                continue
           end
           ss=[];
           %make sure identified signal not within 1 second of
           %already stored TS
           if numDA>1
               ss=samplesTemp.maxTS(samplesTemp.maxTS>0);
               %or if already within 1 s of identified maxTS
               if ~isempty(intersect(maxTS-10:maxTS+10,ss))
                   continue
               end
               ss_min=samplesTemp.minTS(samplesTemp.minTS>0);
               if ~isempty(intersect(minTSPos-10:minTSPos+10,ss_min))
                   continue
               end
           end
           if maxTS==0
               continue
           end
           %check if peak correlated and find best BG sub CV
            maxCV=CVIsubR(:,maxTS-CVavg:maxTS+CVavg); maxCVmean=mean(maxCV,2);
            Dprojmax=pcr.Vc'*maxCVmean;
            Emax=maxCVmean-(pcr.Vc*Dprojmax);
            Qmax=diag(Emax'*Emax)';
            KscaledDA=pcr.K(:,1)*Dprojmax(1,:);     %ideal CV template to correlate for current projected conc
            KscaledPH=pcr.K(:,2)*Dprojmax(2,:);        %ideal CV template to correlate to for pH
            KscaledBG=pcr.K(:,3)*Dprojmax(3,:);        %ideal CV template to correlate to for BG drift
            KscaledM=KscaledPH;                      %FOR anesthetised CFMEA ONLY NO MOVEMENT TEMPLATE FOR ANESTHETISED
            try
                KscaledM=pcr.K(:,4)*Dprojmax(4,:);         %ideal CV template to correlate to for movement
            catch
            end
            Rtemp1=corr2(KscaledDA,maxCVmean);      %correlation to DA
            RRpeak = abs(Rtemp1);                   %correlation to DA
            RM=corr2(KscaledM,maxCVmean);           %correlation to movement
            RRoutliers=abs(RM);
            diffCV=diff(maxCVmean,3);               %digital artifact (ie. Spike)
            if max(abs(diffCV))>dParam.glitchThres
                RRoutliers=1;
            end
            if ~(RRpeak>dParam.RThres) || abs(Qmax) > pcr.QaDAPH || abs(RRoutliers)>dParam.RThresOut
                %11/08/2017 added parenthesis for ~(xx>xx), loop never
                %worked properly before
                %if uncorrelated, try to scan around peak,
                %may have been disturbed by sudden artifact
                for scanPeakCor=maxTS-dParam.scanWindowMaxDetected:maxTS+dParam.scanWindowMaxDetected
                    if scanPeakCor-CVavg<1 || scanPeakCor+CVavg>sizeData(2)
                        continue
                    end
                    maxCV=CVIsubR(:,scanPeakCor-CVavg:scanPeakCor+CVavg); maxCVmean=mean(maxCV,2);
                    Dprojmax=pcr.Vc'*maxCVmean; KscaledDA=pcr.K(:,1)*Dprojmax(1,:);
                    Rtemp1=corr2(KscaledDA,maxCVmean);
                    RM=corr2(KscaledM,maxCVmean);
                    RRpeak = abs(Rtemp1);
                    if RRpeak>dParam.RThres && abs(RM)<dParam.RThresOut
                        maxTS=scanPeakCor;     %now max TS equals this nearby TS away from artifact
                        break
                        %changed to break from continue 11/08/2017, this
                        %loop never worked properly in past
                       
                    end
                end
            end

            peakRR=RRpeak>=dParam.RThres;           %logical 
            if isempty(peakRR) || sum(peakRR)==0
                %if still uncorrelated, skip, go to next loop
                continue
            end

           maxI=mean(CVIsubR(dParam.Voxid, maxTS-2:maxTS+2));      %get raw current at peak
           maxDA=mean(DprojDA(maxTS-2:maxTS+2));            %get proj da at peak
           halfTS=find(DprojDA<=maxDA/2);                   %get half life
           if isempty(halfTS)
               halfTS=sizeData(2);          %default end of sample is end of recording
           end
           halfTSafterMaxID=find(halfTS>maxTS);
           if isempty(halfTSafterMaxID)
               halfTSafterMax=sizeData(2);
           else
                halfTSafterMax=halfTS(halfTSafterMaxID);
                halfTSafterMax=halfTSafterMax(1);
           end
           signalLength=halfTSafterMax-maxTS;       %in samples how long the signal is
           %store max TS and duration
           numPeak=numPeak+1;
           samplesTemp.maxTS(numPeak)=maxTS;                  %TS of correlated peak
           samplesTemp.minTS(numPeak)=minTSPos;               %TS of peak decrease to threshold
           samplesTemp.duration(numPeak)=signalLength;        %in samples, half life
           samplesTemp.maxI(numPeak)=maxI;
           samplesTemp.maxDA(numPeak)=maxDA;                  %projected da at peak
           samplesTemp.maxCV(:,numPeak)=maxCVmean;          %cv at peak
           samplesTemp.R(numPeak)=RRpeak;
           
            noiseI=CVIsubR(:,t_point-2:t_point+2);                    %current noise around BG cursor (subtracted out current that is) ie across various potentials
            noiseI=abs(noiseI(dParam.Voxid,:));             %current noise around BG cursor (time) at oxidation potential
            noiseI=rms(noiseI);
            noiseDA=DprojDA(t_point-2:t_point+2);
            noiseDA=abs(noiseDA);
            noiseDA=rms(noiseDA);
            samplesTemp.noiseI(numPeak)=noiseI;
            samplesTemp.noiseDA(numPeak)=noiseDA;
        end

    end
end

%clean up\
    actualsamples=find(~isnan(samplesTemp.cv(1,:)));
    samplesTemp.cv=samplesTemp.cv(:,actualsamples);
    actualmaxCV=find(~isnan(samplesTemp.maxCV(1,:)));
    samplesTemp.maxCV=samplesTemp.maxCV(:,actualmaxCV);
end