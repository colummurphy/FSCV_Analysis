function getdasignals(timewin,varargin)
%input argument is fscv timestamp id's  to process
%updates processed.glitchids & processed.detected

global processed parameters plotParam

detectionscheme=0;      %correlated search
argnum=1;
t1=timewin(1);      %time window to analyze in fscv samples
t2=timewin(2);
selch=plotParam.selch;
allchdata={};
badids=[];
%processed.glitchids={};
%processed.detected={};

for ii=1:length(selch)
    allchdata{selch(ii)}=processed.Iread(selch(ii)).data(:,t1:t2);
end

%find periods of recording where large & long-lasting movement
%artifacts & shared across channels
%badids=findsigartifacts(allchdata,parameters,'timestamp',...
 %       processed.LFPread.LFPts(1)+(t1-1)/parameters.samplerate);
if isfield(processed,'badids')
    if ~isempty(processed.badids)
    %ask user if replace currently detected bad ids with stored bad ids or
    %to merge
    if strcmp(varargin{argnum},'xcov')
        %xcov implies already did da detect and have bad ids for this sessi
        %use stored
        badids=processed.badids-t1+1;
    else
        do=input(['found already stored bad ids \n' ...
            '1. replace stored ids with detected bad ids \n'...
            '2. merge stored ids with detected bad ids \n'...
            '3. use stored ids \n'],'s');
        switch do
            case '1'
                %replace
                badids=findsigartifacts(allchdata,parameters,'artifactduration',2);  %timestamps
                badids=round(badids.*parameters.samplerate)+1;  %conver to samples
                processed.badids=badids+t1-1;   %store absolute ts ids from start of rec
            case '2'
                %merge with saved bad ids from user relative to current t1 window
                badids=findsigartifacts(allchdata,parameters,'artifactduration',2);  %timestamps
                badids=round(badids.*parameters.samplerate)+1;  %conver to samples
                badstored=processed.badids;
                badstored=badstored(badstored>=t1 & badstored <=t2);
                badstored=badstored-t1+1;
                badids=[badids; badstored];
                badids=unique(sort(badids));
                processed.badids=badids+t1-1;        %store absolute ts ids from start of rec
            case '3'
                %just use alaready stored ids ignore newly detected
                badids=processed.badids-t1+1;
            otherwise
                error('unknown command..');
        end
    end
    else
        %no stored bad ids, find new & create
        badids=findsigartifacts(allchdata,parameters,'artifactduration',2);  %timestamps
        badids=round(badids.*parameters.samplerate)+1;  %conver to samples
        processed.badids=badids+t1-1;   %store absolute ts ids from start of rec
    end
else
    %no stored bad ids, find new & create
    badids=findsigartifacts(allchdata,parameters,'artifactduration',2);  %timestamps
    badids=round(badids.*parameters.samplerate)+1;  %conver to samples
    processed.badids=badids+t1-1;   %store absolute ts ids from start of rec
end

for ii=1:length(selch)
    
    %find ids where HF glitches, to skip these data points in further
    %processing, 
    glitchids=findglitches(...
        processed.Iread(selch(ii)).rawdata(:,t1:t2),...
        parameters.glitchThres);
    %combine hf glitches w/ movement artifact periods & store
    if size(badids,2)<size(badids,1)
        %transpose so row array
        badids=badids';
    end
    processed.glitchids{selch(ii)}=sort([glitchids'; badids'])';

    %find local "transient' increases in DA & their parameters
    switch varargin{argnum}  
        case 'corr'
            disp(['detecting corr ch' num2str(selch(ii))])
        processed.detected{selch(ii)}=detectdatransientschunks(...
        processed.Iread(selch(ii)).data(:,t1:t2),parameters,...
        processed.glitchids{selch(ii)},plotParam,'timestamp',...
        processed.LFPread.LFPts(1)+(t1-1)/parameters.samplerate); 
    
        case 'uncorr'
            %must initialize variables since global, otherwise retained
            %over next loaded files
            processed.resampled=[];
            processed.resampledb=[];
            processed.nlxtsdown=[];
            processed.betadown=[];
            processed.sitelfp={};
            badzone=badids;
            parameters.badzone=badzone+t1-1;        %store non-window values in parameters
             disp(['detecting uncorr ch' num2str(selch(ii))])            
        if isfield(processed.Iread(selch(ii)),'events')
            %use timestamps stored in fscv.events(:,4) more accurate
            processed.detected{selch(ii)}=detectdauncorr(...
            processed.Iread(selch(ii)).data(:,t1:t2),parameters,...
            processed.glitchids{selch(ii)},'ts',...
            processed.Iread(selch(ii)).events(4,t1:t2),...
            'badzone',badids); 
        else
            processed.detected{selch(ii)}=detectdauncorr(...
            processed.Iread(selch(ii)).data(:,t1:t2),parameters,...
            processed.glitchids{selch(ii)},'timestamp',...
            processed.LFPread.LFPts(1)+(t1-1)/parameters.samplerate,...
            'badzone',badids);  
        end
        
        case 'xcov'
            disp(['detecting xcov ch' num2str(selch(ii))])        
            %first filter/down all lfp signals if first sel fscv ch, ii
            if ii==1
                resamplech=1;
                %first find # of chs/bands to resample
                 for ilfp=1:length(processed.LFPread.LFPchNames)
                     if ismember(ilfp,plotParam.lfpid)
                     sitelfp=plotParam.cscNames{plotParam.lfpid(ilfp)};
                     bandid=find(strcmp(sitelfp,plotParam.fbands)==1);
                     for iband=1:length(bandid)
                        bandfilt=plotParam.fbands{bandid(iband)+1};
                        resamplech=resamplech+1;
                     end
                     else
                         resamplech=resamplech+1;
                     end
                 end
                 iid=0;
            for ilfp=1:length(processed.LFPread.LFPchNames)
                sitelfp=[];
                samplesncs=processed.samplesNCS(:,ilfp)';
                if ismember(ilfp,plotParam.lfpid)
                    sitelfp=plotParam.cscNames{plotParam.lfpid(ilfp)};
                    bandid=find(strcmp(sitelfp,plotParam.fbands)==1);
                    for iband=1:length(bandid)
                        bandfilt=plotParam.fbands{bandid(iband)+1};
                        iid=iid+1;
                        %filter at "beta-band" as defined in fbands     
                        processed.resampled(iid,:)=samplesncs;   %store unfiltered data
                        %if another filter band specified 07/07/18
                        processed.resampledb(iid,:)=filterLFP(samplesncs,plotParam.ratelfp,bandfilt);
                        %square & envelope signal
                        winlength=round(plotParam.ratelfp*.5/mean(bandfilt));
                        if isfield(plotParam,'winlength')
                            %if window already defined in parameters load value
                            winlength=round(plotParam.ratelfp*plotParam.winlength);
                        end
                        processed.resampledb(iid,:)=processed.resampledb(iid,:).^2;   %get power V^2
                        processed.resampledb(iid,:)=smoothwin(processed.resampledb(iid,:),winlength);   %smoothing
                        processed.sitelfp{iid}=sitelfp;
                        processed.bandfilt{iid}=bandfilt;
                                         %downsample signal to same rate as fscv rate                
                processed.betadown(iid,:) = downsampleSingleBatch(processed.resampledb(iid,:),round(plotParam.ratelfp/parameters.samplerate));

                    end
               elseif ismember(ilfp,plotParam.lickid)
                   iid=iid+1;
                    sitelfp=plotParam.cscNames{plotParam.lickid};
                    processed.resampled(iid,:)=samplesncs;
                    %lick filter
                    processed.resampledb(iid,:)=filterLFP(samplesncs,plotParam.ratelfp,plotParam.filtlick);
                    %square & envelope signal
                    winlength=round(plotParam.ratelfp*.5/mean(plotParam.filtlick));
                    if isfield(plotParam,'winlengthphys')
                        %if window already defined in parameters load value
                        winlength=round(plotParam.ratelfp*plotParam.winlengthphys);
                    end
                    processed.resampledb(iid,:)=smoothwin(processed.resampledb(iid,:),winlength);   %smoothing
                    processed.sitelfp{iid}=sitelfp;
                        processed.bandfilt{iid}=plotParam.filtlick;
                                         %downsample signal to same rate as fscv rate                
                processed.betadown(iid,:) = downsampleSingleBatch(processed.resampledb(iid,:),round(plotParam.ratelfp/parameters.samplerate));

                %eye/pulse signal --> dont do anything already filt in reconvertFSCV.m            
                else
                    iid=iid+1;
                    sitelfp=plotParam.cscNames{ilfp};
                    processed.resampled(iid,:)=samplesncs;
                    processed.resampledb(iid,:)=samplesncs;
                    processed.sitelfp{iid}=sitelfp;
                        processed.bandfilt{iid}=[0 0];
                                         %downsample signal to same rate as fscv rate                
                processed.betadown(iid,:) = downsampleSingleBatch(processed.resampledb(iid,:),round(plotParam.ratelfp/parameters.samplerate));

                end
 
                if ilfp==1
                    %downsample timestamps to fscv rate
                    %these should be closely matching to those
                    %fscv timestamps synchronized to behav events
                    % stored in fscv.events(:,4)
                    processed.nlxtsdown=processed.LFPread.LFPts(1:round(plotParam.ratelfp/parameters.samplerate):end);
                    if length(processed.nlxtsdown)~=size(processed.betadown(1,:),2)
                        error('downsampled data and timestamps are different lengths')
                    end
                end
                
            end
            
            %end processing for just timestamps & resampled lfps (only done for first fscv
            %ch)
            end
                        
            tsdown=processed.nlxtsdown;
            tstrace=processed.detected{selch(ii)}.tstrace;  %# 
            numts=size(tstrace,2);      %number of sample points should match for all singals
            peakts=tstrace(:,median(1:numts)); %ts of peak da point
            rawts=processed.LFPread.LFPts;
            numtsraw=numts*round(plotParam.ratelfp/parameters.samplerate)+1;
            targids=[];
            rawtargids=[];
            dadata=processed.detected{selch(ii)}.datrace;
            %get filtered downsampeld signal around each detected da trace
            processed.detected{selch(ii)}.tslfp=[];
            processed.detected{selch(ii)}.tsnlx=[];
            processed.detected{selch(ii)}.lfptrace={};
            processed.detected{selch(ii)}.nlxtrace={};
             processed.detected{selch(ii)}.xcov={};  
             processed.detected{selch(ii)}.toounequalts=[];
             processed.detected{selch(ii)}.unequalts=[];    %store ids' of unequal timestamps between lfps
            if ~isempty(tstrace)
            for iid=1:size(tstrace,1)
                %go through all detected signals for curr fscv ch
                for ilfp=1:size(processed.resampled,1)
                    if ilfp==1
                    %find matching timestamps of nlx down (tsdown) & 
                    %fscv detected da timestamps (tstrace)
                    %common to all lfp channels
                    %get tsdown for peak da stored in peakts from tstrace
                    %with .1 resolution
                    peakid=find(tsdown<=peakts(iid)+.1 & tsdown>=peakts(iid)-0.1);
                    if length(peakid)>1
                         diffts=[];
                         for ip=1:length(peakid)
                             diffts(ip)=tsdown(peakid(ip))-peakts(iid);
                         end
                         selp=find(abs(diffts)==min(abs(diffts)));
                         peakid=peakid(selp);
                     end
                    %target ids should be same amount of samples around
                    %peak as stored tstrace
                    targids=peakid-median(1:numts)+1:peakid+median(1:numts)-1;
                    %targids=find(tsdown>=tstrace(iid,1) & tsdown<=tstrace(iid,end));
                          
                     tstarg=tsdown(targids); 
                     tserror=mean(abs(tstrace(iid,:)-tstarg));  %mean difference between ts's
                     if tserror>.1
                         warning(['# ' num2str(iid) ...
                             ': ts do not match for lfp down & fscv trace']);
                         processed.detected{selch(ii)}.unequalts=[...
                             processed.detected{selch(ii)}.unequalts iid];
                     end                     
                     %get ts ids for raw nlx (non-downsampled signal)
                     peakrawid=find(rawts<=peakts(iid)+.001 & rawts>peakts(iid)-.001);
                     if length(peakrawid)>1
                         diffts=[];
                         for ip=1:length(peakrawid)
                             diffts(ip)=rawts(peakrawid(ip))-peakts(iid);
                         end
                         selp=find(abs(diffts)==min(abs(diffts)));
                         peakrawid=peakrawid(selp);
                     end
                     rawtargids=peakrawid-median(1:numtsraw)+1:peakrawid+median(1:numtsraw)-1;
                     %rawtargids=find(rawts>=tstrace(iid,1) & rawts<=tstrace(iid,end));
                     tsrawtarg=rawts(rawtargids);
                     %check equal ts
                     %figure; plot(tstarg); hold on; plot(tsrawtarg(51:100:end));
                     processed.detected{selch(ii)}.tslfp(iid,:)=tstarg;
                     processed.detected{selch(ii)}.tsnlx(iid,:)=tsrawtarg;
                 end
                 processed.detected{selch(ii)}.lfptrace{ilfp}(iid,:)=processed.betadown(ilfp,targids);
                %get original ts's (not downsampled for raw data)
                 processed.detected{selch(ii)}.nlxtrace{ilfp}(iid,:)=processed.resampled(ilfp,rawtargids);
                 %calculate xcov for group of lfp's and each da for given
                 %check signals down/raw lfp side by side
                 %figure; plot(tstarg,processed.betadown(ilfp,targids))
                   % hold on; yyaxis right; plot(tsrawtarg,processed.resampled(ilfp,rawtargids))
                %trace (ie detected da peak)
                %finish getting all traces
                end
            end
            for ilfp=1:size(processed.resampled,1)
                lfpdata=processed.detected{selch(ii)}.lfptrace{ilfp};
                xdata={lfpdata dadata};        
                sitenames={processed.sitelfp{ilfp} plotParam.sites{selch(ii)}};
                xrates=[parameters.samplerate parameters.samplerate];
                processed.detected{selch(ii)}.xcov{ilfp}=xvardata(xdata,xrates,[],plotParam,[],[],...
                    'splitwin','sitename',sitenames,'fbands',plotParam.fbands,'type','none','noplot');
                %detected id's that are not nan, others invalid
                processed.detected{selch(ii)}.xcov{ilfp}.sel=...
                    find(~isnan(processed.detected{selch(ii)}.xcov{ilfp}.xcovda(:,1)));
                processed.detected{selch(ii)}.xcov{ilfp}.sitelfp=processed.sitelfp{ilfp};   
                processed.detected{selch(ii)}.xcov{ilfp}.bandfilt=processed.bandfilt{ilfp};  
                %if ilfp==1
               % plot(axx,processed.detected{2}.xcov{1}.xcovda);
               % hold on;
                %else
                %plot(axx,processed.detected{2}.xcov{1}.xcovda);
                %end 
            end
            end  
    end

end
        
end