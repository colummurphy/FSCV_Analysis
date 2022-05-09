function nlx=reconvertncs(ids,targetts,padding)
%targetts is targeted time point (nlx time point +/- 0.0001 s accurate)
%padding is # seconds before after to save
%save spectra at targetts only
global processed plotParam parameters
tracewin=padding;
samplerate=round(parameters.sampleratencs);
data=processed.samplesNCS;

argnum = 1;
numfilts=0;     %number of filter bands for processing
filtbands=[];
filtlim=[];     %filt band
winlength=[];            %smoothing window
nlx.resampled=[];
cscids=ids;

if isempty(cscids)
    cscids=1:size(data,2);    
end

nlx.eyeid=find(ismember(cscids,plotParam.eyeid)==1);
nlx.pulseid=find(ismember(cscids,plotParam.pulseid)==1);
nlx.lickid=find(ismember(cscids,plotParam.lickid)==1);
nlx.physid=find(ismember(cscids,plotParam.physid)==1);
nlx.lfpid=find(ismember(cscids,plotParam.lfpid)==1);

%get labels from provided ids
nlx.cscNames=plotParam.cscNames(cscids);

%get filter bands to process for lfps
if isfield(plotParam,'filtbetal')
    numfilts=numfilts+1;
    filtbands(numfilts,:)=plotParam.filtbetal;
    winlength(numfilts)=round(samplerate*.5/mean(filtbands(numfilts,:)));    
    nlx.filtbetal=filtbands(numfilts,:);
end
if isfield(plotParam,'filtbetah')
    numfilts=numfilts+1;
    filtbands(numfilts,:)=plotParam.filtbetah;
    winlength(numfilts)=round(samplerate*.5/mean(filtbands(numfilts,:)));    
    nlx.filtbetah=filtbands(numfilts,:);
end
if isfield(plotParam,'filtgammal')
    numfilts=numfilts+1;
    filtbands(numfilts,:)=plotParam.filtgammal;
    winlength(numfilts)=round(samplerate*.5/mean(filtbands(numfilts,:)));   
    nlx.filtgammal=filtbands(numfilts,:);
end

nlx.filtbands=filtbands;
nlx.winlength=winlength;
TS=processed.LFPread.LFPts;

%scroll targeted timestamps & process & store signals around these
for itarget=1:length(targetts)
    %get nlx sample idx for each nlx referred targetts supplied
    targetidx=find(round(TS.*1000)/1000==round(targetts(itarget)*1000)/1000);
    %get +/- idxes from this central point
    idx1=targetidx+tracewin(1)*samplerate;
    idx2=targetidx+tracewin(2)*samplerate;
    ids=idx1:idx2;
    idsneg=ids(ids<=0);
    idsover=ids(ids>size(data,1));
    idsinwin=ids(ids>0 & ids<=size(data,1));
    %fill negative/over time values with nan
    fillnegids=repmat(nan,length(idsneg),1);
    filloverids=repmat(nan,length(idsover),1);

    rets=[fillnegids' TS(idsinwin) filloverids'];
    firstts=TS(idsinwin(1));
    lastts=TS(idsinwin(end));
    idEvents=find(ismember(processed.LFPread.LFPeventTTL,plotParam.disp_events)==1);
    TSEvents=processed.LFPread.LFPeventTS(idEvents);
    idEventsinwin=idEvents(TSEvents>=firstts & TSEvents<=lastts);
    Eventsinwin=processed.LFPread.LFPeventTTL(idEventsinwin);
    tsEventsinwin=processed.LFPread.LFPeventTS(idEventsinwin);
    nlx(itarget).eventts=tsEventsinwin;
    nlx(itarget).eventcode=Eventsinwin;
    nlx(itarget).ts=rets;
    for ii=1:length(cscids)
        if ismember(cscids(ii),[plotParam.physid])
            %not lfp signal, do nothing
            nlx(itarget).resampled(ii,:)=[fillnegids; ...
                data(idsinwin,cscids(ii));...
                    filloverids];
        else
            %lfp signal, filter, envelope at selected filter bands
            for ifilt=1:numfilts
                filtlim=filtbands(ifilt,:);
                filtdata=filterLFP(data(idsinwin,cscids(ii)),samplerate,filtlim);
                nlx(itarget).filtered{ii}(ifilt,:)=[fillnegids;...
                    filtdata; filloverids];
                filtdatasq=filtdata.^2;   %get power V^2
                envdata=smoothwin(filtdatasq,winlength(ifilt));
                nlx(itarget).enveloped{ii}(ifilt,:)=[fillnegids;...
                    envdata; filloverids]; 
                if ifilt==2
                    %only beta high stored in resampled variable, assume 2nd
                    nlx(itarget).resampled(ii,:)=[fillnegids;...
                        envdata; filloverids]; 
                end
            end
            nlx(itarget).raw(ii,:)=[fillnegids; ...
                data(idsinwin,cscids(ii));...
                    filloverids];
        end
    end
end

end