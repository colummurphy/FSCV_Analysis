function [xcovda,xcovlag,xts]=daxcov(alldata,glitchids,varargin)
%export xcov matrices for each channel for each window in xcovda
%export xts nlx timestamps for start of each matrix
%use data stored in processed.nlx.resampled (ie. default high-beta envelope
%filtered) to calc cross cor with da in alldata;
%need global processed, plotParam, & parameters for reconvertncs called

global processed plotParam parameters
numpcs=3;       %da,ph,bg
if isfield(parameters,'numpcs')
    %how many components to use, just ph & da, or inlcude bg/move?
   numpcs=parameters.numpcs; 
end

samplerate=parameters.samplerate;
sampleratenlx=parameters.sampleratencs;
argnum=1;
relts=0;
%default sliding window sizes & overlap
winb=120;          
winbov=round(0.05*winb);
winp=30;        %window size for sliding window for xcov
winpstep=5;       %step size for sliding window for local peak from pcda
percentout=0.3;     %if >60% of window is bad then skip to next bg win sub
BGavgi=mean(alldata,2); 
ImaxBG=max(BGavgi);
Rthres=parameters.RThres;
KDA=parameters.K(:,1);
KPH=parameters.K(:,2);
KBG=parameters.K(:,3);
KM=parameters.K(:,4);
BGavg=2;
cvavg=1;
ats=parameters.CV;
cts=parameters.Cts;
wda=1;
wm=-1;       %weights for how much these bad signals matter in detecting optimal bg sub
wb=-0.6;
wp=-0.3;

while argnum <= length(varargin)
    switch varargin{argnum}       
        case 'timestamp'
            %output detected timestamps relative to indicated absolute time
            %stamp here (ie. neuralynx saved TS)
            argnum=argnum+1;
            relts=varargin{argnum};
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
        otherwise
            error('badoption');
    end
    argnum = argnum + 1;
end

xcovda={};
xts=[];

%Get info on loaded FSCV data
%initialize holding variable detected for detected parameters
deflectionFlag=0;       %movement occured, wait until stable

%generate iwin intervals
fullsize=size(alldata);
winbs=winb*samplerate;              %windows in samples
winbovs=winbov*samplerate;
winps=winp*samplerate;
winpsteps=winpstep*samplerate;

dacount=0;

%get nlx signals for indicated time window
nlxtss=fullsize(2)/2/samplerate+relts;
tracepad=[-(fullsize(2)-1)/2/samplerate fullsize(2)/2/samplerate];
%need global processed, plotParam, & parameters for reconvertncs
nlx=reconvertncs([],nlxtss,tracepad);   %get enveloped lfp data
covdata=nlx.resampled;
numcovchs=size(covdata,1);
dsn=round(sampleratenlx/samplerate);
tsnlx=nlx.ts;
alreadydsts=0;
tsdown=[];
downnlx=[];
covdata(isnan(covdata))=0;     %need to convert nan otherwise does not filter
for inlx=1:numcovchs
    downnlx(inlx,:)=downsampleSingleBatch(covdata(inlx,:),dsn);
    if alreadydsts==0
        %downsample timestamps once per file
        tsdown = reshape(tsnlx(1:dsn:end), 1, []);
        alreadydsts=1;
    end
end
    
lengthdata=0;
for iwin=winbovs+1:winbs-1:fullsize(2)-10
    %find optimal bg sub id in bg window
    %create window boundaries first
    iwinstartid=iwin-winbovs;
    iwinendid=iwinstartid+winbs-1+2*winbovs;
    idsiwin=iwinstartid:iwinendid;
    idswinbg=idsiwin(idsiwin>0 & idsiwin<fullsize(2));
    data=alldata(:,idswinbg);
    relativets=relts+(idswinbg(1)-1)/samplerate;
    if lengthdata==0
        lengthdata=length(idswinbg);        %store size, make sure same for each window
    end
    tsfscv=(idswinbg-1)./samplerate+relts;        %absolute fscv ts relative nlx ts
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
        'qthres',parameters.QaDAPH,'rthresm',parameters.RThresOut,...
        'rthresb',parameters.RThresOut,'rthresp',parameters.RThresOut,...
        'imaxbg',max(abs(bgmean)),'glitchids',glitchidsrel,'numpcs',numpcs);
    daproj=ipcr.DAiso;      %da projected for current bg window    
    
    lengthwin=0;        %make sure each win same length
    %slide window to get xcor matrices for current bg window
    for iiwin=winpsteps+1:winpsteps:length(idswinbg)
        %make windows for sliding for finding local maxima
        iiwinstart=iiwin;
    	iiwinendid=iiwin+winps-1;        
        idspwin=iiwinstart:iiwinendid;
        idsinpwin=idspwin(idspwin>0 & idspwin<size(data,2));
        daprojwin=daproj(idsinpwin);
        winids=idsinpwin;        
        if lengthwin==0
            lengthwin=length(winids);
        end
        if length(winids)~=lengthwin
            continue
        end
        nonnanids=find(~isnan(daprojwin));
        if length(nonnanids)<=5
            continue
        end
        dacount=dacount+1;
        winfirstts=tsfscv(iiwin);
               
        lfpwin=downnlx(:,idsinpwin);
        
        %need to replace nan values in da proj with i or zero other xcov
        %does not work
        %daprojwin(isnan(daprojwin))=1i;
        
        %cannot use above, basically interpolating nan values with zeros
        %makes correlation inaccurate
        %exclude windows of data with any nan values for now
        %later smooth data points surrounding nan values if sufficiently
        %small
        if any(isnan(daprojwin))
            for inlx=1:numcovchs  
            xcovda{inlx}(dacount,:)=nan(1,length(daprojwin)*2-1);
           xcovlag{inlx}(dacount,:)=nan(1,length(daprojwin)*2-1);
           end
        else
        for inlx=1:numcovchs      
            [xcovtemp,lagt]=xcov(daprojwin,lfpwin(inlx,:),'coeff');
           xcovda{inlx}(dacount,:)=xcovtemp;
           xcovlag{inlx}(dacount,:)=lagt;
        end
        end
        %test plot
        %figure; plot(xcovda(1:9,:)')
        %legend(nlx.cscNames{1:9})        
        xts(dacount)=winfirstts;    %save first timestamp that it is relative to
        
       %finish loop for current sliding window
     end
%finish bg sub loop
end

end