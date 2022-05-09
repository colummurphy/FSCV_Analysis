function [optimalidx, isub]=findoptimalsubidx(data,parameters,glitchids,win,varargin)
%find optimal sub idx that maximizes # of correlated CV's & R values
%based on # of correlated signals found within win targeted
%interval & based on minimizing # of movement corelated signals
sizeData=size(data); 
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
t_start=5; t_end=sizeData(2)-10;
%use supplied time window for BG cursor, with dParam.BGwin param offsets
if ~isempty(win)
    t_start=win(1)*samplerate-parameters.BGWin(1);
    t_end=win(2)*samplerate+parameters.BGWin(2);
    if t_start<1
        t_start=5;
    end
    if t_end>sizeData(2)
        t_end=sizeData(2)-10;
    end
end
Rthres=parameters.RThres;
RthresOut=parameters.RThresOut;
KDA=parameters.K(:,1);
KPH=parameters.K(:,2);
KBG=parameters.K(:,3);
KM=parameters.K(:,4);
bgavg=2;
cvavg=1;
QaDAPH=parameters.Qa;
CV_mat_DA=parameters.DAcvs;
ats=parameters.CV;
cts=parameters.cts;
percentout=0.3;     %if >60% of window is bad then skip to next bg win sub

Rs=[];
Ms=[];
Ts=[];
deflectionFlag=0;       %movement occured
count=1;
for tpoint=t_start:10:t_end 
    if ismember(tpoint-bgavg:tpoint+bgavg,glitchids)
        %if within glitch period skip
        continue
    end
    bg=data(:,tpoint-bgavg:tpoint+bgavg);
    bgstd=std(bg,0,2);%standard deviation of BG current versus time
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
        if max(bgstd)<parameters.stableLevel
            deflectionFlag=0;
        end
        continue
    end     
    bgmean=mean(bg,2);
    refmatrix=repmat(bgmean,1,size(data,2));   %tile bgavg vector for entire meas matrix span
    subdata=data-refmatrix;
   % rr=abs(corr(subdata,KDA));
    rr=corr(subdata,KDA);          %only positive corr's
    rda=find(rr>Rthres);
    idsinwin=[rda>win(1)*samplerate & rda<win(2)*samplerate];
    rda=rda(idsinwin);
    numda=length(rda);      %# of correlated dopamine signals (based on 1 sec sampling interval)
    magda=mean(rr(rda));    %mean of correlation coef for detected dopamine signals 
    rrm=abs(corr(subdata,KM));
    rm=find(rrm>RthresOut);
    idsinwin=[rm>win(1)*samplerate & rm<win(2)*samplerate];
    rm=rm(idsinwin);
    numm=length(rm);    %# of correlated movement signals (based on 1 sec sampling interval)
    magm=mean(rrm(rm)); %mean of correlation coef for detected movement signals 
    Rs(count,1)=numda;      %# of DA correlated signals for this indx
    Rs(count,2)=magda;      %avg R DA for this indx
    Ms(count,1)=numm;       %# of M correlated signals for this indx
    Ms(count,2)=magm;       %avg abs R M  for this indx
    Ts(count,1)=tpoint;     %all good time points (ie. without glitches)
    
    count=count+1;
end
    
Rs(isnan(Rs(:,2)),2)=Rthres;
Rscores=zscore(Rs,1);
Rscorestotal=Rscores(:,1)+Rscores(:,2);     %total zscores of RDA mag & RDA frequency
Ms(isnan(Ms(:,2)),2)=RthresOut;  
Mscores=zscore(Ms,1);
Badts=find(Ms(:,1)>=3);         %any subtraction points that create more than 3 seconds of movement artifacts
Rscorestotalwithoutbad=Rscorestotal;
Rscorestotalwithoutbad(Badts)=[]; 
Goodid=find(Rscorestotal==max(Rscorestotalwithoutbad));
optimalidx=Ts(Goodid);      %optimal idx for bg subtraction
bg=data(:,optimalidx-bgavg:optimalidx+bgavg);
bgmean=mean(bg,2);
refmatrix=repmat(bgmean,1,size(data,2)); 
isub=data-refmatrix;
end