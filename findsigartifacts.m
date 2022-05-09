function artifacttss=findsigartifacts(allchdata,parameters,varargin)
%2/9/2019 update KM is now 3rd compoentn
%important to detect larger movement signals because dopamine 
%can be falsely detected around it, want to blank out entire period
%find idx of larger movements based on:
%existence on all channels
%length of movement signal > 5 seconds or defined period
%change to 3 s
samplerate=parameters.samplerate;
argnum=1;
relativets=0;
glitchids=[];
timewin=[];
durmin=3;       %in seconds, min length of artifact to 
selch=[];
thresfactor=.01;        %percent of recorded time window artifact can exist 
%thresfactor=.03;        %percent of recorded time window artifact can exist 
%thresfactor=.005;       %reduced 2/9/2019 from .01;


padwin=5;       % samples +/- artifacts found padd
for ii=1:size(allchdata,2)
    if ~isempty(allchdata{ii})
        selch=[selch ii];
    end
end
data=allchdata{selch(1)};

while argnum <= length(varargin)
    switch varargin{argnum}  
        case 'glitchids'
            argnum=argnum+1;
            glitchids=varargin{argnum};
        case 'timestamp'
            %output detected timestamps relative to indicated absolute time
            %stamp here (ie. neuralynx saved TS)
            argnum=argnum+1;
            relativets=varargin{argnum};            
        case 'window'
            %use defined time window for detection
            argnum=argnum+1;
            timewin=varargin{argnum};
        case 'artifactduration'
            argnum=argnum+1;
            durmin=varargin{argnum};
        case 'pad'
            argnum=argnum+1;
            padwin=varargin{argnum};
        otherwise
            error('badoption');
    end
    argnum = argnum + 1;
end


sizeData=size(data); 
t_start=5; t_end=sizeData(2)-10;
%use supplied time window for BG cursor, with dParam.BGwin param offsets
if ~isempty(timewin)
    t_start=timewin(1)*samplerate-parameters.BGWin(1);
    t_end=timewin(2)*samplerate+parameters.BGWin(2);
    if t_start<1
        t_start=5;
    end
    if t_end>sizeData(2)
        t_end=sizeData(2)-10;
    end
else
    timewin(1)=(t_start-1)./samplerate;
    timewin(2)=(t_end-1)./samplerate;
end
Rthres=parameters.RThres;
RthresOut=parameters.RThresOut;
BGthres=.85;
KDA=parameters.K(:,1);
KM=parameters.K(:,3);
KBG=parameters.K(:,4);
KPH=parameters.K(:,2);
BGavg=2;
QaDAPH=parameters.Qa;
CV_mat_DA=parameters.DAcvs;
CVs=parameters.CV;
Cts=parameters.Cts;

bgwin=600;       %look at +/- this many samples around BG sub point
Rs=[];
Ms=[];
Ts=[];
deflectionFlag=0;       %movement occured
count=1;
detected={};
totalscans=length(t_start:10:t_end);        %total # of scanning intervals
cutoffrate=totalscans*thresfactor;          %previously used entire time to get thres
%cutoffrate=10;      %if more than 10 times detected conside bad 10/2/2018
cutoffrate=thresfactor*bgwin;
counts=[];
badids={};
for ich=1:length(selch)
    data=allchdata{selch(ich)};
    detected(selch(ich)).MTS=[];
for tpoint=t_start:10:t_end 
    if ismember(tpoint-BGavg:tpoint+BGavg,glitchids)
        %if within glitch period skip
        continue
    end
    bg=data(:,tpoint-BGavg:tpoint+BGavg);
    bgstd=std(bg,0,2);%standard deviation of BG current versus time
    if max(bgstd)>parameters.meanBGThres 
        detected(selch(ich)).MTS=[detected(selch(ich)).MTS; tpoint];
        %check if mean BG is unstable, then go to next BG
        continue
    end
    if max(bgstd)>parameters.deflectionThres
        %if large charging current, maybe movement
        %related, causing time to recover
        deflectionFlag=1;
        detected(selch(ich)).MTS=[detected(selch(ich)).MTS; tpoint];
        continue
    end
    if deflectionFlag==1
        %need to wait for current to stabilize, skip until
        %then & check when BG is stable to below the normal
        %threshold
        detected(selch(ich)).MTS=[detected(selch(ich)).MTS; tpoint];
        if max(bgstd)<parameters.stableLevel
            deflectionFlag=0;
        end
        continue
    end     
    bgmean=mean(bg,2);
    refmatrix=repmat(bgmean,1,size(data,2));   %tile BGavg vector for entire meas matrix span
    subdata=data-refmatrix;
    %look at only discrete interval around BG point, not entire file.. so
    %same proportion cut off as needed, 10/18/2018
    currwin=tpoint-bgwin:tpoint+bgwin;
    currwin=currwin(currwin>0 & currwin<=t_end);
    subdatawin=subdata(:,currwin);
   % rr=abs(corr(subdata,KDA));
   % rr=corr(subdata,KDA);          %only positive corr's
    %rda=find(rr>Rthres);
    %idsinwin=[rda>timewin(1)*samplerate & rda<timewin(2)*samplerate];
    %rda=rda(idsinwin);
   % numda=length(rda);      %# of correlated dopamine signals (based on 1 sec sampling interval)
   % magda=mean(rr(rda));    %mean of correlation coef for detected dopamine signals 
   % rrm=abs(corr(subdata,KM));
    rrm=corr(subdatawin,KM);
    [maxrm,maxr]=max(rrm);
    rrbg=corr(subdatawin,KBG);
    [maxrb, maxrbg]=max(rrbg);
    if maxrm<RthresOut
        if maxrb>=BGthres
            idsinwin=maxrbg+currwin(1)-1;
        else
            %max R not above thres, go to next sub idx
            continue
        end
    end
    %rm=find(rrm>RthresOut);        
    idsinwin=maxr+currwin(1)-1;     %actual idx value
      %only get one movement artifact per cycle, otherwise could get 600 if bg
  %sub point at movement artifact
  
   % idsinwin=[rm>timewin(1)*samplerate & rm<timewin(2)*samplerate];
    %rm=rm(idsinwin);
    %numm=length(rm);    %# of correlated movement signals (based on 1 sec sampling interval)
    %magm=mean(rrm(rm)); %mean of correlation coef for detected movement signals 
    %Rs(count,1)=numda;      %# of DA correlated signals for this indx
    %detected(ich).Rs(count,2)=magda;      %avg R DA for this indx
    detected(selch(ich)).MTS=[detected(selch(ich)).MTS; idsinwin];       %# of M correlated signals for this indx

  
    
    % detected(ich).Ms(count,1)=numm;       %# of M correlated signals for this indx
    %detected(ich).Ms(count,2)=magm;       %avg abs R M  for this indx
   % detected(ich).Ts(count,1)=tpoint;     %all good time points (ie. without glitches)
    
  %  count=count+1;
end
    %if within +/- 1 samples consider same badidx
    badtsrange=[detected(selch(ich)).MTS-1; detected(selch(ich)).MTS; detected(selch(ich)).MTS+1];
    uniquets=unique(detected(selch(ich)).MTS);
  %  counts(selch(ich)).mts=uniquets;
    nummts=histc(badtsrange,uniquets);
   % counts(selch(ich)).nummts=nummts;
    badids{selch(ich)}=uniquets(find(nummts>=cutoffrate));
end

sharedbadids=[];
%get ids when movement artifact seen on at least 2 channels at same time
for ich=1:length(selch)
    otherchs=selch(~ismember(selch,selch(ich)));
    for iich=1:length(otherchs)
        sharedids=intersect(badids{selch(ich)},badids{otherchs(iich)});
        sharedbadids=[sharedbadids; sharedids];
    end
end

if length(selch)==1
    %only one channel cannot compare
    sharedbadids=badids{selch(ich)};
end

sharedbadids=unique(sharedbadids); %remove duplicates
extsharedbadids=bsxfun(@plus,(-padwin:padwin),sharedbadids);  %extend through window
extsharedbadids=unique(extsharedbadids);
extsharedbadids=extsharedbadids(extsharedbadids>0);
artifacttss=(extsharedbadids-1)./samplerate+relativets;


end