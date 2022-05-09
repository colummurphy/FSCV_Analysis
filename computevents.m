function info=computevents(ncsread,targetperiod,event_codes)
%computevents computes behavioral metrics and relative time points 
%of relevant behavioral events
%made for patra, can be used for cleo (1-dr in vcortex)
%04/10/2018
%07/27/2018 fixed so display fix used instead of left/right condition to
%calculate fix cue apperance (~1.2 s difference)
ncsdata=ncsread;            %ephys nlx read parameters & events/ttls/ts's
TS=ncsdata.LFPts;
relTS=TS-TS(1);                     %get ts of all recorded signal
eventTTL=ncsdata.LFPeventTTL;           %get event id's
eventTS=ncsdata.LFPeventTS - TS(1);     %get relative timestamps to zero
if isempty(targetperiod)
    targetperiod=[relTS(1) relTS(end)];     %default entire sample
end
alignTS=mean(targetperiod);             %middle of period is targeted alignment
if length(targetperiod)==1
    %only one number (relts) given to define targeted period to search
    %around for alignment mark
    alignTS=targetperiod;
    targetperiod=[relTS(1) relTS(end)];
end
starttarget=str2num(event_codes{find(strcmpi(event_codes,'start_target')==1)-1});   %saccade to target
ontarget=str2num(event_codes{find(strcmpi(event_codes,'display_target')==1)-1});    %target on screen
startfix=str2num(event_codes{find(strcmpi(event_codes,'start_fix')==1)-1});     %saccade to fix
onfix=str2num(event_codes{find(strcmpi(event_codes,'display_fix')==1)-1});      %fix cue on screen
left=str2num(event_codes{find(strcmpi(event_codes,'left_condition')==1)-1});      %left big reward cond
right=str2num(event_codes{find(strcmpi(event_codes,'right_condition')==1)-1});      %right big reward cond
breakfix=str2num(event_codes{find(strcmpi(event_codes,'break_fix')==1)-1});      %error break fix
breaktarg=str2num(event_codes{find(strcmpi(event_codes,'break_target')==1)-1});      %error break target
rewardbig=str2num(event_codes{find(strcmpi(event_codes,'reward_big')==1)-1});      %big reward
rewardsmall=str2num(event_codes{find(strcmpi(event_codes,'reward_small')==1)-1});      %small reward
altids=breaktarg;       %look for these events when trials not aligned to reward

%get targeted event labels and id's to search relative to
targetID=find(strncmpi(event_codes,'reward',6)==1); %get event labels
targetIDs(1)=str2num(event_codes{targetID(1)-1});
if length(targetID)>1
    for ii=2:length(targetID)
        targetIDs(ii)=str2num(event_codes{targetID(ii)-1});      %get event # (-1 from label)
    end
end
idxtargets=find(ismember(eventTTL,targetIDs)==1);       %find when targeted event occurs
tstargets=eventTS(idxtargets);                  %get ts of targeted events
tstargetsbound=tstargets(tstargets<targetperiod(2) & tstargets>targetperiod(1));    %target within period
if (isempty(tstargetsbound) || prod(abs(tstargetsbound-alignTS)>5))
    %too far from alingment time or
    %maybe trials loaded are for error or something else not aligned to
    %reward
    idsalt=eventTTL(eventTS>targetperiod(1) & eventTS<targetperiod(2));
    idxalt=find(ismember(idsalt,altids)==1);  %check if break target event
    %get event closest to align period
    aa=find(abs(eventTS(idxalt)-alignTS)==min(abs(eventTS(idxalt)-alignTS)));
    idxalt=idxalt(aa);
    %idxalt1=find(eventTTL==idsalt(idxalt));  %get ids for event
   % idsalt2=find(eventTS(idxalt1)>targetperiod(1) & eventTS(idxalt1)<targetperiod(2));
    %tstargetsbound=eventTS(idxalt1(idsalt2));
    tstargetsbound=eventTS(idxalt);
end    
if isempty(tstargetsbound)
    %not for reward or target break, then for fix
    idxalt=find(ismember(idsalt,breakfix)==1);  %check if break target event
    %get event closest to align period
    aa=find(abs(eventTS(idxalt)-alignTS)==min(abs(eventTS(idxalt)-alignTS)));
    idxalt=idxalt(aa);
    tstargetsbound=eventTS(idxalt);
end
if length(tstargetsbound)>1
    %more than 1 targeted alignemnt event found
    %get one closest to alignTS
    aa=find(abs(tstargetsbound-alignTS)==min(abs(tstargetsbound-alignTS)));
    tstargetsbound=tstargetsbound(aa);
end
idxTTLtarget=find(eventTS==tstargetsbound);     %new target idx (instead of reward)
targeteventcode=eventTTL(idxTTLtarget);    %get targeted event code

%find previous trial type (ie. rewarded/error/etc.)
priorevents=eventTTL(1:idxTTLtarget-1); %event TTLs prior to targeted event
navigateevents=find(ismember(priorevents,[rewardbig rewardsmall breakfix breaktarg])==1);
%get all idx's of interesting event TTLs to look for to define previous trial
if ~isempty(navigateevents)
    %trial detected before current targeted trial
    navigateevent=navigateevents(end);      %idx of previous event
    prevevent=eventTTL(navigateevent);
    preveventts=eventTS(navigateevent);     %rel ts of previous event
    info.prevevent=event_codes{find(strcmpi(event_codes,num2str(prevevent))==1)+1};   
        info.preveventts=preveventts;
    if isempty(info.prevevent)
        info.prevevent=nan;
        info.preveventts=nan;
    end
    %determine which side was previous trial
    prevpreveevents=eventTTL(1:navigateevent);
    prevside=find(ismember(prevpreveevents,[left right])==1);
    if ~isempty(prevside)
        prevside=prevside(end);
        info.prevside=event_codes{find(strcmpi(event_codes,num2str(eventTTL(prevside)))==1)+1};
    else
        info.prevside=nan;
    end

else
    info.prevevent=nan;
        info.preveventts=nan;

    info.prevside=nan;
end

%find next trial type (ie. rewarded/error/etc.) OUTCOME
nextevents=eventTTL(idxTTLtarget+1:end); %event TTLs after targeted event
navigateevents=find(ismember(nextevents,[rewardbig rewardsmall breakfix breaktarg])==1);
%get all idx's of interesting event TTLs to look for to define next trial
if ~isempty(navigateevents)
    %1st trial detected after target trial
    navigateevent=navigateevents(1);
    nexteventts=eventTS(navigateevent+idxTTLtarget);
    nextevent=eventTTL(navigateevent+idxTTLtarget);
    info.nextevent=event_codes{find(strcmpi(event_codes,num2str(nextevent))==1)+1};  
    info.nexteventts=nexteventts;    
else
    info.nextevent=nan;
        info.nexteventts=nan;

end

%find next trial first event (not outcome) based on fix appearance &
%left/right condition, usually left/right condition encoded before fix cue
nextevents=eventTTL(idxTTLtarget+1:end); %event TTLs after targeted event
navigateevents=find(ismember(nextevents,[onfix left right ])==1);
%get all idx's of interesting event TTLs to look for to define next trial
if ~isempty(navigateevents)
    %1st trial detected after target trial
    navigateevent=navigateevents(1);
    nexteventts=eventTS(navigateevent+idxTTLtarget);
    nextevent=eventTTL(navigateevent+idxTTLtarget);
    info.nexttrial=event_codes{find(strcmpi(event_codes,num2str(nextevent))==1)+1};  
    info.nexttrialts=nexteventts;
    nextfixcue=find(nextevents==onfix);   %fix cue appears
    if ~isempty(nextfixcue)
        nextfixcuets=eventTS(nextfixcue(1)+idxTTLtarget);
        info.nextfixcuets=nextfixcuets;
    else
       %no next fix, maybe paused task?
        info.nextfixcuets=nan;
    end
else
    info.nexttrial=nan;
    info.nexttrialts=nan;
    info.nextfixcuets=nan;

end

%find if current trial target is left or right 
condition=find(ismember(priorevents,[left right])==1);
if isempty(condition)
    condition=0;        %no relevant task event marks to identify 
    info.currside=nan;
    info.currbigside=nan;
    info.target_rt=[];
    info.fix_rt=[];
    info.prevtargcuets=nan;
    info.prevfixcuets=nan;
else
    condition=condition(end);
    info.currside=event_codes{find(strcmpi(event_codes,num2str(eventTTL(condition)))==1)+1};

    %find which side corresponds to big based on all events (Ie. which side
    %precedes big reward?)
    bigrew=find(eventTTL==rewardbig);
    currbigside=[];
    if ~isempty(bigrew)
        for ii=1:length(bigrew)
            %scroll big reward events until find what side it was on based on
            prev=eventTTL(1:bigrew(ii));
            idsidemarker=find(ismember(prev,[left right])==1);
            if ~isempty(idsidemarker)
                currbigside=eventTTL(idsidemarker(end));
                info.currbigside=event_codes{find(strcmpi(event_codes,num2str(currbigside))==1)+1};
                break;
            end
        end
    end

    %find reaction time to target, if successful trial
    ontime=find(priorevents==starttarget);  
    if ~isempty(ontime)
        ontime=eventTS(ontime(end));
        starttime=find(priorevents==ontarget);  starttime=eventTS(starttime(end));
        info.target_rt=ontime-starttime;
    else
        info.target_rt=[];
    end

    %find time between fix cue appearance and eye to fix cue
    ontime2=find(priorevents==startfix);  
    %fixed 07/28/2018
        %find fix cue appearance time
    fixappear=find(priorevents==onfix);

    if ~isempty(ontime2) && ~isempty(fixappear)
        ontime2=eventTS(ontime2(end));
       % starttime2=find(priorevents==onfix);  starttime2=eventTS(starttime2(end));
       % starttime2=condition;  
       % starttime2=eventTS(starttime2);  %base fix on on condition ID
        starttime2=eventTS(fixappear(end));     %add 07/27/2018 to fix calculation
        info.fix_rt=ontime2-starttime2;
    else
        info.fix_rt=[];
    end
    
    %find previous trial reward ts, fix cue ts, target cuet s
    ontimes=find(priorevents==ontarget);  
    prevtargts=nan;
    if length(ontimes)>1
        %get 2nd to last idx for previous trial
        %prevtargcue=ontimes(end-1);
        prevtargts=eventTS(ontimes(end-1));
    end
    info.prevtargcuets=prevtargts;    
    fixtimes=find(priorevents==onfix);  
    prevfixcuets=nan;
    if length(fixtimes)>1
        %get 2nd to last idx for previous trial
        prevfixcue=fixtimes(end-1);
        prevfixcuets=eventTS(fixtimes(end-1));
    end
    info.prevfixcuets=prevfixcuets;
    %find previous outcome
    %already stored in info.prevevent & info.preveventts
   
    
    
end

end
