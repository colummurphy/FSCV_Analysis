function infobehav=calcBehav(processed,targetperiod,event_codes)
%calculate behavior parameters like reaction time, previous trial type & result
%focusing on target period
%07/27/2018 fixed so display fix used instead of left/right condition to
%calculate fix cue apperance (~1.2 s difference)
ncsdata=processed.LFPread;
TS=ncsdata.LFPts;
relTS=TS-TS(1);                     %get ts of all recorded signal
eventTTL=ncsdata.LFPeventTTL;           %get event id's
eventTS=ncsdata.LFPeventTS - TS(1);     %get relative timestamps to zero
%get other relevant event labels & id's for parameters
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

infobehav.fix_rt=[];
infobehav.target_rt=[];
infobehav.currbigside=[];
infobehav.currside=[];
infobehav.prevevent=[];
infobehav.prevside=[];
infobehav.nextevent=[];

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
if isempty(tstargetsbound)
    %maybe trials loaded are for error or something else not aligned to
    %reward
    idsalt=eventTTL(eventTS>targetperiod(1) & eventTS<targetperiod(2));
    idxalt=find(ismember(idsalt,altids)==1);  %check if break target event
    if ~isempty(idxalt)
        idxalt1=find(eventTTL==idsalt(idxalt(1)));  %get ids for event
        idsalt2=find(eventTS(idxalt1)>targetperiod(1) & eventTS(idxalt1)<targetperiod(2));
        tstargetsbound=eventTS(idxalt1(idsalt2));
    else
        tstargetsbound=[];
    end
end    
if length(tstargetsbound)>1
    tstargetsbound=tstargetsbound(1);
end
if isempty(tstargetsbound)
    %no alignment during center of zoomTS
    return
end
idxTTLtarget=find(eventTS==tstargetsbound);
targeteventcode=eventTTL(idxTTLtarget);    %get targeted event code

%find previous trial type (ie. rewarded/error/etc.)
priorevents=eventTTL(1:idxTTLtarget-1); %event TTLs prior to targeted event
navigateevents=find(ismember(priorevents,[rewardbig rewardsmall breakfix breaktarg])==1);
%get all idx's of interesting event TTLs to look for to define previous trial
if ~isempty(navigateevents)
    %trial detected before target trial
    navigateevent=navigateevents(end);
    prevevent=eventTTL(navigateevent);
    infobehav.prevevent=event_codes{find(strcmpi(event_codes,num2str(prevevent))==1)+1};   
    %determine which side was previous trial
    prevpreveevents=eventTTL(1:navigateevent);
    prevside=find(ismember(prevpreveevents,[left right])==1);
    if ~isempty(prevside)
        prevside=prevside(end);
        infobehav.prevside=event_codes{find(strcmpi(event_codes,num2str(eventTTL(prevside)))==1)+1};

    end
end

%find next trial type (ie. rewarded/error/etc.)
nextevents=eventTTL(idxTTLtarget+1:end); %event TTLs after targeted event
navigateevents=find(ismember(nextevents,[rewardbig rewardsmall breakfix breaktarg])==1);
%get all idx's of interesting event TTLs to look for to define next trial
if ~isempty(navigateevents)
    %1st trial detected after target trial
    navigateevent=navigateevents(1);
    nextevent=eventTTL(navigateevent+idxTTLtarget);
    infobehav.nextevent=event_codes{find(strcmpi(event_codes,num2str(nextevent))==1)+1};   

end

%find if current trial target is left or right 
%condition=find(ismember(priorevents,[left right])==1);
condition=find(ismember(priorevents,[left right])==1);
if isempty(condition)
    condition=0;        %no relevant task event marks to identify 
    infobehav.currside=-1;
    infobehav.currbigside=-1;
    infobehav.target_rt=[];
    infobehav.fix_rt=[];
else
    condition=condition(end);
    infobehav.currside=event_codes{find(strcmpi(event_codes,num2str(eventTTL(condition)))==1)+1};

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
                infobehav.currbigside=event_codes{find(strcmpi(event_codes,num2str(currbigside))==1)+1};
                break;
            end
        end
    end

    %find previous successful trial target side & reward size
    navigateevents=find(ismember(priorevents,[rewardbig rewardsmall])==1);
    %get all idx's of interesting event TTLs to look for to define previous trial
    if ~isempty(navigateevents)
        %trial detected before target trial
        navigateevent=navigateevents(end);
        prevevent=eventTTL(navigateevent);
        infobehav.prevrewtype=event_codes{find(strcmpi(event_codes,num2str(prevevent))==1)+1};   
        %determine which side was previous trial
        prevpreveevents=eventTTL(1:navigateevent);
        prevside=find(ismember(prevpreveevents,[left right])==1);
        if ~isempty(prevside)
            prevside=prevside(end);
            infobehav.prevrewside=event_codes{find(strcmpi(event_codes,num2str(eventTTL(prevside)))==1)+1};

        end
    end

    %find reaction time to target, if successful trial
    ontime=find(priorevents==starttarget);  
    if ~isempty(ontime)
        ontime=eventTS(ontime(end));
        starttime=find(priorevents==ontarget);  starttime=eventTS(starttime(end));
        infobehav.target_rt=ontime-starttime;
    end

    %find time between fix cue appearance and eye to fix cue
    ontime2=find(priorevents==startfix);  

    %find fix cue appearance time
    fixappear=find(priorevents==onfix);

    if ~isempty(ontime2) && ~isempty(fixappear)
        ontime2=eventTS(ontime2(end));
       % starttime2=find(priorevents==onfix);  starttime2=eventTS(starttime2(end));
       % starttime2=condition;  
       % starttime2=eventTS(starttime2);  %base fix on on condition ID
        starttime2=eventTS(fixappear(end));     %add 07/27/2018 to fix calculation
        infobehav.fix_rt=ontime2-starttime2;
    else
        infobehav.fix_rt=[];
    end

end

end
