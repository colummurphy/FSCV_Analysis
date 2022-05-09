function [frt,trt,side,ats]=gettrialrts(ttls,tsttls,alnts,outcome)
%calculate behavior parameters like reaction time, previous trial type & result
%function called from gettrialorder.m
frt=nan;
trt=nan;
side=nan;
alnttl=[];
ats=alnts;
event_codes={
    '4'     'display_fix' ...
    '5'     'start_fix' ...
    '6'     'break_fix' ...
    '10'    'display_target' ...
    '11'    'start_target'  ...
    '12'    'break_target'  ...
    '14'    'error' ...
    '29'    'left_condition'    ...
    '30'    'right_condition'   ...
    '45'    'reward_big'    ...
    '46'    'reward_small'  ...
    };

reltsttl=tsttls-alnts;                     %get ts relative to alignment (outcome)
starttarget=str2num(event_codes{find(strcmpi(event_codes,'start_target')==1)-1});   %saccade to target
disptarget=str2num(event_codes{find(strcmpi(event_codes,'display_target')==1)-1});    %target on screen
startfix=str2num(event_codes{find(strcmpi(event_codes,'start_fix')==1)-1});     %saccade to fix
dispfix=str2num(event_codes{find(strcmpi(event_codes,'display_fix')==1)-1});      %fix cue on screen
left=str2num(event_codes{find(strcmpi(event_codes,'left_condition')==1)-1});      %left big reward cond
right=str2num(event_codes{find(strcmpi(event_codes,'right_condition')==1)-1});      %right big reward cond
breakfix=str2num(event_codes{find(strcmpi(event_codes,'break_fix')==1)-1});      %error break fix
breaktarg=str2num(event_codes{find(strcmpi(event_codes,'break_target')==1)-1});      %error break target
rewardbig=str2num(event_codes{find(strcmpi(event_codes,'reward_big')==1)-1});      %big reward
rewardsmall=str2num(event_codes{find(strcmpi(event_codes,'reward_small')==1)-1});      %small reward
if contains(outcome,'big')
    alnttl=rewardbig;
elseif contains(outcome,'small')
    alnttl=rewardsmall;
elseif contains(outcome,'fix')
    alnttl=breakfix;
elseif contains(outcome,'target')
    alnttl=breaktarg;
end

%get actual aln ts (may be off a little because of FSCV interpolation)
ttlout=find(ttls==alnttl);
if ~isempty(ttlout)
    tsout=tsttls(ttlout);
    tsfromout=alnts-tsout;
    [mints,closest]=min(abs(tsfromout)); %only get those TS's closest to outcome time
    if ~isempty(closest)
        ats=tsout(closest);
    end
end
%get target reaction time first, based on time from outcome alnts
idstarg=find(ttls==starttarget);
tstargon=nan;
if ~isempty(idstarg)
    tstarg=tsttls(idstarg);
    tsfromout=alnts-tstarg;
    [mints,closest]=min(tsfromout(tsfromout>0)); %only get those TS's before outcome, get closest to outcome time
    if ~isempty(closest) 
        tstargon=tstarg(closest);
        %find ts for target display previous to this
        idstargdisp=find(ttls==disptarget);
        if ~isempty(idstargdisp)
            tstargdisp=tsttls(idstargdisp);
            tsfromout=tstargon-tstargdisp;
            [mints,closest]=min(tsfromout(tsfromout>0));
            if ~isempty(closest)
                trt=tstargon-tstargdisp(closest);
            end
        end
    end
end
%get fix reaction time 
idsfix=find(ttls==startfix);
tsfixdisp=nan;
tsfixon=nan;
if ~isempty(idsfix)
    tsfix=tsttls(idsfix);
    tsfromout=tstargon-tsfix;       %closest to targ onset
    [mints,closest]=min(tsfromout(tsfromout>0)); %only get those TS's before outcome, get closest to outcome time
    if ~isempty(closest) 
        tsfixon=tsfix(closest);
        %find ts for target display previous to this
        idsfixdisp=find(ttls==dispfix);
        if ~isempty(idsfixdisp)
            tsfixdisp=tsttls(idsfixdisp);
            tsfromout=tsfixon-tsfixdisp;
            [mints,closest]=min(tsfromout(tsfromout>0));
            if ~isempty(closest)
                frt=tsfixon-tsfixdisp(closest);
                tsfixdisp=tsfixdisp(closest);
                if length(closest)>1
                    tsfixdisp=tsfixdisp(closest(1));
                end
            else
                tsfixdisp=nan;      %not found, maybe fixbreak
            end
        end
    end
end   

if isnan(tsfixdisp) && contains(outcome,'fix')
    tsfixdisp=tsfixon;      %for fixbreak trials
end

%get side
idsleft=find(ttls==left);
idsright=find(ttls==right);
tsright=inf;
tsleft=inf;
if ~isnan(tsfixdisp)
if ~isempty(idsright) 
    tsr=tsttls(idsright);
    tsfromout=tsfixdisp-tsr;       %closest to fix onset
    [mints,closest]=min(tsfromout(tsfromout>0)); %only get those TS's before outcome, get closest to outcome time
    if ~isempty(closest) 
            tsright=tsfixdisp-tsr(closest);   
    end
end
if ~isempty(idsleft) 
    tsl=tsttls(idsleft);
    tsfromout=tsfixdisp-tsl;       %closest to fix onset
    [mints,closest]=min(tsfromout(tsfromout>0)); %only get those TS's before outcome, get closest to outcome time
    if ~isempty(closest) 
            tsleft=tsfixdisp-tsl(closest);   
    end
end           
if tsright<tsleft
    side=1;         %if right condition appears closer to trial start cue, then right condition, 1
else
    side=0;
end
end

if isinf(tsright) && isinf(tsleft) && contains(outcome,'fix')
    %fix break trials, just take whatever side actually in trial
    if length(idsleft)>length(idsright)
        side=0;
    else
        side=1;
    end
end

end
