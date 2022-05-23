function [evmat,evmatfids,evmatfscv,evmatfscvnlx] = getEvtMat(plotevents,trlist,eventcodes,varargin)
%05/2022
%Get event matrix containing time stamps (Nlx TSs) for each selected plotting event
%Only get events closest to outcome
%(col) and trial (row) from trlist (variable in trlists containing all trs)
%evmatfscvnlx = Also export Nlx TS's as stored in events FSCV file (ie. after syncsigs
%function has interpolated NLx vents into FSCV events file at 100 ms res
%evmatfids = FSCV sample ids (in 601 trial sample file, aligned to outcome
%at sample 301 = 30s) for the events
%evmatfscv = FSCV time stamps (i.e FSCV timestamp for the sample ids)
%Will choose the earliest time stamp provided for a given event (e.g. for
%18,19, 45 all associated with reward_big, the earliest timestamp is
%output, e.g. in chronic 83, when used just 45, it was correct for all
%trials, but for one trial 18 showed up 100 ms earlier according to fscv
%timestamps).

outcomeFSCVsample=301;  %Sample for all trial based FSCV data at which outcome is delivered as set in syncsigs
trids=1:length(trlist);%Default all trials selected
argnum=1;
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'trialids'
            %User supplied trial id's selected
            argnum=argnum+1;
            trids=varargin{argnum};
    end
    argnum=argnum+1;
end

evmatrix=zeros(length(trids),length(plotevents));
evmatFSCV=zeros(length(trids),length(plotevents));  %From evFSCV, check, fscv TS
evmatFSCVids=zeros(length(trids),length(plotevents));  %From evFSCV, but sample IDs rather than fscv TS
evmatFSCVnlx=zeros(length(trids),length(plotevents));  %From evFSCV, check -- does not match nlx events (e.g. NlxEventTS, most likely during syncsigs did not explicitly write these in for each ev but relied on interp)
for it=1:length(trids)
    outTS=trlist(trids(it)).ts; %Outcome TS, any event within current trial should occur before or at this TS (e.g. reward TS would be at this TS)
    evFSCV=trlist(trids(it)).eventsfscv;            %Events embedded in FSCV sample TS's
    evNlx=trlist(trids(it)).NlxEventTTL;            %Original Nlx events 
    evNlxTS=trlist(trids(it)).NlxEventTS;       %Nlx event TSs
    for ee=1:length(plotevents)
        curEv=plotevents{ee};   %Task event name to find
        curcodes=getEvtID(curEv,eventcodes);     %Get numeric code for event name
        evIDs=find(ismember(evNlx,curcodes));      %Find event id's
        targTS=nan;%NaN default if not found
        if ~isempty(evIDs)
            evTS=evNlxTS(evIDs);%TS's of ev ids
            diffTS=outTS-evTS;  %Find time difference between ev TS and out TS
            posTSid=find(diffTS>=-0.0001);   %Positive TS's = before outcome, 0.1 ms allowed given rounding errors, Cannot use this for spikes
            if ~isempty(posTSid)
                targTS=evTS(posTSid(end));%Target TS is the last positive (minimum) closest to out TS
            end
        end
        evmatrix(it,ee)=targTS;        

        %Check evFSCV to see if matches nlx
        [evIDsR,evIDsC]=find(ismember(evFSCV(:,5:12),curcodes));%Find event ids in fscv ev data
        %[evIDsR,evIDsC]=find(evFSCV(:,5:12)==curcode);%Find event ids in fscv ev data
        evIDsR=sort(evIDsR);    %Get all rows regardless of column in order for code to work below (i.e. finding idx closest to outcome)
        targTS=nan;%target nlx TS in FSCV
        targTSf=nan;%target TS FSCV
        targIDf=nan;%target ID FSCV
        if ~isempty(evIDsR)
            evTS=evFSCV(evIDsR,1);%FSCV TS
            evTSnlx=evFSCV(evIDsR,4);%Nlx TS
            posTSid=find(evIDsR<=outcomeFSCVsample);%Row 301 aka outcomeFSCVsample is outcome in FSCV ev data
            if ~isempty(posTSid)
                targTSf=evTS(posTSid(end));
                targIDf=evIDsR(posTSid(end));
            end
            posTSid=find((outTS-evTSnlx)>=-0.0001);
            if ~isempty(posTSid)
                targTS=evTSnlx(posTSid(end));
            end
        end
        evmatFSCV(it,ee)=targTSf;
        evmatFSCVids(it,ee)=targIDf;
        evmatFSCVnlx(it,ee)=targTS;
    end
end
evdiff=evmatrix-evmatFSCVnlx;%difference in TS's retrieved from NLx orig and FSCV merged ev's --> Use NLx orig
if ~isempty(find(abs(evdiff)>0.1))
    %if error > 0.1 s, check
    [rows,cols]=find(abs(evdiff)>0.1);
    warning([num2str(length(rows)) ' differences between event TS from Nlx and FSCV > 0.1 s check e.g. TRID = ' num2str(rows(1))]);
end

evmat=evmatrix;
evmatfscvnlx=evmatFSCVnlx;
evmatfscv=evmatFSCV;
evmatfids=evmatFSCVids;

