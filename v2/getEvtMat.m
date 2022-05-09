function evmat = getEvtMat(plotevents,trlist,eventcodes,varargin)
%05/2022
%Get event matrix containing time stamps for each selected plotting event
%(col) and trial (row) from trlist (variable in trlists containing all trs)

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
evmatFSCV=zeros(length(trids),length(plotevents));  %From evFSCV, check
evmatFSCVnlx=zeros(length(trids),length(plotevents));  %From evFSCV, check
for it=1:length(trids)
    outTS=trlist(trids(it)).ts; %Outcome TS, any event within current trial should occur before or at this TS (e.g. reward TS would be at this TS)
    evFSCV=trlist(trids(it)).eventsfscv;            %Events embedded in FSCV sample TS's
    evNlx=trlist(trids(it)).NlxEventTTL;            %Original Nlx events 
    evNlxTS=trlist(trids(it)).NlxEventTS;       %Nlx event TSs
    for ee=1:length(plotevents)
        curEv=plotevents{ee};   %Task event name to find
        curcode=getEvtID(curEv,eventcodes);     %Get numeric code for event name
        evIDs=find(evNlx==curcode);      %Find event id's
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
        [evIDsR,evIDsC]=find(evFSCV(:,5:12)==curcode);%Find event ids in fscv ev data
        evIDsR=sort(evIDsR);    %Get all rows regardless of column in order for code to work below (i.e. finding idx closest to outcome)
        targTS=nan;
        targTSf=nan;
        if ~isempty(evIDsR)
            evTS=evFSCV(evIDsR,1);%FSCV TS
            evTSnlx=evFSCV(evIDsR,4);%Nlx TS
            posTSid=find(evIDsR<=301);%Row 301 is outcome in FSCV ev data
            if ~isempty(posTSid)
                targTSf=evTS(posTSid(end));
            end
            posTSid=find((outTS-evTSnlx)>=-0.0001);
            if ~isempty(posTSid)
                targTS=evTSnlx(posTSid(end));
            end
        end
        evmatFSCV(it,ee)=targTSf;
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

