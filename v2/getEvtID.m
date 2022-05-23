function evtcode = getEvtID(eventname,eventcodes)
%05/2022
%From eventcodes list (nums on left column and label on right column)
%output the eventcode num for the given label (event)

evtcodeid=find(strcmp(eventcodes(:,2),eventname));
evtcode=[];
if ~isempty(evtcodeid)
    %Cycle through found event id's since may be multiple
    for ii=1:length(evtcodeid)
        evtcode(ii)=str2num(eventcodes{evtcodeid(ii),1});
    end
else
    error('Event code invalid')
end