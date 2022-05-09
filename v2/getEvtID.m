function evtcode = getEvtID(eventname,eventcodes)
%05/2022
%From eventcodes list (nums on left column and label on right column)
%output the eventcode num for the given label (event)

evtcodeid=find(strcmp(eventcodes(:,2),eventname));
evtcode=[];
if ~isempty(evtcodeid)
    evtcode=str2num(eventcodes{evtcodeid,1});
else
    error('Event code invalid')
end