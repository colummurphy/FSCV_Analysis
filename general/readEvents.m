function [nlx_events,events]=readEvents(eventsTimeStamps, eventsToDisplay)
events=[];
nlx_events=eventsTimeStamps(5:end,:);
aa=1;
    aaa=0;
    for ii=1:length(eventsToDisplay)
        events(ii,1)=eventsToDisplay(ii);
        for aa=1:size(nlx_events,1)
            aaa=find(nlx_events(aa,:)==eventsToDisplay(ii));
            if aa==1
                if length(aaa)~=0
                    events(ii,2:length(aaa)+1)=aaa;
                end
            end
            if aa~=1
                if length(aaa)~=0 && sum(events(ii,:))~=0
                    events(ii,length(events(ii,:))+2:length(events(ii,:))+length(aaa)+1)=aaa;
                end 
            end
        end
    end
end