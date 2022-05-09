function finishsel(hobject,event)
global hgui
%finished selecting channels button
%now update / output

%get check box values from hgui struct
numplots=length(hgui.check1);
checks=zeros(2,numplots);
for id=1:numplots  
    checks(1,id)=hgui.check1{id}.Value;
    checks(2,id)=hgui.check2{id}.Value;  
end

%save variables
hgui.checks=checks;

end