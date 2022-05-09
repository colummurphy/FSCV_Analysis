function checkbox(source,event)
global plotParam processed hgui parameters
val = source.Value;     %checked?  
digids=regexp(source.String,'\d');
chnumsel=str2num(source.String(digids));
if ~ismember(chnumsel,plotParam.selch) && val==1
    %not yet in channel database, but now selected, add 
    plotParam.selch=sort([chnumsel plotParam.selch]);
elseif ismember(chnumsel,plotParam.selch) && val==0
    %member of channel database, deselected, remove
    plotParam.selch(plotParam.selch==chnumsel)=[];
end
[processed.Iread,processed.LFPread,processed.samplesNCS]=...
    loadall(hgui.PathName, hgui.FileName,parameters,plotParam.selch);
compileloaded(hgui);    %update processed, parameters, plotParam globals

end
