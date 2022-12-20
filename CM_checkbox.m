function CM_checkbox(chanNumSelected, checkVal, ...
                    appDesignColorFig, dopaminePlot)

global plotParam processed hgui parameters

%{
  Function updates the selected channels array (selch) in 
  the global plotParam.
  
  The data is read in from the file.
  The color and dopamine plots are updated.
%}

% if the channel is not in selch array AND the check is selected
if ~ismember(chanNumSelected, plotParam.selch) && checkVal == 1
    % add to selch array 
    plotParam.selch=sort([chanNumSelected plotParam.selch]);

% if channel is in the selch array AND the check is deselected    
elseif ismember(chanNumSelected, plotParam.selch) && checkVal == 0
    % remove from selch array
    plotParam.selch(plotParam.selch == chanNumSelected)=[];
end

[processed.Iread, processed.LFPread, processed.samplesNCS]=...
    loadall(hgui.PathName, hgui.FileName, parameters, plotParam.selch);

CM_compileloaded(hgui, appDesignColorFig, dopaminePlot);

end
