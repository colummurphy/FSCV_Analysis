function refreshbutton(hObject,~)
global plotParam processed hgui parameters settings
%refresh button pressed
%set new time scale & prepare fscv data in this selected time range

parameters.sampleratencs=processed.LFPread.LFPsamplingfreq;
%update plotParam (incase channels changed update)
%update t_start/end
plotParam.t_start = sscanf(get(hgui.tstart, 'string'), '%f').*parameters.samplerate+1; 
plotParam.t_end = sscanf(get(hgui.tend, 'string'), '%f').*parameters.samplerate+1; 
if isvalid(hgui.zoomts(1))
plotParam.zoomTS(1) = sscanf(get(hgui.zoomts(1), 'string'), '%f'); 
if plotParam.zoomTS(1)<(plotParam.t_start-1)/parameters.samplerate ||...
        plotParam.zoomTS(1)>(plotParam.t_end-1)/parameters.samplerate
    plotParam.zoomTS(1)=round((plotParam.t_start-1)/parameters.samplerate);
    set(hgui.zoomts(1),'string',num2str(plotParam.zoomTS(1)));
end
plotParam.zoomTS(2) = sscanf(get(hgui.zoomts(2), 'string'), '%f'); 
if plotParam.zoomTS(2)>(plotParam.t_end-1)/parameters.samplerate || ...
    plotParam.zoomTS(2)<=plotParam.zoomTS(1)
    plotParam.zoomTS(2)=round((plotParam.t_end-1)/parameters.samplerate);
    set(hgui.zoomts(2),'string',num2str(plotParam.zoomTS(2)));
end
end
plotParam.cmax = sscanf(get(hgui.cmax, 'string'), '%f'); 
[plotParam.map, plotParam.cmin]=setColorScale(plotParam.cmax);
plotParam.cscale=[plotParam.cmin plotParam.cmax];
if isvalid(hgui.blsub)
plotParam.buttonm =    hgui.checkfft.Value;
plotParam.blsub =    hgui.blsub.Value;
plotParam.mcsc = hgui.mcsc;
plotParam.fftclim(1) = sscanf(get(hgui.fftlimsdown, 'string'), '%f'); 
plotParam.fftclim(2) = sscanf(get(hgui.fftlimsup, 'string'), '%f');  
plotParam.powerscale(1) = sscanf(get(hgui.powerscale(1), 'string'), '%f'); 
plotParam.powerscale(2) = sscanf(get(hgui.powerscale(2), 'string'), '%f');  
plotParam.LFPscale(1) = sscanf(get(hgui.LFPscale(1), 'string'), '%f'); 
plotParam.LFPscale(2) = sscanf(get(hgui.LFPscale(2), 'string'), '%f');  
plotParam.filtlfp(1)=sscanf(get(hgui.bpf(1), 'string'), '%f');  
plotParam.filtlfp(2)=sscanf(get(hgui.bpf(2), 'string'), '%f');  
end
getplotsettings(plotParam.filtlfp,plotParam.cscNames,...
    plotParam.event_codes,parameters.sampleratencs,settings);   
%getplotsettings(plotParam.filtlfp,processed.LFPread.LFPchNames,...
 %   plotParam.event_codes,parameters.sampleratencs,settings);   

 plotParam.refresh=1;       %flag to refresh
%compile/organize data for plotting as loaded above in hgui
if isempty(hgui.ephysid)
compileloaded(hgui);    %update processed, parameters, plotParam globals
end
%plot ephys signals if exist
if ~isempty(processed.samplesNCS) && ~isempty(plotParam.cscNames)
    if isfield(plotParam,'event_codes')
        processed.behav=calcBehav(processed,plotParam.zoomTS,plotParam.event_codes);        %calculate RT, etc.
    end
    guiprolfp('getbursts');  %process ncs signals for plotting
    setguicloseup(hgui.closeup{2},'ncsids',plotParam.lfpid,...
                'sqenv','useprocessed','bursts','scale',plotParam.powerscale);   
    setguicloseup(hgui.closeup{3},'ncsids',[plotParam.eyeid...
        plotParam.pulseid],...
        'norm',1,'xlab'); 
    setguicloseup(hgui.closeup{3},'ncsids',plotParam.lickid,...
        'norm',1,'xlab','env','winlength',plotParam.winlengthphys,'hold');  
    setguicloseup(hgui.closeup{1},'ncsids',plotParam.lfpid,...
        'scale',plotParam.LFPscale,  'bursts');
    if  plotParam.buttonm
        setspectrum(processed,plotParam,hgui.fftplot,plotParam.mcsc); 
    end
end
        
end