function CM_menuselect(appDesColorPlots, appDesDopaminePlot, appDesTitleText)
global hgui plotParam parameters processed csc_map settings

%{
This function loads a new file into the analysis program.
%}

%select file
[hgui.FileName,hgui.PathName] = uigetfile('*.*','Select file');
cd(hgui.PathName)
%get file info, global hgui to get info on file/path, update hgui
        
% CM mod here, new file    
CM_getfileinfo; 
        
% existing session
cursess=hgui.cname;
        
% new session - hgui.cname = 'patra125' for 125, fscv_multi_100
hgui.cname=[hgui.subject num2str(hgui.sessionnum) ];
        
% if different session
if ~isequal(hgui.cname,cursess)
    % if ephys session
    if ~isempty(hgui.ephysid)
        hgui.cname=[hgui.subject num2str(hgui.sessionnum) hgui.ephysid];
    end

    ncschannels={};
    if strcmp(hgui.subject,'cleo')
        % CM - new file
        [parameters,csc_map,eventcodes]=CM_getparams(hgui.cname,'cleo',ncschannels);
    else
        % CM - new file
        [parameters,csc_map,eventcodes]=CM_getparams('patrabipolar','default',ncschannels,'sessnum',hgui.sessionnum);
    end
end

% Function gets the site names from a file in the 1dr folder.
% filename = '1dr_cl6_pl1_xx_cl5_100' sites = {'cl6', 'pl1', 'xx' 'cl5'}
plotParam.sites=getfscvsitenames(hgui.PathName);
        
% initialize
processed={}; parameters.badzone=[];
        
% get listing of every file in current directory
plotParam.dir=dir(hgui.PathName);
        
% load all ephys and fscv data
[processed.Iread,processed.LFPread,processed.samplesNCS]=...
   loadall(hgui.PathName, hgui.FileName,parameters,plotParam.selch,'dir',plotParam.dir);

% parameters.samplingratencs = sampling rate from fscv file
parameters.sampleratencs=processed.LFPread.LFPsamplingfreq;
        
%update plotParam (incase channels changed update)
getplotsettings(plotParam.filtlfp,plotParam.cscNames,...
       plotParam.event_codes,parameters.sampleratencs,settings);       
        
disp(['loading: ' hgui.FileName]);
        
%compile/organize data for plotting as loaded above in hgui
if isempty(hgui.ephysid)
   % CM - new file
   CM_compileloaded(hgui, appDesColorPlots, appDesDopaminePlot);
end

        %{
        % CM - Commented out ephys code.
        % plot ephys signals if exist
        if ~isempty(processed.samplesNCS) && ~isempty(plotParam.cscNames)
            if isfield(plotParam,'event_codes')
                processed.behav=calcBehav(processed,plotParam.zoomTS,plotParam.event_codes);        %calculate RT, etc.
            end
            guiprolfp('getbursts');  %process ncs signals for plotting
            if ~isempty(plotParam.lfpid)
            setguicloseup(hgui.closeup{2},'ncsids',plotParam.lfpid,...
                        'sqenv','useprocessed','bursts');  
                    setguicloseup(hgui.closeup{1},'ncsids',plotParam.lfpid,...
                'scale',plotParam.LFPscale,'bursts');
            end
            setguicloseup(hgui.closeup{3},'ncsids',[plotParam.eyeid...
                plotParam.pulseid],...
                'norm',1,'xlab'); 
            setguicloseup(hgui.closeup{3},'ncsids',plotParam.lickid,...
                'norm',1,'xlab','env','winlength',plotParam.winlengthphys,'hold');  
            
            if plotParam.buttonm==1
                setspectrum(processed,plotParam,hgui.fftplot,plotParam.mcsc); 
            end
        end
        %}

appDesTitleText.Value = [hgui.trialtype ' | trial #: '...
      hgui.trialnum ' | session #: ' num2str(hgui.sessionnum)];            

end

