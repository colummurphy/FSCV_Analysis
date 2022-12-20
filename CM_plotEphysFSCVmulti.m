function CM_plotEphysFSCVmulti(appFigure1, selectedFileName, ...
                selectedPathName, appDesColorPlots, appDesDopaminePlot, ...
                appDesTitleText, checkBoxes, checkNums)

global plotParam processed hgui parameters csc_map settings

% initialize global variables
plotParam={};   processed={}; hgui={}; parameters={}; csc_map={}; settings=[];

%selch=[1:3, 9:12];            
selch=[1 2 3 4];                %fscv channels to plot
differentialch=[];

ncschannels={};
settings=[];                    %load settings file (plotParam)
filtLFP=[13 28];  

t_start=0;  t_end=100;           %defaults in seconds
csc_map={};

% subset of events, should be displayed
disp_events=[ 4, 5, 9, 10,6,12, 14, 15, 29,30 45, 46]; 

%%
%read file

% CM - Modification here 
% File + Path are now set in App Designer startScript method
% [hgui.FileName,hgui.PathName] = uigetfile('*.*','Select file');
hgui.FileName = selectedFileName;
hgui.PathName = selectedPathName;


% get listing of every file in current directory
currentdir=dir(hgui.PathName);

% CM - returns starting char index of each filename in the current dir 
% that matches 'settings', if no match empty array. 
% There is no settings file in dir, returns a cell of empty arrays.
getset=regexpi({currentdir.name},'settings');  %load settings file if in directory

% CM - cellfun applies isempty to every array in getset
% find returns an vector of indexes of each non-zero element (non-empty array)
% note: no settings file, setfile=[]
setfile=find(~cellfun(@isempty,getset));


% CM - If the settings file is in the directory
% load the settings into the workspace variables
if ~isempty(setfile)
    load([hgui.PathName currentdir(setfile).name]); %load settings file if in directory
end

% CM - Add variables to global
% plotParam.selch = [1, 2, 3, 4]
% Flag to refresh
plotParam.selch=selch;      
plotParam.refresh=1;        


% CM - Check if 't_start' is a field in plotParam global
% load variables into plotParam if not in stored settings
% plotParam.t_start = 0;
% plotParam.t_end = 100;
if ~isfield(plotParam,'t_start')
    plotParam.t_start=t_start;
    plotParam.t_end=t_end;
end

%get file info, global hgui to get info on file/path, update hgui
CM_getfileinfo;      

% hgui.cname = 'patra125' for 125, fscv_multi_100
hgui.cname=[hgui.subject num2str(hgui.sessionnum) ];

% CM - Not run
if ~isempty(hgui.ephysid)
    %ephys session
    hgui.cname=[hgui.subject num2str(hgui.sessionnum) hgui.ephysid];
end

% get parameters
if strcmp(hgui.subject,'cleo')
    [parameters,csc_map,eventcodes]=CM_getparams(hgui.cname,'cleo',ncschannels);
elseif strcmp(hgui.subject,'patra')
    [parameters,csc_map,eventcodes]=CM_getparams('patrabipolar','default',ncschannels,'sessnum',hgui.sessionnum);
elseif strcmp(hgui.subject,'cfmea')
    [parameters,csc_map,eventcodes]=CM_getparams('cfmea','rodent',ncschannels,'sessnum',hgui.sessionnum);
end

% if there are no event_codes in plotParam, add event_codes + display_events
% plotParam.event_codes = { '4' 'display fix' '5' 'start fix' ... }
% plotParam.disp_events = [ 4, 5, 9, 10, ... ]
if ~isfield(plotParam,'event_codes')
    plotParam.event_codes=eventcodes;
    plotParam.disp_events=disp_events; 
end

% Function gets the site names from a file in the 1dr folder.
% filename = '1dr_cl6_pl1_xx_cl5_100' sites = {'cl6', 'pl1', 'xx' 'cl5'}
plotParam.sites=getfscvsitenames(hgui.PathName);

% get listing of every file in current directory
plotParam.dir=dir(hgui.PathName);

% load all ephys and fscv data
[processed.Iread,processed.LFPread,processed.samplesNCS]=...
    loadall(hgui.PathName, hgui.FileName,parameters,plotParam.selch,'dir',plotParam.dir);

% parameters.samplingratencs = sampling rate from fscv file
parameters.sampleratencs=processed.LFPread.LFPsamplingfreq;

    
% get settings, update plotParam global, depends on subject not recording 
% need nlx file in loaded file to open to get ch names
getplotsettings(filtLFP,ncschannels,plotParam.event_codes,parameters.sampleratencs,settings);       


% plotParam.t_start = 1
plotParam.t_start=round(plotParam.t_start*parameters.samplerate)+1;       %converted to sample point

% plotParam.t_end = 1001
plotParam.t_end=round(plotParam.t_end*parameters.samplerate)+1;


% set-up check boxes
CM_setUpCheckBox(checkBoxes, checkNums);


% CM - Start of figure code
plotParam.firstplot=1;


% Need color map for the App Designer figure
colormap(appFigure1, plotParam.map); 


%compile/organize data for plotting as loaded above in hgui
% CM - true for patra 125
if isempty(hgui.ephysid)
    %update processed, parameters, plotParam globals
    CM_compileloaded(hgui, appDesColorPlots, appDesDopaminePlot);    
end


% initialize processed.badids
processed.badids=[];


% Assign Title Text in App Designer
appDesTitleText.Value = [ hgui.trialtype ' | trial #: '...
    hgui.trialnum ' | session #: ' num2str(hgui.sessionnum) ];




%{
% CM - Commented out ephys code
% plot ephys signals if exist
% CM - cscNames = {} for patra 125
if ~isempty(processed.samplesNCS) && ~isempty(plotParam.cscNames)
    if isfield(plotParam,'event_codes')
        processed.behav=calcBehav(processed,plotParam.zoomTS,plotParam.event_codes);        %calculate RT, etc.
    end
           % setguicloseup(hgui.closeup{2},'ncsids',plotParam.lfpid,...
             %   'filt',plotParam.filtlfp,'sqenv','units','env-filt z score',...
             %   'scale',plotParam.powerscale,'norm',1,'winlength',plotParam.winlength);     
    guiprolfp('getbursts');  %process ncs signals for plotting
    if ~isempty(plotParam.lfpid)
    setguicloseup(hgui.closeup{2},'ncsids',plotParam.lfpid,...
                'sqenv','useprocessed','bursts','scale',plotParam.powerscale);  
            setguicloseup(hgui.closeup{1},'ncsids',plotParam.lfpid,...
        'scale',plotParam.LFPscale,'bursts');
    end
    setguicloseup(hgui.closeup{3},'ncsids',[plotParam.eyeid...
        plotParam.pulseid],...
        'norm',1,'xlab'); 
    setguicloseup(hgui.closeup{3},'ncsids',plotParam.lickid,...
        'norm',1,'xlab','env','winlength',plotParam.winlengthphys,'hold');  
    
end

%}
   
end