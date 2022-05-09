%plotEphysFSCV.m
%05/11/2018, plot multiple fscv channels
%set     hgui.subject='patra';   %07/02/2021 in getfileinfo
global plotParam processed hgui parameters csc_map settings
plotParam={};   processed={}; hgui={}; parameters={}; csc_map={}; settings=[];
%selch=[1:3, 9:12];            %fscv channels to plot
selch=[1 2 3 4];
differentialch=[];

ncschannels={};
settings=[];                    %load settings file (plotParam)
filtLFP=[13 28];  

t_start=0;  t_end=100;           %defaults in seconds
csc_map={};
disp_events=[ 4, 5, 9, 10,6,12, 14, 15, 29,30 45, 46]; 

%%
%read file 
[hgui.FileName,hgui.PathName] = uigetfile('*.*','Select file');
currentdir=dir(hgui.PathName);
getset=regexpi({currentdir.name},'settings');  %load settings file if in directory
setfile=find(~cellfun(@isempty,getset));
if ~isempty(setfile)
    load([hgui.PathName currentdir(setfile).name]); %load settings file if in directory
end
plotParam.selch=selch;
 plotParam.refresh=1;       %flag to refresh
%load above defined variables into plotParam if not in stored settings
if ~isfield(plotParam,'t_start')
    plotParam.t_start=t_start;
    plotParam.t_end=t_end;
end


%get file info, global hgui to get info on file/path, update hgui
getfileinfo;      
hgui.cname=[hgui.subject num2str(hgui.sessionnum) ];
if ~isempty(hgui.ephysid)
    %ephys session
    hgui.cname=[hgui.subject num2str(hgui.sessionnum) hgui.ephysid];
end
if strcmp(hgui.subject,'cleo')
    [parameters,csc_map,eventcodes]=getparams(hgui.cname,'cleo',ncschannels);
elseif strcmp(hgui.subject,'patra')
    [parameters,csc_map,eventcodes]=getparams('patrabipolar','default',ncschannels,'sessnum',hgui.sessionnum);
elseif strcmp(hgui.subject,'cfmea')
    [parameters,csc_map,eventcodes]=getparams('cfmea','rodent',ncschannels,'sessnum',hgui.sessionnum);
end
if ~isfield(plotParam,'event_codes')
    plotParam.event_codes=eventcodes;
    plotParam.disp_events=disp_events; 
end
plotParam.sites=getfscvsitenames(hgui.PathName);
plotParam.dir=dir(hgui.PathName);
[processed.Iread,processed.LFPread,processed.samplesNCS]=...
    loadall(hgui.PathName, hgui.FileName,parameters,plotParam.selch,'dir',plotParam.dir);
parameters.sampleratencs=processed.LFPread.LFPsamplingfreq;

    
%get settings, update plotParam global, depends on subject not recording 
%need nlx file in loaded file to open to get ch names
%getplotsettings(filtLFP,ncschannels,plotParam.event_codes,parameters.sampleratencs,settings);       

getplotsettings(filtLFP,ncschannels,plotParam.event_codes,parameters.sampleratencs,settings);       
plotParam.t_start=round(plotParam.t_start*parameters.samplerate)+1;       %converted to sample point
plotParam.t_end=round(plotParam.t_end*parameters.samplerate)+1;

%load gui figure
fig1=figure;
numcolor=length(plotParam.selch);     %# color plots
numcolor=16;
numothers=3;     %% behavior plots
%figpos=[300,50,1000+plotParam.widen,900];
figpos=[50,50,1850,950];
set(fig1, 'position',figpos,'color',[1 1 1])
%setup gui & callbacks for gui figure selection
setfscvguilarge(fig1, numcolor);    %set up hgui variable figure    
if strcmp(hgui.subject,'cfmea')
    %different gui
   plotParam.colorplotsize=[200 125];
   setfscvguiarray(fig1, numcolor);    %set up hgui variable figure     
end
plotParam.firstplot=1;
%compile/organize data for plotting as loaded above in hgui
if isempty(hgui.ephysid)
    compileloaded(hgui);    %update processed, parameters, plotParam globals
end
%initialize processed.badids
processed.badids=[];
%update csc_map based on session id
%{
if strcmp(hgui.subject,'patra')
    [~,csc_map,~]=getparams('patrabipolar','default',ncschannels,'sessnum',hgui.sessionnum);
end
%}
%%
%setup click & button press callbacks
set(fig1,'windowbuttondownfcn',@mouseclick);
set(fig1,'WindowKeyPressFcn',@buttonpress);
set(fig1,'DeleteFcn',@closegui);
hgui.title=text(hgui.titletext,0, 0,...
    [hgui.trialtype ' | trial #: ' hgui.trialnum ' | session #: '...
    num2str(hgui.sessionnum)],...
    'interpreter','none'); 
        
%plot ephys signals if exist
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


            