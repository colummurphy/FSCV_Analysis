%plotEphysFSCV.m
%05/11/2018, plot multiple fscv channels
global plotParam processed hgui parameters csc_map settings
plotParam={};
plotParam.selch=[ 1 2 3 4];            %fscv channels to plot
%ncschannels={'cl1','cl1-cl4','lickx','eyed','pulse'};          %plot ncs channels
ncschannels={'p1','p3','p1-p3','pl2','pl3','cl1','cl1-cl4','cl4',...
    'cl5','cl4-cl5','s6','s5','s6-s5','s4','s5-s4','s3','s2','s3-s2',...
    's1','s2-s1','eyed','lickx','pulse'};          %plot ncs channels chronic58 all
ncschannels={'p3','p1-p3','pl3','cl1','cl1-cl4','cl5','cl4-cl5',...
    's1','s2-s1','eyed','lickx','pulse'};          %plot ncs channels chronic58 all

chronic58chconfig
%chronic38chconfig
%chronic67chconfig
chconfigsimple
cleo_simple17
cleo_simple16
%cleo_simple25
%cleo_simple14
settings=[];                    %load settings file (plotParam)
filtLFP=[17 25];  
plotParam.t_start=0;  plotParam.t_end=300;           %in seconds
csc_map={};
%[dParam,csc_map,event_codes]=getparams('patra');
%[parameters,csc_map,plotParam.event_codes]=getparams('patrabipolar','default',ncschannels);
[parameters,csc_map,plotParam.event_codes]=getparams('cleo16','noisycleo',ncschannels);
%[parameters,csc_map,plotParam.event_codes]=getparams('cleo17','default',ncschannels);
%[parameters,csc_map,plotParam.event_codes]=getparams('cleo25','noisycleo',ncschannels);
%[parameters,csc_map,plotParam.event_codes]=getparams('cleo14','noisycleo',ncschannels);

plotParam.disp_events=[ 5, 9, 10,6,12, 14, 15,  29, 30, 45, 46]; 

%%
%read file 
[hgui.FileName,hgui.PathName] = uigetfile('*.*','Select file');

%get file info, global hgui to get info on file/path, update hgui
getfileinfo;        

[processed.Iread,processed.LFPread,processed.samplesNCS]=...
    loadall(hgui.PathName, hgui.FileName,parameters,plotParam.selch);
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
numothers=3;     %% behavior plots
figpos=[300,50,1000+plotParam.widen,900];
set(fig1, 'position',figpos,'color',[1 1 1])
%setup gui & callbacks for gui figure selection
hgui=setfscvgui(fig1, numcolor,hgui);        

%compile/organize data for plotting as loaded above in hgui
compileloaded(hgui);    %update processed, parameters, plotParam globals

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
if ~isempty(processed.samplesNCS)
    if isfield(plotParam,'event_codes')
        processed.behav=calcBehav(processed,plotParam.zoomTS,plotParam.event_codes);        %calculate RT, etc.
    end
            setguicloseup(hgui.closeup{2},'ncsids',plotParam.lfpid,...
                'filt',[16 33],'sqenv','units','env-filt z score',...
                'scale',plotParam.powerscale,'norm',1,'winlength',plotParam.winlength);     
    setguicloseup(hgui.closeup{3},'ncsids',[plotParam.eyeid...
        plotParam.pulseid],...
        'norm',1,'xlab'); 
    setguicloseup(hgui.closeup{3},'ncsids',plotParam.lickid,...
        'norm',1,'xlab','env','winlength',plotParam.winlengthphys,'hold');  
    setguicloseup(hgui.closeup{1},'ncsids',plotParam.lfpid,...
        'scale',plotParam.LFPscale,'fscvch',plotParam.selch(1));
end


            