function []=setfscvguichannelselection(cscmap,varargin)
%set up figure & individual plots
global hgui hsel plotParam parameters processed 
argnum=1;
chmap=[];
chmap={cscmap{2:2:end}};   %just ch names
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'selectchs'
            %provided chmap for current session, does not have ch# between
            chmap=cscmap;
    end
    argnum=argnum+1;
end
hsel.figsel=figure;
ncschannels=plotParam.cscNames;
numsel=length(chmap);     %# channels to select from
selnames={};
figsize=[700 900];
ofigpos=get(hgui.hf,'position');
figpos=[ofigpos(1),ofigpos(2)+50,figsize(1),figsize(2)];
set(hsel.figsel, 'position',figpos,'color',[1 1 1])
set(0,'CurrentFigure',hsel.figsel);    %set figure handle to current figure
posy=figsize(2)/(numsel+1);
posy=posy:posy:figsize(2);
posy=posy(1:end-1);
posx=50;
if numsel>10 && numsel<=20
    %split 2 columns
    posy=figsize(2)/10;
    posy=25:posy:figsize(2);
    posx=figsize(1)/3;
    posx=25:posx:figsize(1); 
    posx=posx(1:end-1);
elseif numsel>20 && numsel<=30
    %split 3 columns
    posy=figsize(2)/10;
    posy=25:posy:figsize(2);
    posx=figsize(1)/4;
    posx=25:posx:figsize(1); 
    posx=posx(1:end-1);
elseif numsel>30 && numsel<=40
    %split 4 columns
    posy=figsize(2)/10;
    posy=25:posy:figsize(2);
    posx=figsize(1)/5;
    posx=25:posx:figsize(1); 
    posx=posx(1:end-1);
elseif numsel>40 && numsel<=50
    %split 5 columns
    posy=figsize(2)/10;
    posy=25:posy:figsize(2);
    posx=figsize(1)/6;
    posx=25:posx:figsize(1); 
    posx=posx(1:end-1);
elseif numsel>50 && numsel<=60
    %split 6 columns
    posy=figsize(2)/10;
    posy=25:posy:figsize(2);
    posx=figsize(1)/7;
    posx=25:posx:figsize(1); 
    posx=posx(1:end-1);
    elseif numsel>60
            posy=figsize(2)/10;
    posy=25:posy:figsize(2);
    posx=figsize(1)/8;
    posx=25:posx:figsize(1); 
    posx=posx(1:end-1);
end
posy=fliplr(posy);
plotvalue=[];
selnames={};
ich=1;
for ii=1:numsel
   selnames{ii}=chmap{ii};
   if ismember(selnames{ii},ncschannels)
       %already selected, set value to 1 to check
       plotvalue(ii)=1;
   else
       plotvalue(ii)=0;
   end   
end
%set up plotting axes
column=1;
row=1;
for ich=1:length(selnames)
    %ch selections
    hsel.chsel{ich}= uicontrol('Style', 'checkbox','value',plotvalue(ich),... 
        'string',selnames{ich},'Position', [posx(column) ...
        posy(row) 15 15]); 
    hsel.chseltxt{ich}= uicontrol('Style', 'text',... 
        'string',selnames{ich},'Position', [posx(column)+20 ...
        posy(row)-10 50 25],'backgroundcolor',[1 1 1],'foregroundcolor',...
        [0 0 0],'fontsize',11);     
    if rem(ich,10)==0
        column=column+1;
        row=0;
    end
    row=row+1;
end
%set up finish button
hgui.menusel = uicontrol('Style', 'pushbutton',...
   'String', 'finish',...
   'Position', [posx(end)+75  posy(1) 100 25],...
   'Callback', @finishsel);  

end


function finishsel(hobject,event)
global hsel plotParam parameters csc_map processed hgui
%finished selecting channels button
%now update plotParam.cscNames list global variable
for ii=1:length(hsel.chsel)
    if isvalid(hsel.chsel{ii})
    if hsel.chsel{ii}.Value==1 && ~ismember(hsel.chsel{ii}.String,plotParam.cscNames)
    %not yet in channel database, but now selected, add 
    plotParam.cscNames=[hsel.chsel{ii}.String plotParam.cscNames];
    elseif ismember(hsel.chsel{ii}.String,plotParam.cscNames) && hsel.chsel{ii}.Value==0
    %member of channel database, deselected, remove
    plotParam.cscNames(find(ismember(plotParam.cscNames,hsel.chsel{ii}.String)))=[];
    end
    end
end
paramtemp={};
if strcmp(hgui.subject,'patra')
    [paramtemp,~,plotParam.event_codes]=getparams('patrabipolar','default',plotParam.cscNames,'sessnum',hgui.sessionnum);
else
    [paramtemp,~,plotParam.event_codes]=getparams(hgui.cname,'cleo',plotParam.cscNames);
end
parameters.NCSchannels=paramtemp.NCSchannels;
%getplotsettings(plotParam.filtlfp,plotParam.cscNames,plotParam.event_codes,parameters.sampleratencs,[]);       
[processed.Iread,processed.LFPread,processed.samplesNCS]=...
    loadall(hgui.PathName, hgui.FileName,parameters,plotParam.selch);
parameters.sampleratencs=processed.LFPread.LFPsamplingfreq;
%update plotParam (incase channels changed update)
getplotsettings(plotParam.filtlfp,plotParam.cscNames,...
    plotParam.event_codes,parameters.sampleratencs,[]);   
disp(['loading: ' plotParam.cscNames]);
close(hsel.figsel); %close figure popup, must close first to allow plotting on main fig
%compile/organize data for plotting as loaded above in hgui
if isempty(hgui.ephysid)
compileloaded(hgui);    %update processed, parameters, plotParam globals
end
if ~isempty(processed.samplesNCS)
    if isfield(plotParam,'event_codes')
        processed.behav=calcBehav(processed,plotParam.zoomTS,plotParam.event_codes);        %calculate RT, etc.
    end
   setguicloseup(hgui.closeup{2},'ncsids',plotParam.lfpid,...
        'filt',plotParam.filtlfp,'sqenv',...
        'scale',plotParam.powerscale,'winlength',plotParam.winlength);  
    setguicloseup(hgui.closeup{3},'ncsids',[plotParam.eyeid...
        plotParam.pulseid],'norm',1,'xlab'); 
setguicloseup(hgui.closeup{3},'ncsids',plotParam.lickid,...
    'norm',1,'xlab','env','winlength',plotParam.winlengthphys,'hold');  
        setguicloseup(hgui.closeup{1},'ncsids',plotParam.lfpid,...
            'scale',plotParam.LFPscale);
    if plotParam.buttonm==1
        setspectrum(processed,plotParam,hgui.fftplot,plotParam.mcsc); 
    end
end
updatelegends;
hgui.loadedcsc=plotParam.cscNames;
%set up menu for fft channel to plot
pos=get(hgui.cscmenu,'position');
hgui.cscmenu = uicontrol('Style', 'popup',...
   'String', hgui.loadedcsc,...
   'Position', pos,...
   'Callback', @cscselect);  
disp('updated parameters');
%delete(hsel);

end