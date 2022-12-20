function [] = CM_compileloaded(hgui, appDesignColorFig, dopaminePlot)

%organize/prepare all loaded data from directory for plotting/analysis
%also load events associated with file
%hgui is struct that holds indidivudal figure/axes 
%called from 

global plotParam processed parameters

% selch = [1,2,3,4]
selch=plotParam.selch;

% sizeData = [175, 601]
sizeData=size(processed.Iread(selch(1)).data); 

% parameters.samplesperscan = 175
parameters.samplesperscan=sizeData(1); 

% lengthData = 601
lengthData=sizeData(2);

% 601 < 1001,   plotParam.t_end = 601 
if lengthData<plotParam.t_end
    plotParam.t_end=lengthData;
end

% display file length in seconds (60.1 s)
disp(['file length = ' num2str(lengthData./parameters.samplerate) ' s'])


% initialize events
plotParam.events=[];

% True
% 12 > 0 && not empty && sum(sum of col) > 0 && subject isn't cfmea
% .events rows = 12, .events is not empty, sum(sum cols)>0, hgui.subject = patra 
if size(processed.Iread(selch(1)).events,1)>0 ...
    && ~isempty(processed.Iread(selch(1)).events) ... 
    && sum(sum(processed.Iread(selch(1)).events~=0)) ...
    && ~contains(hgui.subject,'cfmea')

    %Add processed.Iread(selch(1)).events~=0 above 07/01/2021
    %get events that occured in selected time period
    %Add sum(processed.Iread(selch(1)).events~=0) since an matrix for all channels
    
    % nlx_events not used
    % plotparam.events
    % First column contains the eventsToDisplay (vertical).
    % Remaining columns contain the time indexes (1-601) of those events
    [nlx_events,plotParam.events]=readEvents(...
        processed.Iread(selch(1)).events(:,plotParam.t_start:plotParam.t_end),...
        plotParam.disp_events);
end


%%
%initialize/filter loaded data & plot
% for each 1:4
% for each channel, filter the channel, then subtract the background
for ii=1:length(selch)
    %filter all fscv data
    
    % processed.Iread(ii).data = (ii).rawdata after being passed through a Gauss filter.
    processed.Iread(selch(ii)).data = filterDisplay([],[],processed.Iread(selch(ii)).rawdata);
     
    % processed.Isub(chX).data (175 * 601)
    % Resulting matrix after subtracting the background matrix 
    % from the gauss filtered channel data.

    % processed.BG(chX).data (175 * 1)
    % The mean of each row of a vertical window of filtered channel data. 
    [processed.Isub(selch(ii)).data, processed.BG(selch(ii)).data]=...
        refreshfscv(processed.Iread(selch(ii)).data(:,plotParam.t_start:plotParam.t_end),...
        plotParam); %based on updated xsel & tstart/tend range

end


% clear the dopamine plot
cla(dopaminePlot);

%CM - Comment out
%{
% clear the plots
if length(selch)>4
    cla(hgui.itplotx{1})
     cla(hgui.itplotx{2})
      cla(hgui.itplotx{3})
end
%}

  % for each channel 1:4
  for ii=1:length(selch)

    CM_setguicolorplot( appDesignColorFig(selch(ii)),...
            processed.Isub(selch(ii)).data, ii ); 
    
    %compute pca - hfwidth = []
    hfwidth=[];
    if isfield(parameters,'hfwidth')
        hfwidth=parameters.hfwidth;
    end
    

    % plot pca computed concentrations  

    processed.Ipcr{selch(ii)} = ...
            getpct(processed.Iread(selch(ii)).rawdata(:,plotParam.t_start:plotParam.t_end),...
            processed.Isub(selch(ii)).data,...
            processed.BG(selch(ii)).data,parameters,selch(ii),...
            'removebgph','nanwidth',8,'glitchwidth',hfwidth);
    
        % CM - Comment Out
        %{
        plotax=hgui.itplot;
        
        if selch(ii)<=8 && selch(ii)>4  
            plotax=hgui.itplotx{1};
        elseif selch(ii)<=12 && selch(ii)>8
            plotax=hgui.itplotx{2};
        elseif selch(ii)<=16 && selch(ii)>12
            plotax=hgui.itplotx{3};
        end
        %}
        
     CM_setguipct(dopaminePlot, processed.Ipcr{selch(ii)},plotParam, parameters,...
                plotParam.colorFSCV((selch(ii)-4*(ceil(selch(ii)/4)-1)),:),'plotnum',ii); 
    
  end
end