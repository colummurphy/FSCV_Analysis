function []=compileloaded(hgui)
%organize/prepare all loaded data from directory for plotting/analysis
%also load events associated with file
%hgui is struct that holds indidivudal figure/axes 
%called from 
global plotParam processed parameters
selch=plotParam.selch;
sizeData=size(processed.Iread(selch(1)).data); 
parameters.samplesperscan=sizeData(1); lengthData=sizeData(2);
if lengthData<plotParam.t_end
    plotParam.t_end=lengthData;
end
disp(['file length = ' num2str(lengthData./parameters.samplerate) ' s'])
plotParam.events=[];
if size(processed.Iread(selch(1)).events,1)>0 ...
    && ~isempty(processed.Iread(selch(1)).events) && sum(sum(processed.Iread(selch(1)).events~=0)) && ~contains(hgui.subject,'cfmea')
    %Add processed.Iread(selch(1)).events~=0 above 07/01/2021
%get events that occured in selected time period
%Add sum(processed.Iread(selch(1)).events~=0) since an matrix for all channels
    [nlx_events,plotParam.events]=readEvents(...
        processed.Iread(selch(1)).events(:,plotParam.t_start:plotParam.t_end),...
        plotParam.disp_events);
end

%%
%initialize/filter loaded data & plot
for ii=1:length(selch)
    %filter all fscv data
    processed.Iread(selch(ii)).data = filterDisplay([],[],processed.Iread(selch(ii)).rawdata);
    %prepare process fscv data in time window selected
    [processed.Isub(selch(ii)).data, processed.BG(selch(ii)).data]=...
        refreshfscv(processed.Iread(selch(ii)).data(:,plotParam.t_start:plotParam.t_end),...
        plotParam); %based on updated xsel & tstart/tend range
end
%[processed.Isub, processed.BG]=refreshfscv(processed.Iread,plotParam); %based on updated xsel
cla(hgui.itplot)
if length(selch)>4
    cla(hgui.itplotx{1})
     cla(hgui.itplotx{2})
      cla(hgui.itplotx{3})
end
for ii=1:length(selch)
    setguicolorplot(hgui.hfig{selch(ii)},...
            processed.Isub(selch(ii)).data,ii); 
    %compute pca
    hfwidth=[];
    if isfield(parameters,'hfwidth')
        hfwidth=parameters.hfwidth;
    end
     processed.Ipcr{selch(ii)} = ...
            getpct(processed.Iread(selch(ii)).rawdata(:,plotParam.t_start:plotParam.t_end),...
            processed.Isub(selch(ii)).data,...
            processed.BG(selch(ii)).data,parameters,selch(ii),...
            'removebgph','nanwidth',8,'glitchwidth',hfwidth);
    %plot pca computed concentrations   
        plotax=hgui.itplot;
        if selch(ii)<=8 && selch(ii)>4  
            plotax=hgui.itplotx{1};
        elseif selch(ii)<=12 && selch(ii)>8
            plotax=hgui.itplotx{2};
        elseif selch(ii)<=16 && selch(ii)>12
            plotax=hgui.itplotx{3};
        end
            setguipct(plotax,processed.Ipcr{selch(ii)},plotParam, parameters,...
                plotParam.colorFSCV((selch(ii)-4*(ceil(selch(ii)/4)-1)),:),'plotnum',ii); 
    %{
    if selch(ii)<=8 && selch(ii)>4        
    setguipct(hgui.itplotx{1},processed.Ipcr{selch(ii)},plotParam, parameters,...
        plotParam.colorFSCV((ii-4*(ceil(ii/4)-1)),:),'plotnum',ii);  
    end
    if selch(ii)<=12 && selch(ii)>8
    setguipct(hgui.itplotx{2},processed.Ipcr{selch(ii)},plotParam, parameters,...
        plotParam.colorFSCV((ii-4*(ceil(ii/4)-1)),:),'plotnum',ii);  
    end
    if selch(ii)<=16 && selch(ii)>12
    setguipct(hgui.itplotx{3},processed.Ipcr{selch(ii)},plotParam, parameters,...
        plotParam.colorFSCV((ii-4*(ceil(ii/4)-1)),:),'plotnum',ii);  
    end
        %}
end


end