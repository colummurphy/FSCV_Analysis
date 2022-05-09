function []= mouseclick(hObject,~)
global plotParam processed hgui parameters
%functions to run with mouse clicked on figure
hf = get(hObject,'parent');
b = get(hObject,'selectiontype');
cp = get(hObject,'CurrentPoint');
button=get(hObject,'CurrentKey');  
P=get(gca,'position');
%P = varargin{3}; 
%cp = get(varargin{1},'currentpoint');
if strcmpi(b,'normal')
    clktype=0;
else
    clktype=1;
end

xsel0=cp(1)-P(1);       %position on plot dimensions
xsel1=xsel0/P(3)*(plotParam.t_end-plotParam.t_start);   %relative position mapped to data domain
xsel=round(xsel1+plotParam.t_start);
ysel=cp(2)-P(2);
disp(['x: ' num2str(xsel) ' | y: ' num2str(ysel) ' | clk: ' num2str(clktype)])
xsel=round(xsel1);
plotParam.xsel=xsel;
plotParam.ysel=ysel;
plotParam.clktype=clktype;
%tstart=plotParam.t_start+1;
%tend=plotParam.t_end;

selch=plotParam.selch;
if clktype==0
    %left click, background subtract, update plots & output just
    %selected time period
    %processed.Iread=Iread;
    for ii=1:length(selch)
        [processed.Isub(selch(ii)).data, processed.BG(selch(ii)).data]=...
            refreshfscv(processed.Iread(selch(ii)).data(:,plotParam.t_start:plotParam.t_end),...
            plotParam); %based on updated xsel & tstart/tend range
        setguicolorplot(hgui.hfig{selch(ii)},...
            processed.Isub(selch(ii)).data,ii);      %set colorplots
    end
else
    %right click, plot cv, pca-da vs t, metrics 
    cla(hgui.cv)
    cla(hgui.itplot)
    cla(hgui.data);
if length(selch)>4
    cla(hgui.itplotx{1})
     cla(hgui.itplotx{2})
      cla(hgui.itplotx{3})
end
    %only checked color plots
    for ii=1:length(selch)
        %check if already detected signals for plotting
        detected=[];
        if isfield(processed,'detected')    
            if ~isempty(processed.detected)
                if ii<=length(processed.detected)
            if ~isempty(processed.detected{selch(ii)})
                if isfield(processed.detected{selch(ii)},'maxTS')
                detected=round((processed.detected{selch(ii)}.maxTS-...
                processed.LFPread.LFPts(1)).*parameters.samplerate...
                -plotParam.t_start+2);
                end
            end
                end
            end
        end
        processed.cv{selch(ii)}=getcv(...
            processed.Isub(selch(ii)).data,...
            parameters,plotParam);
        %plot cv
        if ~plotParam.buttonm
            %if not plotting morletgram
            if selch(ii)<=4
        setguicv(hgui.cv,processed.cv{selch(ii)},parameters,plotParam.colorFSCV(selch(ii),:));
            end
        end
        %compute pca    
        hfwidth=[];
        %processed.Ipcr{selch(ii)}=[];
            if isfield(parameters,'hfwidth')
                hfwidth=parameters.hfwidth;
            end
        processed.Ipcr{selch(ii)} = ...
            getpct(processed.Iread(selch(ii)).rawdata(:,plotParam.t_start:plotParam.t_end),...
            processed.Isub(selch(ii)).data,...
            processed.BG(selch(ii)).data,parameters,selch(ii),...
            'removebgph','nanwidth',8,'glitchwidth',hfwidth);
        %plot pca computed concentrations of indicated window
        plotax=hgui.itplot;
        if selch(ii)<=8 && selch(ii)>4  
            plotax=hgui.itplotx{1};
        elseif selch(ii)<=12 && selch(ii)>8
            plotax=hgui.itplotx{2};
        elseif selch(ii)<=16 && selch(ii)>12
            plotax=hgui.itplotx{3};
        end
        %setguipct(plotax,processed.Ipcr{selch(ii)},plotParam, parameters,...
        %    plotParam.colorFSCV((ii-4*(ceil(ii/4)-1)),:),'plotnum',ii); 
        if isfield(processed,'detected')
            if ~isempty(processed.detected{selch(ii)}.maxTS)
                setguipct(plotax,processed.Ipcr{selch(ii)},plotParam, ...
                parameters,plotParam.colorFSCV((selch(ii)-4*(ceil(selch(ii)/4)-1)),:),...
                'detected',round((processed.detected{selch(ii)}.maxTS-...
                processed.LFPread.LFPts(1)).*parameters.samplerate...
                -plotParam.t_start+2),'plotnum',ii);            
            else
            setguipct(plotax,processed.Ipcr{selch(ii)},plotParam, parameters,...
                plotParam.colorFSCV((selch(ii)-4*(ceil(selch(ii)/4)-1)),:),'plotnum',ii); 
            end
        else
            setguipct(plotax,processed.Ipcr{selch(ii)},plotParam, parameters,...
                plotParam.colorFSCV((selch(ii)-4*(ceil(selch(ii)/4)-1)),:),'plotnum',ii); 
        end
            
    end
    calcdataf=1;
    if isfield(parameters,'calcdata')
        calcdataf=parameters.calcdata;
    end
    if calcdataf
    processed.info=calcdata(processed.Iread, processed.Isub, processed.cv,selch);   
    if ~plotParam.buttonm
                    %if not plotting morletgram
        setguidata(hgui.data, processed.info, [parameters.Vrange parameters.Vrange_cathodal])      
    end
    end
end


end