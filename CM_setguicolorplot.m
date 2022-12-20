function [] = CM_setguicolorplot(axPlot,cData,plotnum)
global plotParam parameters


%{
  Input Parameters

  axPlot = hgui.hfig{chX}               Color plot handle
  cData = processed.Isub(chX).data      Channel data after filtering
                                        and background subtraction
  plotnum = chX                         Channel number 
  Note: chX = 1 or 2 or 3 or 4.  

  Function is run for each channel (1 to 4);   



  plotParam.map (233 * 3)
  Create a custom colormap by defining a three-column 
  matrix of values between 0.0 and 1.0. 
  Each row defines a three-element RGB triplet.

  Color map (plotParam.map) is assigned in getplotsettings.m
  Calculated in setColorScale.m
  [plotParam.map, plotParam.cmin]=setColorScale(plotParam.cmax);

%}

% default font size
fontsize=12;

% firstplot = 1
firstplot=plotParam.firstplot;

% if channel data is available
if ~isempty(cData)
    cmax=1;     %default max color
    
    % fontsize = 8
    fontsize=plotParam.fontsize;   
    
    % events = events table
    % First column contains the eventsToDisplay (vertical columns).
    % Remaining columns contain the time indexes (1-601) of those events
    events=plotParam.events;
    
    % False, if there is no figure create a figure
    if isempty(axPlot)
        fign=figure;    
        set(fign,'position',[100 100 900 500],'color',[1 1 1]);
        axPlot=axes;
    end
    
    % Set plotsize = [750, 150]
    plotsize=[500 150];
    if isfield(plotParam,'colorplotsize')
        plotsize=plotParam.colorplotsize;
    end
    
    % Set color max = 1
    if isfield(plotParam,'cmax')
        cmax=plotParam.cmax;        
    end
    
    % default colormap settings, cmin not used, 
    % cScale + cMap overritten below
    [map, cmin]=setColorScale(cmax);        
    cScale=[cmin cmax];
    cMap=map;
    
    % Set color scale = [-0.6667, 1]
    if isfield(plotParam,'cscale')
        cScale=plotParam.cscale;
    end
    
    % Set color map = 233 * 3 vector
    if isfield(plotParam,'map')
        cMap=plotParam.map;
    end
    


    % default sampling_freq = 10
    sampling_freq=10;
    
    % parameters has a 'samplingrate' = 10, 
    if isfield(parameters,'samplingFreq')
        sampling_freq=parameters.samplingFreq;
    end
    
    % t_start = 1, t_end = 601, overritten below
    t_end=plotParam.t_end;
    t_start=plotParam.t_start;
    
    %ignore plotParam, define data size by t_start t_end before
    %inputting to this script now
    
    % t_start = 1, t_end = 601 (no change)
    t_start=1;
    t_end=size(cData,2);


    % axPlot - handle to our current plot
    % Set the colormap for the axPlot figure to cMap
    if firstplot
        colormap(axPlot,cMap)      %defined color map
        
        % toggle plotParam.firstplot = 0 after it is run with the first
        % channel
        plotParam.firstplot=0;
    end
    
    % 1 < 1, False, t_start=1 (no change)
    if t_start<1
        t_start=1;
    end


    % 601 > 601, False
    if t_end>size(cData,2)
        image(axPlot,cData(:,(t_start:size(cData,2))),'cdatamapping','scaled')

    % True, run this one    
    else
        % Displays the data in cData as an image.
        % Each element of cData specifies the color for 1 pixel of the image.
        % image is an m-by-n grid of pixels,(m = num rows, n = num cols of cData)
        % 
        % Scale the values to the full range of the current colormap
        % i.e. Use the full range of colors
        
        % image(fig/subplot, dataMatrix, 'name', 'value')
        image(axPlot,cData(:,(t_start:t_end)),'cdatamapping','scaled')
    end
    
    
    % Set the color limits (clim) for the plot 
    % cScale = [-0.6667, 1]
    % ax.CLim(1) - color bar bottom value. 
    % ax.Clim(2) - color bar top value.
    % 'clim' automatically computed from the range of the data being plotted.
    % Setting 'clim' changes the way the color is scaled from the data values.
    set(axPlot,'clim', cScale)


    % True for channel 1
    % plot the events onto the color plot
    if plotnum==1
        
        % if events table is not empty, True
        if ~isempty(events)
            % Copy of events
            % First column contains the eventsToDisplay (vertical column).
            % Remaining columns contain the time indexes (1-601) of those
            % events.
            events2=events;
            
            % Subtract 1 from the time indices of the events
            events2(:,2:end)=events(:,2:end)-t_start;
            
            % Event time indices, with 1 subtracted, 
            % eventsToDisplay column removed.
            events3=events2(:,2:end);
            
            % row + column indices of non-zero values
            % searches col 1 first then, col 2, ...
            [rowEvent,colEvent]=find(events3 >= 0);
            

            % if there are events(events are not empty)
            if ~isempty(rowEvent)
                
                % YLim = ?
                YLim=get(axPlot,'ylim');
                
                % for each event                
                for ii=1:size(rowEvent,1)
                    text(axPlot, ...
                        events3(rowEvent(ii), colEvent(ii)), ...% x position (event time index)
                        YLim(1)-10, ...                         % y position 
                        num2str(events(rowEvent(ii),1)), ...    % dat from eventsToDisplay column 
                        'color', [0 0 0], ...
                        'fontsize', fontsize);
                end
            end
        end
    end
 
    % time_range = 601 - 1 = 600
    time_range=t_end-t_start;
    
    % ticks = 150
    length_xtick=4;
    ticks=time_range/length_xtick;
    
    % time_scale = 1:150:601 = [1, 151, 301, 451, 601]
    time_scale=t_start:ticks:t_end;

    % time_scale = [0, 15, 30, 45, 60]
    time_scale=round(round(time_scale)./sampling_freq);
    
    % b1 = ['0', '15', '30', '45', '60'] (column array)
    b1=num2str(time_scale');


    set(axPlot,'xlim',[0 time_range])
    set(axPlot,'XTick',[0:ticks:time_range])
     
    % if this is the last color plot - show the X labels
    if plotnum==length(plotParam.selch)
        set(axPlot,'xticklabel',b1)
    else
        set(axPlot,'xticklabel',[])
    end

    
    % color plots 1, 2, 3, 4
    % Add yticks + yticklabels
    if plotnum<=4
        % Add yticks to color plots
        ysize=size(cData,1);
        set(axPlot,'ytick',[1; floor(ysize/2); ysize],'tickdir','out');
        
        % Add y labels to color plots
        %x1=[char(0150) '0.4'; '+1.3'; char(0150) '0.4'];
        x1=['-0.4'; '+1.3'; '-0.4'];
        set(axPlot,'yticklabel',x1);
    end

    % plotnum >= 5, then no yticklabels
    if plotnum>4
        set(axPlot,'yticklabel',[]);
    end

    %{ 
    CM - Commented Out    
    
    if plotParam.refresh
        position=getpixelposition(axPlot);
        set(axPlot, 'Units','Pixels','Position',  [position(1) position(2) plotsize(1) plotsize(2)]);
        %axis ij  %square

        % if its either channel 1 or the first channel to be plotted.
        if firstplot || plotnum==1
            ylabelV=ylabel('Voltage (V)');
            
            % set h as a colorbar attached to the color plot
            h=colorbar(axPlot);
            
            cpos = getpixelposition(h);
            cpos(3) = 10;
            % position the color bar
            set(h,'Units','Pixels','Position', [position(1)+plotsize(1)+10  position(2) cpos(3) plotsize(2)]);
            
            
            
            clabel2=ylabel(h,'Current (nA)','rotation',270,'fontsize',fontsize);
            set(axPlot, 'Units','Pixels','Position',  [position(1) position(2) plotsize(1) plotsize(2)]); 
            
            posclabel=get(clabel2,'position');
            set(clabel2,'Position',[posclabel(1)+1 posclabel(2)]);
            set(findall(h,'-property','FontSize'),'FontSize',fontsize)
        end

        plotParam.refresh=0;        %reset flag
        % end
        set(axPlot, 'Units','Pixels','Position',  [position(1) position(2) plotsize(1) plotsize(2)]);
    end
    %}

end
    
    
   % hAxes=gca;
   % hAxes.YRuler.Axle.LineStyle = 'none';      %remove yaxis
    set(axPlot,'box','off')
    set(axPlot,'xtick',[],'xcolor',[1 1 1])
    set(axPlot,'ycolor',[0 0 0]);
   % hAxes.TickLength=[0,0];
    set(findall(axPlot,'-property','FontSize'),'FontSize',fontsize)


end
