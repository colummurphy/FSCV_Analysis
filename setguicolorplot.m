function [] = setguicolorplot(axPlot,cData,plotnum)
global plotParam parameters
fontsize=12;
firstplot=plotParam.firstplot;
if ~isempty(cData)
    cmax=1;     %default max color
    fontsize=plotParam.fontsize;   
    events=plotParam.events;
    if isempty(axPlot)
        fign=figure;    set(fign,'position',[100 100 900 500],'color',[1 1 1]);
        axPlot=axes;
    end
    plotsize=[500 150];
    if isfield(plotParam,'colorplotsize')
        plotsize=plotParam.colorplotsize;
    end
    if isfield(plotParam,'cmax')
        cmax=plotParam.cmax;        
    end
    [map, cmin]=setColorScale(cmax);        %default colormap settings
    cScale=[cmin cmax];
    cMap=map;
    if isfield(plotParam,'cscale')
        cScale=plotParam.cscale;
    end
    if isfield(plotParam,'map')
        cMap=plotParam.map;
    end
    sampling_freq=10;
    if isfield(parameters,'samplingFreq')
        sampling_freq=parameters.samplingFreq;
    end
    t_end=plotParam.t_end;
    t_start=plotParam.t_start;
    %ignore plotParam, define data size by t_start t_end before
    %inputting to this script now
    %
    t_start=1;
    t_end=size(cData,2);
  %  subplot(axPlot)
    if firstplot
        colormap(axPlot,cMap)      %defined color map
        plotParam.firstplot=0;
    end
    if t_start<1
        t_start=1;
    end
    if t_end>size(cData,2)
        image(axPlot,cData(:,(t_start:size(cData,2))),'cdatamapping','scaled')
    else
        image(axPlot,cData(:,(t_start:t_end)),'cdatamapping','scaled')
    end
    set(axPlot,'clim', cScale)

    if plotnum==1
        if ~isempty(events)
            events2=events;
            events2(:,2:end)=events(:,2:end)-t_start;
            events3=events2(:,2:end);
            [rowEvent,colEvent]=find(events3>=0);
            if ~isempty(rowEvent)
                YLim=get(axPlot,'ylim');
                for ii=1:size(rowEvent,1)
                    text(axPlot,events3(rowEvent(ii),colEvent(ii)),YLim(1)-10,num2str(events(rowEvent(ii),1)),'color',[0 0 0],'fontsize',fontsize);
                end
            end
        end
    end
    
 time_range=t_end-t_start;
    length_xtick=4;
    ticks=time_range/length_xtick;
    time_scale=t_start:ticks:t_end;
    time_scale=round(round(time_scale)./sampling_freq);
    b1=num2str(time_scale');

    set(axPlot,'xlim',[0 time_range])
    set(axPlot,'XTick',[0:ticks:time_range])
    if plotnum==length(plotParam.selch)
        set(axPlot,'xticklabel',b1)
    else
        set(axPlot,'xticklabel',[])
    end
    if plotnum<=4
    ysize=size(cData,1);
    set(axPlot,'ytick',[1; floor(ysize/2); ysize],'tickdir','out');
    %x1=num2str([-0.4 ;1.3; -0.4]);
    x1=[char(0150) '0.4'; '+1.3'; char(0150) '0.4'];
    set(axPlot,'yticklabel',x1);
    end
    if plotnum>4
        set(axPlot,'yticklabel',[]);
    end
    if plotParam.refresh
            position=getpixelposition(axPlot);
    set(axPlot, 'Units','Pixels','Position',  [position(1) position(2) plotsize(1) plotsize(2)]);
    %axis ij  %square

    if firstplot || plotnum==1
            ylabelV=ylabel('Voltage (V)');
        h=colorbar(axPlot);
        cpos = getpixelposition(h);
        cpos(3) = 10;
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

end
    
    
   % hAxes=gca;
   % hAxes.YRuler.Axle.LineStyle = 'none';      %remove yaxis
    set(axPlot,'box','off')
    set(axPlot,'xtick',[],'xcolor',[1 1 1])
    set(axPlot,'ycolor',[0 0 0]);
   % hAxes.TickLength=[0,0];
    set(findall(axPlot,'-property','FontSize'),'FontSize',fontsize)


end
