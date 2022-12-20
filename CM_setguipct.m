function Ipcr = CM_setguipct(axPlot, pct, plotParam,parameters,colorfscv,varargin)
%call getpct to compute PCA components (output pct) based on input params first
%input pct here to plot


%{  
    Input Parameters

    axPlot = dopaminePlot
    pct = processed.Ipcr{selch(ii)} (for this channel)
    plotParam = plotParam
    parameters = parameters
    colorfscv = plotParam.colorFSCV((selch(ii)-4*(ceil(selch(ii)/4)-1)),:)
    varargin = 'plotnum',ii 
%}


argnum=1;
detected=[];
plotnum=1;

while argnum <= length(varargin)
    switch varargin{argnum}       
        case 'detected'
            %output detected timestamps 
            argnum=argnum+1;
            detected=varargin{argnum};
        case 'plotnum'
            %trace #, if more than 1, may not want to repeat event marks
            argnum=argnum+1;
            plotnum=varargin{argnum};
        otherwise
            error('badoption');
    end
    argnum = argnum + 1;
end


Ipcr={};
fontsize=plotParam.fontsize;
events=plotParam.events;        %events already shifted rel t_start
plotsize=plotParam.colorplotsize;
scaleBar=plotParam.scaleBar;
yscale=plotParam.fscvscale;
scaleBarConc=plotParam.scaleBarConc;
sampling_freq=parameters.samplerate;
t_start=plotParam.t_start;
t_end=plotParam.t_end;
time_window=(t_start:t_end);
%%%
time_window=1:size(pct.DAiso,2);
dispOthers=plotParam.displaypH;
badids=[];

if isfield(parameters,'badzone')
    %plot badids
badids=parameters.badzone;
end

%set up plots

%{
CM - Commented Out
  subplot(axPlot)
%}

hold(axPlot,'on');
%set(axPlot,'XTick',[])
%set(axPlot,'XTicklabel',[])

%{
CM - Commented Out
  Position set in App Designer
  position=getpixelposition(axPlot);
%}

set(axPlot,'FontName','Arial','fontsize',fontsize)
time_range=t_end-t_start;
length_xtick=4;
ticks=time_range/length_xtick;
%reltimetick=time_window/40;
reltimewin=1:ticks:time_window(end);
time_scale=t_start:ticks:t_end;
time_scale=round(round(time_scale)./sampling_freq);
b1=num2str(time_scale');
set(axPlot,'xlim',[time_window(1) time_window(end)])
set(axPlot,'XTick',[reltimewin(1):ticks:reltimewin(end)])
set(axPlot,'xticklabel',b1)
xlabel(axPlot, 'time (s)');

%if not only displaying da set up right axis
if sum(dispOthers)>0
    yyaxis('right');
    cla(axPlot)
    yyaxis('left');
end

if sum(pct.pHplot>0) && dispOthers(1)~=0
    yyaxis('right');
    set(axPlot,'ycolor',[ 0 0 0]) 
    plot(axPlot, time_window,pct.pHplot(time_window),'color',colorfscv,'linestyle',':')
    ylabel('pH','rotation',270);
    yyaxis('left');
    hAxes=gca;
    hAxes.YRuler.Axle.LineStyle = 'none';      %remove yaxis
end

if sum(pct.Mplot>0) && dispOthers(3)~=0
    yyaxis('right');
    hold(axPlot,'on')
    set(axPlot,'ycolor',[ 0 0 0]) 
    plot(axPlot, time_window,abs(pct.Mplot(time_window)),'color',colorfscv,'linestyle','--','linewidth',0.5)
    yyaxis('left');
end    

if sum(pct.BGplot>0) && dispOthers(2)~=0
    yyaxis('right');
    hold(axPlot,'on')
    set(axPlot,'ycolor',[ 0 0 0]) 
    plot(axPlot, time_window,pct.BGplot(time_window),'color',colorfscv,'linestyle','-.','linewidth',0.5)
    yyaxis('left');
end
if sum(pct.Iox>0) && dispOthers(4)~=0
    yyaxis('right');
    hold(axPlot,'on')
    set(axPlot,'ycolor',[ 0 0 0]) 
    plot(axPlot, time_window,pct.Iox(time_window),'color',colorfscv,'linestyle','-.','linewidth',0.5)
    yyaxis('left');
end

plot(axPlot, time_window,pct.DAiso(time_window),'color',colorfscv,'linewidth',1,'marker','none')

if plotnum<=4
    ylabel(axPlot, '[DA] (nM)');
end
set(axPlot,'ycolor',[ 0 0 0]); 

%{
% CM - Commented Out
% Position set in App Designer
set(axPlot, 'Units','Pixels','Position',  [position(1) position(2) plotsize(1) plotsize(2)]);
%}

if ~isempty(yscale)
    ylim(axPlot, yscale);
end

if scaleBar~=0 && scaleBarConc~=0
    set(axPlot,'YColor','none');    
    %plot conc scale bar
    plot(axPlot, [time_window(5); time_window(5)], [scaleBarConc+5; 5], '-k',  'LineWidth', 1);
end
YLim=get(axPlot,'ylim');
if ~isempty(events) && plotnum==1
    eventbreaks=find(events(:,1)==6);       %plot fix break (not engaged int ask)
    eventsbre=events(eventbreaks,2:end);
    eventsbre=eventsbre((eventsbre)~=0);
    eventsbig=find(events(:,1)==45);
    eventsbig=events(eventsbig,2:end);
    eventsbig=eventsbig((eventsbig)~=0);
    eventssmall=find(events(:,1)==46);
    eventssmall=events(eventssmall,2:end);
    eventssmall=eventssmall((eventssmall)~=0);
    
    %color coded
    if ~isempty(eventsbre)
        for ii=1:size(eventsbre,2)
            x1=eventsbre(ii);
            plot(axPlot, [x1 x1],YLim,'color',[1 0 0],'linestyle','-','linewidth',...
                2,'marker','none');
        end
    end
    if ~isempty(eventsbig)
        for ii=1:size(eventsbig,2)
            x1=eventsbig(ii);
            plot(axPlot, [x1 x1],YLim,'color',[0 .7 .4],'linestyle','--','marker','none');
        end
    end
    if ~isempty(eventssmall)
        for ii=1:size(eventssmall,2)
            x1=eventssmall(ii);
            plot(axPlot, [x1 x1],YLim,'color',[0 0 .5],'linestyle','--','marker','none');
        end
    end
    %plot targeted event markers only for first plot
    events2=reshape(events(:,2:end),[],1);
    for ii=1:size(events2,1)
    %txt1=num2str(events2(ii,1));
        if events2(ii)~=0
            x1=events2(ii);
            if ~ismember(x1,[eventsbre eventssmall eventsbig])
            plot(axPlot, [x1 x1],YLim,'color',[.6 .6 .6],'linestyle','--','marker','none');
            end
        end
    end

end

if ~isempty(detected)
    %plot detected da marks
    %these are already shifted in time, just need to get in window
    dinwin=detected(detected>=1 & detected<=plotParam.t_end-plotParam.t_start+1);
    for ii=1:length(dinwin)
        text(dinwin(ii),YLim(2),'*','color',colorfscv,'fontsize',12,...
            'horizontalalignment','center');
    end        
end
if plotParam.plotbad && plotnum==1
    %plot bad ids
    badinwin=badids(badids>=plotParam.t_start & badids<=plotParam.t_end);
    for ii=1:length(badinwin)
        %{
    pp=badids(ii)./(plotParam.t_end-plotParam.t_start)*P(3);
    xp=pp+P(1);
    yp=P(4);
    yp=get(hgui.itplot,'ylim');
    hgui.bad{ii}=text(hgui.itplot,badids(ii)-plotParam.t_start,yp(2),'x','fontsize',12,'horizontalalignment',...
        'center','color',[1 0 0]);
        %}
    %text(badinwin(ii)-plotParam.t_start,YLim(2),'x','fontsize',10,'horizontalalignment',...
     %   'center','color',[1 0 0]);
    plot(axPlot, [badinwin(ii)-plotParam.t_start badinwin(ii)-plotParam.t_start],YLim,'color',[.6 .6 .6],'linestyle',':','marker','none');
    end
end
hold(axPlot,'off')

end
