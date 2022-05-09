function togglem(source,event)
global plotParam processed hgui parameters
pressed = source.Value;     %checked?  
badids=[];
selch=plotParam.selch;
if pressed==1 
    if isfield(processed,'badids')
        badids=processed.badids;
    end
    if isfield(hgui,'badids')
        badids=unique([badids hgui.badids]);
    end
    if isempty(badids)
        %if still empty bad ids
        return
    end
  %  badids=processed.badids;
 % badids=hgui.badids;
    %P=get(hgui.itplot,'position');
    if isfield(hgui,'bad')
        if ~isempty(hgui.bad)
            for ii=1:length(hgui.bad)
                delete(hgui.bad{ii})
            end
            hgui=rmfield(hgui,'bad');
        end
    end
    plotParam.plotbad=1;        %flag to plot bad ids
    parameters.badzone=badids;      %temp store for bad ids for plotting
    %plot pca computed concentrations
    for ii=1:length(selch)
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
    end
    %rmfield(hgui,'bad');    %clear existing variable/text first
    %{
        for ii=1:length(badids)
            pp=badids(ii)./(plotParam.t_end-plotParam.t_start)*P(3);
            xp=pp+P(1);
            yp=P(4);
            yp=get(hgui.itplot,'ylim');
            hgui.bad{ii}=text(hgui.itplot,badids(ii)-plotParam.t_start,yp(2),'x','fontsize',12,'horizontalalignment',...
                'center','color',[1 0 0]);
        end
    %}
    else
    %not pressed
        if isfield(hgui,'bad')
        if ~isempty(hgui.bad)
            %for ii=1:length(hgui.bad)
            %delete(hgui.bad{ii})
           % end
            hgui=rmfield(hgui,'bad');
             %rmfield(hgui,'bad');    %clear existing variable/text first
        end
        end
        plotParam.plotbad=0;
        for ii=1:length(selch)
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
        end
end

end
