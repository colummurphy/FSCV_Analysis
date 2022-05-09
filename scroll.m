function scroll(hObject,~)
global hgui processed plotParam parameters
%scroll time window left/right, set up bottons in setfscvgui.m
    twin=[sscanf(get(hgui.tstart, 'string'), '%f') ...
        sscanf(get(hgui.tend, 'string'), '%f')];
    if isvalid(hgui.zoomts(1))
        %not deleted uicontrol
    tzoom=[sscanf(get(hgui.zoomts(1), 'string'), '%f') ...
        sscanf(get(hgui.zoomts(2), 'string'), '%f')];
    end
    winwidth=twin(2)-twin(1);
    overlap=round(.15*winwidth);
    %backward
    twin2=twin-winwidth+overlap;
    if strcmp(hObject.String,'>>')
        %forward        
        twin2=twin+winwidth-overlap;   
        if twin2(2)>plotParam.timeWin(2)
            twin2(2)=plotParam.timeWin(2);
            if twin2(1)>=twin2(2)
                twin2(1)=twin(1);
            end
        end
        disp('forward scroll');
    else
        if twin2(1)<0
            twin2(1)=0;
            if twin2(2)<=0
                twin2(2)=twin(2);
            end
        end
        disp('backward scroll');
    end
    set(hgui.tstart,'string',num2str(twin2(1)));
    
    set(hgui.tend,'string',num2str(twin2(2)));
    if isvalid(hgui.zoomts(1))
    set(hgui.zoomts(1),'string',num2str(twin2(1)));
    set(hgui.zoomts(2),'string',num2str(twin2(2)));    
    end
    refreshbutton;


end

