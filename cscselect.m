function cscselect(source,event)
%select channel for morlet fft from dropdown menu created in setfscvgui
global hgui
val = source.Value;     %selected # setting on drop down list
listoptions = source.String;   %all strings on drop down list    
% For R2014a and earlier: 
% val = get(source,'Value');
% maps = get(source,'String'); 
disp(listoptions{val});
if val>0
hgui.mcsc=hgui.loadedcsc{val};
pos=get(hgui.dispch,'position');
hgui.dispch = uicontrol('Style', 'text',...
   'String', hgui.mcsc,...
   'Position', pos,...
   'BackgroundColor', [ 1 1 1]);
end
end
