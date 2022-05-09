function updatelegends
%update text legends, called by setfscvgui and setfscvguichannelselection
global plotParam hgui

hposcloseup{1}=getpixelposition(hgui.closeup{1});
txtpad=hgui.txtpad;
figpos=get(hgui.hf,'position');
cscnames=plotParam.cscNames;
lfpid=plotParam.lfpid;
physid=plotParam.physid;
colortab=plotParam.colormaptab;
%delete current labels
if isfield(hgui,'labelslfp')
for ii=1:length(hgui.labelslfp)
    delete(hgui.labelslfp{ii});
    %hgui.labelslfp{ii}.String='';
end
end
for ii=1:length(lfpid)
hgui.labelslfp{ii}= uicontrol('Style', 'text',...
   'String', cscnames{lfpid(ii)},'fontangle','italic',...
   'Position', [figpos(3)-60 ...
   hposcloseup{1}(2)-120-txtpad*(ii-1) 60 25],...
   'BackgroundColor', [ 1 1 1],'foregroundcolor',colortab(ii,:));
end
for ii=1:length(physid)
hgui.labelslfp{ii}= uicontrol('Style', 'text',...
   'String', cscnames{physid(ii)},'fontangle','italic',...
   'Position', [figpos(3)-60 ...
   hposcloseup{1}(2)-320-txtpad*(ii-1) 60 25],...
   'BackgroundColor', [ 1 1 1],'foregroundcolor',colortab(ii,:));
end



end