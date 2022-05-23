function [figsess,axa]=setupFig(figpos,numplots)

figsess=figure('position',figpos,'color',[1 1 1]);
set(0,'CurrentFigure',figsess);    %set figure handle to current figure
axa={};
for ip=1:numplots       
    axa{ip}=subplot(1,numplots,ip);
    hold(axa{ip},'on');
    set(axa{ip},'units','pixels');
    axpos{ip}=get(axa{ip},'position');
end