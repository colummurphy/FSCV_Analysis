function [] = CM_setguicv(axPlot,cvdata,parameters,colorCV)

%{
  Comments here
%} 

% if color is not set use black 
if isempty(colorCV)
    colorCV=[0 0 0];        % black
end

% get variables from parameters
Vrange_anodal=parameters.Vrange;
Vrange_cathodal=parameters.Vrange_cathodal;

% plotting twice
hold(axPlot,'on');

% axis square
axPlot.PlotBoxAspectRatio = [1 1 1];
axPlot.PlotBoxAspectRatioMode = 'manual';

% split the data in two
Irange_anodal=cvdata(1:length(Vrange_anodal));
Irange_cathodal=cvdata(length(Vrange_anodal)+1:end);

% plot data
plot(axPlot, Vrange_anodal,Irange_anodal,'color',colorCV)
plot(axPlot, Vrange_cathodal,Irange_cathodal,'color',colorCV)

% change to dynamic
% xlim([min(Vrange_anodal) max(Vrange_anodal)])
axPlot.XLim = [min(Vrange_anodal) max(Vrange_anodal)];

% plotting completed
hold(axPlot,'off')

end
