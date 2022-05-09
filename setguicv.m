function [] = setguicv(axPlot,cvdata,parameters,colorCV)
if isempty(colorCV)
    colorCV=[0 0 0];
end
Vrange_anodal=parameters.Vrange;
Vrange_cathodal=parameters.Vrange_cathodal;
subplot(axPlot)
hold(axPlot,'on');
plotsize=[200 150];
position=getpixelposition(axPlot);
set(axPlot, 'Units','Pixels','Position',  [position(1) position(2) plotsize(1) plotsize(2)]);
axis  square
Irange_anodal=cvdata(1:length(Vrange_anodal));
Irange_cathodal=cvdata(length(Vrange_anodal)+1:end);
plot(Vrange_anodal,Irange_anodal,'color',colorCV)
plot(Vrange_cathodal,Irange_cathodal,'color',colorCV)
xlim([min(Vrange_anodal) max(Vrange_anodal)])
set(axPlot,'XTick',[-0.4 -0.2 0 0.2 0.4 0.6 0.8 1 1.3]);
set(axPlot,'FontName','Arial','fontsize',8)
xlabel('E_a (V)','fontsize',8)
ylabel('I (nA)','fontsize',8)
nn=title('CV at selected time point','fontsize',8);
%nnpos=get(nn, 'Position');
hold(axPlot,'off')

end
