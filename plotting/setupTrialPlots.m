function hfig=setupTrialPlots(fh,plotparam,numcolor,numbehav)
set(0,'CurrentFigure',fh);    %set figure handle to current figure
numtrials=plotparam.numtrials;
hfig={};
hpos={};
%figpos=get(fh,'position');
figpos=plotparam.figpos;
set(fh,'position',figpos);

for ifig=1:numcolor+numbehav
    hfig{ifig}=subplot(1,numcolor+numbehav,ifig); 
    hpos{ifig}=getpixelposition(hfig{ifig});
end
%{
if numtrials<50
figpos=[200 200 1800 600]; plotparam.colorsize(2)=400;
if ~ispc
    figpos=[0 0 1000 600]; plotparam.colorsize(2)=400;
end
end
if numtrials<10
figpos=[200 200 1800 300]; plotparam.colorsize(2)=200;
if ~ispc
    figpos=[0 0 1000 300]; plotparam.colorsize(2)=200;
end
end
%}


set(fh, 'position',figpos,'color',[1 1 1])

    %resize color plots stretch
    margins=25;
    margins=plotparam.margins;

    xoff=70;
    yoff=30;

    firstbehavwidth=plotparam.vertplotwidth;
if ~ispc
margins=10;
xoff=30;
yoff=15;
firstbehavwidth=30;
end    
    totalcolorpos=xoff+margins*numcolor+(plotparam.colorsize(1))*numcolor;
    countb=0;
    for ifig=1:numcolor+numbehav
        set(hfig{ifig},'Units','Pixel','Position',...
            [xoff+margins*(ifig-1)+(plotparam.colorsize(1))*(ifig-1) 100+yoff ...
            plotparam.colorsize(1) plotparam.colorsize(2)]);
        if ifig==numcolor+1
            %first behav plot is time reaction time
            set(hfig{ifig},'Units','Pixel','Position',...
            [totalcolorpos 100-plotparam.colorsize(2)/numtrials/3+yoff ...
            firstbehavwidth plotparam.colorsize(2)]);      
        end
        if ifig>numcolor+1
            set(hfig{ifig},'Units','Pixel','Position',...
            [totalcolorpos+firstbehavwidth+...
            plotparam.vertplotwidth2*countb 100-plotparam.colorsize(2)/numtrials/3+yoff  ...
            plotparam.vertplotwidth2 plotparam.colorsize(2)]);
            countb=countb+1;
        end
    end


end