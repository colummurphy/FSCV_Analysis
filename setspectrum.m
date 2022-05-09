function hI=setspectrum(processed,plotParam,fighandle,cscch,varargin)
axPlot=fighandle;
subplot(axPlot)
hold(axPlot,'on')
cla(axPlot);
cscale=plotParam.fftclim;
bl=plotParam.blsub;     %subtract baseline;
plotsize=[650 300];        %plot size for close up
plotsize=[plotParam.ncsplotsize(1) 300];        %plot size for close up
argnum = 1;
dalabel='[DA] (nM)';
norm=0;
units=[];
scalebar=[];
labels={};
xlab=0;         %disable xlabel
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'xlab'
            %put xlabel or not
            xlab=1;           
        otherwise
            error('badoption');
    end
    argnum = argnum + 1;
end

legend(axPlot,'off')
yyaxis right
cla(axPlot)        %delete da trace overlaid
set(axPlot,'ytick',[]);
aa=get(axPlot,'ylabel');
    delete(aa);      
    
yyaxis left
%ids=find(ismember(plotParam.NCSchannels,cscch)==1);  %channel ids to be processed and plotted here 
%cscid=find(ismember(data.cscNames,varargin{argnum}));
ids=find(ismember(plotParam.cscNames,cscch)==1);
if isempty(ids)
    error(['channel ' cscch ' not loaded']);
end
position=getpixelposition(axPlot);
TS=processed.LFPread.LFPts;
zoomTS=plotParam.zoomTS;
yscale=plotParam.freqlim;
idEvents=find(ismember(processed.LFPread.LFPeventTTL,plotParam.disp_events)==1);
relTSEvents=processed.LFPread.LFPeventTS(idEvents)-TS(1);
relTS=TS-TS(1);
idxStart=find(relTS<=zoomTS(1)); idxStart=idxStart(end);
idxEnd=find(relTS>=zoomTS(2)); 
if isempty(idxEnd)
    idxEnd=length(relTS);
end
idxEnd=idxEnd(1);
relTS=TS(idxStart:idxEnd)-TS(1);

tdata=processed.samplesNCS(idxStart:idxEnd,ids);
[morletdata, allP]=morletgram(tdata,processed.LFPread.LFPsamplingfreq,[2 150]);     %default log

%lin scale does not work after helen modifieid morletgram code, check with
%dan
%[morletdata, allP]=morletgram(tdata,processed.LFPread.LFPsamplingfreq,[2 150],'lin');   %does not output freq
imagedata = 10*log10(morletdata.P);
%imagedata=morletdata.P;        %linear power
offset=size(tdata,1)-size(imagedata,2);
if size(tdata,2)>size(tdata,1)
offset=size(tdata,2)-size(imagedata,2);
end
offset1=floor(offset/2); offset2=ceil(offset/2);
%pad morlet spectrum at edges where fft was not calculated by morletgram
morletdata.Presize=[repmat(nan,length(morletdata.freqs),offset1) imagedata repmat(nan,length(morletdata.freqs),offset1)];
plotdata=morletdata.Presize;
if bl
    %baseline sub
    avgfft=nanmean(plotdata,2);
    plotdata=plotdata-avgfft;
end
%relt=[relTS(1)+offset/2/processed.LFPread.LFPsamplingfreq relTS(end)-offset/2/processed.LFPread.LFPsamplingfreq];
hI = imagesc(relTS, morletdata.freqs, plotdata,'Parent', axPlot,'cdatamapping','scaled');

%image(Isub_display(:,(t_start:t_end)),'cdatamapping','scaled')

for ii=1:length(relTSEvents)
   if relTSEvents(ii)~=0
        YLim=get(axPlot,'ylim');
        x1=relTSEvents(ii);
        aa=plot(axPlot,[x1 x1],YLim,'color',[1 0 .75],'linestyle','--','marker','none');
        aa.Color(4)=.75;
   end
end

if xlab==1
xlabel(axPlot,'time (s)')
end
xlim=set(axPlot,'xlim',zoomTS);
colormap(axPlot,parula)      %defined color map

set(axPlot, 'visible','on', 'box','off');  %make invisble until prompted

set(axPlot, 'YScale', 'log')
%set(axPlot, 'YScale', 'linear')    lin scale does not wokr

%set(axPlot,'ytick',round(morletdata.freqs(1:10:end)))
set(axPlot,'ytick',round(yscale(1):10:yscale(2)),'color',[0 0 0])
set(axPlot,'ycolor',[0 0 0])

ylabel(axPlot,'frequency (Hz)','color',[0 0 0])
set(axPlot, 'Units','Pixels','Position',  [position(1) position(2) plotsize(1) plotsize(2)]);

ylim=set(axPlot,'ylim',yscale);
if isempty(get(axPlot,'colorbar'))
    h1=colorbar(axPlot,'eastoutside');
    cpos = getpixelposition(h1);
    aa=get(h1,'ylabel');
     cpos(3) = 15; 
    set(h1,'Units','Pixels','Position', [position(1)+plotsize(1)+5  cpos(2) cpos(3) cpos(4)]);
end
%{
if isempty(aa.String)
    delete(aa);            
    y_h=ylabel(h1,'log V^2');  
    ylp = get(y_h, 'Position');
    ext=get(y_h,'Extent');
    %pause(.05);
    set(y_h, 'Rotation',270, 'Position',ylp+[ext(3) 0 0])
    ypos = getpixelposition(y_h);
    cpos(3) = 15; 
set(h1,'Units','Pixels','Position', [cpos(1)  cpos(2) cpos(3) cpos(4)]);

end
%ylabelbar=ylabel(h1,'log V^2','rotation', 270); 
%}
%set(h1,'Units','Pixels','Position', [cpos(1)+50  cpos(2) cpos(3) cpos(4)]);
%set(h1, 'ylim', [cminshow cscale(2)])
if ~isempty(cscale)
    caxis(axPlot,cscale); 
end
if bl
    stdfft=nanstd(plotdata,[],2);
    meanstd=nanmean(stdfft);
    meanfft=nanmean(plotdata,2);
    meanfft2=nanmean(meanfft);
    cmax=meanfft2+meanstd*2;
    cmin=meanfft2-meanstd;
    caxis(axPlot,[cmin cmax]);
end
%colormap(axPlot,'hot');
aa=get(axPlot,'title');
%{
if isempty(aa.String)
    titlelabel=({['csc ' num2str(cscch) ]});
    title1=title(axPlot,titlelabel,'Interpreter','none','color',[0 0 0]);
elseif ~isequal(aa.String{:}, ['csc ' num2str(cscch) ])
titlelabel=({['csc ' num2str(cscch) ]});
title1=title(axPlot,titlelabel,'Interpreter','none','color',[0 0 0]);
end
  %}

hold(axPlot,'off')


end