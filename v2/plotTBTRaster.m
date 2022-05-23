function ax=plotTBTRaster(TS_plot,tbt_winData,varargin)
%Plot TBT trial by trial data

fontsize=10;
figpos=[50,50,700,800];
numplots=1;
plotname='';
cscale=[];
markers=[];
tbt_evmatTSplot=[];
normscale=1;
axa={};
figsess={};

argnum=1;
while length(varargin)>argnum
    switch varargin{argnum}
        case 'axes'
            argnum=argnum+1;
            axa{1}=varargin{argnum};
        case 'title'
            argnum=argnum+1;
            plotname=varargin{argnum};
        case 'cscale'
            argnum=argnum+1;
            cscale=varargin{argnum};
        case 'markers'
            %Event markers {{name of events}, [timestamps on each column for
            %each event relative to the alignment event at 0s]}
            argnum=argnum+1;
            getevents=varargin{argnum};
            markers=getevents{1};
            tbt_evmatTSplot=getevents{2};
    end
    argnum=argnum+1;
end

if isempty(axa)
    %No figure axes supplie
    [figsess,axa]=setupFig(figpos,numplots);
end


pEvtSize=15;    %Marker size for events
pEvtColors=cool(length(markers));  %Plot markers colors

cla(axa{1});
set(axa{1},'Units','normalized');   %relative position for axes so depends on position of fig set above
set(axa{1},'Position',[ 0.15 0.1 .7 .75]);
set(axa{1},'Fontsize',fontsize);

title(axa{1},plotname,'interpreter','none');
imagetrials=image(axa{1},TS_plot,1:size(tbt_winData,1), tbt_winData,'cdatamapping','scaled');
axcb=colorbar(axa{1},'location','eastoutside','Fontsize',10);
set(axa{1},'YDir','reverse');        %flip y-axis values so first trial on top
artTime=isnan(tbt_winData);   %find artifact points (nan periods)
artTime=abs(artTime-1);         %make alpha data mask by inverting 1/0's
artTime2=artTime;
maskGray=artTime2==0;             %find Zero indices representing artifact mask
maskGray=maskGray*.15;            %make gray rather than white default by making non-zero
artTime=artTime+maskGray;
set(imagetrials, 'AlphaData', artTime);    
set(axa{1},'tickdir','out','box','off')
set(axa{1},'xlim',[min(TS_plot) max(TS_plot)]);
set(axa{1},'ylim',[1 size(tbt_winData,1)]);

%Set colorscale
if ~isempty(cscale)
    caxis(axa{1},cscale)
else
    if normscale
        %normalize to mean / +2 * std and -1 * std, positive bias
        meanSig=nanmean(nanmean(tbt_winData,1));
        stdSig=nanmean(nanstd(tbt_winData,1));
        cscale=[meanSig-1*stdSig meanSig+2*stdSig];
        caxis(axa{1},cscale);
    end
end

%Set markers for all markers
set(axa{1},'clipping','off');   %To display text/scatter outside of axes limits for legend
ylims=max(ylim(axa{1}));    %Y boundary for axis for scaling points at which to put legend
yinc=ylims/8;%Use this range above the axis to plot marker legend 
xstart=min(xlim(axa{1}));
xscale=abs(min(xlim(axa{1}))-max(xlim(axa{1})));
xinc=xscale/4.5;
deltaX=xstart;
for ic=1:length(markers)
    evtid=find(contains(markers,markers{ic}));     %Get numeric code for specified event name
    for i3=1:size(tbt_evmatTSplot,3)
        %scroll through eavh event found within period of trial in 3rd
        %dimension
        scatter(axa{1},tbt_evmatTSplot(:,evtid,i3),1:size(tbt_winData,1),pEvtSize,pEvtColors(ic,:),'MarkerEdgeAlpha',0.7) %Overlay event markers  
    end
        text(axa{1},deltaX,-yinc,markers{ic},'color',pEvtColors(ic,:),'interpreter','none');%Text legend for event markers   
        scatter(axa{1},deltaX-xinc/15,-yinc,pEvtSize*2,pEvtColors(ic,:),'MarkerEdgeAlpha',1)    %Symbol legend event markers    
        deltaX=deltaX+xinc;
end

ax=axa{1};