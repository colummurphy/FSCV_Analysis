function setguicloseup(ax,varargin)
global processed plotParam parameters
%varargin: 'units' (string array), scalebar [min max] lfp, norm (ie. zscore)
%plot idata (da) overlaid
axPlot=ax;
%reset(axPlot);      %clear all parameters, ie ylabel positions
plotsize=[650 150];        %plot size for close up
plotsize=plotParam.ncsplotsize;        %plot size for close up
fontsize=plotParam.fontsize;
argnum = 1;
dalabel='[DA] (nM)';
norm=0;     %flag if want z-score normalized signal
units=[];
scalebar=[];
xlab=0;         %disable xlabeling
filtlim=[];     %filt band
cscids=[];
ids=[];
labels={};
sqenv=0;        %flag to square/envelope data
env=0;          %flag to just env data
winlength=0;            %smoothing window
samplerate=parameters.sampleratencs;
data=processed.samplesNCS;
idata=[];       %fscv data, may not be generated yet from pca
fscvch=plotParam.selch;
holdax=0;
if isfield(processed,'Ipcr')
    for ich=1:length(fscvch)
        idata(ich,:)=processed.Ipcr{fscvch(ich)}.DAiso;
    end
end
usep=0;
bursts=0;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'norm'
            argnum = argnum + 1;
            norm=1; %normalize signals (no extra arg needed)
        case 'ncs'
            %plot ncs channels overlay, indicate csc #s
            argnum = argnum + 1;
            cscids = varargin{argnum};
        case 'ncsids'
            %plot ncs channels overlay, here give the id's of loaded
            %channels directly instead of csc # above
            argnum = argnum + 1;
            ids = varargin{argnum};      
        case 'fscvch'
            %ch #s to plot for fscv
            argnum = argnum + 1;
            fscvch = varargin{argnum};    
            if isfield(processed,'Ipcr')
                if isfield(processed.Ipcr,'DAiso')
                    idata=processed.Ipcr{fscvch}.DAiso ;
                else
                    idata=[];
                end
            end
        case 'units'
            argnum = argnum + 1;
            units = varargin{argnum};
        case 'scale'
            argnum = argnum + 1;
            scalebar = varargin{argnum};
        case 'filt'
            argnum=argnum+1;
            filtlim=varargin{argnum};
        case 'sqenv'
            %square and envelope data
            sqenv=1;  
            winlength=round(samplerate*.5/mean(filtlim));       %envelop over 1/2 cycle of targetd osc 
        case 'env'
            %env data
            env=1; sqenv=1;
            winlength=round(samplerate*.5/mean(filtlim)); 
        case 'winlength'
            %give default winlength if provided
            argnum=argnum+1;
            if ~isempty(varargin{argnum})
                winlength=round(varargin{argnum}*samplerate); %in seconds
            end
        case 'xlab'
            %put xlabel or not
            xlab=1;      
        case 'hold'
            %hold plot, do not clear previous axes
            holdax=1;
        case 'useprocessed'
            %use already processed data from guiprolfp in
            %processed.resampled & processed.resampledb
            usep=1;
        case 'bursts'
            bursts=1;
            
        otherwise
            error('badoption');
    end
    argnum = argnum + 1;
end

if isempty(ids) && ~isempty(cscids)
    ids=find(ismember(plotParam.NCSchannels,cscids)==1);  %channel ids to be processed and plotted here 
end

position=getpixelposition(axPlot);
set(axPlot, 'Units','Pixels','Position',  [position(1) position(2) plotsize(1) plotsize(2)]);
subplot(axPlot)
TS=processed.LFPread.LFPts;
zoomTS=plotParam.zoomTS;
idEvents=find(ismember(processed.LFPread.LFPeventTTL,plotParam.disp_events)==1);
relTSEvents=processed.LFPread.LFPeventTS(idEvents)-TS(1);
relTS=TS-TS(1);
idxStart=find(relTS<=zoomTS(1)); idxStart=idxStart(end);
idxEnd=find(relTS>=zoomTS(2)); 
eventcodes=processed.LFPread.LFPeventTTL(idEvents);
if isempty(idxEnd)
    idxEnd=length(TS);
else
    idxEnd=idxEnd(1);
end
relTS=TS(idxStart:idxEnd)-TS(1);
fscv_freq=parameters.samplerate;
time_window_fscv=(round(relTS(1)):1/fscv_freq:round(relTS(end)));
idx_fscv=round(time_window_fscv.*fscv_freq)-plotParam.t_start+2;
blinew=1;
hold(axPlot,'on')        
if ~isempty(ids)
    labels=plotParam.cscNames(ids);     %function automatically sorts ids
   % [~,sortid]=sort(ids);
   % labels=labels(sortid);      %recover original unsorted list
    yyaxis left
    if holdax==0
        cla(axPlot)
    end
    colormap=plotParam.colormaptab;
    if holdax==1
        colormap=plotParam.colormaptab(2,:);        %default 2nd color???
    end
    plotdata=[];
    ylims=[];
    for ii=1:length(ids)        
        if norm==1 && sqenv==0 && usep==0
            plotdata=zscore(data(idxStart:idxEnd,ids(ii)));
            ylims=[-3 4];
            units=plotParam.unitsz;
        elseif ~isempty(filtlim) || env==1 && usep==0
            if ~isempty(filtlim)
                plotdata=filterLFP(data(idxStart:idxEnd,ids(ii)).*1e6,samplerate,filtlim);
            else
                plotdata=data(idxStart:idxEnd,ids(ii));
            end
            units=plotParam.unitslfpfilt;
            ylims=scalebar;
            if sqenv==1
                %if want envelop for this channel, envelope
                if env==0
                    %not just enveloping, but need square
                    plotdata=plotdata.^2;   %get power V^2
                end
               plotdata=envwave(plotdata);     % 1/2019
                plotdata=smoothwin(plotdata,winlength);   %smoothing
                %  samples(ids(ii),:)=waveenv(samples(ids(ii),:));
                units=plotParam.unitslfpfiltenv;
                if norm==1 && env==0
                    plotdata(1:100)=0;
                    plotdata(end-100:end)=0;
                    plotdata=zscore(plotdata);
                    ylims=[-0.5 3];
                    %ylims=[-0.5 1.5];
                    units=plotParam.unitsz;
                elseif norm==1 && env==1
                    plotdata=zscore(plotdata);
                    ylims=[-3 4];
                    units=plotParam.unitsz;
                end
            end
        elseif usep==1
            %use processed data accoridng to choices
            if sqenv==1
                plotdata=processed.resampledb(ids(ii),idxStart:idxEnd).*1e12;
                units=plotParam.unitslfpfiltenv;
                ylims=scalebar;
            elseif sqenv==0
                plotdata=processed.resampled(ids(ii),idxStart:idxEnd).*1e6;
                units=plotParam.unitslfpfilt;
                ylims=scalebar;
            elseif norm==1 && sqenv==0
                plotdata=zscore(processed.resampled(ids(ii),idxStart:idxEnd));
                units=plotParam.unitsz;
                ylims=[-3 4];
            elseif norm==1 && sqenv==1
                plotdata=zscore(processed.resampledb(ids(ii),idxStart:idxEnd));
                units=plotParam.unitsz;
                ylims=[-3 4];
            end            
        else
            %just plot raw data filtered down to 300 hz
            plotdata=filterLFP(data(idxStart:idxEnd,ids(ii)).*1e6,samplerate,[.5 100]);
            ylims=scalebar;
            units=plotParam.unitslfp;
            blinew=3;
        end
        
        plot(axPlot,relTS,plotdata,'color',colormap(ii,:),'linestyle','-','marker','none')
        if bursts==1
            %plot bursts overlay
            idsb=processed.betaids{ids(ii)};
            colorburst=colormap(ii,:)-.3;
            colorburst(colorburst<=0)=0;
            idsplot=[];
            if ~isempty(idsb)
                for ib=1:length(idsb)
                    if iscell(idsb)
                    idsplot=idsb{ib}(idsb{ib}>=idxStart & idsb{ib}<=idxEnd)-idxStart+1;
                        plot(axPlot,relTS(idsplot),plotdata(idsplot),...
                        'color',colorburst,'linestyle','-','marker','none',...
                        'linewidth',blinew);
                    else
                        %single cell loaded
                      idsplot=idsb(idsb>=idxStart & idsb<=idxEnd)-idxStart+1;
                       plot(axPlot,relTS(idsplot),plotdata(idsplot),...
                        'color',colorburst,'linestyle','-','marker','none',...
                        'linewidth',blinew);
                    break;
                    end

                end
            end
            
        end
    end
    if ~isempty(ylims)
        ylim = set(axPlot,'ylim',ylims);
    end
end
if isempty(ids)
    dalabel=units;
end
yyaxis(axPlot,'right')
%Plot Dopamine signal on right axis (idata)
cla(axPlot)
if ~isempty(idata)
    if length(idx_fscv)>size(idata,2)
        idx_fscv=1:length(idata);
    end
plot(axPlot,time_window_fscv,idata(:,idx_fscv),'k')
plt = gca;
if bursts==1
    %plot bursts overlay
    idsb=processed.betaids{ids(1)};
    if ~iscell(processed.betaids{ids(1)})
        %single id case not cell, just array
        idsb={};
        idsb{1}=processed.betaids{ids(1)};
    end
    if ~isempty(idsb)
        iwindata=idata(:,idx_fscv);
        for ib=1:length(idsb)
            idsplot=idsb{ib}(idsb{ib}>=idxStart & idsb{ib}<=idxEnd)-idxStart+1;
            tsb=relTS(idsplot);
            if ~isempty(tsb)
            ftsb=find(time_window_fscv>=tsb(1) & ...
                time_window_fscv<=tsb(end));
            if ~isempty(ftsb)
            if length(ftsb)==1
                %extend window if only one fscv sample for visualization
                ftsb=ftsb-1:ftsb+1;
                ftsb=ftsb(ftsb>1 & ftsb<size(iwindata,2));
            end
            plotchs=find(~isnan(mean(iwindata(:,ftsb),2)));
            if ~isempty(plotchs)
                for ix=1:length(plotchs)
                plot(axPlot,time_window_fscv(ftsb),iwindata(plotchs(ix),ftsb),...
                    'color',[1 0 0],'linestyle','-','marker','none',...
                    'linewidth',2);
                end
            end
            end
            end
        end
    end
end


limdamax=nanmean(nanmean(idata(:,idx_fscv),2))+nanmean(nanstd(idata(:,idx_fscv),[],2))*2.5;
limdamin=nanmean(nanmean(idata(:,idx_fscv),2))-nanmean(nanstd(idata(:,idx_fscv),[],2))*2.5;
if ~any(isnan([limdamin limdamax]))
set(axPlot,'ylim',[limdamin limdamax]);
end
end
plt.YAxis(1).Color='r';
plt.YAxis(2).Color='k';
events=[];
if isfield(plotParam,'events')
    %combined fscv /ephys
events=plotParam.events;
else
    events=processed.LFPread.LFPeventTS;
end
YLim=get(axPlot,'ylim');
if ~isempty(events)
    relTSinwin=relTSEvents(relTSEvents>=relTS(1) &relTSEvents<=relTS(end));
    eventsinwin=eventcodes(relTSEvents>=relTS(1) &relTSEvents<=relTS(end));
    eventsbig=find(eventsinwin==45);
    eventssmall=find(eventsinwin==46);
    %color coded    
    if ~isempty(eventsbig)
        for ii=1:size(eventsbig,2)
            x1=relTSinwin(eventsbig(ii));
            plot(axPlot,[x1 x1],YLim,'color',[0 .7 .4],'linestyle','--','marker','none');
        end
    end
    if ~isempty(eventssmall)
        for ii=1:size(eventssmall,2)
            x1=relTSinwin(eventssmall(ii));
            plot(axPlot,[x1 x1],YLim,'color',[0 0 .5],'linestyle','--','marker','none');
        end
    end
    for ii=1:length(relTSinwin)
       if relTSinwin(ii)~=0
            x1=relTSinwin(ii);
            if ~ismember(x1,relTSinwin([eventssmall eventsbig]))
                plot(axPlot,[x1 x1],YLim,'color',[.5 .5 .5],'linestyle','-','marker','none');
            end
       end
    end
end
if xlab==1
    xlabel(axPlot,'time (s)','fontsize',fontsize)
end
aa=get(axPlot,'ylabel');
delete(aa);
y_h=ylabel(axPlot,dalabel,'fontsize',fontsize);  
ylp = get(y_h, 'Position');
ext=get(y_h,'Extent');
set(y_h, 'Rotation',270, 'Position',ylp+[ext(3) 0 0])

yyaxis(axPlot,'left')
set(axPlot, 'YScale', 'linear')
aa=get(axPlot,'ylabel');
delete(aa);
%if isempty(aa.String)
ylabel(axPlot,units,'fontsize',fontsize);
%end
set(axPlot,'fontsize',fontsize);
set(axPlot,'xlim',zoomTS);

hold(axPlot,'off')
end