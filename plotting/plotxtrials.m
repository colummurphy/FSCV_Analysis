function plotxtrials(hax,xdata,varargin)
%min info nneeded is xdata.tsx, rate, xcovda, sitename, eventnames,
%freqband
tsx=xdata.tsx;
rate=xdata.rate;
plotdata=xdata.xcovda;
if isempty(plotdata) || size(plotdata,1)<5
    return;
end
plotwin=0;
sitename=xdata.sitename;
eventnames=xdata.eventnames;
freqband=xdata.freqband;
fontsize=10;
plotleft=0;
argnum=1;
sortpos=0;
sortneg=0;
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'plotwin'
            argnum=argnum+1;
            plotwin=varargin{argnum};
        case 'leftlabel'
            plotleft=1;
        case 'sortpos'
            sortpos=1;
        case 'sortneg'
            sortneg=1;
        case 'groupxcov'
            %get grouped xcov data instead of all trials
            %according to supplied name in subseq arg
            argnum=argnum+1;
            groupname=varargin{argnum};
            plotdata=getfield(xdata.groupxcov,groupname);
            eventnames=[eventnames ' | ' groupname];
    end        
    argnum=argnum+1;
end
cla(hax);
interval=1;       %in seconds
xticklabels=min(tsx):interval:max(tsx);
xticklabels=round(xticklabels.*10)./10;
xticklabels=num2str(xticklabels');
xticks=1:round(interval*rate):length(tsx);
if plotwin~=0
    %plot select interval
    tplot=tsx(find(tsx>=-plotwin & tsx<=plotwin));
    xticklabels=min(tplot):interval:max(tplot);
    xticklabels=round(xticklabels.*10)./10;
    xticklabels=num2str(xticklabels');
    xticks=1:round(interval*rate):length(tplot);
    plotdata=plotdata(:,(tsx>=-plotwin & tsx<=plotwin));
end
if sortpos
[~, sortid]=sort(xdata.xcovlag(:,1));
plotdata=plotdata(sortid,:);
end
if sortneg
    [~,sortid]=sort(xdata.xcovlaganti(:,1));
    plotdata=plotdata(sortid,:);
end
imagetrials=image(hax, plotdata,'cdatamapping','scaled');
hold(hax,'on')
set(hax,'YDir','reverse')        %flip y-axis values so first trial on top
%organize plot
yints=round(size(plotdata,1)/10);      %row intervals
trialsort=sort(1:size(plotdata,1));
yticks=trialsort(1:yints:end);
set(hax,'XTick',xticks)
set(hax,'xticklabel',xticklabels)
set(hax,'tickdir','out','box','off')
xlabel(hax,'time (s)')
set(hax,'ylim',[1 size(plotdata,1)]);
set(hax,'xlim',[1 size(plotdata,2)]);
origpos=getpixelposition(hax);      %get original position 
if isempty(sitename)
    sitename{1}='';
    sitename{2}='';
end
if plotleft
set(hax,'ytick',yticks);
ylabel(hax,'trial #')
title(hax,[eventnames ' | xcov(lfp [' num2str(freqband(1))...
    '-' num2str(freqband(2)) ']' sitename{1} ' , DA ' sitename{2} ')'])
else
set(hax,'Ytick',[])
title(hax,['xcov(lfp [' num2str(freqband(1))...
    '-' num2str(freqband(2)) ']' sitename{1} ' , DA ' sitename{2} ')'])
end
clabel='[DA]-beta xvar';


matlabver=version('-release');
matlabver=str2num(matlabver(regexp(matlabver,'\d')));

if matlabver>2013
    %colorbar function introduced in matlab 2014
    h1=colorbar(hax,'southoutside');
    cpos = getpixelposition(h1);
    %ylabelbar=ylabel(h1,clabel); ypos = getpixelposition(ylabelbar);
    cpos(4) = 15; 
    set(hax,'Units','Pixels','Position', [origpos(1) origpos(2) origpos(3) origpos(4)]);

    set(h1,'Units','Pixels','Position', [cpos(1) origpos(2)-75 cpos(3)*2/3 cpos(4)]);
else
    set(hax,'Units','Pixels','Position', [origpos(1) origpos(2) origpos(3) origpos(4)]);
  h1=colorbar(hax,'eastoutside');  
end

%set c scale
climsmax=max(nanmean(plotdata,1))+median(nanstd(plotdata,[],1))*3;
climsmin=min(nanmean(plotdata,1))-median(nanstd(plotdata,[],1))*3;

if ~isempty(climsmin) && ~isempty(climsmax) && ~isnan(climsmax) && ~isnan(climsmin)
    set(hax,'clim',[climsmin climsmax])
end

set(findall(hax,'-property','FontSize'),'FontSize',fontsize)

hold(hax,'off')

end
