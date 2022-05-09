function plotxtraces(xinfo,trialinfo,plotparam,varargin)
%trialinfo must be given in same order as sesstypes
%plot bar scatter plots of properties of cross-covariance lags/waveforms
%12/31/2018 updated for timing characteristics from new xclust 
%same as plotx, but plot big/small/targ on same plots
%for plotting multiple sessions
%1/3/2018, udpates as in plotx/plotxallsess
%PLOT TRACES for DA GROUPS
rate=10;      %10hz default for downsampled da/lfp from xcov
win=[-2 2];     %+/-2 s from aln idx
interval=1;       %in seconds

argnum=1;
datype='dapos';
sortda=0;
plotz=0;
tstype='dapos';
tstraces={};
conditions={'all'};
sesstypes={'big','small'};
evt={};
lfpsites={};
sorttype='damax';
plotmarks=0;
sortrts=0;
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'daneg'
            datype='daneg';
        case 'dareb'
            datype='dareb';
        case 'daall'
            datype='daall';
            tstype='daall';
        case 'sort'
            %sort da increasing
            sortda=1;
            argnum=argnum+1;
            sorttype=varargin{argnum};
            if contains(sorttype,'rts')
                sortrts=1;
            end
        case 'win'
            %user specified window +/- s
            argnum=argnum+1;
            win=varargin{argnum};
        case 'plotz'
            plotz=1;           
        case 'plotmarks'
            plotmarks=1;
            tstraces{1}.tstrace={'damaxts'};
            tstraces{2}.tstrace={'lfppostmaxts','lfpmints'};
        case 'lfpmints'
            tstraces{2}.tstrace={'','lfpmints'};
        case 'condition'
            argnum=argnum+1;
            conditions=varargin{argnum};
        case 'sesstype'
            argnum=argnum+1;
            sesstypes=varargin{argnum};
        case 'evtype'
            argnum=argnum+1;
            evt=varargin{argnum};   %user provides aln event  
        case 'lfps'
            argnum=argnum+1;
            lfpsites=varargin{argnum};      %user provides lfp sites
    end
    argnum=argnum+1;
end
fontsize=12;
savepath=plotparam.savepath;
if sortda
    savepath=[savepath 'sorted' filesep];
end
lfpchs=plotparam.lfpchs;
if ~isempty(lfpsites)
    lfpchs=lfpsites;
end
%get lfp ch groups p & c
pgroup=find(contains(lfpchs,'p')==1);
cgroup=find(contains(lfpchs,'c')==1);
dasites={xinfo(1:end).siteda};
lfprows={xinfo(1:end).sitelfp};
%lfpchs=unique(lfprows);
targdasites=unique(dasites);
figpos=[50,50,1000,800];
figsess=figure('visible','off');     %figure for each channel
if ispc
figsess=figure('visible','on');     %figure for each channel
end
set(figsess,'position',figpos,'color',[1 1 1]);

set(0,'CurrentFigure',figsess);    %set figure handle to current figure
axa={};

eventrows={xinfo(1:end).event};
eventtypes=unique(eventrows);

sessionrows={xinfo(1:end).sessiontype};
sessiontypes=unique({xinfo(1:end).sessiontype});
axpos={};
axsiz=[750,320];
axoff=100;
sessiontypes=sesstypes;
if ~isempty(evt)
    eventtypes=evt;
end
sesslabel=sessiontypes;

for ida=1:length(targdasites)
    %each da separate subplot  
    lfpsites={};
    daregion=contains(targdasites(ida),'c');     %1=='c'
    if daregion==1
        lfpsites=lfpchs(cgroup);
    else
        lfpsites=lfpchs(pgroup);
    end
for ilfp=1:length(lfpsites)    
for ievent=1:length(eventtypes)
    scalesda=[];
    scaleslfp=[];
for itype=1:length(sessiontypes)
    %for each pair, event, plot da & lfp signals side by side trial by
    %trial for selected pos/neg/reb groups
     targrow=find((contains({xinfo.siteda},targdasites(ida)) & ...
        contains({xinfo.sitelfp},lfpsites(ilfp)) & ...
        strcmp({xinfo.event},eventtypes(ievent))) & ...
        contains({xinfo.sessiontype},sessiontypes(itype))==1);
    if ~isempty(targrow)
    %curdata=xinfo(targrow).dapos;
    curdata=getfield(xinfo(targrow),datype);
    alnidx=curdata.mididx;
    wins=[alnidx+win(1)*rate:alnidx+win(2)*rate];       %time window to plot
    da=curdata.datracesaln(:,wins);                             %signals aligned to aln idx
    lfp=curdata.lfptracesaln(:,wins);                             %signals aligned to aln idx
    scaleslfp(itype)=nanmean(nanstd(curdata.lfptracesaln(:,wins),[],2));
    scalesda(itype)=nanmean(nanstd(curdata.datracesaln(:,wins),[],2));
    if plotz
        scaleslfp(itype)=nan;
        scalesda(itype)=nan;
    end
    end
end   
for itype=1:length(sessiontypes)
    %for each pair, event, plot da & lfp signals side by side trial by
    %trial for selected pos/neg/reb groups
     targrow=find((contains({xinfo.siteda},targdasites(ida)) & ...
        contains({xinfo.sitelfp},lfpsites(ilfp)) & ...
        strcmp({xinfo.event},eventtypes(ievent))) & ...
        contains({xinfo.sessiontype},sessiontypes(itype))==1);
if ~isempty(targrow)
for icond=1:length(conditions)
    %curdata=xinfo(targrow).dapos;
    da=[];
    lfp=[];
    sortdata=[];
    tsdatat={};
    tsdata{1}=[];
    tsdata{2}=[];
    for irow=1:length(targrow)
        curdata=getfield(xinfo(targrow(irow)),datype);
        trialtypes=trialinfo(itype);
        targch=find(strcmp({trialtypes(1:end).site},targdasites(ida)));
        trialid=find(contains(trialtypes(targch).names,conditions{icond}));
        origtrialids=trialtypes(targch).nums{trialid};
        seltrials=find(ismember(curdata.trials,origtrialids));
        alnidx=curdata.mididx;
        wins=[alnidx+win(1)*rate:alnidx+win(2)*rate];       %time window to plot
        dat=curdata.datracesaln(seltrials,wins);                             %signals aligned to aln idx
        lfpt=curdata.lfptracesaln(seltrials,wins); 
        if sortrts
            sortdatat=getfield(xinfo(targrow(irow)),'trt');
            if strcmp(sorttype,'fixrts')
            sortdatat=getfield(xinfo(targrow(irow)),'frt');
            end
                
            origseltrialids=intersect(curdata.trials,origtrialids);
            rttrials=intersect(xinfo(targrow(irow)).goodtrials,origseltrialids);
            sortdatat=sortdatat(rttrials);
        else
        sortdatat=getfield(curdata,sorttype);
        sortdatat=sortdatat(seltrials);
        end
        if ~isempty(tstraces)
            for id=1:2
                %da & lfp
                tstrace=tstraces{id}.tstrace;
                tsdatat={};
                %plot overlying ts parameter
                for it=1:length(tstrace)
                    %max first row, min 2nd row
                    if ~isempty(tstrace{it})
                    dataforts=getfield(xinfo(targrow(irow)),tstype);
                    tsdatatemp=getfield(dataforts,tstrace{it})-alnidx-win(1)*rate+1;
                    tsdatat{id}(it,:)=tsdatatemp(seltrials);
                    end
                end
                tsdata{id}=[tsdata{id} tsdatat{id}];
            end
        end
        sortdata=[sortdata sortdatat];
        da=[da; dat];
        lfp=[lfp; lfpt];
    end
    if sortda
        [maxvals,sortids]=sort(sortdata);
        sortids=sortids(~isnan(maxvals));       %remove nan values (usually placed as max)
        da=da(sortids,:);
        lfp=lfp(sortids,:);
                tsdatatemp=tsdata;
        tsdata={};
        if ~isempty(tstraces)
            for id=1:2
                tstrace=tstraces{id}.tstrace;
                for it=1:length(tstrace)
                    %max first row, min 2nd row
                    tsdata{id}(it,:)=tsdatatemp{id}(it,sortids);
                end
            end
        end
    end
    if plotz
        %da=zscore(da,[],2);            %ZSCORE ONLY NON-NAN IDS
        for it=1:size(da,1)
            %one by one to ignore nan's
            nonnan=~isnan(da(it,:));
            daz=zscore(da(it,nonnan));
            da(it,nonnan)=daz;
        end
        lfp=zscore(lfp,[],2);
    end
    freqband=getfield(xinfo(targrow(1)),'freq');
    clf(figsess,'reset');
    set(figsess,'color',[1 1 1]);
    axpos={};
    pdata=[];
    clabel=[];
    for ip=1:2        
        axa{ip}=subplot(1,2,ip);   hold(axa{ip},'on');
        set(axa{ip},'units','pixels');
        axpos{ip}=get(axa{ip},'position');
        if size(da,1)>1
        if ip==1
            pdata=da;
            title(axa{ip},[datype ' | ' sesslabel{itype} ' | ' ...
                conditions{icond} ' | ' eventtypes{ievent} ' | da ' targdasites{ida}]);
            clabel='\DeltaDA (nM)';
            if plotz
                clabel='z-score DA';
            end
            ylabel(axa{ip},'trial #')
            if sortda
                ylabel(axa{ip}, [sorttype ' sorted trials']);
            end
        else
            pdata=lfp;
             title(axa{ip},[datype ' | ' sesslabel{itype} ' | ' ...
                conditions{icond} ' | ' eventtypes{ievent} ' | \beta-lfp ' lfpsites{ilfp}]);
             clabel=['\beta-LFP ' num2str(freqband(1)) ...
                 '-' num2str(freqband(2)) ' hz'];
        end  
        imagetrials=image(axa{ip}, pdata,'cdatamapping','scaled');
        set(axa{ip},'YDir','reverse')        %flip y-axis values so first trial on top
        
        xticklabels=min(win):interval:max(win);
        xticklabels=round(xticklabels.*rate)./rate;
        xticklabels=num2str(xticklabels');
        xticks=1:round(interval*rate):size(pdata,2);
        set(axa{ip},'xtick',xticks,'xticklabel',xticklabels);
        yints=round(size(pdata,1)/10);      %row intervals
        trialsort=sort(1:size(pdata,1));
        yticks=trialsort(1:yints:end);
        set(axa{ip},'tickdir','out','box','off')
        xlabel(axa{ip},'time (s)')
        set(axa{ip},'xlim',[0 size(pdata,2)]);
        set(axa{ip},'ylim',[1 size(pdata,1)]);
        origpos=getpixelposition(axa{ip});      %get original position 
        set(axa{ip},'ytick',yticks);
        h1=colorbar(axa{ip},'southoutside');  
        set(h1,'units','pixels');
        set(axa{ip},'Units','Pixels','Position', [axpos{ip}(1) axpos{ip}(2) axpos{ip}(3) axpos{ip}(4)]);
        set(h1,'position',[axpos{ip}(1) 30 100 10]);
        %climsmax=nanmean(pdata(mididx,:))+nanmean(nanstd(pdata,[],2))*3;
        %climsmin=nanmean(pdata(mididx,:))-nanmean(nanstd(pdata,[],2))*1.5;
        climsmax=nanmean(scalesda)*4;
        climsmin=-nanmean(scalesda)*1.5;
        if ip==2
           % climsmax=nanmean(pdata(mididx,:))+nanmean(nanstd(pdata,[],2))*4;
            %climsmin=nanmean(pdata(mididx,:))-nanmean(nanstd(pdata,[],2))*1;
            climsmax=nanmean(scaleslfp)*5;
            climsmin=nanmean(scaleslfp)*1;
            set(axa{ip},'ytick',[]);
        end
        %climsmax=max(nanmean(pdata,1))+median(nanstd(pdata,[],1))*3;
        %climsmin=min(nanmean(pdata,1))-median(nanstd(pdata,[],1))*3;
        if ~plotz
        if ~isempty(climsmin) && ~isempty(climsmax) && ~isnan(climsmax) && ~isnan(climsmin)
            set(axa{ip},'clim',[climsmin climsmax])
        end
        end
        if ~isempty(tstraces) && plotmarks
            tstrace=tstraces{ip}.tstrace;
            markc=[.75 .15 .75; 0 0 0];
            markc=[1 .75 1; 0 0 0];
            marka=[1 .75];
            markt={'s','o'};
            marks=[20,15];
            markl=[2 1.5];
            %plot overlying ts parameter
            for it=1:length(tstrace)
                if ~isempty(tstrace{it})
                dataforts=getfield(xinfo(targrow),tstype);
                tsdatatemp=getfield(dataforts,tstrace{it})-alnidx-win(1)*rate+1;
                tsdataplot=tsdata{ip}(it,:);
                meants=nanmean(tsdataplot);
                stdts=nanstd(tsdataplot);
                tsci=stdts./sqrt(length(tsdataplot))*1.96; 
                for jj=1:size(da,1)
                    %ama=plot(axa{ip},[tsdata(jj)],jj,'Color',markc(it,:),'LineWidth',2,'LineStyle',':','marker','s','markersize',4);           
                    scatter(axa{ip},tsdataplot(jj),jj,marks(it),markt{it},'markeredgecolor',markc(it,:),'MarkerEdgeAlpha',marka(it),'linewidth',markl(it));      
                end
                end
                
            end
        end
        set(findall(axa{ip},'-property','FontSize'),'FontSize',fontsize)
        set(h1,'fontsize',8);
        end
    end
savedir=[savepath 'allsess_traces'  filesep ...
    targdasites{ida} 'x' lfpsites{ilfp} filesep];
if ~isdir(savedir)
    mkdir(savedir);
end
savename=[savedir  datype '_' sesslabel{itype} '_' eventtypes{ievent} ...
    '_' conditions{icond} '_' targdasites{ida} 'x' lfpsites{ilfp}];
if plotz
    savename=[savedir  datype '_' sesslabel{itype} '_' eventtypes{ievent} ...
    '_' conditions{icond} '_' targdasites{ida} 'x' lfpsites{ilfp} '_' sorttype '_z'];
end
if sortda    
    savename=[savedir  datype '_' sesslabel{itype} '_' eventtypes{ievent} ...
    '_' conditions{icond} '_' targdasites{ida} 'x' lfpsites{ilfp} '_' sorttype];
end
if plotmarks
        savename=[savedir  datype '_' sesslabel{itype} '_' eventtypes{ievent} ...
    '_' conditions{icond} '_' targdasites{ida} 'x' lfpsites{ilfp} '_' sorttype '_marks'];
end
savefig(figsess,savename);
saveas(figsess,savename,'tif')
delete(findall(figsess,'type','text')) 
end
end
end
end
end
end


end


%set c scale


